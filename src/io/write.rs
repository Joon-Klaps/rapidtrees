//! Writing a distance matrix as TSV, optionally gzipped.

use flate2::{Compression, write::GzEncoder};
use std::io::{self, Write};
use std::path::Path;

/// Write a labeled square matrix as TSV to a file or stdout.
pub fn write_matrix_tsv<P: AsRef<Path>, T: std::fmt::Display + Sync>(
    path: P,
    names: &[String],
    mat: &[T],
    n_trees: usize,
) -> io::Result<()> {
    use std::fs::File;
    use std::io::BufWriter;

    let p = path.as_ref();
    if p.as_os_str() == "-" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "writing to stdout is not supported by write_matrix_tsv",
        ));
    }

    let is_gz = p.to_string_lossy().ends_with(".gz");

    let mut out: Box<dyn Write> = if is_gz {
        let f = File::create(p)?;
        let enc = GzEncoder::new(f, Compression::default());
        Box::new(BufWriter::new(enc))
    } else {
        Box::new(BufWriter::new(File::create(p)?))
    };

    write!(&mut out, "\t")?;
    for (k, name) in names.iter().enumerate() {
        if k > 0 {
            write!(&mut out, "\t")?;
        }
        write!(&mut out, "{}", name)?;
    }
    writeln!(&mut out)?;

    let rows_per_block = (BLOCK_BYTES / (n_trees.max(1) * CELL_BYTES)).max(1);
    write_rows(&mut out, names, mat, n_trees, rows_per_block)?;

    out.flush()?;
    Ok(())
}

/// Rough size of one formatted cell, tab included. It only sizes the blocks in
/// [`write_rows`], so a wrong guess costs memory or parallelism, never output.
const CELL_BYTES: usize = 16;

/// Formatted text [`write_matrix_tsv`] holds at once: about 400 rows of a
/// 10 000-tree matrix.
const BLOCK_BYTES: usize = 64 << 20;

/// Write the matrix body, `rows_per_block` rows at a time: each block's rows
/// are formatted in parallel, then written in order.
fn write_rows<T: std::fmt::Display + Sync>(
    out: &mut dyn Write,
    names: &[String],
    mat: &[T],
    n_trees: usize,
    rows_per_block: usize,
) -> io::Result<()> {
    use crate::par::*;

    if n_trees == 0 {
        return Ok(());
    }
    for (block, rows) in mat.chunks(n_trees * rows_per_block).enumerate() {
        let first = block * rows_per_block;
        let texts = rows
            .par_chunks(n_trees)
            .enumerate()
            .map(|(k, row)| -> io::Result<Vec<u8>> {
                let mut text = Vec::with_capacity(row.len() * CELL_BYTES);
                write!(text, "{}", names[first + k])?;
                for val in row {
                    write!(text, "\t{val}")?;
                }
                text.push(b'\n');
                Ok(text)
            })
            .collect::<io::Result<Vec<_>>>()?;
        for text in &texts {
            out.write_all(text)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{write_matrix_tsv, write_rows};
    use std::io::Read;

    /// The matrix body as the writer produced it before it went parallel: one
    /// row after another, one cell at a time.
    fn serial_body<T: std::fmt::Display>(names: &[String], mat: &[T], n: usize) -> String {
        let mut text = String::new();
        for (i, row) in mat.chunks(n).enumerate() {
            text.push_str(&names[i]);
            for val in row {
                text.push_str(&format!("\t{val}"));
            }
            text.push('\n');
        }
        text
    }

    fn tree_names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("tree_{i}")).collect()
    }

    #[test]
    fn rows_match_serial_formatting_across_block_boundaries() {
        let n = 7;
        let names = tree_names(n);
        let weighted: Vec<f64> = (0..n * n).map(|k| k as f64 / 3.0).collect();
        let rf: Vec<u32> = (0..n * n).map(|k| k as u32 * 2).collect();
        let (want_weighted, want_rf) = (
            serial_body(&names, &weighted, n),
            serial_body(&names, &rf, n),
        );

        for rows_per_block in [1, 2, 3, n, n + 5] {
            let mut got = Vec::new();
            write_rows(&mut got, &names, &weighted, n, rows_per_block).unwrap();
            assert_eq!(
                String::from_utf8(got).unwrap(),
                want_weighted,
                "f64, rows_per_block={rows_per_block}"
            );

            let mut got = Vec::new();
            write_rows(&mut got, &names, &rf, n, rows_per_block).unwrap();
            assert_eq!(
                String::from_utf8(got).unwrap(),
                want_rf,
                "u32, rows_per_block={rows_per_block}"
            );
        }
    }

    #[test]
    fn writes_header_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("matrix.tsv");
        let mat: Vec<u32> = vec![0, 4, 2, 4, 0, 6, 2, 6, 0];
        write_matrix_tsv(&path, &tree_names(3), &mat, 3).unwrap();
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "\ttree_0\ttree_1\ttree_2\ntree_0\t0\t4\t2\ntree_1\t4\t0\t6\ntree_2\t2\t6\t0\n"
        );
    }

    #[test]
    fn gz_output_decompresses_to_the_plain_text() {
        let dir = tempfile::tempdir().unwrap();
        let (plain, gz) = (dir.path().join("m.tsv"), dir.path().join("m.tsv.gz"));
        let names = tree_names(4);
        let mat: Vec<f64> = (0..16).map(|k| k as f64 * 0.125).collect();
        write_matrix_tsv(&plain, &names, &mat, 4).unwrap();
        write_matrix_tsv(&gz, &names, &mat, 4).unwrap();

        let mut decoded = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&gz).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, std::fs::read_to_string(&plain).unwrap());
    }
}
