// Vanilla driver for the rapidtrees WebAssembly module.
//
// Flow: read the uploaded BEAST `.trees` file as text → hand it to the wasm
// `compute_distances` → receive the matrix as a flat little-endian byte buffer
// → stream it to disk as CSV, one row at a time, so we never build a second
// multi-gigabyte copy in JS.

import init, { compute_distances } from "./pkg/rapidtrees.js";

const $ = (id) => document.getElementById(id);
const form = $("form");
const runBtn = $("run");
const statusEl = $("status");

function setStatus(msg, kind = "") {
  statusEl.textContent = msg;
  statusEl.className = `status ${kind}`;
}

// Initialise the wasm module up front so the first compute isn't slow.
init()
  .then(() => {
    runBtn.disabled = false;
    runBtn.textContent = "Compute & Save";
    setStatus("Ready.");
  })
  .catch((err) => {
    setStatus(`Failed to load WebAssembly: ${err}`, "error");
  });

// Quote a CSV field only when it contains a comma, quote, or newline.
function csvField(value) {
  const s = String(value);
  if (/[",\n\r]/.test(s)) {
    return `"${s.replace(/"/g, '""')}"`;
  }
  return s;
}

// Read one matrix value at flat index `k` from the byte buffer.
function valueReader(view, dtype) {
  if (dtype === "u32") {
    return (k) => view.getUint32(k * 4, true);
  }
  // f64
  return (k) => view.getFloat64(k * 8, true);
}

// Yield CSV lines: a header row of tree names, then one row per tree.
function* csvLines(names, view, dtype, n) {
  const at = valueReader(view, dtype);
  yield ["", ...names].map(csvField).join(",") + "\n";
  for (let i = 0; i < n; i++) {
    const cells = new Array(n + 1);
    cells[0] = csvField(names[i]);
    const base = i * n;
    for (let j = 0; j < n; j++) {
      cells[j + 1] = at(base + j);
    }
    yield cells.join(",") + "\n";
  }
}

// Stream CSV directly to a file handle (Chrome/Edge), gzipping on the fly.
async function streamToFileHandle(handle, lines, gzip) {
  const writable = await handle.createWritable();
  const encoder = new TextEncoder();

  let writer;
  let pipePromise;
  if (gzip) {
    const cs = new CompressionStream("gzip");
    pipePromise = cs.readable.pipeTo(writable);
    writer = cs.writable.getWriter();
  } else {
    writer = writable.getWriter();
  }

  try {
    for (const line of lines) {
      await writer.write(encoder.encode(line));
    }
    await writer.close();
  } catch (err) {
    await writer.abort?.(err);
    throw err;
  }
  if (pipePromise) await pipePromise;
}

// Fallback for browsers without showSaveFilePicker: build a Blob and download.
// This buffers the whole file in memory.
async function downloadViaAnchor(filename, lines, gzip) {
  const encoder = new TextEncoder();
  const chunks = [];
  for (const line of lines) chunks.push(encoder.encode(line));
  let blob = new Blob(chunks, { type: "text/csv" });
  if (gzip) {
    const stream = blob.stream().pipeThrough(new CompressionStream("gzip"));
    blob = await new Response(stream).blob();
  }
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

form.addEventListener("submit", async (event) => {
  event.preventDefault();

  const file = $("file").files[0];
  if (!file) {
    setStatus("Please choose a .trees file first.", "error");
    return;
  }

  const method = $("method").value;
  const rooted = $("rooted").checked;
  const burninTrees = parseInt($("burnin_trees").value, 10) || 0;
  const burninStates = parseInt($("burnin_states").value, 10) || 0;
  const useRealTaxa = $("use_real_taxa").checked;
  const gzip = $("gzip").checked;
  const suggestedName = `distances.csv${gzip ? ".gz" : ""}`;

  // The save dialog MUST be opened synchronously inside the user gesture
  // (this submit handler), before any `await`. If we computed first, the
  // gesture would have expired and showSaveFilePicker would throw. So we
  // acquire the file handle now; the file is only written to after the matrix
  // is computed below. Browsers without the API fall back to a download.
  let handle = null;
  if (typeof window.showSaveFilePicker === "function") {
    try {
      handle = await window.showSaveFilePicker({ suggestedName });
    } catch (err) {
      if (err && err.name === "AbortError") {
        setStatus("Save cancelled.");
        return;
      }
      setStatus(`Error: ${err && err.message ? err.message : err}`, "error");
      return;
    }
  }

  runBtn.disabled = true;
  try {
    setStatus("Reading file…");
    const content = await file.text();

    setStatus("Computing distance matrix… (this can take a while)");
    // Yield to the event loop so the status text paints before we block.
    await new Promise((r) => setTimeout(r, 0));

    const result = compute_distances(
      content,
      method,
      rooted,
      burninTrees,
      burninStates,
      useRealTaxa,
    );

    const names = result.names;
    const n = result.n;
    const dtype = result.dtype;
    const matrix = result.take_matrix();
    const view = new DataView(matrix.buffer, matrix.byteOffset, matrix.byteLength);

    setStatus(`Computed ${n}×${n} matrix. Saving ${suggestedName}…`);

    const lines = csvLines(names, view, dtype, n);

    if (handle) {
      await streamToFileHandle(handle, lines, gzip);
    } else {
      await downloadViaAnchor(suggestedName, lines, gzip);
    }

    setStatus(`Done — saved ${n}×${n} ${dtype} matrix to ${suggestedName}.`, "ok");
  } catch (err) {
    setStatus(`Error: ${err && err.message ? err.message : err}`, "error");
  } finally {
    runBtn.disabled = false;
  }
});
