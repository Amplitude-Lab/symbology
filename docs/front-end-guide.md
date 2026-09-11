# Operating the front-end — an illustrated guide

Symbology Studio is the browser interface for building and running symbol
computations: you describe a flow as a graph of operations, the server
compiles it into concrete `bootstrap` / `tensor_ops` / `tensor_add` /
Wolfram steps, runs them with verification and caching, and you inspect the
exact tensors that come out. This guide walks through the whole loop with
short screen recordings. Setup instructions (build the C++ core, start the
server) are in the [main README](../README.md#the-visual-front-end-web-editor);
everything below assumes the app is open at `http://127.0.0.1:8321`.

> **The five-minute tour:** create a project from a built-in example →
> open its flow → **Compile** → **▶ Run this plan** → watch the steps go
> green → **Results** → **Inspect** a tensor → (later) **⇪ Export
> standalone script** to run the same flow on a cluster.

## 1. Create a project from a built-in example

![Creating a project from the E6 template](gifs/create-project.gif)

Click **Create your first project** (or **+ New** in the top bar at any
time), give the project a name, and pick a template from the dropdown:

- **4-Point Form Factor (93 letters)** — ships with a ready-made alphabet
  and precomputed condition tensors;
- **Pentagon** / **E6** — the example bootstrap workflows (the recording
  uses E6);
- *Empty project* — start from scratch.

After creation you land on **Materials**, the project's alphabet library:
letters, computed properties (integrability, first/last entry, symmetry
transformations), and the tensors shipped or computed for them. Properties
show a **Compute** action while pending; everything else in the app reads
from here.

## 2. Compile and run a flow

![Compiling and running the E6 example flow](gifs/compile-and-run.gif)

Open **Flows** and click the example flow. The editor has three parts: a
palette of node types on the left, the graph canvas in the middle, and an
inspector for the selected node on the right. Connect nodes by dragging
from an output port to a compatible input port.

Then:

1. Click **Compile** in the toolbar. The server validates the graph and
   shows the concrete step list in the **Plan** panel on the right — each
   step names the binary it will run and the outputs it will produce.
2. Click **▶ Run this plan**. You are taken to the run page, where steps
   turn green as they finish and the full log streams below. The log is
   also saved to the project as `runs/<run_id>.log`.
3. Done — every output the flow declares now exists in the project's
   `output/` directory.

A step is only skipped when its outputs **and** a stored fingerprint of
every input match — editing anything upstream automatically re-computes
the affected chain on the next run. One run may be active per project;
concurrent attempts get a clear message instead of a race.

## 3. Inspect results

![Inspecting a produced tensor](gifs/inspect-results.gif)

The **Results** page lists every tensor file in `output/` (and `data/`).
Click **Inspect** on any file to load its dimensions, non-zero count, and
a sample of entries. The summary calls Wolfram, so it can take a few
seconds for large tensors. This is the quickest way to sanity-check a
shape (e.g. the `28 × 7 × 42` FEC tensor above) before building on it.

## 4. Export a standalone script for a cluster

![Exporting the flow as a standalone script](gifs/export-script.gif)

Design on your laptop, compute anywhere the C++ core builds: with the
plan compiled, click **⇪ Export standalone script** in the Plan panel.
The server writes `exported/<flow>.sh` into the project directory and
shows the invocation to copy. On the target machine:

```bash
SYMBOLOGY_ROOT=/path/to/symbology \
PROJ_DIR=/path/to/synced/project-dir \
bash /path/to/project-dir/exported/<flow>.sh
```

The exported script embeds the same verification and caching as the
local runner, needs Bash and Python but no front-end dependencies, and
needs no Mathematica when the Wolfram-dependent tensors are already in
the project. Details, environment variables (`STEPS`, `STEP_TIMEOUT`,
`WOLFRAM_MODE`, …) and cluster prerequisites are in the
[README's export section](../README.md#design-locally-run-on-the-cluster-standalone-script-export).

## Everyday notes

- **Saving.** The editor autosaves ~0.8 s after your last edit (toggle in
  the toolbar); the indicator next to it shows *saved / unsaved / saving /
  error*. On a save conflict (the flow changed elsewhere, e.g. another
  tab) the server refuses the overwrite and offers **Download draft** /
  **Reload saved flow** — your draft is never silently discarded.
- **Wolfram.** Mathematica (`wolframscript`) is needed only for computing
  alphabet properties and for the Results tensor summaries. Drawing,
  compiling, running native steps, and exporting scripts work without it.
- **Run history.** The **Runs** tab lists every run with its status;
  finished runs survive server restarts. Cancellation is immediate and
  kills the whole process group of the running step.
- **Where files live.** A project is a directory under
  `front-end/projects/<id>/` — `data/` (inputs and shipped tensors),
  `output/` (results), `runs/` (logs and run records), `exported/`
  (standalone scripts). You can copy or sync the whole directory to move
  work between machines.
