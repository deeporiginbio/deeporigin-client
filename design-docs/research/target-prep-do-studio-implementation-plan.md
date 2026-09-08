# Target Preparation in DO Studio — implementation plan

**PRD:** [Target Preparation Tool Workflow and Application](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1063289062/Target+Preparation+Tool+Workflow+and+Application) (epic DDOS-7091)
**Child PRD:** [Structure Report Tool](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1066631169/Structure+Report+Tool)

**Primary sources**

- `do-dd-client` — source of truth for the tool contracts:
  `src/drug_discovery/protein_prep.py`, `src/drug_discovery/structure_report.py`,
  `src/drug_discovery/pocket_finder.py`, `src/platform/constants.py`,
  `src/platform/executions.py`, `tests/fixtures/executions/protein-prep-*.json`
- `platform-ui` — the App Engine that Studio apps are built on:
  `apps/uui/src/app-engine/**`, `apps/uui/src/app-schemas/**`,
  `packages/molstar/**`

Scope of this document is **UI**: what DO Studio has to build, which tools-API calls it
makes at each step, and which data-platform records get read or written. Backend asks are
called out explicitly where the UI cannot proceed without them.

---

## 1. What already exists

### 1.1 `deeporigin.protein-prep` (version `latest`, major `2`)

One tool, two operations discriminated by `inputs.action`. Not billable — there is no
quote path (`ProteinPrep` has no `run(quote=True)`).

**Execution endpoint** (same for every tool):

```
POST /tools/{orgKey}/tools/deeporigin.protein-prep/{major}/executions
body: { inputs, outputs: {}, metadata: {}, sync: <bool>, clusterId, projectId?, name? }
```

Note `sync` is a **top-level body key**, not an input — see `_make_protein_prep_payload`
(`src/drug_discovery/protein_prep.py:1030`). The protein-prep JSON Schema has no `sync`
property, so it must not be sent inside `inputs`.

**`action: "recommend"`** — inventories the structure. Always run with `sync: true`.

```jsonc
// inputs
{ "action": "recommend", "protein": { "file_path": "entities/proteins/<hash>.pdb", "id": "<protein_id>" } }
```

```jsonc
// jobOutputs.recommendation
{
  "source_sha256": "aaaa…",           // digest of the uploaded structure bytes
  "analyzer_version": "1.0.0",
  "chain_id_mapping": {},
  "components": [
    { "id": "chain:A",           "kind": "chain",  "subtype": "protein",
      "label": "Chain A", "recommendation": "keep",
      "reason": "Ordinary protein chain", "reason_code": "ordinary_protein_chain",
      "author": { "chain_id": "A" } },
    { "id": "ligand:LIG:A:100",  "kind": "ligand", "subtype": "small_molecule",
      "label": "LIG", "recommendation": "review",
      "reason": "Ambiguous ligand", "reason_code": "ambiguous_ligand",
      "author": { "chain_id": "A", "resname": "LIG", "resseq": 100 } },
    { "id": "water:A:HOH:310:",  "kind": "water",  "subtype": "crystal",
      "recommendation": "skip", "evidence": {} },
    { "id": "water:A:HOH:311:",  "kind": "water",  "subtype": "coordinating",
      "recommendation": "keep", "evidence": { "n_coord": 1 } }
  ]
}
```

- `kind` ∈ `chain | ligand | cofactor | water` (`PROTEIN_PREP_COMPONENT_KINDS`)
- `recommendation` ∈ `keep | review | skip` — the analyzer's frozen tag
- `author` carries the PDB addressing (`chain_id`, `resname`, `resseq`) — **this is what
  the Mol\* keep/exclude renderer selects on.**

**`action: "prepare"`** — applies the caller's resolved decisions, models loops, protonates.

```jsonc
// inputs
{
  "action": "prepare",
  "protein": { "file_path": "…", "id": "…" },
  "selection": {
    "source_sha256": "aaaa…",        // must match the recommend digest
    "analyzer_version": "1.0.0",
    "decisions": { "chain:A": "keep", "ligand:LIG:A:100": "skip", … }
  },
  "model_missing_loops": false,      // omitted when true (the default)
  "pdb_id": "1EBY"                   // REQUIRED when model_missing_loops is true
}
```

```jsonc
// jobOutputs.protein
{ "protein_pdb_file_path": "…", "pdb_id": "1EBY", "protein_id": "brd",
  "ph": 7.4, "force_field": "amber14", "protonate_protein": true }
```

Also written to result-explorer as `result_type: "preparedprotein"` with
`data.protein_pdb_file_path` (`protein_prep.py:1318`).

Hard rules the UI has to honour:

- Every `review` decision must be resolved to `keep` or `skip` before prepare, or the
  tool rejects the payload (`_validate_for_submit`).
- The `selection` is **digest-bound**: if the protein file changes, `source_sha256` no
  longer matches and prepare fails. Re-running recommend is the only fix.
- `model_missing_loops: true` requires a 4-character `pdb_id`. Proteins in the project
  without a `pdb_id` cannot use loop modelling.
- Loops-off prepare may be synchronous (`sync: true`); loops-on must be async.

### 1.2 `deeporigin.structure-report` (version `latest`)

Synchronous, billable (`approveAmount` supported). Takes a protein, a PDB ID, or both.

```jsonc
// inputs (at least one of the two)
{ "protein": { "file_path": "…", "id": "…" }, "pdb_id": "6GOG" }
```

```jsonc
// jobOutputs.structure_reports[0]
{
  "grade": "A", "weighted_score": 0.83,
  "resolution_score": …, "coverage_score": …, "rfree_score": …,
  "inhibitor_score": …, "method_score": …, "organism_score": …,
  "metadata_source": "rcsb" | "file_header" | "file_header+rcsb",
  "field_status": { "<field>": "value" | "not_applicable" | "unknown" },
  "coverage": 0.826, "has_ligand": true,
  "method": "X-RAY DIFFRACTION", "method_class": "x-ray",
  "organism": "Homo sapiens", "organism_class": "human",
  "resolution": 2.31, "rfree": 0.24,
  "pdb_id": "6GOG", "protein_id": "…", "source_sha256": "…"
}
```

This maps 1:1 onto the PRD's report card: grade badge `A`, the pills
`X-ray` / `Human` / `2.31Å` / `Ligand`, and the `6GOG · 1,036 Residues · 82.6% Coverage`
line. `field_status` says which pills to render as unknown rather than fabricating a value.
**Residue count is not in the payload** — it comes from `proteins.protein_length` on the
entity row, or has to be added to the tool output (see §7).

### 1.3 `deeporigin.pocket-finder` (version `2`)

Billable, async. Auto-find mode:

```jsonc
{ "protein": { "file_path": "…", "id": "…" }, "pocket_count": 5, "pocket_min_size": 30 }
```

Results land as `result_type: "pocket"` rows and are already rendered by the existing
`renderStructureAndPockets` Mol\* renderer (see `pocket-finder.json`).

### 1.4 The prepared-protein stamp

`src/drug_discovery/structures/prepared_protein_stamp.py` writes
`REMARK  99 DO_PREPARED` (PDB) / `_deeporigin.prepared` (mmCIF). Downstream Pocket Finder,
Docking and System Prep skip their AUTO cleanup when they see it. The Target Prep output
must carry the stamp, and the Mol\* viewer already surfaces it
(`packages/molstar/src/components/prepared-badge.tsx`).

---

## 2. The gap between the PRD and the shipped tools

The PRD asks for one button that runs, *in this order*:

1. Add missing atoms & residues
2. Loop Builder
3. Protonation
4. Pocket Finder
5. Structure Report

Steps 1–3 are all inside `protein-prep`'s `prepare` action today. Steps 4 and 5 are
separate tools. So the PRD workflow is **three tool executions**, and the PRD's four
"Advanced Parameters" checkboxes do not all have inputs behind them:

| PRD checkbox | Backing input today | Status |
| --- | --- | --- |
| Add missing atoms & residues | none — always on inside `prepare` | **needs a tool input** |
| Add missing loops | `model_missing_loops` | exists |
| Protonate (pH 7.4) | none — always on; `jobOutputs.protein.protonate_protein` only echoes it | **needs a tool input** |
| Find Pockets | n/a — separate tool | **needs workflow step gating** |

---

## 3. Architecture decision: one workflow tool, one app

The App Engine is deliberately **1 app ↔ 1 backend tool**
(`apps/uui/CLAUDE.md`, "Backend Tool Integration"). `manifest.steps` exists for
*pre-submit, synchronous* intermediate tools only (`useRunStep`) — it cannot chain
long-running executions after submit. There is no post-submit orchestration in the engine,
and building one would mean the UI owns retry, partial failure and three Activity rows for
one user action.

**Recommendation: a backend workflow tool `deeporigin.target-prep`**, modelled exactly on
`deeporigin.abfe-end-to-end` (which already chains `system-prep` → `abfe` via a `steps`
array; see `apps/uui/src/app-schemas/abfe.json`). One execution, one Activity row, one
billing transaction, one results view.

```jsonc
// deeporigin.target-prep inputs (proposed)
{
  "protein":   { "file_path": "…", "id": "…" },
  "selection": { "source_sha256": "…", "analyzer_version": "…", "decisions": { … } },
  "pdb_id": "6GOG",
  "add_missing_atoms": true,
  "model_missing_loops": true,
  "protonate": true,
  "find_pockets": true,
  "output_name": "Default_TargetPrep_ProteinId"
}
```

```jsonc
// jobOutputs (proposed)
{
  "protein": { "protein_pdb_file_path": "…", "protein_id": "<new proteins row id>", … },
  "pockets": [ … ],
  "structure_reports": [ { "grade": "A", … } ]
}
```

The `selection` is produced in the UI from a **`recommend` step run before submit** — that
part *does* fit `manifest.steps`, because recommend is synchronous and cheap.

Everything below assumes this shape. If the workflow tool slips, §10 has a phase-1 fallback
that ships the app against `deeporigin.protein-prep` alone.

---

## 4. End-to-end flow, with the exact calls

### Step 0 — route and shell

`/target-prep` renders `AppPage` → `AppEngineProvider` → mosaic layout + form sidebar.
No API calls beyond the manifest fetch (`/app-manifests/deeporigin.target-prep/{major}.json`
in staging/prod; the bundled `src/app-schemas/target-prep.json` locally and on PR previews).

### Step 1 — proteins table

`TableWrapper` with `config.entity: "proteins"`, `singleSelection: true`.
Reads `POST /data-platform/{orgKey}/entities/proteins_with_results/search` through the
existing server-side row model (`use-server-datasource.ts`). **No new data-platform work**
— the same table pocket-finder and ABFE use.

The PRD screenshot shows a `Docking Score` column; that is the existing
`showsResults: true` results-column machinery, nothing new.

### Step 2 — row selected → two synchronous steps auto-run

`TableWrapper.handleSelectionChanged` → `engine.setSelectedProteins([row])` → Zustand.
Two `manifest.steps` entries fire via hidden `ActionField`s with `autoRun: true`
(`apps/uui/src/app-engine/form-engine/fields/action-field.tsx`). `useRunStep` handles both.

**2a. Structure Report (initial)**

```
POST /tools/{orgKey}/tools/deeporigin.structure-report/{major}/executions
body: { inputs: { protein: {id, file_path}, pdb_id? }, outputs: {}, clusterId,
        visibility: "hidden" }
bind: { "structure_reports": "structureReport" }
```

**2b. Protein Prep recommend**

```
POST /tools/{orgKey}/tools/deeporigin.protein-prep/{major}/executions
body: { inputs: { action: "recommend", protein: {id, file_path} }, outputs: {},
        clusterId, sync: true, visibility: "hidden" }
bind: { "recommendation": "recommendation" }
```

Both land in `useAppStore.stepOutputs`. `useSelectionResets` already clears `stepOutputs`
whenever the entity selection changes, and `ActionField` re-arms per selection — so picking
a different protein re-runs both with no extra wiring.

`visibility: "hidden"` keeps these out of the user's Activity page
(`AppStep.visibility`, already supported).

> **Engine gap #1** — `useRunStep` builds its body as `{ inputs, outputs, clusterId, … }`
> and never sends a top-level `sync`. Protein-prep recommend needs `sync: true` at the top
> level. Add `AppStep.sync?: boolean` and spread it into the body
> (`apps/uui/src/app-engine/hooks/use-run-step.tsx:135`). Precedent for a top-level `sync`
> already exists in `renderer/table-wrapper/use-export-dataset.ts:109`.

### Step 3 — Structure Report tile renders

New tile reads `stepOutputs.structureReport[0]` reactively and renders the grade card.
No API call.

### Step 4 — Structure Filtering tile renders

New tile reads `stepOutputs.recommendation`, seeds a local decision map from each
component's `recommendation`, and renders one Keep/Skip row per component ordered
chain → ligand → water → cofactor (matching the PRD screenshot). Every edit writes the
full Selection object back to the store:

```ts
engine.setStateValue('selection', {
  source_sha256:    recommendation.source_sha256,
  analyzer_version: recommendation.analyzer_version,
  decisions:        { 'chain:A': 'keep', 'ligand:XO4:A:1': 'keep', … },
});
```

> **Engine gap #2** — `stepOutputs` is currently written only by `useRunStep`. Expose a
> generic `engine.setStateValue(key, value)` (thin wrapper over the store's existing
> `setStepOutputs`) on the `AppEngine` interface so a tile can contribute state that
> `x-from-state` picks up. This stays tool-agnostic — the engine never learns what a
> "selection" is.

`review` components must render in a visually distinct state and must **not** be silently
coerced. The Run button stays disabled while any decision is `review`.

> **Engine gap #3** — `manifest.parameters.requiredEntities` can only gate on entity
> selections. Add a `SubmitValidation` kind for "no unresolved value in a state key", e.g.
> `{ type: 'state-resolved', stateKey: 'selection', path: 'decisions', disallow: 'review',
> message: 'Resolve every component before running Target Preparation' }`, evaluated in
> `engine.submit` alongside the existing `disjoint` rule.

### Step 5 — Mol\* viewer colours keep vs exclude

`ProteinViewerWrapper` with a new renderer. Structure file is fetched from UFA the same way
pocket-finder does it (`use-structure-file.ts`). Component → Mol\* selection uses the
`author` block from each recommendation component:

```ts
MS.struct.generator.atomGroups({
  'chain-test':   MS.core.rel.eq([MS.ammp('auth_asym_id'), author.chain_id]),
  'residue-test': MS.core.rel.eq([MS.ammp('auth_seq_id'), author.resseq]),   // ligand/water/cofactor
})
```

then `plugin.builders.structure.tryCreateComponentFromExpression(...)` +
`buildUniformColor(value)`. Both primitives already exist
(`packages/molstar/src/api/loaders.ts:448`, `packages/molstar/src/utils/color-themes.ts`).

PRD palette, hex → the `number` that `buildUniformColor` takes:

| Component | Keep | Exclude |
| --- | --- | --- |
| Protein chain | `#2563eb` | `#bfdbfe` |
| Ligand | `#f97316` | `#fed7aa` |
| Co-factor | `#9333ea` | `#e9d5ff` |
| Water | `#06b6d4` | `#cffafe` |

The viewer re-renders on every decision toggle, so the renderer must be cheap to re-apply —
build the components once per structure load and only swap the colour theme on toggle.

### Step 6 — sidebar

Driven by `manifest.parameters.sections` (`AppFormEngine`). Three sections, matching the
screenshot:

1. **Structure Components** — read-only chips of the currently-kept components
   (`Chain A`, `Ligand XO4`, `Mn2+`), or `None` before a selection. Reads
   `stepOutputs.selection` + `stepOutputs.recommendation` for labels.
2. **Advanced Parameters** (collapsible, open) — four `boolean` fields:
   `add_missing_atoms`, `model_missing_loops`, `protonate`, `find_pockets`. `helpText`
   carries the sub-labels ("Fills incomplete residues", "Assign states @ pH 7.4", …).
   `model_missing_loops` must be disabled with an explanatory hint when the selected
   protein has no `pdb_id` — the tool rejects loops-on without one.
3. **Run Details** — `output_name` string field (the screenshot's "Output Property Name",
   placeholder `Default_TargetPrep_ProteinId`) plus the standard `run_name` with
   `x-body-key: "name"`.

> **Engine gap #4** — a read-only "chips from state" field type does not exist. Add a
> `component-summary` field (or generalise `entity-ref` to a read-only multi-value display).
> Small; alternative is to drop the chips from the sidebar and show them in the filtering
> tile only, which loses PRD parity.

### Step 7 — Run Target Preparation

`engine.submit()` → `useSubmitExecution`. `buildToolPayload` walks the manifest
`inputSchema`; `stateValues` already includes everything in `stepOutputs`, so `selection`
resolves through `x-from-state` with no engine change.

```
POST /tools/{orgKey}/tools/deeporigin.target-prep/{major}/executions
body: {
  inputs: {
    protein:  { id, file_path },              // x-data-type: Protein
    selection: { … },                          // x-from-state: "selection"
    pdb_id, add_missing_atoms, model_missing_loops, protonate, find_pockets, output_name,
  },
  outputs: {}, clusterId, projectId, name: "<run_name>"
}
```

Billing: the existing quote → `price-confirmation-panel.tsx` →
`insufficient-funds-modal.tsx` path applies unchanged, assuming the workflow tool quotes
like ABFE does.

### Step 8 — results mode

Route `/activity/{executionId}` → `appMode: 'results'`. Use a dedicated `resultsLayout` +
`resultsComponents` (supported today, see `resolveManifestView` in `init-app.tsx`):

- Structure Report tile — final report from the run's results
- Mol\* viewer — prepared structure + pockets, renderer `renderStructureAndPockets`
- Pockets table

---

## 5. Data platform: reads and writes

| Record | Direction | Who | New? |
| --- | --- | --- | --- |
| `proteins` / `proteins_with_results` rows | read | Table tile, existing datasource | no |
| Structure file bytes (UFA `file_path`) | read | Mol\* viewer | no |
| `results__pocket` rows (`result_type: "pocket"`) | write | pocket-finder step of the workflow tool | no |
| `results__preparedprotein` rows | write | protein-prep step | no (already emitted) |
| **`proteins` row for the prepared structure** | write | **needs an owner — see below** | **yes** |
| **`results__structurereport` rows** | write | structure-report tool | **yes** |
| `executions` row | write | tools-service, automatic | no |

Two genuine data-platform asks:

**(a) The prepared protein must become a `proteins` entity row.** Today the CLI is explicit
that it does not create one — `ProteinPrep.get_results()` returns an in-memory `Protein`
with `id is None` and only `remote_path` set; the caller has to `sync()` or `update()`
(`src/drug_discovery/protein_prep.py:1280`). But the whole point of Target Prep is to hand
a cleaned structure to Docking / ABFE / HTVS, and those apps' tables read `proteins`. So
one of:

- **Preferred** — the `target-prep` workflow tool registers the prepared protein as a new
  `proteins` row (carrying `pdb_id`, `protein_name` from `output_name`, `project_id`,
  and the `DO_PREPARED` stamp in the file), mirroring how docking indexes poses. The UI
  then needs nothing.
- Fallback — the UI POSTs `/data-platform/{orgKey}/proteins` after the run completes.
  Rejected: it puts a write on a client that may be closed before the async run finishes.

**(b) Structure Report needs indexed result rows.** Today the report exists only in
`jobOutputs` — `StructureReport` in the client parses `jobOutputs.structure_reports`
directly and there is no `_RESULT_TYPE_*` for it (compare `pose`, `pocket`,
`preparedsystem`, `abferesult`, `preparedprotein`). That is fine for the in-app card, but
the PRD's first requirement — *"a Structure Report should be generated and displayed
automatically"* for each structure in the table, and the child PRD's framing of the report
as a per-protein statistic — implies grade/resolution/coverage as **columns on the proteins
table**. That needs a `result_type: "structurereport"` joined to `protein_id`, at which
point the manifest's `results.detailColumns` / `results.aggregates` machinery renders them
for free (exactly as `pocket-finder.json` does with `min_volume` / `max_druggability`).

Without (b), each row's report can only be fetched on selection, one protein at a time —
which is what §4 step 2a does, and is enough for v1.

---

## 6. New and changed UI components

### New tiles (`apps/uui/src/app-engine/renderer/`)

| Component | Purpose | Registration |
| --- | --- | --- |
| `StructureFiltering` | Keep/Skip list per component, writes `selection` to store | `renderer/registry.ts`, `types/tiles.ts` `LeafComponent`, `renderer/types.ts` config type |
| `StructureReportCard` | Grade badge + pills + score breakdown | same three files |

Both follow the `ComponentWrapperProps` contract (`tileId`, `config`, `title`; state via
`useAppStore` / `useEngine`, never props). `AdmetScoreCard` is the closest existing model
for a read-only card fed by tool results.

### New Mol\* renderer (`packages/molstar/`)

- `src/api/components.ts` — `renderStructureComponents(plugin, proteinContent, format, components)`
- add `'renderStructureComponents'` to `RendererName` (`src/types/index.ts:243`)
- add the case to the `render()` dispatcher (`src/api/index.ts:526`)
- export the keep/exclude palette as constants so the tile and the viewer agree

### New form field

- `component-summary` — read-only chips (engine gap #4)

### Engine changes

| Gap | File | Change |
| --- | --- | --- |
| #1 | `hooks/use-run-step.tsx` | `AppStep.sync?: boolean` → top-level body key |
| #2 | `engine-provider.tsx`, `engine-context.ts` | expose `setStateValue(key, value)` |
| #3 | `engine-provider.tsx` submit path, `types/app.ts` | `SubmitValidation` kind for unresolved state |
| #4 | `form-engine/form-field.tsx`, `types/form.ts` | `component-summary` field type |

None of these introduce tool-specific logic into the engine — that invariant holds.

### Manifest and registration

- `apps/uui/src/app-schemas/target-prep.json` (sketch in §8)
- register in `app-schemas/index.ts` under key `target-prep`
- nav entry in `packages/global-provider/src/containers/subscription.container.tsx`
- app card in `apps/uui/src/pages/applications/constants.ts` + card image asset
- pick an `identityHue` ≥10° from existing tools (current cluster ~197–323; pocket-finder
  is 185)
- tool display name in `components/data-platform-tables/job-manager-table/index.tsx`
- add `deeporigin.target-prep` to `ALL_RESULTS_TOOL_KEYS` (`hooks/use-manifest.tsx:14`)
  once (b) above lands, so the column manager offers its results
- publish the manifest to S3 via `.github/workflows/register-tool-manifest.yml`
  (major version only; `latest` for the first release)

---

## 7. Loose ends the UI can't answer alone

1. **Residue count / coverage denominator.** The PRD card shows `1,036 Residues`.
   `structure_reports` returns `coverage` (a fraction) but no residue count. Either add it
   to the tool output or read `proteins.protein_length` off the entity row.
2. **`pdb_id`-less proteins.** Loop modelling is impossible without one. Decide whether the
   UI hard-disables the checkbox (proposed) or the workflow tool resolves a template
   another way.
3. **Cofactor labelling.** The screenshot names `Mn2+` and `Zn2+` chips. The recommend
   payload gives `kind: "cofactor"` and a `label`; confirm the label is already the ion
   name with charge, or format it client-side from `author.resname`.
4. **Re-prep of an already-prepared protein.** The `DO_PREPARED` stamp means downstream
   tools skip cleanup; running Target Prep on a prepared protein should probably warn.
5. **Billing.** protein-prep is not billable; pocket-finder and structure-report are.
   Confirm the workflow tool quotes as a single line item.

---

## 8. Manifest sketch

```jsonc
{
  "id": "target-prep",
  "toolKey": "deeporigin.target-prep",
  "toolVersion": "1",
  "identityHue": 210,
  "url": "/target-prep",
  "name": "Target Preparation",

  "steps": [
    {
      "id": "structure-report",
      "toolKey": "deeporigin.structure-report",
      "toolVersion": "latest",
      "visibility": "hidden",
      "inputSchema": { "type": "object", "properties": {
        "protein": { "type": "object", "x-data-type": "Protein", "properties": {
          "id":        { "type": "string", "x-data-type": "Protein.id" },
          "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } },
        "pdb_id":  { "type": "string", "x-data-type": "Protein.pdb_id" } } },
      "bind": { "structure_reports": "structureReport" }
    },
    {
      "id": "recommend",
      "toolKey": "deeporigin.protein-prep",
      "toolVersion": "2",
      "sync": true,                       // engine gap #1
      "visibility": "hidden",
      "inputSchema": { "type": "object", "properties": {
        "action":  { "type": "string", "x-from-form": { "field": "__const_recommend" } },
        "protein": { "type": "object", "x-data-type": "Protein", "properties": {
          "id":        { "type": "string", "x-data-type": "Protein.id" },
          "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } } } },
      "bind": { "recommendation": "recommendation" }
    }
  ],

  "layout": { "tiles": [
    { "id": "left", "component": "VerticalMosaic", "splitPercentage": 35, "tiles": [
      { "id": "protein-table",     "component": "Table" },
      { "id": "structure-filter",  "component": "StructureFiltering" } ] },
    { "id": "protein-viewer", "component": "ProteinViewer" } ] },

  "components": {
    "protein-table": { "component": "Table", "title": "Proteins", "config": {
      "entity": "proteins", "singleSelection": true, "isEditable": false,
      "previousRunsMarker": { "entityType": "proteins" },
      "columns": [ { "key": "display_id" }, { "key": "protein_name" },
                   { "key": "gene_symbol" }, { "key": "file_path" } ],
      "showsResults": true } },
    "structure-filter": { "component": "StructureFiltering", "title": "Structure Filtering",
      "config": { "recommendationKey": "recommendation", "selectionKey": "selection",
                  "reportKey": "structureReport" } },
    "protein-viewer": { "component": "ProteinViewer", "title": "Mol* Viewer", "config": {
      "molstar": { "edit":    { "renderer": "renderStructureComponents" },
                   "results": { "renderer": "renderStructureAndPockets" } } } }
  },

  "parameters": {
    "resultLabel": { "singular": "prepared structure", "plural": "prepared structures" },
    "requiredEntities": [ { "type": "Protein",
      "message": "Select a protein before running Target Preparation" } ],
    "validations": [ { "type": "state-resolved", "stateKey": "selection",
      "path": "decisions", "disallow": "review",
      "message": "Resolve every reviewed component before running" } ],
    "sections": [
      { "id": "structure-components", "label": "Structure Components",
        "fields": [ { "id": "components_summary", "type": "component-summary",
                      "optionsFrom": "selection", "labelsFrom": "recommendation" } ] },
      { "id": "advanced", "label": "Advanced Parameters", "collapsible": true,
        "fields": [
          { "id": "add_missing_atoms",   "type": "boolean", "default": true,
            "label": "Add missing atoms & residues", "helpText": "Fills incomplete residues" },
          { "id": "model_missing_loops", "type": "boolean", "default": true,
            "label": "Add missing loops", "helpText": "Models unresolved regions" },
          { "id": "protonate",           "type": "boolean", "default": true,
            "label": "Protonate", "helpText": "Assign states @ pH 7.4" },
          { "id": "find_pockets",        "type": "boolean", "default": true,
            "label": "Find Pockets", "helpText": "Detect druggable pockets" } ] },
      { "id": "run-details", "label": "Run Details", "fields": [
          { "id": "output_name", "type": "string", "label": "Output Property Name",
            "placeholder": "Default_TargetPrep_ProteinId" },
          { "id": "run_name", "type": "string", "label": "Run Name", "required": true,
            "x-body-key": "name", "x-run-name-default": "TargetPrep" } ] }
    ]
  },

  "inputSchema": { "type": "object",
    "required": ["protein", "selection"],
    "properties": {
      "protein":   { "type": "object", "x-data-type": "Protein", "properties": {
        "id":        { "type": "string", "x-data-type": "Protein.id" },
        "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } },
      "selection": { "type": "object", "x-from-state": "selection" },
      "pdb_id":              { "type": "string",  "x-data-type": "Protein.pdb_id" },
      "add_missing_atoms":   { "type": "boolean", "x-user-input": true },
      "model_missing_loops": { "type": "boolean", "x-user-input": true },
      "protonate":           { "type": "boolean", "x-user-input": true },
      "find_pockets":        { "type": "boolean", "x-user-input": true },
      "output_name":         { "type": "string",  "x-user-input": true } } }
}
```

The `action: "recommend"` constant is awkward to express with today's annotations — the
cleanest fix is `ActionField.inputs` (`{ "action": "recommend" }`), which `useRunStep`
already merges over the schema-built inputs as `inputOverrides`. Use that rather than the
`x-from-form` placeholder shown above.

---

## 9. Test surface

- **Unit** — `build-tool-payload` with `x-from-state` carrying a full `selection`; the
  decision reducer in the filtering tile (`review` never auto-resolves); the
  hex→Mol\*-colour mapping.
- **Storybook** — `StructureReportCard` across grades A–D and with
  `field_status: "unknown"` pills; `StructureFiltering` with a mixed keep/review/skip
  inventory.
- **E2E** (`apps/platform-e2e`) — select a protein, assert both steps fire and both tiles
  populate; toggle a component, assert the viewer recolours; assert Run stays disabled
  while a `review` remains.
- **Contract** — a fixture check that the manifest's `inputSchema` matches the published
  `deeporigin.target-prep` tool definition, so a backend schema change fails CI rather
  than production.

---

## 10. Phasing

**Phase 1 — app shell against `protein-prep` alone (no new backend tool).**
Manifest `toolKey: "deeporigin.protein-prep"`, submit sends `action: "prepare"`. Ships the
proteins table, both auto-run steps, both new tiles, the Mol\* renderer, and all four
engine gaps. Output is a prepared structure; no pockets, no final report. Everything built
here is reused verbatim in phase 2 — only `toolKey` and `inputSchema` change.

**Phase 2 — `deeporigin.target-prep` workflow tool.** Swap `toolKey`, add the four
parameter inputs, add `resultsLayout` with the pockets table and the final report.

**Phase 3 — report as table columns.** Depends on the `structurereport` result type
(§5b). Adds `results.detailColumns` / `aggregates` to the manifest; the table renders them
with no new component work.

**Backend critical path:** the four missing tool inputs (§2) and the workflow tool (§3)
gate phase 2; the `structurereport` result type gates phase 3; the prepared-protein
`proteins` row (§5a) gates Target Prep being useful to Docking/ABFE at all.
