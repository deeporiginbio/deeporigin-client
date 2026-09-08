# UniProt discovery in DO Studio — integration research

**Goal.** Ship the PRD
[UNIPROT Addition via Protein Table](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1067679762/UNIPROT+Addition+via+Protein+Table)
in DO Studio (`platform-ui/apps/uui`): an "Import from UNIPROT" item in the
proteins table's **+ Data** menu, backed by the same platform tools this SDK
already drives.

Status: research + implementation plan. No code in either repo yet.

Related work items:

| Item | What it is | Status |
|---|---|---|
| [DDOS-7091](https://deeporigin.atlassian.net/browse/DDOS-7091) | Target Preparation Tool (parent epic; the PRD is one of its inputs) | In Progress |
| [DDOS-7094](https://deeporigin.atlassian.net/browse/DDOS-7094) | Ship UniProt-to-PDB discovery served tool | Done |
| [DDOS-7380](https://deeporigin.atlassian.net/browse/DDOS-7380) | `Protein.from_uniprot` + register selected structures (this repo) | In review |
| [DDOS-7555](https://deeporigin.atlassian.net/browse/DDOS-7555) | Build `deeporigin.pdb-import` for registered PDB coordinate acquisition | Ready for Release |
| [DDOS-7039](https://deeporigin.atlassian.net/browse/DDOS-7039) | Bug: duplicate proteins can be added to the Proteins table | To Do |

No PUI ticket exists for the UI work yet.

---

## 1. What the PRD asks for

From the requirements table (priority 0):

1. **+ Data → "Import from UNIPROT"** on the Protein table opens a modal with a
   UniProt code input and a **Find** action (Enter also submits).
2. On submit, **all related PDBs are listed**, one row per structure, in
   **tabular Structure Report format**, sorted **best to worst by Grade**, with
   the top row highlighted as the **"Recommended Structure"**.
3. Each row has a **preview / eye button** that opens a **Mol\*** window showing
   that structure.
4. The user can select **any number (≥1)** of structures and **"Import to
   Project"**, which brings them into that project's protein table.

The PRD embeds Figma screenshots (modal, results table, Mol\* preview) that are
Confluence media blobs — **the exact column list, copy and layout must be read
off those designs before building**; the column set below is derived from the
tool contract, not from the mock.

### Row columns — Structure Report output

The results table shows the
[Structure Report Tool](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1066631169/Structure+Report+Tool)
fields, which `deeporigin.uniprot-discovery` returns per candidate:

| Column | Field | Notes |
|---|---|---|
| PDB | `pdb_id` | 4-character code |
| Grade | `grade` | `A`–`D`; sort key (desc) |
| Score | `weighted_score` | weighted sum behind the grade |
| Organism | `organism` / `organism_class` | human, mammal, vertebrate, other |
| Method | `method` / `method_class` | X-ray, cryo-EM, NMR, other |
| Resolution | `resolution` | Å |
| Coverage | `coverage` | 0–1 fraction; PRD shows a percentage |
| Rfree | `rfree` | X-ray only |
| Inhibitor | `has_ligand` | non-water, non-ion ligand present |
| — | `recommended` | exactly one `true` → "Recommended Structure" badge |
| — | `field_status` | per field: `value` / `not_applicable` / `unknown` — render "n/a" vs "unknown" instead of a blank cell |

Grade weights and the A/B/C/D thresholds are documented in the Structure Report
PRD but are **tool-owned**: the client must display, never recompute. The grade
is informational — epic decision #5: "No grade threshold accepts or rejects a
structure."

---

## 2. What the CLI does today

Source: [`src/drug_discovery/uniprot_discovery.py`](../../src/drug_discovery/uniprot_discovery.py),
docs: [`docs/dd/tools/uniprot-discovery.md`](../../docs/dd/tools/uniprot-discovery.md).

### The discovery tool

| | |
|---|---|
| Tool key | `deeporigin.uniprot-discovery` |
| Tool version | `latest` (`src/platform/constants.py:161`) |
| Endpoint | `POST /tools/{orgKey}/tools/{toolKey}/{toolVersion}/executions` (`src/platform/executions.py:56`) |
| Inputs | `{"uniprot_accession": "P00533"}` — the whole input schema |
| Body | `{inputs, outputs: {}, metadata: {}, sync: true, name?, approveAmount?}` (`_default_execution_payload`, `src/drug_discovery/execution.py:54`) |
| Outputs | `jobOutputs.candidates[]`; an empty list is a valid success |

It is a **tool**, not a tools-service *function* (`/tools/{orgKey}/functions/{key}`),
so DO Studio must use `fetchToolsExecuteTool`, not `useToolsRunFunction`.

Required per candidate: `pdb_id`, `grade`, `recommended`, `weighted_score`, the
six component scores (`coverage_score`, `inhibitor_score`, `method_score`,
`organism_score`, `resolution_score`, `rfree_score`) and `field_status`.
Optional: `coverage`, `has_ligand`, `method`, `method_class`, `organism`,
`organism_class`, `resolution`, `rfree`. A representative payload lives in the
mock server (`tests/mock_server/routers/tools.py:167`); accession `P99999`
returns an empty list.

Client-side validation before any call: accession is 6 or 10 characters
(`_UNIPROT_ACCESSION_RE`; isoform suffixes rejected), PDB IDs are 4 alphanumeric
characters, uppercased.

### Import is a separate step

`run()` only ranks — **the discovery tool never writes protein rows**. The CLI
composes the import itself (`import_proteins`): `Protein.from_pdb_id(pdb_id)`
downloads from RCSB, sets `uniprot_accession` + `project_id`, then `sync()`
uploads to UFA and creates the `proteins` row (`file_path`, `pdb_id`,
`uniprot_accession`, `protein_name`, `protein_length`) —
`src/drug_discovery/structures/protein.py:1741`.

**The platform has since moved that composition server-side.**

### `deeporigin.pdb-import` — the import step DO Studio should use

DDOS-7555 (Ready for Release) adds a served tool that closes the gap after
discovery:

| | |
|---|---|
| Tool key | `deeporigin.pdb-import` |
| Inputs | `pdb_id` (required); optional `uniprot_accession`, `project_id` |
| Behavior | download from `files.rcsb.org` → upload to UFA → `create_protein` |
| Outputs | `protein_id`, `file_path`, `pdb_id`; persists `uniprot_accession` on the entity |
| Cardinality | **exactly one PDB per execution** — callers must not batch |

The ticket is explicit that the CLI's `from_pdb_id()` + `sync()` pair is only a
reference for the contract, and that "Target Preparation **and the Protein
table** must not require caller-side composition". Epic decision #33 says the
same. So the DO Studio import must call this tool once per selected candidate —
**not** re-implement fetch + upload + `entity.create` in the browser the way the
RCSB action does.

This tool is not yet in this repo's `TOOL_KEYS_AND_VERSIONS`; `import_proteins`
still composes client-side. Worth its own SDK ticket so CLI and UI converge.

---

## 3. What DO Studio has today

DO Studio is `platform-ui/apps/uui` (README: "DO Studio (formerly known as
UUI)"). The proteins table is `ChemicalTable` with `entity: "proteins"`, used by
the app-engine `Table` tile (`app-engine/renderer/table-wrapper/table-wrapper.tsx:1033`)
and the dashboard (`app-engine/available-components.ts:4`).

### The + Data menu

`apps/uui/src/components/data-platform-tables/chemical-table.tsx`:

- `isProteinTable` — `:350`
- `extraToolbarItems` builds the `ActionMenu title="Data"` — `:1061`, menu at
  `:1135` (`testId="add-data-menu"`)
- The one protein-only plugin today: `importFromRcsbAction` — `:470`
- The same actions also appear as empty-state cards — `emptyTableActions`, `:1034`

An `ActionMenuPlugin` needs only `menuItems`
(`components/data-platform-tables/action-menu.tsx`), and `chemical-table.tsx`
already builds one inline (`combinedUploadPlugin`, `:1121`) — so a menu entry
does **not** have to be a `@platform-ui/data-table` hook.

### Import from RCSB — the UX precedent

`packages/data-table/src/hooks/use-action-import-from-rcsb.tsx`: modal, one
input, browser `fetch('https://files.rcsb.org/download/{ID}.pdb')`, then
`onImport(file)` → `handleRcsbImport` (`chemical-table.tsx:451`) →
`handleProteinFiles` (`:374`) → `uploadFileToFileService` +
`serverInterface.entity.create` + `applyServerSideTransaction`.

Take the **shape** of this flow (menu item, modal, empty-state card, grid
refresh, `DataUploaded` capture) but not its import mechanics — `pdb-import`
owns those now.

### Running a tool from a modal

`components/datasets/import-modal/use-import-dataset.ts` is the precedent:
`useToolsListClusters<any>` (prefer `us-west-2`, else first) →
`fetchToolsExecuteTool({ pathParams: { orgKey, toolKey }, body })` →
`usePollToolExecution(orgKey)` (`use-poll-tool-execution.ts`, 3 s poll to a
terminal status, hands back `jobOutputs`).

`ExecuteToolSchemaDto` (`packages/api-gateway/openapi.yaml:11244`) is
`{ inputs, outputs?, clusterId, name?, projectId?, metadata?, approveAmount? }`
— **no `sync` field**, so DO Studio polls rather than relying on the SDK's
`sync: true`.

### Mol\* in a modal

`packages/chat-engine/src/viewers/structure-viewer.tsx` is the right precedent
for the per-row preview: a few dozen lines around `MolstarViewer` from
`@platform-ui/molstar`, taking structure text + extension, built for "one file
and a modal" rather than the app-engine's manifest-driven
`protein-viewer-wrapper` (which resolves `file_path` → blob URL via
`use-structure-files.ts` and expects overlays, pockets and docking boxes).

A candidate is **not imported yet**, so the preview has no `file_path`: fetch the
PDB text straight from `files.rcsb.org` (same URL the RCSB action uses) and hand
it to that viewer. Cache per PDB ID so re-opening a row is free.

---

## 4. Proposed implementation

### 4.1 Where the code lives

Put the whole feature in **`apps/uui`**, not `packages/data-table`:

```
apps/uui/src/components/data-platform-tables/uniprot-import/
  uniprot-import-modal.tsx      // accession → Find → candidate table → Import
  candidate-table.tsx           // Structure Report columns, grade sort, recommended badge
  structure-preview-modal.tsx   // eye button → Mol* (chat-engine StructureViewer pattern)
  use-uniprot-discovery.ts      // execute deeporigin.uniprot-discovery + poll
  use-pdb-import.ts             // execute deeporigin.pdb-import per selected candidate
  index.ts
```

Reason: the modal needs the tools-service client, `useToolExecutionContext`, and
Mol\*. `packages/data-table` is a generic grid package with none of those
dependencies, and the RCSB hook only lives there because it is pure fetch. Wire
it into `chemical-table.tsx` as an inline `ActionMenuPlugin` next to
`combinedUploadPlugin`, gated on `isProteinTable`, plus one `emptyTableActions`
card.

If a later product decision wants the same modal on other protein surfaces,
promote it to a package then — not before.

### 4.2 Flow

```
+ Data → "Import from UNIPROT"
   │
   ├─ accession input (validate 6/10-char UniProtKB shape client-side) → Find / Enter
   │     └─ execute deeporigin.uniprot-discovery { uniprot_accession }  → poll → candidates[]
   │
   ├─ candidate table: sort by grade desc (tie-break weighted_score desc),
   │     recommended row badged + pre-checked, checkbox multi-select,
   │     eye button per row → Mol* preview from files.rcsb.org
   │
   └─ "Import to Project" (n ≥ 1)
         └─ for each checked pdb_id: execute deeporigin.pdb-import
              { pdb_id, uniprot_accession, project_id } → poll → { protein_id, file_path }
         └─ refresh the grid, notify with a per-PDB success/failure summary
```

Details that matter:

- **One execution per PDB** (pdb-import cardinality). Run them with bounded
  concurrency, report per-PDB outcomes, and do not abort the batch on the first
  failure — unlike the CLI, which raises.
- **Grid refresh.** The row is created server-side, so there is no response row
  to hand `applyServerSideTransaction`. Either re-read the created rows by
  `protein_id` through `serverInterface.entity` and add them, or invalidate the
  datasource and `invalidateRowCount()` / `markNotEmpty()`. Prefer the read —
  it keeps the new rows visible without a full refetch.
- **Provenance.** Pass `projectId` and spread `useToolExecutionContext()` into
  both execution bodies, as every other execution in the app does. Rows written
  by the tool get the backend's own tags, so the manual `provenanceTags`
  workaround (PUI-2146) is not needed on this path.
- **Empty result** (`candidates: []`) is a success: "No experimental PDB
  structures for {accession}", not an error toast.
- **Analytics.** `EventName.DataUploaded` with `uploadType: 'uniprot'`,
  `recordCount`, `projectId`.

### 4.3 Column addition

`uniprot_accession` is not among the dashboard proteins-table default columns
(`available-components.ts:12`). Add it so an imported row visibly carries its
accession.

---

## 5. Relation to the app engine

The App Engine is manifest-driven: one JSON per app in `apps/uui/src/app-schemas/`,
mapped 1:1 to a backend tool, with a sidebar form (`inputSchema` +
`x-user-input` / `x-data-type`) and `engine.submit()` → `buildToolPayload` →
`executeToolWithVersion` (`apps/uui/CLAUDE.md`).

**This should not be an app.** The engine's contract is *select entities in the
project → run a tool on them*; UniProt discovery takes a free-text accession,
selects nothing, and exists to *create* the proteins an app later consumes. The
PRD places it in + Data for exactly that reason.

Three real connections remain:

1. **Same execution plumbing** — same endpoint, same `clusterId` resolution,
   same `useToolExecutionContext()` provenance, same `projectId` scoping.
2. **Activity page.** Both executions show up in `JobManagerTable`, which
   resolves tool identity through `manifestByToolKey` (`app-schemas/index.ts:33`,
   `job-manager-table/index.tsx:66`); without an entry they render a raw tool key
   and no `identityHue`. `system-prep.json` is the precedent for a
   **metadata-only manifest** (in `manifestByToolKey`, absent from
   `appManifests`, so no route). Ship one each for `uniprot-discovery` and
   `pdb-import`. Note both will appear in Activity per import — worth a product
   decision on whether `pdb-import` runs should be filtered out of the user's
   activity list (`staticFilterProps` on `tool_key` is the existing mechanism).
3. **Target Preparation will reuse this table.** Epic decision #29: for a UniProt
   input, AUTO takes the recommended candidate, while Expert mode returns
   `status: needs_structure_selection` with ranked candidates and resumes after
   the caller picks one. That is the same candidate list with the same columns
   and the same preview — so build `candidate-table.tsx` as a presentational
   component over `candidates[]`, independent of the modal's execution
   machinery, and the Target Preparation app can mount it unchanged.

`ligand-search` is the nearest existing app-engine cousin (external query →
`CsvTable` hits → `addToProject` writes selected rows into an entity;
`app-engine/renderer/csv-table/use-add-to-project.ts`), but `addToProject` does a
plain `entity.create` and cannot acquire coordinates — another reason the
protein path goes through `pdb-import`.

---

## 6. CLI ↔ DO Studio parity

| CLI | DO Studio |
|---|---|
| `UniprotDiscovery(uniprot_accession=...)` | accession input + Find |
| `run()` | execute `deeporigin.uniprot-discovery` + poll → `candidates` |
| `candidates`, `recommended` | candidate table, recommended badged and pre-checked |
| — | eye button → Mol\* preview from RCSB (no CLI equivalent) |
| `import_proteins()` (no args) | Import with only the recommended row checked |
| `import_proteins([...])` | Import with a multi-row selection |
| `Protein.from_pdb_id` + `sync()` (client-side) | `deeporigin.pdb-import` (server-side) |
| `uniprot_accession` / `pdb_id` on the row | same fields, written by the tool |
| `run(quote=True)` / `approve_amount` | not exposed (see open questions) |

---

## 7. Open questions / risks

- **Is `deeporigin.pdb-import` deployed** in the environments DO Studio targets,
  and what is its registered version? DDOS-7555 is "Ready for Release", which is
  not the same as available. Confirm before building against it; the browser
  fetch + upload path is the fallback, and it is the thing the epic says not to
  build.
- **Figma fidelity.** Column list, ordering, empty/loading states and copy must
  come from the PRD screenshots; this doc derives columns from the tool contract
  only.
- **Duplicates.** DDOS-7039 is open against the proteins table, and importing the
  same PDB twice is the obvious way to hit it. Decide whether the modal
  pre-checks "already in this project" per candidate (an `entity.search` on
  `pdb_id` for the project) and disables those rows.
- **Cost.** Two executions per import (discovery + one pdb-import per structure).
  The CLI supports `quote=True` / `approveAmount`; the PRD says nothing about
  cost. Decide whether the modal quotes or just runs.
- **Latency.** Discovery is a blocking tool run behind a 3 s poll; the modal needs
  a real loading state, and Find should be cancellable / re-runnable.
- **Preview fetch.** Browser fetch from `files.rcsb.org` already works in the RCSB
  action, so CORS is not a new risk, but large structures should stream into the
  viewer with a spinner and a per-PDB cache.
- **Accession validation.** Mirror `_UNIPROT_ACCESSION_RE` so a typo never costs
  an execution, and reject isoform suffixes (`P00533-2`) with a clear message.
- **Tool version pinning.** The CLI pins `latest`; decide whether the UI pins a
  version for reproducibility (the versioned execute variant exists).

---

## 8. Implementation checklist

### platform-ui (`apps/uui`)

1. `components/data-platform-tables/uniprot-import/` — modal, candidate table,
   preview modal, the two execution hooks (§4.1).
2. Move `use-poll-tool-execution.ts` out of `components/datasets/import-modal/`
   into a shared hooks folder and reuse it; don't copy it.
3. `chemical-table.tsx` — inline `ActionMenuPlugin` + `emptyTableActions` card,
   both gated on `isProteinTable`; portal the modal like the RCSB one.
4. Grid refresh after import (read created rows by `protein_id`, else invalidate).
5. `available-components.ts` — add the `uniprot_accession` column.
6. `app-schemas/uniprot-discovery.json` and `app-schemas/pdb-import.json` —
   metadata-only manifests registered in `manifestByToolKey` only.
7. Tests: extend
   `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx`
   (it already asserts the + Data menu composition); unit-test the candidate
   sort/recommended logic and the empty-candidates state; add an e2e page-object
   entry beside `addDataButton` in `apps/platform-e2e/src/pages/`.
8. Analytics: `uploadType: 'uniprot'`.

### do-dd-client (this repo)

9. Add `deeporigin.pdb-import` to `TOOL_KEYS_AND_VERSIONS` and switch
   `UniprotDiscovery.import_proteins` (and `Protein.from_uniprot`) to it, so the
   CLI and DO Studio import through the same contract instead of two
   implementations of download-upload-register. Needs its own ticket.
