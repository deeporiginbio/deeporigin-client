# Import from UniProt — DO Studio implementation plan

**Feature.** A new **Import from UNIPROT** entry in the proteins table's
**+ Data** menu in DO Studio (`platform-ui/apps/uui`): enter a UniProtKB
accession → see every experimental PDB for it, graded and sorted → preview any
of them in Mol\* → import the selected ones into the project's proteins table.

PRD:
[UNIPROT Addition via Protein Table](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1067679762/UNIPROT+Addition+via+Protein+Table).
This document is the **UI plan only**. `do-dd-client` is the source of truth for
what the tools accept and return and for what a correct import writes.

> **The PRD's three screenshots could not be read from this environment**
> (Confluence media blobs; direct HTTPS to `atlassian.net` is blocked by the
> egress policy and the Atlassian MCP exposes no attachment reader). Everything
> visual below is derived from the tool contract, the PRD text, and the ADF
> media dimensions — modal `509×155`, results table `749×608`, Mol\* preview
> `759×660`. Every place the mock would settle a detail is marked
> **[confirm against Figma]**.

---

## 1. What the user does

| Step | UI | Backend effect |
|---|---|---|
| 1 | + Data → **Import from UNIPROT** | none |
| 2 | Type accession, press **Find** / Enter | 1 tools execution (`uniprot-discovery`) |
| 3 | Results table: every candidate PDB, graded, sorted best→worst, top row badged **Recommended Structure** | 1 data-platform read (which PDBs are already in the project) |
| 4 | Eye button on a row | RCSB fetch in the browser, no platform call |
| 5 | Check ≥1 rows → **Import to Project** | 1 tools execution **per selected PDB** (`pdb-import`), each of which writes a `proteins` row |
| 6 | Modal closes, rows appear in the table | 1 data-platform read-back per created row |

---

## 2. Tool contracts (from `do-dd-client`)

### 2.1 `deeporigin.uniprot-discovery` — ranking

Source: `src/drug_discovery/uniprot_discovery.py`, registry
`src/platform/constants.py:161` (`tool_version: "latest"`).

```
POST /tools/{orgKey}/tools/deeporigin.uniprot-discovery/executions
body: { inputs: { uniprot_accession: "P00533" }, outputs: {}, clusterId, projectId?, name?, metadata? }
```

`inputs` has exactly one property. The SDK additionally sends `sync: true`
(`_default_execution_payload`, `src/drug_discovery/execution.py:54`) and reads
`jobOutputs` off the create response; **`sync` is not in the UI's
`ExecuteToolSchemaDto`** (`packages/api-gateway/openapi.yaml:11244`), so DO
Studio polls instead (§4.2).

Output — `jobOutputs.candidates: []`, one row per experimental PDB:

| Field | Type | Required | Notes |
|---|---|---|---|
| `pdb_id` | string | ✔ | 4 alphanumeric chars; SDK uppercases |
| `grade` | `"A"｜"B"｜"C"｜"D"` | ✔ | primary sort key |
| `weighted_score` | number | ✔ | tie-break |
| `recommended` | boolean | ✔ | exactly one `true` unless the list is empty |
| `coverage_score`, `inhibitor_score`, `method_score`, `organism_score`, `resolution_score`, `rfree_score` | number | ✔ | component scores behind the grade |
| `field_status` | `{ [field]: "value"｜"not_applicable"｜"unknown" }` | ✔ | why a value is missing |
| `coverage` | number | – | 0–1 fraction |
| `has_ligand` | boolean | – | non-water, non-ion ligand |
| `method`, `method_class` | string | – | e.g. `X-RAY DIFFRACTION` / `x-ray` |
| `organism`, `organism_class` | string | – | e.g. `Homo sapiens` / `human` |
| `resolution` | number | – | Å |
| `rfree` | number | – | X-ray only |

`jobOutputs` may arrive as a **dict or as a list whose first element is the
dict** — `_execution_outputs_dict` (`execution.py:79`) handles both, and the UI
must too (the existing `normalizeJobOutputs` in
`use-poll-tool-execution.ts` accepts only a plain object).

An **empty `candidates` list is a success**, not an error: unknown accession or
no experimental structures. Missing `candidates` altogether *is* an error.

Reference payload for mocks/tests:
`tests/mock_server/routers/tools.py:167` (`P99999` → empty list).

### 2.2 `deeporigin.pdb-import` — coordinate acquisition + registration

**Not in `do-dd-client`'s `TOOL_KEYS_AND_VERSIONS`** — the CLI still composes the
import client-side. The contract below comes from the platform ticket that
shipped it; `platform-toolbox` is not reachable from this session, so
**the exact input/output key names must be confirmed against
`tools/pdb-import/tool-definition.json` before coding**.

```
POST /tools/{orgKey}/tools/deeporigin.pdb-import/executions
body: { inputs: { pdb_id, uniprot_accession?, project_id? }, outputs: {}, clusterId, projectId, name? }
```

- Behaviour: download from `files.rcsb.org` → upload to UFA → `create_protein`.
- Output: `protein_id`, `file_path`, `pdb_id`; persists `uniprot_accession` on
  the entity when supplied.
- **One PDB per execution** — no batching. N selected candidates = N executions.

This is the step that writes to the data platform. The UI does **not** create
the protein row itself on this path.

### 2.3 What a correct protein row looks like (`do-dd-client`)

`Protein.register` (`src/drug_discovery/structures/protein.py:1741`) writes:

| Column | Value |
|---|---|
| `file_path` | UFA path of the uploaded structure (`entities/proteins/…`) |
| `pdb_id` | the PDB ID |
| `uniprot_accession` | the accession the discovery ran on |
| `protein_name` | the protein's name (`from_pdb_id` sets it to the PDB ID) |
| `protein_length` | residue count, when a local file was parsed |
| `project_id` | resolved project |

`Protein.sync` (`:1795`) is register plus a **dedupe**: it searches
`proteins` by `file_path` (+ `project_id`) first and reuses the existing row
instead of creating a second one (`client.entities.search_proteins`,
`src/platform/entities.py:684`, filter `{deleted: false, file_path, project_id}`).
`from_pdb_id` (`:280`) downloads `https://files.rcsb.org/download/{id}.pdb`.

All of these columns exist on the entity —
`PROTEIN_RETURNING_FIELDS`, `src/platform/entities.py:59`.

**Path B (fallback).** If `pdb-import` turns out not to be deployed where DO
Studio runs, the UI can mirror the SDK exactly with what the proteins table
already does: browser `fetch` from RCSB → `uploadFileToFileService({ orgKey,
file, subdir: 'structure-imports' })` → `serverInterface.entity.create(orgKey,
'proteins', { set: { file_path, pdb_id, uniprot_accession, protein_name,
project_id, tags: { app } }, returning })`. That is `handleProteinFiles`
(`chemical-table.tsx:374`) plus two fields. Keep the import step behind one
interface (`importCandidate(pdbId, accession)`) so A and B are swappable and
the rest of the feature does not care. **Decide A vs B before building the
import step** — everything else is identical.

---

## 3. Where each call goes

```
                          ┌─────────────────────────── DO Studio (apps/uui) ───────────────────────────┐

 Find ─────────────────►  useToolsListClusters            GET  /tools/{orgKey}/clusters
                          fetchToolsExecuteTool           POST /tools/{orgKey}/tools/
                                                               deeporigin.uniprot-discovery/executions
                          poll fetchToolsGetToolExecution GET  /tools/{orgKey}/tools/executions/{id}
                          → jobOutputs.candidates[]

 render rows ──────────►  serverInterface.search          POST /{orgKey}/proteins/search
                          (which pdb_ids are already                filter: project_id eq + pdb_id in [...]
                           in this project)

 eye ──────────────────►  browser fetch                   GET  https://files.rcsb.org/download/{ID}.pdb
                          MolstarViewer                   (no platform call, no import)

 Import to Project ────►  for each selected pdb_id:
                          fetchToolsExecuteTool           POST /tools/{orgKey}/tools/
                                                               deeporigin.pdb-import/executions
                          poll fetchToolsGetToolExecution GET  /tools/{orgKey}/tools/executions/{id}
                          → { protein_id, file_path }     ← the tool created the proteins row

 refresh ──────────────►  serverInterface.entity.get      GET  /{orgKey}/proteins/{protein_id}
                          grid.applyServerSideTransaction
                          tableRef.invalidateRowCount()
```

Every execution body also carries the standard execution context —
`useToolExecutionContext()` → `{ app, session, projectId }`
(`packages/global-provider/src/hooks/use-tool-execution-context.ts`) — plus
`clusterId` resolved the way `use-import-dataset.ts` does it (prefer a
`us-west-2` cluster, else the first).

**Data-platform writes: exactly one, and the tool makes it.** On path A the UI
never calls `entity.create` — it only *reads* (dedupe + read-back). On path B
the UI makes the write and must reproduce §2.3.

---

## 4. New code

### 4.1 File layout

Everything lives in `apps/uui`, not `packages/data-table`: the feature needs the
tools-service client, the data-platform provider and Mol\*, none of which the
generic grid package depends on. `useActionImportFromRcsb` only lives there
because it is a pure `fetch`.

```
apps/uui/src/components/data-platform-tables/uniprot-import/
├── index.ts
├── uniprot-import-modal.tsx        NEW  the whole flow, two steps in one Modal
├── candidate-table.tsx             NEW  presentational: DataTable over candidates[]
├── candidate-columns.tsx           NEW  column defs + grade/status/boolean cells
├── grade-badge.tsx                 NEW  A–D pill  [confirm against Figma]
├── structure-preview-modal.tsx     NEW  eye → Mol* over an RCSB blob
├── use-uniprot-discovery.ts        NEW  execute + poll uniprot-discovery
├── use-pdb-import.ts               NEW  execute + poll pdb-import, once per selection
├── use-existing-pdb-ids.ts         NEW  which candidates are already in the project
└── rcsb.ts                         NEW  shared RCSB URL + fetch-to-File/text helper

apps/uui/src/utils/
└── run-tool-and-wait.ts            NEW  promise-based execute+poll (see 4.2)
```

Changed files:

| File | Change |
|---|---|
| `components/data-platform-tables/chemical-table.tsx` | mount the modal (portaled) + add an inline `ActionMenuPlugin` and an `emptyTableActions` card, both gated on `isProteinTable` |
| `app-engine/available-components.ts` | add `uniprot_accession` to the proteins table's default columns |
| `packages/data-table/src/hooks/use-action-import-from-rcsb.tsx` | export its RCSB URL builder (or move it to the new `rcsb.ts` and import from there) so both import paths agree |
| `app-schemas/uniprot-discovery.json`, `app-schemas/pdb-import.json` | metadata-only manifests (`id`, `toolKey`, `name`, `identityHue`) registered in `manifestByToolKey` only, so Activity rows render with an identity instead of a raw tool key — mirrors `system-prep.json`, which is in `manifestByToolKey` but not in `appManifests` |

### 4.2 `run-tool-and-wait.ts` — why a new helper

`usePollToolExecution` (`components/datasets/import-modal/use-poll-tool-execution.ts`)
is a react-query hook holding **one** execution id in state. The import step
runs N executions in a loop, so it needs a promise:

```ts
runToolAndWait({ orgKey, toolKey, body, signal, intervalMs = 3000, timeoutMs })
  → { executionId, jobOutputs }   // rejects on Failed/Cancelled/InsufficientFunds/timeout
```

- `fetchToolsExecuteTool` → read `executionId ?? id` off the response.
- Poll `fetchToolsGetToolExecution` until `status` is terminal. Terminal set and
  success semantics are already settled in `do-dd-client`
  (`src/platform/constants.py`): success is `Completed` **or** legacy
  `Succeeded`; `Failed`, `Cancelled`, `InsufficientFunds`, `FailedQuotation`,
  `Quoted` are the other terminals. Reuse the same list as
  `use-poll-tool-execution.ts` plus `FailedQuotation`.
- Normalize `jobOutputs` for the dict-or-list shape (§2.1).
- On failure surface `statusReason`.

Use it for both tool calls; the discovery step can keep a thin hook around it
for loading state. Leave `use-poll-tool-execution.ts` where it is — the datasets
import is not part of this change.

### 4.3 `uniprot-import-modal.tsx`

One Mantine `Modal`, two visual states in the same shell so Find can be re-run
without losing context. **[confirm against Figma]** whether the modal grows
(509×155 → 749×608) or the results open in a second, larger modal.

State machine:

```
idle ──Find──► discovering ──ok──► results(candidates[])
                    │                   │
                    │ error             ├─ empty  → "No experimental structures for {ACC}"
                    ▼                   ├─ preview(pdbId) → structure-preview-modal
                  error                 └─ Import → importing(n) ──► done → close + toast
                 (retry)                                        └─► partial → keep open, mark failures
```

Rules:

- **Accession validation before the call.** Mirror `_UNIPROT_ACCESSION_RE`
  (`uniprot_discovery.py`): `^(?:[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z0-9]{2}\d){1,2})$`,
  case-insensitive, uppercase before sending. Reject isoform suffixes
  (`P00533-2`) with a specific message — a typo must never cost an execution.
- Enter in the input = Find (`onKeyDown`), matching the RCSB modal.
- Find is disabled while `discovering`; the request is abortable and re-runnable
  (new accession replaces the result set).
- Closing the modal mid-import is allowed: imports continue, the toast reports
  the outcome. Closing mid-discovery aborts the poll (the execution itself is
  already running and will land in Activity).

### 4.4 `candidate-table.tsx` + `candidate-columns.tsx`

Client-side `DataTable` from `@platform-ui/data-table` (same entry point
`CsvTable` uses), so selection, sorting and styling match the rest of DO Studio.

Columns **[confirm against Figma]**:

| Header | Field | Cell |
|---|---|---|
| _(checkbox)_ | – | `rowSelection: multiRow`, `headerCheckbox`, `enableClickSelection: false` |
| PDB | `pdb_id` | mono; the recommended row also carries the **Recommended Structure** badge |
| Grade | `grade` | `grade-badge.tsx` (A–D pill) |
| Score | `weighted_score` | 2 dp |
| Method | `method` | title-case `method_class` as tooltip |
| Resolution | `resolution` | `x.xx Å` |
| Coverage | `coverage` | `Math.round(v * 100)%` — tool emits 0–1, PRD shows a percentage |
| Rfree | `rfree` | 3 dp |
| Organism | `organism` | italic species, `organism_class` as tooltip |
| Ligand | `has_ligand` | check / dash |
| _(actions)_ | – | eye button → preview |

- **Sort**: `grade` A→D, then `weighted_score` desc. Do it in the component, on
  the array — the tool already returns ranked rows, but the table must not
  depend on server order for the "Recommended" row to sit on top.
- **Never recompute the grade.** The weights and thresholds live in the
  Structure Report tool; the client displays what it is given.
- **`field_status` drives empty cells**: `not_applicable` → `—` with tooltip
  "not applicable for this method"; `unknown` → `?` with tooltip "not reported".
  A blank cell that could mean either is the thing to avoid (Rfree on a cryo-EM
  entry is `not_applicable`, not missing data).
- **Already-in-project rows**: pattern already exists — copy
  `in-project-cell.tsx` / `inProjectRowSelection` from
  `app-engine/renderer/csv-table/`: unselectable row, dimmed, "Already in this
  project" tooltip.
- The component takes `candidates`, `existingPdbIds`, `selected`,
  `onSelectionChange`, `onPreview` and nothing else, so the Target Preparation
  app can mount the same table when it needs a structure chosen.

### 4.5 `structure-preview-modal.tsx`

Nested modal (or side panel — **[confirm against Figma]**; the 759×660 mock
looks like a full-width panel under the row) showing one candidate in Mol\*.

Follow `packages/chat-engine/src/viewers/structure-viewer.tsx`, not the
app-engine `protein-viewer-wrapper`: that one is manifest-driven, resolves a
`file_path` through `use-structure-files.ts`, and expects overlays, pockets and
docking boxes. Here there is one structure and no entity yet:

```ts
const text = await fetchRcsbPdbText(pdbId);   // https://files.rcsb.org/download/{ID}.pdb
new MolstarViewer(containerId) → init() → load(text, 'pdb')
```

- Cache per PDB ID for the modal's lifetime; re-opening a row must not refetch.
- Own container id per instance (`useId`), and a `live` flag so a late load never
  draws into a torn-down container — both mistakes the chat-engine viewer already
  documents.
- Header: PDB ID, grade badge, title/organism; a **Select this structure**
  action that checks the row and returns is worth having **[confirm against
  Figma]**.

### 4.6 `use-existing-pdb-ids.ts`

Direct analogue of `use-existing-smiles.ts` (`app-engine/renderer/csv-table/`):

```ts
serverInterface.search(orgKey, 'proteins', {
  select: ['pdb_id'],
  filter: { props: [
    { column: 'project_id', op: 'eq', value: projectId },
    { column: 'pdb_id',     op: 'in', value: candidatePdbIds },
  ]},
  limit: candidatePdbIds.length,
})
```

Cursor-paginate with the same `collectMatches` helper. This is what keeps a
second import of the same structure from happening — the CLI gets that for free
because `Protein.sync` dedupes on `file_path`; on path A the UI is the only
place that can.

### 4.7 `chemical-table.tsx` wiring

```tsx
// alongside importFromRcsbAction (:470)
const uniprotImport = useUniprotImportAction({ orgKey, projectId, baseEntity,
                                               agentIdPrefix: agentTableId,
                                               onImported: handleImportedProteins });

// extraToolbarItems (:1061) — portal the modal like the RCSB one
if (isProteinTable) items.push(<Fragment key="uniprot-modal-inline">{uniprotImport.modal}</Fragment>);

// the Data menu (:1135) — an inline plugin, same shape as combinedUploadPlugin (:1121)
plugins={[ …, ...(isProteinTable ? [importFromRcsbAction, uniprotImport.plugin] : []) ]}

// emptyTableActions (:1034)
...(isProteinTable ? [uniprotImport.emptyTableAction] : [])
```

`ActionMenuPlugin` only requires `menuItems`
(`components/data-platform-tables/action-menu.tsx`), and `EmptyTableAction`
requires `{ title, subtitle, onClick, icon, agentId? }`
(`packages/data-table/src/hooks/types.ts`).

Menu label **"Import from UniProt"**, empty-state card "Import from UniProt" /
"Rank experimental PDB structures for a UniProtKB accession", `icon: 'download'`
**[confirm against Figma]** (the PRD writes "UNIPROT"; the rest of DO Studio
uses sentence case).

### 4.8 Grid refresh after import

The row is created by the tool, so there is no response row to insert. For each
successful import:

```ts
const row = await serverInterface.entity.get(orgKey, 'proteins', protein_id);
gridApiRef.current?.applyServerSideTransaction({ add: [row.data] });
tableRef.current?.invalidateRowCount();
tableRef.current?.markNotEmpty();
```

and once for the batch:

```ts
queryClient.invalidateQueries({ queryKey: buildDataPlatformSearchQueryKey(orgKey, 'proteins') });
```

so other tiles (protein viewer, other tables) see the new rows. If `entity.get`
returns fewer columns than the grid shows, fall back to one
`serverInterface.search` with `filter: { props: [{ column: 'id', op: 'in',
value: createdIds }] }` and `select: activeColumnsFromGrid()`.

---

## 5. Errors, edge cases, limits

| Case | Behaviour |
|---|---|
| Invalid accession shape | inline field error, no execution |
| Discovery fails / times out | inline error in the modal + retry; keep the typed accession |
| `candidates: []` | empty state, not an error: "No experimental PDB structures for {ACC}" |
| `candidates` missing from `jobOutputs` | treat as a tool failure (matches the SDK, which raises) |
| Some imports fail | do **not** abort the batch — the CLI raises on the first failure, but a user who checked six rows should keep the five that worked. Per-PDB status in the table, one summary toast |
| PDB already in project | row unselectable + "Already in this project" |
| Many selections | N executions; run with bounded concurrency (3–4) and show `n/m imported`. **[confirm against Figma]** whether to cap the selection |
| Modal closed mid-import | imports continue; toast reports the result |
| Cost | each import is a billable execution. The SDK exposes `quote=True` / `approveAmount`; the PRD is silent. **Product decision** |

Activity page: both tools appear there per import. If `pdb-import` runs are
noise, filter them out with the existing mechanism —
`staticFilterProps={[{ column: 'tool_key', op: 'neq', value: 'deeporigin.pdb-import' }]}`
on the activity `JobManagerTable`.

---

## 6. Analytics

`useCaptureEvent` (`@platform-ui/global-provider`), `EventName` in
`providers/analytics-provider.tsx`:

- `ToolExecutionStarted` / `ToolExecutionResult` for each execution
  (`toolKey`, `status`, `projectId`).
- `DataUploaded` on a successful import: `uploadType: 'uniprot'`,
  `recordCount`, `status`, `projectId` — same shape `handleRcsbImport` uses
  (`chemical-table.tsx:451`).

---

## 7. Tests

| Level | What |
|---|---|
| Unit | accession regex (valid 6/10-char, isoform rejected); grade sort + recommended-first; `field_status` cell rendering; `jobOutputs` dict-vs-list normalization; terminal-status handling in `run-tool-and-wait` |
| Component | modal state machine (idle → discovering → results → importing), empty candidates, partial import failure, already-in-project rows unselectable |
| Table | extend `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx` — it already asserts the + Data menu composition and that protein-only plugins are absent on ligand tables |
| E2E | `apps/platform-e2e/src/pages/` — page-object entry next to `addDataButton`; mock both executions |

Mock payloads: copy `tests/mock_server/routers/tools.py:167` from this repo so
the UI and SDK test against the same candidate shape, including the `P99999`
empty case.

---

## 8. Open items

1. **`pdb-import` availability and exact schema.** `platform-toolbox` is not
   reachable from this session — confirm `tools/pdb-import/tool-definition.json`
   input/output key names (`pdb_id` / `uniprot_accession` / `project_id` →
   `protein_id` / `file_path` / `pdb_id`) and that the tool is deployed in the
   environments DO Studio targets. This is the one thing that decides path A vs
   path B (§2.3).
2. **The three Figma screenshots** — paste them into a comment or share the
   Figma link and every **[confirm against Figma]** above closes.
3. **Cost/quote behaviour** for the import step.
4. **Selection cap** on Import to Project.
5. **`uniprot_accession` column** — worth adding to the proteins table's default
   columns so an imported row shows where it came from.
