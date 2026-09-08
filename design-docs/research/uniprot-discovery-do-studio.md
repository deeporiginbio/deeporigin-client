# Import from UniProt — DO Studio implementation plan

**Feature.** A new entry in the proteins table's **+ Data** menu in DO Studio
(`platform-ui/apps/uui`): enter a UniProtKB accession → see every experimental
PDB for it, graded and sorted → preview any of them in Mol\* → import the
selected ones into the project's proteins table.

PRD:
[UNIPROT Addition via Protein Table](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1067679762/UNIPROT+Addition+via+Protein+Table).
This document is the **UI plan only**. `do-dd-client` is the source of truth for
what the tools accept and return, and for what a correct import writes.

---

## 1. UI spec, from the mocks

The mocks reuse an existing "Import PDB Structure" modal as their frame. **This
flow is UniProt-first**: the user types an accession, and PDB structures are the
*result*, not the input. Copy accordingly:

| Element | Mock | Use |
|---|---|---|
| Modal title | "Import PDB Structure" | **"Import from UniProt"** |
| Input placeholder | `e.g. P01116` | keep — `P01116` is a UniProtKB accession (KRAS) |
| Input label / aria | – | "UniProt accession" |
| Helper text | "Press enter to search" | keep |
| Count line | "6 structures available" | keep — "{n} structures available" |
| Primary button | "Import Selected" | keep |
| Menu item | – | "Import from UniProt" |

### 1.1 Modal — empty state (`509×155`)

```
┌────────────────────────────────────┐
│ Import from UniProt                │
│ ┌────────────────────────────────┐ │
│ │ e.g. P01116                    │ │   ← text input
│ └────────────────────────────────┘ │
│ Press enter to search              │   ← helper text, dimmed, xs
└────────────────────────────────────┘
```

There is **no Find button** in the mock — Enter is the submit. The PRD prose
says "hit *Find*"; the mock wins on layout. Put a search icon button in the
input's `rightSection` so the action is also clickable.

### 1.2 Modal — results state (`749×608`)

Same modal, grown; the input stays at the top so a new accession can be searched
in place.

```
┌──────────────────────────────────────────────────────────────┐
│ Import from UniProt                                        ✕ │
│ ┌──────────────────────────────────────────────────────────┐ │
│ │ e.g. P01116                                              │ │
│ └──────────────────────────────────────────────────────────┘ │
│ Press enter to search                                        │
│ 6 structures available                                       │  ← bold, xs
│ ┌──────────────────────────────────────────────────────────┐ │
│ │    ID    Organism  Method   Resolution Coverage Rfree Grade│ │
│ │ ☆ ☑ 1ON1  Human    X-ray      1.8Å       83%    0.231  A 👁│ │ ← highlighted
│ │   ☐ 1EZY  Rat      Cryo-EM    3.8Å       93%    0.333  A 👁│ │
│ │   ☐ 4M35  Human    X-ray      2.5Å       91%    0.101  B 👁│ │
│ │   ☐ 2D5G  Human    X-ray      2.1Å       80%    0.100  C 👁│ │
│ │   ☐ 9EHU  Rat      Cryo-EM    2.7Å       79%    0.294  C 👁│ │
│ │   ☐ 2DD6  Mouse    X-ray      2.0Å       66%    0.231  D 👁│ │ ← scrolls
│ └──────────────────────────────────────────────────────────┘ │
│                                          [ Import Selected ] │  ← primary
└──────────────────────────────────────────────────────────────┘
```

Settled by the mock:

- **Columns**: `ID · Organism · Method · Resolution · Coverage · Rfree · Grade`,
  plus a leading star column, a leading selector, and a trailing eye button.
  **No score column, no ligand column, no title column** — the tool returns
  `weighted_score` and `has_ligand`, and the mock deliberately leaves them out.
- **Recommended row** = star icon (`☆`) plus a tinted row background,
  pre-selected. No "Recommended Structure" text badge; put the wording in the
  star's tooltip.
- **Grade** is a bare coloured letter, not a pill: A green, B amber, C orange,
  D red.
- Table body **scrolls** at a fixed max height; header stays.
- Footer: single primary button, bottom-right. `✕` top-right.

**Selector: checkboxes, not the mock's radios.** The mock draws radio buttons
(single selection); the requirement text says "select any number (e.g. at least
1)". Nothing in the import mechanism constrains this (§2.2), so follow the
requirement — a checkbox column degrades to the mock's behaviour when one row is
checked, and "Import Selected" reads correctly either way.

### 1.3 Preview window (`759×660`)

The eye raises a **floating window layered over the import modal** — the modal
stays mounted behind it.

```
┌─────────────────────────────────────────┐
│ 1ON1                        ⚙  ⤢   ✕   │
├─────────────────────────────────────────┤
│        [ Mol* canvas — cartoon +        │
│          translucent surface ]          │
└─────────────────────────────────────────┘
```

- Title is the **PDB ID alone**.
- `⚙` is Mol\*'s own settings-panel toggle — already shipped in
  `CustomViewportControls`
  (`packages/molstar/src/components/custom-viewport-controls.tsx`,
  `toggleSettingsPanel`), wired by `createViewerSpec`
  (`config/viewer-config.ts`: `controls.right = SettingsPanel`). It comes free
  with `MolstarViewer` — do not rebuild it. That control group also renders
  screenshot and sequence-viewer buttons the mock doesn't show: accept them or
  hide them in the preview's CSS.
- `⤢` toggles the preview to full screen (Mantine `Modal fullScreen`).
- `✕` closes the preview; the modal's selection is untouched.

---

## 2. The tools

### 2.1 `deeporigin.uniprot-discovery` — the one tool this feature calls

Source: `src/drug_discovery/uniprot_discovery.py`; registry
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

Output — `jobOutputs.candidates[]`, one row per experimental PDB:

| Field | Type | Req | Used by the mock |
|---|---|---|---|
| `pdb_id` | string | ✔ | **ID** |
| `grade` | `A｜B｜C｜D` | ✔ | **Grade** + primary sort |
| `weighted_score` | number | ✔ | sort tie-break only (not shown) |
| `recommended` | boolean | ✔ | star + row highlight; exactly one `true` unless empty |
| `coverage_score`, `inhibitor_score`, `method_score`, `organism_score`, `resolution_score`, `rfree_score` | number | ✔ | not shown |
| `field_status` | `{ [field]: "value"｜"not_applicable"｜"unknown" }` | ✔ | empty-cell semantics |
| `coverage` | number | – | **Coverage** (0–1 → `%`) |
| `has_ligand` | boolean | – | not shown |
| `method`, `method_class` | string | – | **Method** (from `method_class`) |
| `organism`, `organism_class` | string | – | **Organism** (see §4.5) |
| `resolution` | number | – | **Resolution** (Å) |
| `rfree` | number | – | **Rfree** |

`jobOutputs` may arrive as a **dict or as a list whose first element is the
dict** — `_execution_outputs_dict` (`execution.py:79`) handles both, and the UI
must too (the existing `normalizeJobOutputs` in `use-poll-tool-execution.ts`
accepts only a plain object).

An **empty `candidates` list is a success** (unknown accession, or no
experimental structures). A missing `candidates` key *is* a failure.

Reference payload for mocks/tests: `tests/mock_server/routers/tools.py:167`
(`P99999` → empty list).

### 2.2 There is no served tool that imports a PDB — the UI composes the import

`deeporigin.pdb-import` was specified (download from RCSB → upload to UFA →
`create_protein`) but **it is not a registered tool as far as anything reachable
here shows**:

- `TOOL_KEYS_AND_VERSIONS` (`src/platform/constants.py:99`) is this repo's
  single registry of platform tools — 17 entries, no `pdb-import`. Nothing in
  `src/`, `tests/` or `docs/` mentions it.
- `platform-ui` never references it either.

The other registered tools do not do this job:

| Tool | Why not |
|---|---|
| `deeporigin.uniprot-discovery` | ranks candidates; explicitly never downloads or writes |
| `deeporigin.structure-report` | grades a structure; "does not mutate or download a structure" |
| `deeporigin.import-dataset` | bulk catalog import — needs `csv_path`, `database_key`, `database_version`, `dataset_schema`, optionally a pre-uploaded `protein_zip_path` (`apps/uui/.../use-import-dataset.ts`). Not fetch-by-ID |
| `deeporigin.protein-prep` | prepares an already-registered protein |

**So the import step is client-side, exactly as both the SDK and DO Studio
already do it today:**

| | SDK (`do-dd-client`) | DO Studio (existing, Import from RCSB) |
|---|---|---|
| fetch | `Protein.from_pdb_id` → `https://files.rcsb.org/download/{id}.pdb` (`protein.py:280`) | `fetch(rcsbPdbUrl(id))` in `use-action-import-from-rcsb.tsx` |
| upload | `protein.upload()` → UFA | `uploadFileToFileService({ orgKey, file, subdir: 'structure-imports' })` |
| register | `protein.sync()` → `proteins` row (`protein.py:1741`, `:1795`) | `serverInterface.entity.create(orgKey, 'proteins', { set, returning })` (`chemical-table.tsx:374`) |

This is a proven path in this exact table, and it removes the one-PDB-per-
execution constraint the served tool would have imposed — multi-select is free.

**If `pdb-import` does land**, swap it in behind a single
`importCandidates(pdbIds, accession)` interface: one execution per PDB with
`{ pdb_id, uniprot_accession, project_id }`, polled like the discovery run, and
the row read back by `protein_id`. Nothing else in the feature changes. Verify
by listing the org's tools (`GET /tools/{orgKey}/tools`) or checking
`platform-toolbox` — neither is reachable from this session.

### 2.3 What the import must write

`Protein.register` (`src/drug_discovery/structures/protein.py:1741`) is the
reference:

| Column | Value |
|---|---|
| `file_path` | UFA path returned by the upload |
| `pdb_id` | the candidate's PDB ID (uppercase) |
| `uniprot_accession` | the accession the discovery ran on |
| `protein_name` | the PDB ID (what `from_pdb_id` sets) |
| `project_id` | the active project |
| `tags` | `{ app }` from `useToolExecutionContext()` — PUI-2146, so the row matches the dashboard's Source filter like tool-written rows do |
| `protein_length` | residue count — the SDK sets it from the parsed structure; the UI has no PDB parser, so leave it unset (optional field) or count distinct `CA` records in the fetched text |

Every column exists on the entity — `PROTEIN_RETURNING_FIELDS`,
`src/platform/entities.py:59`.

**Dedupe differs from the SDK.** `Protein.sync` (`:1795`) dedupes by searching
`proteins` for the same `file_path` (+ `project_id`) before creating a row —
that works because the SDK's remote path is content-derived. The UI's
`uploadFileToFileService` builds a **timestamped** path
(`structure-imports_{Date.now()}/{name}`), so the same PDB uploaded twice never
collides. The UI must therefore dedupe on **`pdb_id` + `project_id`** *before*
uploading (§4.6).

---

## 3. Where each call goes

```
 Enter in the input ──►  useToolsListClusters            GET  /tools/{orgKey}/clusters
                         fetchToolsExecuteTool           POST /tools/{orgKey}/tools/
                                                              deeporigin.uniprot-discovery/executions
                         poll fetchToolsGetToolExecution GET  /tools/{orgKey}/tools/executions/{id}
                         → jobOutputs.candidates[]

 render rows ─────────►  serverInterface.search          POST /{orgKey}/proteins/search
                         (which pdb_ids are already            filter: project_id eq + pdb_id in [...]
                          in this project — §4.6)

 eye ─────────────────►  browser fetch                   GET  https://files.rcsb.org/download/{ID}.pdb
                         viewer.api.loadFromRawContent   (no platform call, no import)

 Import Selected ────►   per selected pdb_id (bounded concurrency 3–4):
                           browser fetch                 GET  https://files.rcsb.org/download/{ID}.pdb
                           fetchFilesPutObject           PUT  /files/{orgKey}/{filePath}
                                                              (via uploadFileToFileService)
                         then once for the batch:
                           serverInterface.batch.create  POST /{orgKey}/proteins/batch-create
                                                              rows: [{ file_path, pdb_id,
                                                                       uniprot_accession, protein_name,
                                                                       project_id, tags }]
                                                              returning: [...grid columns, 'id']

 refresh ─────────────►  grid.applyServerSideTransaction ← rows come back from batch.create
                         tableRef.invalidateRowCount()
```

The discovery execution body also carries `useToolExecutionContext()` →
`{ app, session, projectId }`
(`packages/global-provider/src/hooks/use-tool-execution-context.ts`), plus
`clusterId` resolved the way `use-import-dataset.ts` does it (prefer a
`us-west-2` cluster, else the first).

**One tool execution per search, zero per import. One data-platform write per
import batch.**

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
├── use-uniprot-import-action.tsx   NEW  menu item + empty-state card + mounted modal
├── uniprot-import-modal.tsx        NEW  input → results → import
├── candidate-table.tsx             NEW  presentational DataTable over candidates[]
├── candidate-columns.tsx           NEW  column defs + grade/eye/star cells
├── display-maps.ts                 NEW  organism + method display names (§4.5)
├── structure-preview-modal.tsx     NEW  eye → floating Mol* window
├── use-uniprot-discovery.ts        NEW  execute + poll uniprot-discovery
├── use-import-candidates.ts        NEW  fetch → upload → batch.create (§4.7)
├── use-existing-pdb-ids.ts         NEW  which candidates are already in the project
└── rcsb.ts                         NEW  shared RCSB URL + fetch-to-text/File helper

apps/uui/src/utils/
└── run-tool-and-wait.ts            NEW  promise-based execute+poll (§4.2)
```

Changed files:

| File | Change |
|---|---|
| `components/data-platform-tables/chemical-table.tsx` | mount the modal (portaled) + an inline `ActionMenuPlugin` and an `emptyTableActions` card, both gated on `isProteinTable` |
| `app-engine/available-components.ts` | add `uniprot_accession` to the proteins table's default columns |
| `packages/data-table/src/hooks/use-action-import-from-rcsb.tsx` | export its RCSB URL builder (or move it into the new `rcsb.ts`) so both import paths agree on one URL |
| `app-schemas/uniprot-discovery.json` | metadata-only manifest (`id`, `toolKey`, `name`, `identityHue`) registered in `manifestByToolKey` only, so Activity rows render an identity instead of a raw tool key — mirrors `system-prep.json`, which is in `manifestByToolKey` but absent from `appManifests` |

### 4.2 `run-tool-and-wait.ts`

`usePollToolExecution` (`components/datasets/import-modal/use-poll-tool-execution.ts`)
is a react-query hook holding one execution id in state. A promise is easier to
drive from the modal's state machine and is reusable:

```ts
runToolAndWait({ orgKey, toolKey, body, signal, intervalMs = 3000, timeoutMs })
  → { executionId, jobOutputs }   // rejects on Failed/Cancelled/InsufficientFunds/timeout
```

- `fetchToolsExecuteTool` → read `executionId ?? id` off the response.
- Poll `fetchToolsGetToolExecution` until terminal. Status semantics are settled
  in `do-dd-client` (`src/platform/constants.py`): success is `Completed` **or**
  legacy `Succeeded`; `Failed`, `Cancelled`, `InsufficientFunds`,
  `FailedQuotation`, `Quoted` are the other terminals.
- Normalize `jobOutputs` for the dict-or-list shape (§2.1).
- Surface `statusReason` on failure.

Leave `use-poll-tool-execution.ts` where it is — the datasets import is not part
of this change.

### 4.3 `uniprot-import-modal.tsx`

One Mantine `Modal`, two visual states in the same shell (§1.1 → §1.2); the
input stays mounted so a second search replaces the result set in place.

```
idle ──Enter──► discovering ──ok──► results(candidates[])
                     │                   │
                     │ error             ├─ empty  → "No experimental structures for {ACC}"
                     ▼                   ├─ preview(pdbId) → structure-preview-modal
                   error                 └─ Import Selected → importing(n) ─► done → close + toast
                  (retry)                                                └─► partial → stay open
```

- **Validate the accession before the call.** Mirror `_UNIPROT_ACCESSION_RE`
  (`uniprot_discovery.py`):
  `^(?:[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z0-9]{2}\d){1,2})$`,
  case-insensitive, uppercased before sending. Reject isoform suffixes
  (`P00533-2`) with a specific message.
- Enter submits (`onKeyDown`); the input is disabled while `discovering`, and a
  new search aborts the previous poll.
- The recommended row starts checked (the mock shows the top row selected).
- **Import Selected** disabled at 0 selected; shows `n/m` while importing.
- Closing mid-import is allowed: imports continue, the toast reports the
  outcome. Closing mid-search aborts the poll (the execution still lands in
  Activity).

### 4.4 `candidate-table.tsx` + `candidate-columns.tsx`

Client-side `DataTable` from `@platform-ui/data-table` — the same entry point
`CsvTable` uses, so selection, sorting and styling match the rest of DO Studio.

| Column | Field | Cell |
|---|---|---|
| ☆ | `recommended` | star on the recommended row only; tooltip "Recommended structure"; row also gets a tinted background |
| ☑ | – | `rowSelection: multiRow`, `headerCheckbox`, `enableClickSelection: false` |
| ID | `pdb_id` | uppercase, monospace |
| Organism | `organism` → common name | §4.5 |
| Method | `method_class` | §4.5 |
| Resolution | `resolution` | `1.8Å` — one decimal, no space (mock) |
| Coverage | `coverage` | `83%` — `Math.round(v * 100)`; tool emits 0–1 |
| Rfree | `rfree` | `0.231` — three decimals |
| Grade | `grade` | bare coloured letter: A green, B amber, C orange, D red |
| 👁 | – | eye `ActionIcon` (`variant="eye"` from `@platform-ui/icons`) → preview |

- **Sort** `grade` A→D, then `weighted_score` desc, in the component. The tool
  already returns ranked rows, but the recommended row must not sit on top only
  because the server happened to order it there.
- **Never recompute the grade** — weights and thresholds belong to the Structure
  Report tool; the client displays what it is given.
- **`field_status` drives empty cells**: `not_applicable` → `—` ("not applicable
  for this method"); `unknown` → `?` ("not reported"). Rfree on a cryo-EM entry
  is `not_applicable`, not missing data — a blank cell that could mean either is
  the thing to avoid.
- Fixed body height with vertical scroll (the mock shows 6 rows and a
  scrollbar).
- Props only — `candidates`, `existingPdbIds`, `selected`, `onSelectionChange`,
  `onPreview` — so Target Preparation's expert mode can mount the same table when
  it needs a structure chosen.

### 4.5 `display-maps.ts`

The mock's **Organism** column reads `Human` / `Rat` / `Mouse`. Neither tool
field gives that: `organism` is the scientific name (`Homo sapiens`) and
`organism_class` is a scoring bucket (`human` / `mammal` / `vertebrate` /
`other`) — rat and mouse are both `mammal`. So a scientific → common name lookup
is needed, falling back to the italic scientific name:

```
Homo sapiens → Human, Rattus norvegicus → Rat, Mus musculus → Mouse,
Bos taurus → Bovine, Sus scrofa → Pig, Oryctolagus cuniculus → Rabbit,
Canis lupus familiaris → Dog, Macaca mulatta → Rhesus, Gallus gallus → Chicken,
Danio rerio → Zebrafish, Drosophila melanogaster → Fruit fly,
Caenorhabditis elegans → C. elegans, Saccharomyces cerevisiae → Yeast,
Escherichia coli → E. coli
```

Full scientific name in the tooltip either way.

**Method** maps off `method_class`: `x-ray → X-ray`, `cryo-em → Cryo-EM`,
`nmr → NMR`, else `Other`; raw `method` (`X-RAY DIFFRACTION`) as the tooltip.

### 4.6 `use-existing-pdb-ids.ts` — now load-bearing

With a client-composed import there is no `file_path` dedupe (§2.3), so this is
the only thing standing between a user and a duplicate protein row. Direct
analogue of `use-existing-smiles.ts` (`app-engine/renderer/csv-table/`):

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

Cursor-paginate with the same `collectMatches` helper, then reuse
`in-project-cell.tsx` / `inProjectRowSelection`: unselectable row, dimmed,
"Already in this project" tooltip. Re-check right before the write as well, so a
row added in another tab doesn't slip through.

### 4.7 `use-import-candidates.ts`

```ts
importCandidates(pdbIds: string[], accession: string) → { created: Row[]; failed: {pdbId, error}[] }
```

1. For each `pdbId`, bounded concurrency 3–4:
   `fetch(rcsbPdbUrl(pdbId))` → `new File([blob], '{PDBID}.pdb', { type: 'chemical/x-pdb' })`
   → `uploadFileToFileService({ orgKey, file, subdir: 'structure-imports' })`.
   A 404 here means an obsolete/withdrawn entry — record it in `failed` and keep
   going.
2. One `serverInterface.batch.create(orgKey, 'proteins', { rows, returning })`
   with the §2.3 columns per row (the `useAddToProject` precedent, one round
   trip instead of N). `returning` = the grid's active columns plus `id`,
   excluding read-only/computed columns (`READ_ONLY_COLUMNS`,
   `chemical-table.tsx:66`) — the backend rejects those in `returning`.
3. Return the created rows so the caller can insert them straight into the grid.

Failures are per-PDB and never abort the batch — the CLI raises on the first
error, but a user who checked six rows should keep the five that worked.

### 4.8 `structure-preview-modal.tsx`

Follow `packages/chat-engine/src/viewers/structure-viewer.tsx`, not the
app-engine `protein-viewer-wrapper` (that one is manifest-driven, resolves a
`file_path` through `use-structure-files.ts`, and expects overlays, pockets and
docking boxes). Here there is one structure and no entity yet:

```ts
const text = await fetchRcsbPdbText(pdbId);            // files.rcsb.org/download/{ID}.pdb
const viewer = new MolstarViewer(containerId);
await viewer.init({ showDockingBoxControls: false });
await viewer.api.loadFromRawContent(text, 'pdb', 'structure');
```

- Own container id per instance (`useId()` with `:` stripped) and a `live` flag,
  so a late load never draws into a torn-down container — both mistakes the
  chat-engine viewer already documents.
- Cache the fetched text per PDB ID for the modal's lifetime, and share that
  cache with the import step: previewing then importing the same structure
  should not fetch twice.
- Mantine `Modal` **stacked over the import modal** (higher `zIndex`, import
  modal stays mounted), title = PDB ID, `⤢` toggles `fullScreen`, `✕` closes.
  `⚙` needs no work — it is Mol\*'s own settings toggle.
- Failure state: "This structure could not be displayed."

### 4.9 `chemical-table.tsx` wiring

```tsx
// alongside importFromRcsbAction (:470)
const uniprotImport = useUniprotImportAction({
  orgKey, projectId, baseEntity, agentIdPrefix: agentTableId,
  onImported: handleImportedProteins,
});

// extraToolbarItems (:1061) — portal the modal like the RCSB one
if (isProteinTable) items.push(<Fragment key="uniprot-modal-inline">{uniprotImport.modal}</Fragment>);

// the Data menu (:1135) — an inline plugin, same shape as combinedUploadPlugin (:1121)
plugins={[ …, ...(isProteinTable ? [importFromRcsbAction, uniprotImport.plugin] : []) ]}

// emptyTableActions (:1034)
...(isProteinTable ? [uniprotImport.emptyTableAction] : [])
```

`ActionMenuPlugin` requires only `menuItems`
(`components/data-platform-tables/action-menu.tsx`); `EmptyTableAction` requires
`{ title, subtitle, onClick, icon, agentId? }`
(`packages/data-table/src/hooks/types.ts`).

Empty-state card: "Import from UniProt" / "Find experimental structures for a
UniProtKB accession", `icon: 'download'`.

### 4.10 Grid refresh after import

`batch.create` returns the created rows, so no read-back is needed:

```ts
gridApiRef.current?.applyServerSideTransaction({ add: created });
tableRef.current?.invalidateRowCount();
tableRef.current?.markNotEmpty();
queryClient.invalidateQueries({ queryKey: buildDataPlatformSearchQueryKey(orgKey, 'proteins') });
```

The last line is for other tiles (protein viewer, other tables). Reuse
`handleProteinFiles`'s notification wording (`chemical-table.tsx:374`) so
UniProt imports and file uploads read the same.

---

## 5. Errors and edge cases

| Case | Behaviour |
|---|---|
| Invalid accession shape | inline field error, no execution |
| Discovery fails / times out | inline error under the input + retry; keep the typed accession |
| `candidates: []` | empty state where the table would be: "No experimental PDB structures for {ACC}" — a success, not an error |
| `candidates` missing from `jobOutputs` | treat as a tool failure (matches the SDK, which raises) |
| RCSB 404 on a candidate | per-row failure, batch continues; row marked, summary toast says how many landed |
| Upload fails | same — per-row failure |
| `batch.create` fails wholesale | keep the modal open, surface `detail`, leave the uploaded files (harmless orphans) |
| PDB already in project | row unselectable + tooltip (§4.6) |
| Many selections | bounded concurrency 3–4 on fetch+upload; `n/m` on the button |
| Modal closed mid-import | imports continue; toast reports the result |
| Cost | none — the flow is free; no quote step, no `approveAmount` |

Activity page: one `uniprot-discovery` row per search. Nothing else shows up
there, since the import no longer runs a tool.

---

## 6. Analytics

`useCaptureEvent` (`@platform-ui/global-provider`), `EventName` in
`providers/analytics-provider.tsx`:

- `ToolExecutionStarted` / `ToolExecutionResult` for the discovery run
  (`toolKey`, `status`, `projectId`).
- `DataUploaded` on a successful import: `uploadType: 'uniprot'`, `recordCount`,
  `status`, `projectId` — the shape `handleRcsbImport` already uses
  (`chemical-table.tsx:451`).

---

## 7. Tests

| Level | What |
|---|---|
| Unit | accession regex (valid 6/10-char, isoform rejected); grade sort + recommended-first; organism/method display maps; `field_status` cell rendering; `jobOutputs` dict-vs-list normalization; terminal-status handling in `run-tool-and-wait`; per-row failure isolation in `importCandidates` |
| Component | modal state machine (idle → discovering → results → importing); empty candidates; partial import failure; preview opens/closes over the modal without losing selection; already-in-project rows unselectable |
| Table | extend `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx` — it already asserts the + Data menu composition and that protein-only plugins are absent on ligand tables |
| E2E | `apps/platform-e2e/src/pages/` — page-object entry next to `addDataButton`; mock the discovery execution and the RCSB fetch |

Mock payloads: copy `tests/mock_server/routers/tools.py:167` from this repo so UI
and SDK test against the same candidate shape, including the `P99999` empty case.

---

## 8. Open items

1. **Confirm no `pdb-import` tool exists** in the target org
   (`GET /tools/{orgKey}/tools`). The plan assumes client composition (§2.2); if
   the tool is live, swap it in behind `importCandidates` and drop the
   fetch/upload/batch-create step.
2. **`protein_length`** — leave unset, or count `CA` records client-side? The SDK
   sets it; nothing in the mock shows it.
3. **Already-in-project gating** (§4.6) is an addition to the mock — confirm the
   treatment (unselectable + tooltip) reads right.
4. **Selection cap** on Import Selected, if any.
