# Import from UniProt — DO Studio implementation plan

**Feature.** A new entry in the proteins table's **+ Data** menu in DO Studio
(`platform-ui/apps/uui`): enter a UniProtKB accession → see every experimental
PDB for it, graded and sorted → preview any of them in Mol\* → import the
selected one(s) into the project's proteins table.

PRD:
[UNIPROT Addition via Protein Table](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1067679762/UNIPROT+Addition+via+Protein+Table).
This document is the **UI plan only**. `do-dd-client` is the source of truth for
what the tools accept and return, and for what a correct import writes.

---

## 1. UI spec, from the mocks

### 1.1 Modal — empty state (`509×155`)

```
┌────────────────────────────────────┐
│ Import PDB Structure               │
│ ┌────────────────────────────────┐ │
│ │ e.g. P01116                    │ │   ← text input, placeholder
│ └────────────────────────────────┘ │
│ Press enter to search              │   ← helper text, dimmed, xs
└────────────────────────────────────┘
```

Title is **"Import PDB Structure"** — not "Import from UniProt". There is **no
Find button**: the mock's affordance is "Press enter to search". The PRD prose
says "hit *Find*"; the mock wins on layout, and Enter is the submit. Add a
subtle submit affordance in the input's `rightSection` (search icon button) so
the action is reachable by click too — it costs nothing and keeps the mock's
shape.

### 1.2 Modal — results state (`749×608`)

Same modal, grown; the input stays at the top so a new accession can be searched
in place.

```
┌──────────────────────────────────────────────────────────────┐
│ Import PDB Structure                                       ✕ │
│ ┌──────────────────────────────────────────────────────────┐ │
│ │ e.g. P01116                                              │ │
│ └──────────────────────────────────────────────────────────┘ │
│ Press enter to search                                        │
│ 6 structures available                                       │  ← bold, xs
│ ┌──────────────────────────────────────────────────────────┐ │
│ │    ID    Organism  Method   Resolution Coverage Rfree Grade│ │
│ │ ☆ ◉ 1ON1  Human    X-ray      1.8Å       83%    0.231  A 👁│ │ ← highlighted
│ │   ○ 1EZY  Rat      Cryo-EM    3.8Å       93%    0.333  A 👁│ │
│ │   ○ 4M35  Human    X-ray      2.5Å       91%    0.101  B 👁│ │
│ │   ○ 2D5G  Human    X-ray      2.1Å       80%    0.100  C 👁│ │
│ │   ○ 9EHU  Rat      Cryo-EM    2.7Å       79%    0.294  C 👁│ │
│ │   ○ 2DD6  Mouse    X-ray      2.0Å       66%    0.231  D 👁│ │ ← scrolls
│ └──────────────────────────────────────────────────────────┘ │
│                                          [ Import Selected ] │  ← primary
└──────────────────────────────────────────────────────────────┘
```

Settled by the mock:

- **Count line** `"{n} structures available"` above the table.
- **Columns**: `ID · Organism · Method · Resolution · Coverage · Rfree · Grade`,
  plus a leading star column, a leading selector column, and a trailing eye
  button. **No score column, no ligand column, no title/description column** —
  the tool returns `weighted_score` and `has_ligand`, and the mock deliberately
  leaves them out. Keep them out.
- **Recommended row** = star icon (`☆`) in the leading column **plus** a tinted
  row background, pre-selected. No "Recommended Structure" text badge. Put the
  wording in the star's tooltip.
- **Grade** is a bare coloured letter, not a pill: A green, B amber, C orange,
  D red.
- Table body **scrolls** at a fixed max height; header stays.
- Footer: single primary button **"Import Selected"**, bottom-right.
- Closing ✕ top-right.

**Open conflict — single vs multi select.** The mock draws **radio buttons**
(one filled, the rest empty) — single selection. The PRD requirement text says
"A user can select any number (e.g. at least 1) of the proteins to *Import to
Project*". The platform's import tool is one-PDB-per-execution, which is
consistent with either (N selections = N executions). Build **checkboxes
(multi-select)**: it satisfies the written requirement, degrades to the mock's
behaviour when one row is checked, and the footer label "Import Selected"
already reads correctly for both. **Confirm with product** — this is the one
place where implementing the mock literally would violate the requirement text.

### 1.3 Preview window (`759×660`)

Pressing the eye raises a **floating window layered above the import modal** —
the modal stays visible behind it, dimmed.

```
┌─────────────────────────────────────────┐
│ 1ON1                        ⚙  ⤢   ✕   │
├─────────────────────────────────────────┤
│                                         │
│        [ Mol* canvas — cartoon +        │
│          translucent surface ]          │
│                                         │
└─────────────────────────────────────────┘
```

- Title is the **PDB ID alone**.
- `⚙` is Mol\*'s own settings panel toggle — it already ships in
  `CustomViewportControls` (`packages/molstar/src/components/custom-viewport-controls.tsx`,
  `toggleSettingsPanel`), wired by `createViewerSpec`
  (`config/viewer-config.ts`: `controls.right = SettingsPanel`). It comes free
  with `MolstarViewer`; do not rebuild it. Note the same control group also
  renders screenshot and sequence-viewer buttons, which the mock does not show —
  either accept them or hide them in the preview's CSS.
- `⤢` toggles the preview window to full screen (Mantine `Modal fullScreen`).
- `✕` closes the preview and returns to the list; the modal's selection is
  untouched.

---

## 2. Tool contracts (from `do-dd-client`)

### 2.1 `deeporigin.uniprot-discovery` — ranking

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

### 2.2 `deeporigin.pdb-import` — coordinate acquisition + registration

**Not in `do-dd-client`'s `TOOL_KEYS_AND_VERSIONS`** — the CLI still composes the
import client-side. The contract below comes from the platform ticket that
shipped the tool; `platform-toolbox` is not reachable from this session, so
**confirm the exact key names against `tools/pdb-import/tool-definition.json`
before coding**.

```
POST /tools/{orgKey}/tools/deeporigin.pdb-import/executions
body: { inputs: { pdb_id, uniprot_accession?, project_id? }, outputs: {}, clusterId, projectId, name? }
```

- Behaviour: download from `files.rcsb.org` → upload to UFA → `create_protein`.
- Output: `protein_id`, `file_path`, `pdb_id`; persists `uniprot_accession` on
  the entity when supplied.
- **One PDB per execution** — no batching. N selected candidates = N executions.

This is the step that writes to the data platform. On this path the UI does not
create the protein row itself.

### 2.3 What a correct protein row looks like (`do-dd-client`)

`Protein.register` (`src/drug_discovery/structures/protein.py:1741`) writes:

| Column | Value |
|---|---|
| `file_path` | UFA path of the uploaded structure |
| `pdb_id` | the PDB ID |
| `uniprot_accession` | the accession the discovery ran on |
| `protein_name` | `from_pdb_id` sets it to the PDB ID |
| `protein_length` | residue count, when a local file was parsed |
| `project_id` | resolved project |

`Protein.sync` (`:1795`) is register plus a **dedupe**: it searches `proteins` by
`file_path` (+ `project_id`) and reuses the existing row rather than creating a
second one (`client.entities.search_proteins`, `src/platform/entities.py:684`,
filter `{deleted: false, file_path, project_id}`). `from_pdb_id` (`:280`)
downloads `https://files.rcsb.org/download/{id}.pdb`. Every column exists on the
entity — `PROTEIN_RETURNING_FIELDS`, `src/platform/entities.py:59`.

**Path B (fallback).** If `pdb-import` is not deployed where DO Studio runs, the
UI can mirror the SDK with what the proteins table already does: browser `fetch`
from RCSB → `uploadFileToFileService({ orgKey, file, subdir: 'structure-imports' })`
→ `serverInterface.entity.create(orgKey, 'proteins', { set: { file_path, pdb_id,
uniprot_accession, protein_name, project_id, tags: { app } }, returning })`.
That is `handleProteinFiles` (`chemical-table.tsx:374`) plus two fields. Keep the
import step behind one interface (`importCandidate(pdbId, accession)`) so A and B
are swappable. **Decide A vs B before building the import step** — nothing else
in the feature changes.

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

 Import Selected ─────►  per selected pdb_id:
                         fetchToolsExecuteTool           POST /tools/{orgKey}/tools/
                                                              deeporigin.pdb-import/executions
                         poll fetchToolsGetToolExecution GET  /tools/{orgKey}/tools/executions/{id}
                         → { protein_id, file_path }     ← the tool created the proteins row

 refresh ─────────────►  serverInterface.entity.get      GET  /{orgKey}/proteins/{protein_id}
                         grid.applyServerSideTransaction
                         tableRef.invalidateRowCount()
```

Every execution body also carries `useToolExecutionContext()` →
`{ app, session, projectId }`
(`packages/global-provider/src/hooks/use-tool-execution-context.ts`), plus
`clusterId` resolved the way `use-import-dataset.ts` does it (prefer a
`us-west-2` cluster, else the first).

**Data-platform writes: exactly one per structure, and the tool makes it.** On
path A the UI never calls `entity.create` — it only reads (dedupe + read-back).
On path B the UI makes the write and must reproduce §2.3.

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
├── uniprot-import-modal.tsx        NEW  "Import PDB Structure": input → results → import
├── candidate-table.tsx             NEW  presentational DataTable over candidates[]
├── candidate-columns.tsx           NEW  column defs + grade/eye/star cells
├── display-maps.ts                 NEW  organism + method display names (§4.5)
├── structure-preview-modal.tsx     NEW  eye → floating Mol* window
├── use-uniprot-discovery.ts        NEW  execute + poll uniprot-discovery
├── use-pdb-import.ts               NEW  execute + poll pdb-import, once per selection
├── use-existing-pdb-ids.ts         NEW  which candidates are already in the project
└── rcsb.ts                         NEW  shared RCSB URL + fetch-to-text/File helper

apps/uui/src/utils/
└── run-tool-and-wait.ts            NEW  promise-based execute+poll (§4.2)
```

Changed files:

| File | Change |
|---|---|
| `components/data-platform-tables/chemical-table.tsx` | mount the modal (portaled) + add an inline `ActionMenuPlugin` and an `emptyTableActions` card, both gated on `isProteinTable` |
| `app-engine/available-components.ts` | add `uniprot_accession` to the proteins table's default columns |
| `packages/data-table/src/hooks/use-action-import-from-rcsb.tsx` | export its RCSB URL builder (or move it into the new `rcsb.ts`) so both import paths agree on one URL |
| `app-schemas/uniprot-discovery.json`, `app-schemas/pdb-import.json` | metadata-only manifests (`id`, `toolKey`, `name`, `identityHue`) registered in `manifestByToolKey` only, so Activity rows render an identity instead of a raw tool key — mirrors `system-prep.json`, which is in `manifestByToolKey` but absent from `appManifests` |

### 4.2 `run-tool-and-wait.ts` — why a new helper

`usePollToolExecution` (`components/datasets/import-modal/use-poll-tool-execution.ts`)
is a react-query hook holding **one** execution id in state. The import step runs
N executions in a loop, so it needs a promise:

```ts
runToolAndWait({ orgKey, toolKey, body, signal, intervalMs = 3000, timeoutMs })
  → { executionId, jobOutputs }   // rejects on Failed/Cancelled/InsufficientFunds/timeout
```

- `fetchToolsExecuteTool` → read `executionId ?? id` off the response.
- Poll `fetchToolsGetToolExecution` until terminal. The status semantics are
  settled in `do-dd-client` (`src/platform/constants.py`): success is
  `Completed` **or** legacy `Succeeded`; `Failed`, `Cancelled`,
  `InsufficientFunds`, `FailedQuotation`, `Quoted` are the other terminals.
- Normalize `jobOutputs` for the dict-or-list shape (§2.1).
- Surface `statusReason` on failure.

Use it for both tool calls. Leave `use-poll-tool-execution.ts` where it is — the
datasets import is not part of this change.

### 4.3 `uniprot-import-modal.tsx`

One Mantine `Modal`, two visual states in the same shell, matching §1.1 → §1.2
(the mock keeps the input mounted, so a second search replaces the result set in
place).

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
  (`P00533-2`) with a specific message — a typo must never cost an execution.
- Enter submits (`onKeyDown`), matching the mock's "Press enter to search"; the
  input is disabled while `discovering`, and a new search aborts the previous
  poll.
- The recommended row starts selected (mock shows the top row pre-selected).
- **"Import Selected"** is disabled with 0 selected; shows `n/m` progress while
  importing.
- Closing the modal mid-import is allowed: the imports continue and the toast
  reports the outcome. Closing mid-search aborts the poll (the execution itself
  still lands in Activity).

### 4.4 `candidate-table.tsx` + `candidate-columns.tsx`

Client-side `DataTable` from `@platform-ui/data-table` — the same entry point
`CsvTable` uses, so selection, sorting and styling match the rest of DO Studio.

| Column | Field | Cell |
|---|---|---|
| ☆ | `recommended` | star icon on the recommended row only; tooltip "Recommended structure"; the row also gets a tinted background |
| ◉/☑ | – | selector: `rowSelection: multiRow`, `headerCheckbox`, `enableClickSelection: false` (see the single-vs-multi conflict in §1.2) |
| ID | `pdb_id` | uppercase, monospace |
| Organism | `organism` → common name | §4.5 |
| Method | `method_class` | §4.5 |
| Resolution | `resolution` | `1.8Å` — one decimal, no space (mock) |
| Coverage | `coverage` | `83%` — `Math.round(v * 100)`, tool emits 0–1 |
| Rfree | `rfree` | `0.231` — three decimals |
| Grade | `grade` | bare coloured letter: A green, B amber, C orange, D red |
| 👁 | – | eye `ActionIcon`, `variant="eye"` from `@platform-ui/icons`; opens the preview |

- **Sort** `grade` A→D, then `weighted_score` desc, in the component. The tool
  already returns ranked rows, but the recommended row must not sit on top only
  because the server happened to order it there.
- **Never recompute the grade** — the weights and thresholds belong to the
  Structure Report tool; the client displays what it is given.
- **`field_status` drives empty cells**: `not_applicable` → `—` with tooltip
  "not applicable for this method"; `unknown` → `?` with tooltip "not reported".
  Rfree on a cryo-EM entry is `not_applicable`, not missing data — a blank cell
  that could mean either is the thing to avoid.
- Fixed body height with vertical scroll (mock shows 6 rows and a scrollbar).
- Props only: `candidates`, `existingPdbIds`, `selected`, `onSelectionChange`,
  `onPreview` — no execution machinery — so Target Preparation's expert mode can
  mount the same table when it needs a structure chosen.

### 4.5 `display-maps.ts`

The mock's **Organism** column reads `Human` / `Rat` / `Mouse`. Neither tool
field gives that directly: `organism` is the scientific name
(`Homo sapiens`) and `organism_class` is a score bucket
(`human` / `mammal` / `vertebrate` / `other`) — rat and mouse both fall in
`mammal`. So the UI needs a small scientific-name → common-name lookup, falling
back to the scientific name (italic) when unmapped:

```
Homo sapiens → Human, Rattus norvegicus → Rat, Mus musculus → Mouse,
Bos taurus → Bovine, Sus scrofa → Pig, Oryctolagus cuniculus → Rabbit,
Canis lupus familiaris → Dog, Macaca mulatta → Rhesus, Gallus gallus → Chicken,
Danio rerio → Zebrafish, Drosophila melanogaster → Fruit fly,
Caenorhabditis elegans → C. elegans, Saccharomyces cerevisiae → Yeast,
Escherichia coli → E. coli
```

Keep the full scientific name in a tooltip either way.

**Method** maps straight off `method_class`: `x-ray → X-ray`,
`cryo-em → Cryo-EM`, `nmr → NMR`, else `Other`; raw `method`
(`X-RAY DIFFRACTION`) as the tooltip.

### 4.6 `use-existing-pdb-ids.ts` — beyond the mock, worth having

The mock has no "already in this project" state, but nothing else prevents a
user importing the same structure twice (the CLI gets that free from
`Protein.sync`'s `file_path` dedupe; on path A the UI is the only place that
can). Direct analogue of `use-existing-smiles.ts`
(`app-engine/renderer/csv-table/`):

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
"Already in this project" tooltip. **Flag for product** — it is an addition to
the mock, not a contradiction of it.

### 4.7 `structure-preview-modal.tsx`

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
- Cache the fetched text per PDB ID for the modal's lifetime; re-opening a row
  must not refetch.
- Rendered as a Mantine `Modal` **stacked over the import modal** (higher
  `zIndex`, import modal stays mounted), title = PDB ID, `⤢` toggles
  `fullScreen`, `✕` closes. `⚙` needs no work — it is Mol\*'s own settings
  toggle from `CustomViewportControls`.
- Failure state: "This structure could not be displayed." — the RCSB fetch can
  404 on obsolete entries.

### 4.8 `chemical-table.tsx` wiring

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

Copy: menu item **"Import from UniProt"** (the modal itself is titled "Import PDB
Structure" per the mock); empty-state card "Import from UniProt" / "Find
experimental structures for a UniProtKB accession", `icon: 'download'`.

### 4.9 Grid refresh after import

The row is created by the tool, so there is no response row to insert. Per
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
| Discovery fails / times out | inline error under the input + retry; keep the typed accession |
| `candidates: []` | empty state where the table would be: "No experimental PDB structures for {ACC}" — a success, not an error |
| `candidates` missing from `jobOutputs` | treat as a tool failure (matches the SDK, which raises) |
| Some imports fail | do **not** abort the batch — the CLI raises on the first failure, but a user who checked six rows should keep the five that worked. Per-row status in the table, one summary toast |
| PDB already in project | row unselectable + tooltip (§4.6) |
| Many selections | N executions, bounded concurrency (3–4), `n/m imported` on the button |
| Modal closed mid-import | imports continue; toast reports the result |
| Cost | each import is a billable execution. The SDK exposes `quote=True` / `approveAmount`; the PRD is silent — **product decision** |

Activity page: both tools appear there per import. If `pdb-import` runs are
noise, filter them with the existing mechanism —
`staticFilterProps={[{ column: 'tool_key', op: 'neq', value: 'deeporigin.pdb-import' }]}`
on the activity `JobManagerTable`.

---

## 6. Analytics

`useCaptureEvent` (`@platform-ui/global-provider`), `EventName` in
`providers/analytics-provider.tsx`:

- `ToolExecutionStarted` / `ToolExecutionResult` per execution (`toolKey`,
  `status`, `projectId`).
- `DataUploaded` on a successful import: `uploadType: 'uniprot'`, `recordCount`,
  `status`, `projectId` — the shape `handleRcsbImport` already uses
  (`chemical-table.tsx:451`).

---

## 7. Tests

| Level | What |
|---|---|
| Unit | accession regex (valid 6/10-char, isoform rejected); grade sort + recommended-first; organism/method display maps; `field_status` cell rendering; `jobOutputs` dict-vs-list normalization; terminal-status handling in `run-tool-and-wait` |
| Component | modal state machine (idle → discovering → results → importing); empty candidates; partial import failure; preview opens/closes over the modal without losing selection |
| Table | extend `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx` — it already asserts the + Data menu composition and that protein-only plugins are absent on ligand tables |
| E2E | `apps/platform-e2e/src/pages/` — page-object entry next to `addDataButton`; mock both executions |

Mock payloads: copy `tests/mock_server/routers/tools.py:167` from this repo so UI
and SDK test against the same candidate shape, including the `P99999` empty case.

---

## 8. Open items

1. **Single vs multi select** (§1.2) — mock draws radios, requirement text says
   "any number". Plan assumes checkboxes.
2. **`pdb-import` availability and exact schema** — `platform-toolbox` is not
   reachable from this session; confirm
   `tools/pdb-import/tool-definition.json` key names and that the tool is
   deployed where DO Studio runs. This decides path A vs path B (§2.3).
3. **Cost/quote behaviour** for the import step.
4. **Already-in-project gating** (§4.6) — an addition to the mock.
5. **Selection cap** on Import Selected, if any.
