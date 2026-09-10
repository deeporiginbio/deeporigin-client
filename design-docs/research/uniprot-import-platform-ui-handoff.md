# Handoff prompt — Import from UniProt (platform-ui)

Paste everything below the line as the first message to an agent working in
`deeporiginbio/platform-ui`. It is self-contained: it does not require access to
`do-dd-client`, Confluence, or Jira.

---

Implement **Import from UniProt** in DO Studio end-to-end.

## What you are building

DO Studio (`apps/uui`) shows a Proteins table with a **+ Data** menu. Add a new
menu item, **Import from UniProt**, that opens a modal where the user types a
UniProtKB accession, sees every experimental PDB structure for that accession
(graded and ranked by a platform tool), previews any of them in Mol\*, and
imports the selected ones as rows in the project's `proteins` table.

Ship it complete: UI, tool call, import, grid refresh, analytics, tests.

## Before you write code

Read these, in order:

1. `AGENTS.md` and `apps/uui/CLAUDE.md` (App Engine, data-platform tables,
   submit flow) and `packages/data-table/CLAUDE.md` (action-hook pattern,
   server filter contract).
2. `apps/uui/src/components/data-platform-tables/chemical-table.tsx` — the
   Proteins/Ligands table. Key anchors: `isProteinTable` (~:350),
   `handleProteinFiles` (~:374), `handleRcsbImport` (~:451),
   `importFromRcsbAction` (~:470), `emptyTableActions` (~:1034),
   `extraToolbarItems` (~:1061), the inline `combinedUploadPlugin` (~:1121),
   the `ActionMenu title="Data"` with `testId="add-data-menu"` (~:1135),
   `READ_ONLY_COLUMNS` (~:66).
3. `packages/data-table/src/hooks/use-action-import-from-rcsb.tsx` — the
   existing "Import from RCSB" modal. Same shape of feature, one input.
4. `apps/uui/src/components/datasets/import-modal/use-import-dataset.ts` and
   `use-poll-tool-execution.ts` — how this app runs a tools-service tool from a
   modal and polls it.
5. `apps/uui/src/app-engine/renderer/csv-table/` — `csv-table.tsx`,
   `use-add-to-project.ts`, `use-existing-smiles.ts`, `in-project-cell.tsx`.
   This is the closest existing "results table with checkboxes that writes rows
   into the project" UI. Reuse its patterns.
6. `packages/chat-engine/src/viewers/structure-viewer.tsx` — the minimal
   Mol\* viewer for "one structure in a modal".

Line numbers are from a recent `main`; find the symbols, don't trust the
numbers.

---

## 1. UI specification

Copy strings below are final unless marked otherwise.

### 1.1 Menu entry

- Item label: **"Import from UniProt"**, in the proteins table's **+ Data**
  (`ActionMenu title="Data"`) menu, gated on `isProteinTable` — it must not
  appear on the ligands table.
- Also add an empty-state card (`emptyTableActions`): title "Import from
  UniProt", subtitle "Find experimental structures for a UniProtKB accession",
  `icon: 'download'`.

### 1.2 Modal — empty state (~509×155 in the design)

```
┌────────────────────────────────────┐
│ Import from UniProt                │
│ ┌────────────────────────────────┐ │
│ │ e.g. P01116                    │ │   ← text input
│ └────────────────────────────────┘ │
│ Press enter to search              │   ← helper text, dimmed, xs
└────────────────────────────────────┘
```

- Title: "Import from UniProt". Input aria-label/label: "UniProt accession".
  Placeholder: `e.g. P01116`. Helper text: "Press enter to search".
- **Enter submits.** Also put a search icon button in the input's
  `rightSection` so it is clickable.
- Uppercase the input as the user types (the RCSB modal does this).

### 1.3 Modal — results state (~749×608)

Same modal, grown. The input stays mounted at the top so a second search
replaces the results in place.

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
│ │ ☆ ☑ 1ON1  Human    X-ray      1.8Å       83%    0.231  A 👁│ │ ← tinted row
│ │   ☐ 1EZY  Rat      Cryo-EM    3.8Å       93%    0.333  A 👁│ │
│ │   ☐ 4M35  Human    X-ray      2.5Å       91%    0.101  B 👁│ │
│ │   ☐ 2D5G  Human    X-ray      2.1Å       80%    0.100  C 👁│ │
│ │   ☐ 9EHU  Rat      Cryo-EM    2.7Å       79%    0.294  C 👁│ │
│ │   ☐ 2DD6  Mouse    X-ray      2.0Å       66%    0.231  D 👁│ │ ← body scrolls
│ └──────────────────────────────────────────────────────────┘ │
│                                          [ Import Selected ] │  ← primary
└──────────────────────────────────────────────────────────────┘
```

- Count line above the table: `"{n} structures available"`.
- Table body scrolls at a fixed max height; the header stays.
- Footer: one primary button, **"Import Selected"**, bottom-right. Disabled at
  zero selected. While importing, show progress (`Importing 2/5…`).
- Use `DataTable` from `@platform-ui/data-table` (client-side row data), the
  same entry point `CsvTable` uses.

Columns:

| Column | Source field | Rendering |
|---|---|---|
| ☆ | `recommended` | star icon on the recommended row only, tooltip "Recommended structure"; that row also gets a tinted background and starts **checked** |
| ☑ | – | `rowSelection: { mode: 'multiRow', headerCheckbox: true, enableClickSelection: false }` |
| ID | `pdb_id` | uppercase, monospace |
| Organism | `organism` | common name, see §1.5 |
| Method | `method_class` | display name, see §1.5 |
| Resolution | `resolution` | `1.8Å` — one decimal, no space before Å |
| Coverage | `coverage` | `83%` — `Math.round(value * 100)`; the tool emits a 0–1 fraction |
| Rfree | `rfree` | `0.231` — three decimals |
| Grade | `grade` | bare coloured letter, no pill: **A** green, **B** amber, **C** orange, **D** red |
| 👁 | – | `ActionIcon` with `<Icon variant="eye" />` from `@platform-ui/icons`; opens the preview |

Do **not** add a score column or a ligand column. The tool returns
`weighted_score` and `has_ligand`; the design deliberately omits them.

**Sorting.** Sort in the component: `grade` A→D, then `weighted_score`
descending. The tool already returns ranked rows, but the recommended row must
not be on top only because the server happened to order it that way.

**Never recompute the grade.** Grade weights and thresholds are owned by the
platform's Structure Report tool. Display what you are given.

**Empty cells are meaningful.** Each candidate carries `field_status`, a map of
field name → `"value" | "not_applicable" | "unknown"`. Render `not_applicable`
as `—` with tooltip "not applicable for this method" and `unknown` as `?` with
tooltip "not reported". Rfree on a cryo-EM entry is *not applicable*, not
missing — a blank cell that could mean either is the bug to avoid.

### 1.4 Mol\* preview (~759×660)

The eye button raises a floating window **layered over the import modal** (the
modal stays mounted behind it):

```
┌─────────────────────────────────────────┐
│ 1ON1                        ⚙  ⤢   ✕   │
├─────────────────────────────────────────┤
│        [ Mol* canvas ]                  │
└─────────────────────────────────────────┘
```

- Title is the bare PDB ID. `⤢` toggles the Mantine `Modal`'s `fullScreen`.
  `✕` closes it and leaves the selection untouched.
- `⚙` is **Mol\*'s own settings-panel toggle** — it already ships in
  `packages/molstar/src/components/custom-viewport-controls.tsx`
  (`toggleSettingsPanel`) and is wired by `createViewerSpec`
  (`packages/molstar/src/config/viewer-config.ts`, `controls.right =
  SettingsPanel`). Do not rebuild it. Note that the same control group also
  renders screenshot and sequence-viewer buttons that the design does not show —
  either accept them or hide them with CSS scoped to this modal.
- Build it like `packages/chat-engine/src/viewers/structure-viewer.tsx`, **not**
  like `apps/uui/src/app-engine/renderer/protein-viewer/protein-viewer-wrapper.tsx`
  (that one is manifest-driven, resolves a `file_path` through
  `use-structure-files.ts`, and expects overlays, pockets and docking boxes;
  none of that applies to a structure that is not imported yet):

```ts
const text = await fetchRcsbPdbText(pdbId);   // https://files.rcsb.org/download/{ID}.pdb
const viewer = new MolstarViewer(containerId);
await viewer.init({ showDockingBoxControls: false });
await viewer.api.loadFromRawContent(text, 'pdb', 'structure');
```

- Give each instance its own container id (`useId()` with `:` stripped) and use
  a `live` flag in the effect so a late load never draws into a torn-down
  container — both mistakes the chat-engine viewer documents.
- Cache fetched PDB text per PDB ID and share that cache with the import step:
  previewing then importing the same structure must not fetch twice.
- Failure copy: "This structure could not be displayed." (obsolete entries 404).

### 1.5 Display maps (new file)

**Organism.** The design shows `Human` / `Rat` / `Mouse`. Neither tool field
gives that directly: `organism` is the scientific name (`Homo sapiens`) and
`organism_class` is a scoring bucket (`human` / `mammal` / `vertebrate` /
`other`) in which rat and mouse are both `mammal`. Add a scientific → common
name lookup, falling back to the scientific name in italics, with the full
scientific name always in the tooltip:

```
Homo sapiens → Human            Rattus norvegicus → Rat
Mus musculus → Mouse            Bos taurus → Bovine
Sus scrofa → Pig                Oryctolagus cuniculus → Rabbit
Canis lupus familiaris → Dog    Macaca mulatta → Rhesus
Gallus gallus → Chicken         Danio rerio → Zebrafish
Drosophila melanogaster → Fruit fly    Caenorhabditis elegans → C. elegans
Saccharomyces cerevisiae → Yeast       Escherichia coli → E. coli
```

**Method.** From `method_class`: `x-ray → X-ray`, `cryo-em → Cryo-EM`,
`nmr → NMR`, anything else → `Other`. Put the raw `method`
(`X-RAY DIFFRACTION`) in the tooltip.

---

## 2. Backend contract

### 2.1 The search — `deeporigin.uniprot-discovery`

One tools-service execution per search.

```
POST /tools/{orgKey}/tools/deeporigin.uniprot-discovery/executions
body: {
  inputs: { uniprot_accession: "P00533" },   // the entire input schema
  outputs: {},
  clusterId,
  projectId,
  name: `UniProt ${accession}`,              // shows in Activity
  ...useToolExecutionContext()               // { app, session, projectId }
}
```

- Use the generated `fetchToolsExecuteTool` from `@platform-ui/api-gateway`.
- Resolve `clusterId` the way `use-import-dataset.ts` does: `useToolsListClusters<any>`,
  prefer a cluster whose name includes `us-west-2`, else the first.
- The response carries `executionId` (older shapes: `id`).

Poll `fetchToolsGetToolExecution` (`GET /tools/{orgKey}/tools/executions/{executionId}`)
every 3s until the status is terminal:

- **Success**: `Completed` **or** the legacy `Succeeded`.
- **Failure**: `Failed`, `Cancelled`, `InsufficientFunds`, `FailedQuotation`,
  `Quoted`. Surface `statusReason`.

Then read `jobOutputs.candidates`.

> **`jobOutputs` can be a dict *or* a single-element array containing the
> dict.** The platform SDK normalizes both; the existing
> `normalizeJobOutputs` in `use-poll-tool-execution.ts` only accepts a plain
> object. Handle both shapes.

Each candidate:

| Field | Type | Always present |
|---|---|---|
| `pdb_id` | string (4 alphanumeric chars) | ✔ |
| `grade` | `"A" \| "B" \| "C" \| "D"` | ✔ |
| `weighted_score` | number | ✔ |
| `recommended` | boolean — exactly one `true` unless the list is empty | ✔ |
| `coverage_score`, `inhibitor_score`, `method_score`, `organism_score`, `resolution_score`, `rfree_score` | number | ✔ |
| `field_status` | `{ [field]: "value" \| "not_applicable" \| "unknown" }` | ✔ |
| `coverage` | number (0–1) | – |
| `has_ligand` | boolean | – |
| `method`, `method_class` | string | – |
| `organism`, `organism_class` | string | – |
| `resolution` | number (Å) | – |
| `rfree` | number | – |

Sample row for fixtures:

```json
{
  "pdb_id": "1M17", "grade": "A", "weighted_score": 0.9, "recommended": true,
  "coverage_score": 0.9, "inhibitor_score": 1.0, "method_score": 0.95,
  "organism_score": 1.0, "resolution_score": 0.75, "rfree_score": 0.8,
  "field_status": { "coverage": "value", "inhibitor": "value", "method": "value",
                    "organism": "value", "resolution": "value", "rfree": "value" },
  "coverage": 0.9, "has_ligand": true, "method": "X-RAY DIFFRACTION",
  "method_class": "x-ray", "organism": "Homo sapiens", "organism_class": "human",
  "resolution": 1.5, "rfree": 0.2
}
```

**`candidates: []` is a success**, not an error — it means the accession is
unknown or has no experimental structures. Show an empty state where the table
would be: `No experimental PDB structures for {ACCESSION}`. A **missing**
`candidates` key *is* a tool failure.

**Validate the accession before spending an execution.** UniProtKB accessions
are 6 or 10 characters:

```
/^(?:[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z0-9]{2}\d){1,2})$/i
```

Uppercase before sending. Reject isoform suffixes (`P00533-2`) with a specific
message rather than letting the tool fail.

### 2.2 The import — no tool; the client composes it

There is **no served tool that downloads a PDB by ID**. `uniprot-discovery`
ranks and never writes; `structure-report` grades and never downloads;
`import-dataset` is the catalog bulk path (needs `csv_path`, `database_key`,
`dataset_schema`); `protein-prep` operates on an already-registered protein. A
`deeporigin.pdb-import` tool was specified but is not registered.

So the import is client-side, exactly the way `handleProteinFiles` +
`useActionImportFromRcsb` already import a PDB in this same table today:

For each selected candidate (bounded concurrency, 3–4 at a time):

1. `fetch('https://files.rcsb.org/download/{PDBID}.pdb')` → `Blob` →
   `new File([blob], '{PDBID}.pdb', { type: 'chemical/x-pdb' })`.
   Reuse the URL builder from `use-action-import-from-rcsb.tsx` (export it, or
   move it into a shared `rcsb.ts` and import it from both) so the two paths
   can never disagree.
2. `uploadFileToFileService({ orgKey, file, subdir: 'structure-imports' })`
   (`apps/uui/src/components/data-platform-tables/structure-file-upload/`) →
   `{ filePath }`.

Then **one** write for the whole batch:

```ts
serverInterface.batch.create(orgKey, 'proteins', {
  rows: selected.map(({ pdbId, filePath }) => ({
    file_path: filePath,
    pdb_id: pdbId,                    // uppercase
    uniprot_accession: accession,     // uppercase
    protein_name: pdbId,
    project_id: projectId,
    tags: { app },                    // from useToolExecutionContext()
  })),
  returning: [...activeGridColumns, 'id'],
});
```

- `tags: { app }` matters — without it the row never matches the dashboard's
  Source filter (this is why the other direct write paths in
  `chemical-table.tsx` set it).
- Exclude read-only/computed columns from `returning` (see `READ_ONLY_COLUMNS`
  in `chemical-table.tsx`); the backend rejects them.
- `batch.create` returns the created rows, so no read-back is needed.

Wrap this as `importCandidates(pdbIds, accession) → { created, failed }` so the
whole mechanism can later be swapped for a served tool without touching the UI.

**Per-candidate failures must not abort the batch.** A 404 from RCSB (obsolete
entry) or a failed upload marks that row failed; every other selection still
lands.

### 2.3 Duplicate protection

The platform SDK dedupes imports by looking up an existing row with the same
`file_path`. That does not work here: `uploadFileToFileService` builds a
timestamped path (`structure-imports_{Date.now()}/{name}`), so the same
structure uploaded twice never collides. **Dedupe on `pdb_id` + `project_id`
instead**, before uploading.

Copy `use-existing-smiles.ts` (`app-engine/renderer/csv-table/`) into a
`use-existing-pdb-ids.ts`:

```ts
serverInterface.search(orgKey, 'proteins', {
  select: ['pdb_id'],
  filter: { props: [
    { column: 'project_id', op: 'eq', value: projectId },
    { column: 'pdb_id',     op: 'in', value: candidatePdbIds },
  ]},
  limit: candidatePdbIds.length,
});
```

Cursor-paginate with the same `collectMatches` helper. Then reuse
`in-project-cell.tsx` / `inProjectRowSelection`: the row is unselectable and
dimmed, with an "Already in this project" tooltip. Re-check immediately before
the write so a row added in another tab can't slip through.

### 2.4 Grid refresh

```ts
gridApiRef.current?.applyServerSideTransaction({ add: created });
tableRef.current?.invalidateRowCount();
tableRef.current?.markNotEmpty();
queryClient.invalidateQueries({ queryKey: buildDataPlatformSearchQueryKey(orgKey, 'proteins') });
```

Match the notification wording `handleProteinFiles` already uses so a UniProt
import and a file upload read the same.

---

## 3. Suggested file layout

```
apps/uui/src/components/data-platform-tables/uniprot-import/
├── index.ts
├── use-uniprot-import-action.tsx   menu item + empty-state card + mounted modal
├── uniprot-import-modal.tsx        input → results → import
├── candidate-table.tsx             presentational DataTable over candidates[]
├── candidate-columns.tsx           column defs + grade / star / eye cells
├── display-maps.ts                 organism + method display names (§1.5)
├── structure-preview-modal.tsx     eye → floating Mol* window
├── use-uniprot-discovery.ts        execute + poll the discovery tool
├── use-import-candidates.ts        fetch → upload → batch.create
├── use-existing-pdb-ids.ts         already-in-project lookup
└── rcsb.ts                         RCSB URL + fetch-to-text/File + per-id cache

apps/uui/src/utils/
└── run-tool-and-wait.ts            promise-based execute + poll
```

`run-tool-and-wait.ts` exists because `usePollToolExecution` is a react-query
hook holding one execution id in state; a promise is easier to drive from the
modal's state machine. Signature:

```ts
runToolAndWait({ orgKey, toolKey, body, signal, intervalMs = 3000, timeoutMs })
  → { executionId, jobOutputs }    // rejects on any non-success terminal state or timeout
```

Leave `use-poll-tool-execution.ts` alone — the datasets import is not part of
this change.

**Put all of this in `apps/uui`, not in `packages/data-table`.** The feature
needs the tools-service client, the data-platform provider and Mol\*, and the
grid package depends on none of those. An `ActionMenuPlugin` only needs
`{ menuItems }` (`components/data-platform-tables/action-menu.tsx`), and
`chemical-table.tsx` already builds one inline (`combinedUploadPlugin`) — follow
that, no new hook in the shared package.

Wiring in `chemical-table.tsx`:

```tsx
const uniprotImport = useUniprotImportAction({
  orgKey, projectId, baseEntity, agentIdPrefix: agentTableId,
  onImported: /* apply the grid transaction */,
});

// extraToolbarItems — portal the modal, like the RCSB one
if (isProteinTable) items.push(<Fragment key="uniprot-modal-inline">{uniprotImport.modal}</Fragment>);

// the Data ActionMenu
plugins={[ …, ...(isProteinTable ? [importFromRcsbAction, uniprotImport.plugin] : []) ]}

// emptyTableActions
...(isProteinTable ? [uniprotImport.emptyTableAction] : [])
```

Also add `uniprot_accession` to the proteins table's default columns in
`apps/uui/src/app-engine/available-components.ts`, so an imported row visibly
carries its accession.

Optional, nice for Activity: a metadata-only manifest
`apps/uui/src/app-schemas/uniprot-discovery.json` (`id`, `toolKey`, `name`,
`identityHue`) registered in `manifestByToolKey` **only** (not `appManifests`,
so it gets no route) — mirror `system-prep.json`. Without it, a UniProt run in
Activity renders with a raw tool key and no colour.

---

## 4. Modal state machine

```
idle ──Enter──► discovering ──ok──► results(candidates[])
                     │                   │
                     │ error             ├─ empty  → "No experimental PDB structures for {ACC}"
                     ▼                   ├─ preview(pdbId) → structure-preview-modal
                   error                 └─ Import Selected → importing(n) ─► done → close + toast
                  (retry)                                                 └─► partial → stay open
```

- Input disabled while `discovering`; a new search aborts the in-flight poll.
- Closing the modal mid-import is allowed: imports continue and the toast
  reports the outcome. Closing mid-search aborts the poll (the execution itself
  still completes and lands in Activity).
- Recommended row starts checked.

Edge cases:

| Case | Behaviour |
|---|---|
| Invalid accession shape | inline field error, no execution |
| Discovery fails / times out | inline error under the input, retry, keep the typed accession |
| `candidates: []` | empty state, not an error toast |
| `candidates` key missing | tool failure |
| RCSB 404 / upload failure for one candidate | that row marked failed, batch continues, summary toast says how many landed |
| `batch.create` fails wholesale | keep the modal open, surface the response `detail`; leave the uploaded files (harmless orphans) |
| Candidate already in project | unselectable, dimmed, tooltip |
| Modal closed mid-import | imports continue, toast reports |

This flow is **free** — no quoting, no `approveAmount`, no cost confirmation.

---

## 5. Analytics

`useCaptureEvent` from `@platform-ui/global-provider`, `EventName` from
`packages/global-provider/src/providers/analytics-provider.tsx`:

- `ToolExecutionStarted` / `ToolExecutionResult` for the discovery run
  (`toolKey`, `status`, `projectId`).
- `DataUploaded` on a successful import — `uploadType: 'uniprot'`,
  `recordCount`, `status: 'succeeded'`, `projectId`. Match the shape
  `handleRcsbImport` already uses.

---

## 6. Tests

Jest tests live in `apps/uui/test/src/...`, mirroring `src`.

- **Unit**: accession regex (valid 6- and 10-char, isoform rejected); grade sort
  with recommended-first; organism/method display maps; `field_status` cell
  rendering (`—` vs `?`); `jobOutputs` dict-vs-array normalization; terminal
  status handling in `run-tool-and-wait` (including legacy `Succeeded`);
  per-candidate failure isolation in `importCandidates`.
- **Component**: modal state machine (idle → discovering → results →
  importing); empty candidates; partial import failure; preview opens over the
  modal and closing it preserves the selection; already-in-project rows
  unselectable.
- **Table**: extend
  `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx` —
  it already asserts the + Data menu composition and that protein-only plugins
  are absent on the ligands table. Add the UniProt entry to those assertions.
- **Storybook**: a story for `candidate-table` with a fixture of ~6 candidates
  covering every grade, a `not_applicable` Rfree (cryo-EM), an `unknown` field,
  and one already-in-project row.
- **E2E**: a page-object entry beside `addDataButton` in
  `apps/platform-e2e/src/pages/`, mocking the discovery execution and the RCSB
  fetch.

Commands:

```bash
pnpm nx test uui
pnpm nx lint uui --fix
pnpm nx serve uui          # localhost:4208
pnpm nx storybook uui
```

Conventional commits (`@commitlint/config-conventional`).

---

## 7. Acceptance criteria

1. **Import from UniProt** appears in the proteins table's + Data menu and as an
   empty-state card; neither appears on the ligands table.
2. Typing a valid accession and pressing Enter runs exactly one
   `deeporigin.uniprot-discovery` execution and renders its candidates.
3. Rows are sorted best→worst by grade; the recommended row is starred, tinted
   and pre-checked.
4. Columns render exactly as §1.3, including `—` / `?` from `field_status`.
5. The eye opens a Mol\* window over the modal with that structure; closing it
   preserves the selection.
6. "Import Selected" imports every checked structure, creates one `proteins` row
   each with `file_path`, `pdb_id`, `uniprot_accession`, `protein_name`,
   `project_id` and `tags`, and the new rows appear in the table without a
   manual refresh.
7. A candidate already in the project cannot be selected.
8. A failing candidate does not prevent the others from importing.
9. An accession with no structures shows an empty state, not an error.
10. `pnpm nx test uui` and `pnpm nx lint uui` pass.

## 8. Constraints

- Do not put this feature in `packages/data-table`.
- Do not recompute or re-threshold the grade.
- Do not build the preview on the app-engine `protein-viewer-wrapper`.
- Do not rebuild Mol\*'s settings panel.
- Do not change `use-poll-tool-execution.ts`'s behaviour for the datasets
  import.
- Do not add a cost/quote step.

## 9. Decide these as you go (don't block)

- `protein_length` — the platform SDK sets it from the parsed structure; there
  is no PDB parser in the UI. Leave it unset unless counting distinct `CA`
  records in the fetched text is trivial in your implementation.
- Whether to cap the number of selectable rows. Default: no cap, bounded
  concurrency instead.
- Exact tint/star treatment for the recommended row — follow the app's Mantine
  theme rather than inventing colours.

Call out anything in the spec that conflicts with what you find in the code, and
say what you did instead.
