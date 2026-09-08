# UniProt discovery in DO Studio — integration research

**Goal.** Expose UniProt discovery in DO Studio (`platform-ui/apps/uui`) as an
entry in the proteins table's **Data** (Add Data) menu, matching what
`deeporigin.drug_discovery.UniprotDiscovery` already does in this SDK.

Status: research + implementation plan. No code in either repo yet.

---

## 1. What the CLI does today

Source: [`src/drug_discovery/uniprot_discovery.py`](../../src/drug_discovery/uniprot_discovery.py),
docs: [`docs/dd/tools/uniprot-discovery.md`](../../docs/dd/tools/uniprot-discovery.md).

### Tool

| | |
|---|---|
| Tool key | `deeporigin.uniprot-discovery` |
| Tool version | `latest` (registry: `src/platform/constants.py:161`, `TOOL_KEYS_AND_VERSIONS["uniprot_discovery"]`) |
| Endpoint | `POST /tools/{orgKey}/tools/{toolKey}/{toolVersion}/executions` (`src/platform/executions.py:56`) |
| Inputs | `{"uniprot_accession": "P00533"}` — that is the whole input schema |
| Body | `{inputs, outputs: {}, metadata: {}, sync: true, name?, approveAmount?}` (`_default_execution_payload`, `src/drug_discovery/execution.py:54`) |
| Outputs | `jobOutputs.candidates: []` — ranked PDB rows; empty list is valid |

Note the tool is a **tool**, not a tools-service *function*
(`/tools/{orgKey}/functions/{key}`), so DO Studio must go through
`fetchToolsExecuteTool`, not `useToolsRunFunction`.

### Candidate row (`UniprotDiscoveryCandidate`)

Required on every row: `pdb_id`, `grade` (`A`–`D`), `recommended` (exactly one
`true` unless the list is empty), `weighted_score`, and the component scores
`coverage_score`, `inhibitor_score`, `method_score`, `organism_score`,
`resolution_score`, `rfree_score`, plus `field_status` (per-field
`value` / `not_applicable` / `unknown`).

Optional descriptive fields: `coverage`, `has_ligand`, `method`,
`method_class`, `organism`, `organism_class`, `resolution`, `rfree`.

A representative payload is in the mock server
(`tests/mock_server/routers/tools.py:167`); accession `P99999` yields an empty
candidate list.

### Two phases — discover, then import

`run()` only ranks. **The tool never writes protein rows.** Import is entirely
client-side (`import_proteins`, `uniprot_discovery.py`):

1. Resolve selection — `pdb_ids` argument, or the single `recommended` row.
   Anything not in the candidate list raises.
2. `Protein.from_pdb_id(pdb_id)` → downloads the PDB from RCSB.
3. Set `uniprot_accession` (from the job) and `project_id`.
4. `protein.sync()` → uploads the file to UFA, then looks for an existing row
   by `file_path`; otherwise `register()` creates the entity row with
   `file_path`, `pdb_id`, `uniprot_accession`, `protein_name`,
   `protein_length`, `project_id`
   (`src/drug_discovery/structures/protein.py:1741`).

`Protein.from_uniprot(...)` is sugar for `import_proteins()[0]`.

Validation the client enforces before any network call: UniProtKB accession is
6 or 10 chars (`_UNIPROT_ACCESSION_RE`, isoform suffixes rejected), PDB IDs are
4 alphanumeric chars, uppercased.

---

## 2. What DO Studio has today

DO Studio is `platform-ui/apps/uui` (README: "DO Studio (formerly known as
UUI)"). The proteins table is `ChemicalTable` with `entity: "proteins"` —
rendered by the app-engine `Table` tile
(`app-engine/renderer/table-wrapper/table-wrapper.tsx:1033` passes
`allowAddData`) and by the dashboard's `proteins-table`
(`app-engine/available-components.ts:4`).

### The Add Data menu

`apps/uui/src/components/data-platform-tables/chemical-table.tsx`:

- `isProteinTable` — `:350`
- `extraToolbarItems` builds the `ActionMenu title="Data"` — `:1061`,
  menu at `:1135` (`testId="add-data-menu"`)
- Protein-only plugin today: `importFromRcsbAction` — `:470`
- Same actions also surface as empty-state cards — `emptyTableActions`, `:1034`

Menu items come from **action hooks** in `@platform-ui/data-table`
(`packages/data-table/src/hooks/`), each returning
`{ context, gridProps, menuItems, components, emptyTableActions }`
(contract: `packages/data-table/CLAUDE.md` §"Action Hook Pattern").

### The closest existing analogue — Import from RCSB

`packages/data-table/src/hooks/use-action-import-from-rcsb.tsx`: a modal with
one text input, `fetch('https://files.rcsb.org/download/{ID}.pdb')` in the
browser, then `onImport(file)`.

The consumer side is `handleRcsbImport` (`chemical-table.tsx:451`) →
`handleProteinFiles` (`:374`), which for each file:

1. `uploadFileToFileService({ orgKey, file, subdir: 'structure-imports' })`
2. `serverInterface.entity.create(orgKey, 'proteins', { set: { file_path, protein_name, project_id?, tags: provenanceTags }, returning })`
3. `applyServerSideTransaction({ add })` + `invalidateRowCount()` + a
   notification, and a `DataUploaded` PostHog event.

**This is the exact browser equivalent of `Protein.sync()`** — same upload,
same entity write. UniProt import needs only two extra `set` fields:
`uniprot_accession` and `pdb_id` (both exist on the entity —
`PROTEIN_RETURNING_FIELDS`, `src/platform/entities.py:59`).

### How DO Studio runs a tool outside the app engine

`apps/uui/src/components/datasets/import-modal/use-import-dataset.ts` is the
precedent for "run a tool from a modal, not from an app manifest":

- `useToolsListClusters<any>` → prefer a `us-west-2` cluster, else the first
- `fetchToolsExecuteTool({ pathParams: { orgKey, toolKey }, body })` with
  `ExecuteToolSchemaDto = { inputs, outputs, clusterId, name?, projectId? }`
- `usePollToolExecution(orgKey)` (`use-poll-tool-execution.ts`) polls
  `useToolsGetToolExecution` every 3 s to a terminal status and hands back
  `jobOutputs`.

`ExecuteToolSchemaDto` (generated from `packages/api-gateway/openapi.yaml:11244`)
has **no `sync` field**, so DO Studio should poll rather than rely on the
SDK's `sync: true`. Sending `sync: true` anyway is harmless if the gateway
forwards it, but the polling state machine has to exist regardless.

---

## 3. Proposed integration

### 3.1 New action hook — `useActionImportFromUniprot`

`packages/data-table/src/hooks/use-action-import-from-uniprot.tsx`, exported
from `hooks/index.ts` next to `useActionImportFromRcsb`.

Keep the hook **transport-free**, like the RCSB one: it owns the modal, the
accession input, candidate selection and validation; the caller supplies the
two side effects.

```ts
useActionImportFromUniprot({
  agentIdPrefix,                              // `table-proteins`
  onDiscover: (accession: string) => Promise<UniprotCandidate[]>,
  onImport: (rows: { pdbId: string; accession: string; file: File }[]) => Promise<void>,
})
```

Modal, two steps in one surface:

1. **Accession input** — mirror `_UNIPROT_ACCESSION_RE` client-side so a typo
   never costs an execution; `Discover` runs the tool with a loading state.
2. **Candidate table** — rows sorted by `weighted_score`, `recommended` row
   pre-checked and badged; columns `pdb_id`, `grade`, `weighted_score`,
   `resolution`, `method`, `organism`, `has_ligand`. Multi-select checkbox,
   `Import` imports every checked row. Empty `candidates` → "No experimental
   PDB structures for {accession}".

`menuItems`: `{ text: 'Import from UniProt', action: openModal }`.
`emptyTableActions`: one card, `icon: 'download'`, subtitle "Rank PDB
structures for a UniProtKB accession".

### 3.2 Wiring in `chemical-table.tsx`

- Gate on `isProteinTable`, exactly like `importFromRcsbAction`: add to the
  `ActionMenu` `plugins` array (`:1135`), to `extraToolbarItems`' portaled
  modal list, and to `emptyTableActions` (`:1034`).
- `onDiscover` — new hook `use-uniprot-discovery.ts` under
  `components/data-platform-tables/hooks/`: `useToolsListClusters` +
  `fetchToolsExecuteTool({ toolKey: 'deeporigin.uniprot-discovery' })` +
  the existing polling state machine, returning `jobOutputs.candidates`.
  `use-poll-tool-execution.ts` is generic (`orgKey` only) but currently lives
  under `components/datasets/import-modal/` — move it to a shared
  `components/hooks/` and re-export, rather than copying it.
- `onImport` — a `handleUniprotImport` next to `handleRcsbImport` (`:451`):
  browser-fetch each selected PDB from `files.rcsb.org` (reuse the RCSB hook's
  URL builder — worth lifting to a shared util so both call sites agree), then
  a `handleProteinFiles` variant that accepts extra `set` fields so the row
  carries `uniprot_accession` and `pdb_id`. Add both to the `returning` array.
  Keep `tags: provenanceTags` (PUI-2146) and the `DataUploaded` capture with
  `uploadType: 'uniprot'`.

### 3.3 Suggested column addition

`uniprot_accession` is not in the dashboard proteins-table default columns
(`available-components.ts:12`). Add it there (and to the app manifests'
protein tables that show `gene_symbol`) so an imported row visibly carries its
accession.

---

## 4. Relation to the app engine

The App Engine is manifest-driven: each app is one JSON in
`apps/uui/src/app-schemas/`, mapped **1:1 to a backend tool**, with a sidebar
form (`inputSchema` + `x-user-input` / `x-data-type` annotations), tiles, and
`engine.submit()` → `buildToolPayload` → `executeToolWithVersion`
(`apps/uui/CLAUDE.md` §"Backend Tool Integration", §"Submit Flow").

**UniProt discovery does not fit that shape, and should not be an app.** The
engine's contract is *select entities in the project → run a tool on them*.
UniProt discovery runs on a free-text accession with **no** entity selection,
and its purpose is to *create* the proteins an app would later consume. The
Add Data menu is the right surface, exactly as Import from RCSB is.

Three things still connect it to the engine:

1. **Same execution plumbing.** The modal uses the same tools-service endpoint,
   the same `clusterId` resolution, and the same
   `useToolExecutionContext()` provenance (`app` / `session`) the engine spreads
   into every execution body. Pass `projectId` so the run is scoped like an
   app run.
2. **Activity page.** Every execution shows up on the Activity page's
   `JobManagerTable`, which resolves tool identity through
   `manifestByToolKey` (`app-schemas/index.ts:33`,
   `job-manager-table/index.tsx:66`). Without an entry, a UniProt run renders
   with a raw tool key and no `identityHue`. Precedent exists for a
   **metadata-only manifest**: `system-prep.json` is in `manifestByToolKey`
   but not in `appManifests`, so it has no route. Ship a small
   `uniprot-discovery.json` the same way (`id`, `toolKey`, `name`,
   `identityHue`) — identity without a page.
3. **`ligand-search` is the nearest app-engine cousin** if a full page is ever
   wanted: it queries an external source and renders hits in a `CsvTable` whose
   `addToProject` config writes selected rows into the `ligands` entity
   (`app-schemas/ligand-search.json`,
   `app-engine/renderer/csv-table/use-add-to-project.ts`). The equivalent for
   proteins would need `addToProject` to handle a structure download +
   file-service upload, which it does not do today — another reason to start
   with the menu item.

---

## 5. CLI ↔ DO Studio parity

| CLI | DO Studio |
|---|---|
| `UniprotDiscovery(uniprot_accession=...)` | accession input in the modal |
| `run()` | `fetchToolsExecuteTool` + poll → `jobOutputs.candidates` |
| `candidates`, `recommended` | candidate table, recommended row pre-checked |
| `import_proteins()` (no args) | Import with only the recommended row checked |
| `import_proteins([...])` | Import with a multi-row selection |
| `Protein.from_pdb_id` | browser `fetch(files.rcsb.org/download/{ID}.pdb)` |
| `protein.sync()` | `uploadFileToFileService` + `entity.create` |
| `uniprot_accession` / `pdb_id` on the row | same two fields in the `set` payload |
| `run(quote=True)` / `approve_amount` | not exposed (see open questions) |

---

## 6. Open questions / risks

- **Cost & quoting.** The CLI supports `quote=True` (`approveAmount: 0`). The
  RCSB menu item costs nothing; a UniProt item costs an execution. Decide
  whether the modal shows a quote or just runs.
- **Duplicate proteins.** `Protein.sync()` dedupes on `file_path`;
  `handleProteinFiles` does not — importing the same PDB twice will create two
  rows. Match the CLI by searching `proteins` for the resolved `file_path`
  first, or accept the duplicate.
- **`sync` on the execution body.** Not in the published DTO; polling is the
  contract to build against.
- **RCSB fetch from the browser.** Already proven by the RCSB action, so CORS
  is not a new risk — but a per-candidate multi-import means N sequential
  fetch+upload round trips; cap the selection or run them with limited
  concurrency and report partial failures per PDB (the CLI fails the whole
  call on the first error).
- **Tool version.** The CLI pins `latest`. `fetchToolsExecuteTool` targets the
  latest enabled version by tool key; if the UI needs a pin, use the
  versioned execute variant.
- **Empty candidate list is a success**, not an error — the modal must say so
  rather than showing a failure.

---

## 7. Implementation checklist (platform-ui)

1. `packages/data-table/src/hooks/use-action-import-from-uniprot.tsx` + export
   in `hooks/index.ts`; unit test alongside the RCSB hook's.
2. Lift `rcsbPdbUrl` into a shared util used by both protein import hooks.
3. Move `use-poll-tool-execution.ts` out of `components/datasets/import-modal/`
   into a shared hooks folder; keep the datasets import working.
4. `components/data-platform-tables/hooks/use-uniprot-discovery.ts` — cluster
   resolution, execute, poll, parse `candidates`.
5. `chemical-table.tsx` — `handleUniprotImport`, extra `set` fields on the
   protein create, plugin registration in `emptyTableActions`,
   `extraToolbarItems`, and the `add-data-menu` `plugins` array.
6. `available-components.ts` — add the `uniprot_accession` column to the
   proteins table.
7. `app-schemas/uniprot-discovery.json` — metadata-only manifest registered in
   `manifestByToolKey` only (mirror `system-prep.json`).
8. Tests: `apps/uui/test/src/components/data-platform-tables/chemical-table.test.tsx`
   already asserts the Add Data menu composition — extend it; add an e2e page
   object entry next to `addDataButton` in `apps/platform-e2e/src/pages/`.
9. Analytics: `uploadType: 'uniprot'` on `EventName.DataUploaded`.
