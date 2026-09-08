# Target Preparation in DO Studio — implementation plan

**PRD:** [Target Preparation Tool Workflow and Application](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1063289062/Target+Preparation+Tool+Workflow+and+Application) (epic DDOS-7091)
**Child PRD:** [Structure Report Tool](https://deeporigin.atlassian.net/wiki/spaces/PR/pages/1066631169/Structure+Report+Tool)

**Primary sources**

- **[`target-preparation-backend-report.md`](./target-preparation-backend-report.md)** — verified
  contract for `deeporigin.target-preparation` v1.0.0, read out of
  `deeporiginbio/platform-toolbox`. Every backend claim in this plan is cited to it.
- `do-dd-client` — the tools `target-preparation` wraps:
  `src/drug_discovery/protein_prep.py`, `structure_report.py`, `pocket_finder.py`,
  `src/platform/constants.py`, `src/platform/executions.py`
- `platform-ui` — the App Engine Studio apps are built on:
  `apps/uui/src/app-engine/**`, `apps/uui/src/app-schemas/**`, `packages/molstar/**`

Scope is **UI**: what DO Studio builds, which API calls it makes at each step, and which
data-platform records get read or written.

> **Revision note.** This document was first written against the PRD plus the shipped
> Python client, with the `target-preparation` schema unread. The backend report has since
> landed and **contradicted several assumptions**. §2 lists every conflict. The three that
> change the most work: the four "Advanced Parameters" checkboxes in the mockup do not
> exist (§2.1); `keep` on a ligand means *extract to a separate file*, not *retain in the
> receptor* (§2.3); and the keep/skip map must be **exhaustive over every component,
> including every water** (§2.4). Two earlier claims of my own were wrong and are corrected
> in §2.7.

---

## 1. The verified contract

### 1.1 `deeporigin.target-preparation` v1.0.0

| | |
| --- | --- |
| Tool key | `deeporigin.target-preparation` (dotted form; the Argo name `deeporigin-target-preparation` is not the key) |
| Version to pin | major **`1`** |
| Billing code | `DO_POCKET_FINDER` |
| Execution mode | **Argo workflow, always** — never direct, never synchronous |
| `jobOutputs` | **none** — the executions API omits it for workflow mode by design |
| Root `required` | `["action", "protein"]` |
| `additionalProperties` | **`false`** at the root and on `pocket`, `protein`, `selection`, `pocket.selections[]`, `pocket.selections[].author` |

Two mutually exclusive branches on `action`:

**`action: "recommend"`** — runs source Structure Report + Protein Prep recommend.
`selection`, `model_missing_loops` and `pocket` are **forbidden**.

```jsonc
{ "action": "recommend",
  "protein": { "id": "<protein entity id>", "file_path": "<UFA path>" },
  "pdb_id": "6GOG" }                       // optional, for RCSB metadata
```

**`action: "prepare"`** — applies the Selection, prepares, reports, finds pockets.
`selection`, `model_missing_loops` and `pocket` are all **required**.

```jsonc
{ "action": "prepare",
  "protein": { "id": "…", "file_path": "…" },
  "selection": {
    "source_sha256":    "<64 hex>",        // must match the analyzed bytes exactly
    "analyzer_version": "1.1.0",           // must match the analyzer that ran
    "decisions": { "chain:A": "keep", "water:A:HOH:310:": "skip", … }
  },
  "model_missing_loops": true,             // no default — must be sent
  "pdb_id": "6GOG",                        // required unless model_missing_loops === false
  "pocket": { "mode": "auto-find", "pocket_count": 5, "pocket_min_size": 30 } }
```

`protein` requires **both** `id` and `file_path`. `id` is *"required for pocket indexing"* —
an unregistered structure cannot be run.

**Outputs** are read from indexed result tables, addressed by `x-result-group`:

| Output key | `x-data-type` | `x-result-group` | Notes |
| --- | --- | --- | --- |
| `protein` | `PreparedProtein` | `preparedproteins` | `.protein_pdb_file_path` is the prepared PDB, stamped `REMARK  99 DO_PREPARED` |
| `pockets[]` | `Pocket` | `pockets` | carries the PCA-aligned `box` for docking |
| `structure_reports[]` | `StructureReport` | `structurereports` | `report_role: "source"` on recommend, `"prepared"` on prepare |
| `extracted_ligands[]` | `ExtractedLigand` | `extractedligands` | kept ligands, pulled out as separate PDBs |
| `recommendation` | **none** | **none** | the component inventory — see §3.1 |
| `selection_file_path` | — | — | UFA path to the frozen Selection JSON |
| `audit_file_path` | — | — | UFA path to the transformation audit |

Workflow order on `prepare`: prepare (SystemPrep) → prepared Structure Report →
Pocket Finder. The last two have no `when:` guard — **both always run**.

### 1.2 `deeporigin.protein-prep` — still needed for the inventory

Same analyzer, same component ids, same Selection shape. The difference that matters:
its `recommend` action quotes **`"direct"`**, so it executes synchronously and returns
`jobOutputs.recommendation` in the create response. See §3.1.

### 1.3 `deeporigin.structure-report` and `deeporigin.pocket-finder`

Both run inside `target-preparation`. The UI calls neither directly on submit; it may call
structure-report on selection (§3.2).

---

## 2. Where the mockup and the tool disagree

Every row is a product decision that has to be made before the sidebar is built.

### 2.1 Three of the four Advanced Parameters checkboxes do not exist

| Mockup checkbox | Reality |
| --- | --- |
| Add missing atoms & residues | **No input.** Always on. |
| Add missing loops | `model_missing_loops` — exists, **required on prepare, no default** |
| Protonate (pH 7.4) | **No input.** Protonation always runs. `protonate_protein` exists only as an *output* echo. |
| Find Pockets | **No input.** Pocket Finder always runs on prepare. Only the *mode* is choosable. |

Root is `additionalProperties: false`, so inventing the three missing keys is a hard
schema rejection at execution-create, not a warning.

**Recommendation:** render only `model_missing_loops` as a control. Show the other three as
static "always applied" text if the reassurance is wanted, but not as checkboxes — a toggle
that silently does nothing is worse than no toggle. Take this back to the PRD owner.

### 2.2 There is no run-name or output-name input

No `name`, `output_name`, `run_name` or label input exists.

**But the mockup's "Run Name" field is still fine**, and this is worth being precise about
because the backend report's advice ("do not send `name`") is about `inputs`. The engine's
`x-body-key` annotation puts a field at the **execution body** level —
`useSubmitExecution` builds `{ inputs, outputs, clusterId, ...bodyFields }`, so `name`
lands as a sibling of `inputs`, never inside it, and `additionalProperties: false` does not
see it. **`run_name` with `x-body-key: "name"` stays.**

What dies is the mockup's **"Output Property Name"** field — there is nothing to bind it
to. Drop it.

### 2.3 `keep` on a ligand means *extract*, not *retain*

A kept ligand is routed to `extract_ligand_keys`, written out as its own PDB, and surfaces
in `extracted_ligands[]`. It is **not** in the prepared receptor. The analyzer's own reason
string says so: *"Recognized ligand {resname}; keep means extract"*.

The PRD's two-tone palette (saturated = keep, desaturated = exclude) is therefore **wrong
for ligands**. A user who sees a bright orange ligand labelled "Keep" will reasonably expect
it in the output structure, and it will not be there.

**Recommendation: three states, not two** — Keep (in the receptor) / Extract (separate
file) / Skip (discarded) — with "Extract" being what `keep` means on a `ligand` component.
Chains, cofactors and waters keep the two-state meaning. This needs a palette decision and
a PRD conversation; it is the single most user-visible correction in this document.

### 2.4 The decisions map must be exhaustive — every water included

A partial map is rejected (`selection.decisions is missing components: …`), an unknown id
is rejected, and a lingering `review` is rejected. Every water in the file is inventoried as
its own component, and non-coordinating waters default to `review`.

For a typical crystal structure that is **hundreds of rows the user must resolve** before
the Run button can enable. A flat list of per-row toggles, as the mockup shows, does not
survive contact with a real PDB.

**Recommendation:** group the panel by `kind`, collapse waters into a single summarised
group with bulk "keep all / skip all" plus a count, and virtualize the list. The mockup's
one-row-per-component layout works for chains, ligands and cofactors — which is what the
screenshot happens to show — and needs a different treatment for waters.

### 2.5 Some Keep choices are illegal

Both are 400s from the prepare child, so both must be blocked client-side:

- Keeping a chain whose `subtype` is `nucleic`, `mixed` or `other` →
  *"Cannot keep unsupported polymer chain"*. The inventory flags these already:
  `recommendation: "review"`, `reason_code: "non_protein_chain"`.
- Keeping **zero** protein chains → *"selection must keep at least one supported protein
  chain"*.

### 2.6 `prepare` requires a whole pocket configuration the mockup has no room for

`pocket` is required, with three modes and mutually exclusive per-mode requirements:

| `mode` | Required | Forbidden |
| --- | --- | --- |
| `auto-find` (default) | `pocket_count`, `pocket_min_size` | `crystal_ligand` |
| `define-by-selection` | `selections[]` | `crystal_ligand` |
| `from-crystal-ligand` | `crystal_ligand` | `selections`, `pocket_count`, `pocket_min_size` |

`pocket_count` and `pocket_min_size` have **no schema defaults**, so the UI must supply
numbers. Cheapest viable default: `{ "mode": "auto-find", "pocket_count": 5,
"pocket_min_size": 30 }` — matching the shipped `pocket-finder.json` manifest defaults.

Worth noting as a v2 opportunity: `from-crystal-ligand` takes a `component_id` of a ligand
this very run extracted. A user who keeps `Ligand XO4` (as the mockup screenshot does) could
have the pocket defined from it in one gesture. Out of scope for v1, but the schema is
already there.

### 2.7 Two corrections to earlier drafts of this plan

- **Component id grammar.** An earlier draft claimed ligand and water ids use different
  field orders. That was wrong. There is **one** residue id builder and the grammar is
  uniform: `kind:chain:resname:resseq:icode`. The shape `ligand:LIG:A:100` in
  `do-dd-client`'s fixtures cannot be produced by the analyzer — it is a bad fixture (as are
  `ligand:A:120:IBP` in several `platform-toolbox` tests). The practical advice was right
  for the wrong reason: **do not parse ids** — echo them back verbatim and use `author` for
  addressing. See §9.3.
- **The preflight route.** An earlier draft floated calling it for a cheaper inventory. It
  is a create-time billing/validation gate returning `{method, counts, validations?}`,
  cluster-internal, not browser-reachable. That option is void. Its seven validations are
  still worth mirroring client-side (§7).

---

## 3. Architecture: two executions, one button press

`target-preparation` cannot be driven by a single execution — `prepare` needs a
`source_sha256`, an `analyzer_version` and an exact component id set that only a
`recommend` run produces.

That is **not** a problem for the App Engine, and does not mean a two-step wizard. The
engine already has the shape: `manifest.steps` runs pre-submit tools, and an `ActionField`
with `autoRun: true` fires one the moment its `requires` selection is satisfied. So the
inventory is fetched **when the user clicks a protein row**, not when they press Run. By the
time Run is pressed the Selection exists, and the user presses one button.

```
row select ──► [step] inventory ─────► filtering panel + viewer
           └─► [step] structure report ─► report card
                                              │
                          user toggles ───────┤
                                              ▼
                                     [submit] target-preparation action=prepare
```

### 3.1 Which tool serves the inventory — the one open runtime question

`target-preparation` `action: "recommend"` is the obvious candidate and is probably
**unusable from the UI**:

- it runs in **workflow mode**, so the create response carries no `jobOutputs`;
- its `recommendation` output has **no `x-result-group`**, so it has no indexed result
  table either.

Between those two, there may be nowhere to read it from. `deeporigin.protein-prep`
`action: "recommend"` quotes **`"direct"`**, executes synchronously, and returns
`jobOutputs.recommendation` in the create response — which is exactly what
`useRunStep`'s `bind` reads.

**Plan of record: get the inventory from `deeporigin.protein-prep` `action: "recommend"`.**

This is safe on identity grounds — both tools drive the same `toolbox_core` analyzer, the
same `ANALYZER_VERSION`, the same component ids and the same Selection schema, and
`target-preparation`'s prepare child *is* the sysprep service that `protein-prep` calls. But
it is a cross-tool assumption, so **verify it once**: run `protein-prep recommend`, feed its
Selection to `target-preparation prepare`, confirm the digest and analyzer checks pass.

If instead the platform *does* surface a workflow-mode `recommendation`, switch the step's
`toolKey` to `target-preparation` and add `"action": "recommend"` to the step's `inputs`.
Nothing else changes — the tile's state contract is identical. Build the tile against a
fixture and this stays a one-line manifest change.

### 3.2 The source structure report

`target-preparation recommend` emits a `report_role: "source"` row, and `prepare` emits a
`"prepared"` row. Both are indexed (`x-result-group: structurereports`), so both are
readable even in workflow mode.

But if the inventory comes from `protein-prep` (§3.1), no `target-preparation recommend`
execution exists, and with it no source report. So the app calls
**`deeporigin.structure-report` directly** as a second pre-submit step — synchronous, with
`jobOutputs.structure_reports`, exactly as the original plan had it.

### 3.3 What the app calls, end to end

| Call | Tool | Action | When | Mode |
| --- | --- | --- | --- | --- |
| Grade the structure | `deeporigin.structure-report` | — | on row select | sync, `jobOutputs` |
| Inventory components | `deeporigin.protein-prep` | `recommend` | on row select | direct/sync, `jobOutputs`, `__SKIP`-billed |
| The run | `deeporigin.target-preparation` | `prepare` | on submit | workflow/async, indexed results, 1 × `DO_POCKET_FINDER` |

`deeporigin.pocket-finder` and `protein-prep`'s `prepare` action are **never** called
directly — both are inside `target-preparation`.

---

## 4. End-to-end flow, with the exact calls

### Step 0 — route and shell

`/target-prep` renders `AppPage` → `AppEngineProvider` → mosaic layout + form sidebar.
Manifest from `/app-manifests/deeporigin.target-preparation/1.json` in staging/prod, from
the bundled `src/app-schemas/target-prep.json` locally and on PR previews.

### Step 1 — proteins table

`TableWrapper`, `config.entity: "proteins"`, `singleSelection: true`. Existing server-side
row model, no new data-platform work.

**One gate the other apps do not have:** `target-preparation` requires `protein.id` **and**
`protein.file_path`, both non-blank after trimming (whitespace-only passes schema but fails
preflight). Rows missing either must be unselectable or the Run button disabled with a
reason.

### Step 2 — row selected → two synchronous steps auto-run

`engine.setSelectedProteins([row])` → Zustand. Two hidden `ActionField`s with
`autoRun: true` fire through `useRunStep`.

**2a. Structure Report**

```
POST /tools/{orgKey}/tools/deeporigin.structure-report/{major}/executions
body: { inputs: { protein: {id, file_path}, pdb_id? }, outputs: {}, clusterId,
        visibility: "hidden" }
bind: { "structure_reports": "structureReport" }
```

**2b. Component inventory**

```
POST /tools/{orgKey}/tools/deeporigin.protein-prep/{major}/executions
body: { inputs: { action: "recommend", protein: {id, file_path} }, outputs: {},
        clusterId, sync: true, visibility: "hidden" }
bind: { "recommendation": "recommendation" }
```

`action: "recommend"` comes from `ActionField.inputs` (merged over the schema-built inputs
as `inputOverrides`), not from a schema annotation.

`useSelectionResets` already clears `stepOutputs` when the selection changes, and
`ActionField` re-arms per selection — changing protein re-runs both.

> **Engine gap #1** — `useRunStep` never sends a top-level `sync`. Add
> `AppStep.sync?: boolean` and spread it into the body
> (`apps/uui/src/app-engine/hooks/use-run-step.tsx:135`). Precedent:
> `renderer/table-wrapper/use-export-dataset.ts:109`.

Both steps must handle the **HTTP 200 + `status: "Failed"`** pattern `useRunStep` already
guards, and the filtering tile must render the step's `stepErrors` entry rather than an
empty list — a failed inventory and an empty inventory look identical otherwise.

### Step 3 — Structure Report card

Reads `stepOutputs.structureReport[0]`. Renders grade badge, the
`X-ray` / `Human` / `2.31Å` / `Ligand` pills, and the coverage line. `field_status` (six
fixed keys) says which pills are `unknown` or `not_applicable` rather than absent — render
those states, do not blank them. Never recompute the grade.

Residue count (`1,036 Residues` in the mockup) is **not** in the report; read
`proteins.protein_length` off the entity row.

### Step 4 — Structure Filtering tile

Reads `stepOutputs.recommendation`, seeds a decision map from each component's
`recommendation`, and writes the full Selection back on every edit:

```ts
engine.setStateValue('selection', {
  source_sha256:    recommendation.source_sha256,
  analyzer_version: recommendation.analyzer_version,
  decisions:        { 'chain:A': 'keep', 'water:A:HOH:310:': 'skip', … },  // exhaustive
});
```

Requirements this tile carries, beyond the mockup:

- **Exhaustive** — one entry per inventoried component, no exceptions (§2.4)
- **Grouped, with bulk actions and virtualization** for waters (§2.4)
- **Three states on ligands** — Keep / Extract / Skip (§2.3)
- **Illegal keeps blocked** — non-protein chains, and the last protein chain (§2.5)
- **`review` blocks submit** until resolved

> **Engine gap #2** — `stepOutputs` is written only by `useRunStep`. Expose
> `engine.setStateValue(key, value)` over the store's existing `setStepOutputs` so a tile
> can contribute state that `x-from-state` reads. Stays tool-agnostic.

> **Engine gap #3** — `requiredEntities` only gates on entity selections. Add a
> `SubmitValidation` kind for state-shaped rules, evaluated in `engine.submit`:
> no `review` remaining, decisions exhaustive against the inventory, ≥1 protein chain kept,
> no unsupported chain kept.

### Step 5 — Mol\* viewer

Per-component colouring. Build components once per structure load, swap only the colour
theme on toggle.

```ts
MS.struct.generator.atomGroups({
  'chain-test':   MS.core.rel.eq([MS.ammp('auth_asym_id'), author.chain_id]),
  'residue-test': MS.core.rel.eq([MS.ammp('auth_seq_id'), author.resseq]),
})
```

`tryCreateComponentFromExpression` + `buildUniformColor` already exist
(`packages/molstar/src/api/loaders.ts:448`, `utils/color-themes.ts`).

**Address from `author`, never from the id.** And apply `recommendation.chain_id_mapping`
first — for mmCIF input, author chain ids are remapped to single-character PDB chain ids,
so the inventory's `chain_id` may not match the file the viewer loaded.

PRD palette, plus the third state §2.3 requires:

| Component | Keep | Exclude |
| --- | --- | --- |
| Protein chain | `#2563eb` | `#bfdbfe` |
| Ligand — **needs a third "extract" treatment** | `#f97316` | `#fed7aa` |
| Co-factor | `#9333ea` | `#e9d5ff` |
| Water | `#06b6d4` | `#cffafe` |

### Step 6 — sidebar

| Section | Contents |
| --- | --- |
| Structure Components | read-only chips of kept components (engine gap #4) |
| Preparation | `model_missing_loops` boolean; `pdb_id` string, **required and validated `^[A-Za-z0-9]{4}$` whenever loops are on** |
| Pocket Detection | `pocket.mode` enum (default `auto-find`), `pocket_count`, `pocket_min_size` |
| Run Details | `run_name` (`x-body-key: "name"`) |

**The `model_missing_loops` / `pdb_id` coupling is the sharpest UX edge in the app.** Loops
default on everywhere else, and loops-on makes `pdb_id` mandatory. A structure the user
uploaded that has no PDB entry therefore **cannot run with loop modelling at all** — they
must turn it off. Surface that as an explanatory state on the checkbox, not as a
post-submit error.

> **Engine gap #4** — no read-only "chips from state" field type. Add `component-summary`.

### Step 7 — Run

```
POST /tools/{orgKey}/tools/deeporigin.target-preparation/1/executions
body: { inputs: { action: "prepare", protein, selection, model_missing_loops, pdb_id?, pocket },
        outputs: {}, clusterId, projectId, name: "<run_name>" }
```

`buildToolPayload` emits exactly the declared keys. `action: "prepare"` is a constant —
supply it via an `x-from-form` with a hidden field, or a hidden `enum` field defaulted to
`prepare`.

Billing: one `DO_POCKET_FINDER` unit. The existing quote →
`price-confirmation-panel` → `insufficient-funds-modal` path works unchanged.

**Errors worth handling explicitly**, because retrying makes them worse:

| Error | Cause | Right response |
| --- | --- | --- |
| `selection.source_sha256 does not match…` | the structure file changed since the inventory | re-run the inventory step; do not retry submit |
| `selection.analyzer_version does not match…` | the sysprep serving was upgraded mid-session | re-run the inventory and **diff** — the user's decisions may not all map |
| `decisions is missing/has unknown components` | a UI bug — the map drifted from the inventory | never reachable if the tile is exhaustive |

### Step 8 — results mode

`/activity/{executionId}` → `appMode: 'results'`, using `resultsLayout` +
`resultsComponents`. **There are no `jobOutputs`** — every tile reads indexed result rows
via the result-explorer, filtered by `compute_job_id` (the pattern in
`components/pipelines/hooks/use-step-results.ts`).

- Structure Report tile — the `report_role: "prepared"` row
- Mol\* viewer — prepared structure + pockets, renderer `renderStructureAndPockets`
- Pockets table
- Extracted ligands list, if any

**A FAILED execution can still have produced everything except pockets.** The prepared
protein is written and stamped before Pocket Finder runs, and nothing rolls it back, so a
pocket-finder failure reads FAILED with a perfectly good prepared structure in
`results__preparedproteins`. The results view must render retained artifacts on a FAILED
run rather than showing an empty error state. Zero pockets, separately, is a **success**.

---

## 5. Data platform: reads and writes

| Record | Direction | Who | New? |
| --- | --- | --- | --- |
| `proteins` / `proteins_with_results` | read | table tile | no |
| Structure file bytes (UFA) | read | Mol\* viewer | no |
| `results__preparedproteins` | write | `target-preparation` | no |
| `results__pockets` | write | `target-preparation` | no |
| `results__structurereports` | write | `target-preparation`, `structure-report` | no |
| `results__extractedligands` | write | `target-preparation` | no |
| **`proteins` row for the prepared structure** | write | **nobody** | **missing — see below** |

### 5a. The prepared protein never becomes a selectable Protein — backend gap

Confirmed negative. The prepared structure is emitted with
`x-result-group: "preparedproteins"`, a dynamic `results__*` table. The repo's authoring
rules explicitly forbid tool rows from being typed `Protein`, and the only
`entities.create_protein` call site in the whole toolbox belongs to a different tool
(`deeporigin.pdb-import`). `target-preparation`'s children only call `send_result`.

So after a successful prepare, the structure exists as a UFA file and a
`results__preparedproteins` row — and **a protein picker backed by the `proteins` table
will not offer it to Docking or ABFE.**

If the point of Target Prep is a prepared receptor you can then dock against, that wiring
does not exist. **This is a backend change, not a UI workaround** — a client-side
`POST /data-platform/{orgKey}/proteins` would have to run after an async workflow the user
may have navigated away from. Raise it with the tool owner; it is the largest gap between
what ships and what the PRD is for.

### 5b. Structure reports are already indexed — better than assumed

An earlier draft asked for a new `structurereport` result type. It exists:
`x-result-group: "structurereports"`, with `protein_id` marked `x-key: true` and
`report_role` distinguishing `source` from `prepared`. Grade/resolution/coverage as columns
on the proteins table is therefore feasible with the manifest's existing
`results.detailColumns` / `aggregates` machinery — no backend work.

One unknown: the literal singular `result_type` string the platform derives from the plural
`x-result-group`. That mapping lives outside the toolbox repo. Confirm against one real
completed execution before writing the result-explorer filters.

---

## 6. New and changed UI components

### New tiles (`apps/uui/src/app-engine/renderer/`)

| Component | Purpose | Registration |
| --- | --- | --- |
| `StructureFiltering` | grouped, virtualized, bulk-actioned Keep/Extract/Skip list; writes `selection` | `renderer/registry.ts`, `types/tiles.ts` `LeafComponent`, `renderer/types.ts` config |
| `StructureReportCard` | grade badge, pills, score breakdown | same three files |

`StructureFiltering` is materially bigger than the mockup suggests — §2.3, §2.4 and §2.5 are
all its responsibility. Size it accordingly.

### New Mol\* renderer (`packages/molstar/`)

- `src/api/components.ts` — `renderStructureComponents(plugin, content, format, components)`
- add `'renderStructureComponents'` to `RendererName` (`src/types/index.ts:243`)
- add the dispatch case (`src/api/index.ts:526`)
- export the palette, including the third ligand state, as shared constants
- apply `chain_id_mapping` before addressing

### Engine changes

| Gap | File | Change |
| --- | --- | --- |
| #1 | `hooks/use-run-step.tsx` | `AppStep.sync?: boolean` → top-level body key |
| #2 | `engine-provider.tsx`, `engine-context.ts` | `setStateValue(key, value)` |
| #3 | `engine-provider.tsx`, `types/app.ts` | `SubmitValidation` kinds for state-shaped rules |
| #4 | `form-engine/form-field.tsx`, `types/form.ts` | `component-summary` field type |

None of these put tool-specific logic in the engine.

### Manifest and registration

- `apps/uui/src/app-schemas/target-prep.json` (§8)
- register in `app-schemas/index.ts` under key `target-prep`
- nav entry in `packages/global-provider/src/containers/subscription.container.tsx`
- app card in `apps/uui/src/pages/applications/constants.ts` + image asset
- `identityHue` ≥10° from neighbours (current cluster ~185–323)
- tool display name in `components/data-platform-tables/job-manager-table/index.tsx`
- add `deeporigin.target-preparation` to `ALL_RESULTS_TOOL_KEYS`
  (`hooks/use-manifest.tsx:14`)
- publish via `.github/workflows/register-tool-manifest.yml`, major `1`

---

## 7. Client-side validations to implement

The preflight route is not browser-reachable, so every one of its checks surfaces as an
execution-create failure unless the UI gets there first. All seven are trivially enforceable:

| Check | Rule |
| --- | --- |
| `protein` | object present |
| `protein.id` | non-blank **after trimming** — whitespace-only passes schema, fails preflight |
| `protein.file_path` | non-blank after trimming |
| `model_missing_loops` | present on prepare |
| `selection` | object present on prepare |
| `pocket` | object present on prepare |
| `pdb_id` | present when `model_missing_loops` is true |

Plus the schema and runtime rules preflight does *not* cover: `pdb_id` matches
`^[A-Za-z0-9]{4}$`; no unknown keys; pocket mode exclusivity; decisions exhaustive, two-
valued, no `review`; ≥1 protein chain kept; no `nucleic`/`mixed`/`other` chain kept.

---

## 8. Open questions

1. **Is the tool registered on dev / staging / prod?** Not declarable from the repo; the
   toolbox's own E2E test still skips pending publication. Ask platform to run
   `dump-enabled-tools.yml` per environment. **Blocks any integration testing.**
2. **Does a `protein-prep`-produced Selection validate against `target-preparation`
   prepare?** (§3.1) Same analyzer, same ids — but verify once end to end. **Blocks the
   architecture.**
3. **What is the literal `result_type` for `preparedproteins` / `structurereports`?**
   (§5b) Needed for the result-explorer filters.
4. **Who creates the `proteins` row for the prepared structure?** (§5a) Product-blocking.
5. **Do the three missing checkboxes get tool inputs, or does the mockup lose them?**
   (§2.1)
6. **Does Keep/Extract/Skip get a designed three-state treatment?** (§2.3)
7. **How does a user get a `protein.id` for an uploaded, non-PDB structure?** Only
   `deeporigin.pdb-import` creates Protein entities. If the answer is "they cannot", the
   app only works on imported PDB entries.
8. **Cofactor ids** are derived, not observed — no fixture anywhere contains one. Expect
   `cofactor:A:MN:501:` / `cofactor:A:ZN:302:`; confirm on a real metalloprotein.

---

## 9. Field and enum reference

Verified against the tool definition. **Read the Closed/Open column before typing anything
as a union** — an OPEN set is pinned only in producer code and can grow without a schema
change, so accept any string and render a fallback.

### 9.1 `deeporigin.target-preparation` inputs

| Field | Type | Required | Default | Enum |
| --- | --- | --- | --- | --- |
| `action` | string | **always** | none | **Closed:** `recommend` \| `prepare` |
| `protein` | object `{id, file_path}` | **always** | none | both sub-keys required |
| `selection` | object | prepare only | none | forbidden on recommend |
| `model_missing_loops` | boolean | prepare only | **none** | — |
| `pocket` | object | prepare only | none | forbidden on recommend |
| `pdb_id` | string | when loops on | none | `^[A-Za-z0-9]{4}$` |

`pocket`:

| Field | Type | Required | Default | Enum |
| --- | --- | --- | --- | --- |
| `mode` | string | no | `auto-find` | **Closed:** `auto-find` \| `define-by-selection` \| `from-crystal-ligand` |
| `pocket_count` | integer ≥1 | if `auto-find` | **none** | — |
| `pocket_min_size` | number ≥1 | if `auto-find` | **none** | — |
| `selections[]` | array, minItems 1 | if `define-by-selection` | none | items `{kind, author}` |
| `selections[].kind` | string | yes | — | **Closed:** `residue` \| `ligand` \| `cofactor` |
| `crystal_ligand` | object | if `from-crystal-ligand` | none | `{component_id?, file_path?, ligand_id?}` |
| `pocket_radius` | number >0 | no | `10.0` | — |
| `align_to_pocket` | boolean | no | `false` | — |
| `box_geometry` | string | no | runtime `ligand-extents` | **Closed:** `ligand-extents` \| `fixed-radius` |
| `box_padding` | number ≥0 | no | runtime `4.0` | — |

`selection`:

| Field | Type | Required | Enum |
| --- | --- | --- | --- |
| `source_sha256` | string | yes | must equal the analyzed digest |
| `analyzer_version` | string | yes | currently `1.1.0` |
| `decisions` | `Record<componentId, string>` | yes | **Closed:** `keep` \| `skip` — **no `review`**, and exhaustive |

### 9.2 `recommendation.components[]`

| Field | Type | Enum | Closed? |
| --- | --- | --- | --- |
| `id` | string | §9.3 | — |
| `kind` | string | `chain` \| `ligand` \| `cofactor` \| `water` | **Closed** |
| `subtype` | string | chain: `protein`/`nucleic`/`mixed`/`other`/`short_peptide` · water: `coordinating`/`crystal` · cofactor: `metal`/`ion`/`organic`/`other` · ligand: `organic` | **OPEN** — bare `string` in schema. **But semantically load-bearing:** `nucleic`/`mixed`/`other` chains cannot be kept, so string-compare those three. |
| `label` | string | — | display only |
| `recommendation` | string | `keep` \| `skip` \| `review` | **Closed** |
| `reason` | string | — | required; always populated; **render this** |
| `reason_code` | string | 12 values today: `keep_chain`, `duplicate_chain`, `non_protein_chain`, `short_peptide`, `coordinating_water`, `water_review`, `coordinated_metal`, `crystallization_artifact`, `under_coordinated_metal`, `organic_cofactor`, `recognized_ligand`, `ambiguous_ligand` | **OPEN** — iconography at most |
| `author` | object `{chain_id, resname?, resseq?, icode?, label_asym_id?}` | only `chain_id` required | **the Mol\* addressing key** |
| `evidence` | object | free-form | tool-owned |

Top level also carries `source_sha256`, `analyzer_version`, and **`chain_id_mapping`**
(required) — author chain ids remapped to single-character PDB chain ids for mmCIF.

### 9.3 Component id grammar

One builder, uniform shape. `icode` is `""` in the common case, hence the trailing colon.

| Kind | Grammar | Example |
| --- | --- | --- |
| chain | `chain:<chain_id>` | `chain:A` |
| ligand | `<kind>:<chain>:<resname>:<resseq>:<icode>` | `ligand:A:IBP:100:` |
| cofactor | same | `cofactor:A:MN:501:` |
| water | same | `water:A:HOH:310:` |

**Treat ids as opaque.** Echo back exactly what the inventory returned; address the viewer
from `author`. Fixtures in both `do-dd-client` (`ligand:LIG:A:100`) and `platform-toolbox`
(`ligand:A:120:IBP`) use shapes the analyzer cannot emit — do not copy ids out of tests.

### 9.4 `structure_reports[]`

| Field | Type | Enum | Closed? |
| --- | --- | --- | --- |
| `grade` | string | `A` \| `B` \| `C` \| `D` | **Closed** |
| `weighted_score` | number | ≥0.8 A · 0.66–0.79 B · 0.5–0.65 C · <0.5 D | from the tool — never recompute |
| `metadata_source` | string | `rcsb` \| `file_header` \| `file_header+rcsb` | **Closed** |
| `field_status` | object | six fixed keys `organism`/`method`/`resolution`/`coverage`/`rfree`/`inhibitor`, each `value` \| `not_applicable` \| `unknown` | **Closed** — keys *and* values |
| `report_role` | string | `source` \| `prepared` | **Closed** |
| `*_score` (6) | number | — | resolution, coverage, rfree, inhibitor, method, organism |
| `coverage` | number \| null | 0–1 fraction | render as `82.6%` |
| `has_ligand` | boolean \| null | — | the `Ligand` pill |
| `method` / `organism` | string \| null | raw | — |
| `method_class` | string \| null | over the wire: `cryo-em` \| `x-ray` \| `nmr` \| `other` \| **`null`** | **OPEN in schema** — `unknown` is converted to `null` by the producer; handle null |
| `organism_class` | string \| null | over the wire: `human` \| `mammal` \| `vertebrate` \| `other` \| **`null`** | **OPEN in schema**, same null conversion |
| `resolution` / `rfree` | number \| null | — | rfree is X-ray only |
| `pdb_id` / `protein_id` | string \| null | — | `protein_id` is `x-key` |
| `source_sha256` | string | — | absent in PDB-ID-only mode |

### 9.5 `pockets[]`

`protein_id` (`x-key`), `file_path`, `pocket_center` (3-array), `box`
(`{box_size_x, box_size_y, box_size_z, rotation_deg[3]}` — PCA-aligned, prefer over the
deprecated flat `box_size_*`), `volume`, `total_SASA`, `polar_SASA`, `apolar_SASA`,
`polar_apolar_SASA_ratio`, `hydrophobicity`, `drugability_score`, `polarity`,
`pocket_count`, `pocket_min_size`.

Every metric is **`null` for `define-by-selection` and `from-crystal-ligand`** — a pockets
table must render nulls, not zeros.

### 9.6 Platform execution status

**Closed:** `Quoted` · `Created` · `Queued` · `Running` · `Completed` · `Succeeded` ·
`Failed` · `Cancelled` · `InsufficientFunds` · `FailedQuotation`. `Succeeded` is legacy for
`Completed`. Terminal is everything except `Created`/`Queued`/`Running`.

Tools answer **HTTP 200 with `status: "Failed"`** and no outputs;
`statusReason.message` is itself a JSON string. `useRunStep` already handles this.

### 9.7 App form fields

| Section | Field id | Type | Default | Notes |
| --- | --- | --- | --- | --- |
| structure-components | `components_summary` | `component-summary` ⚠ new | — | read-only chips |
| preparation | `model_missing_loops` | `boolean` | `true` | disable + explain when no `pdb_id` |
| preparation | `pdb_id` | `string` | from entity | required when loops on; `^[A-Za-z0-9]{4}$` |
| pocket | `pocket_mode` | `enum` | `auto-find` | drives `showIf` on the next two |
| pocket | `pocket_count` | `number` | `5` | `auto-find` only |
| pocket | `pocket_min_size` | `number` | `30` | `auto-find` only |
| run-details | `run_name` | `string` | — | required; `x-body-key: "name"` |
| *(hidden)* | `action` | `enum` | `prepare` | constant into `inputs` |

Existing engine `FieldType` (**closed**): `string` · `number` · `boolean` · `enum` ·
`multi-enum` · `result` · `segmented` · `radio` · `entity-ref` · `action` · `site-select` ·
`admet-properties` · `structure`. This app adds `component-summary`.

### 9.8 Store state keys

| Key | Written by | Read by |
| --- | --- | --- |
| `selectedProteins` | proteins table | both steps, all tiles, submit |
| `stepOutputs.structureReport` | structure-report step | report card |
| `stepOutputs.recommendation` | inventory step | filtering tile, viewer, sidebar chips |
| `stepOutputs.selection` | **filtering tile** (gap #2) | sidebar chips, submit via `x-from-state` |
| `stepErrors[stepId]` | `useRunStep` | `ActionField`, tile error states |

Cleared by `useSelectionResets` on selection change.

### 9.9 Engine enums (all closed, all existing)

`AppMode` `edit`\|`results` · `ExecutionVisibility` `visible`\|`hidden` ·
`ResultsMergeMode` `flat`\|`per-group`\|`aggregate` · `ResultsAppliesTo`
`proteins`\|`ligands` · `LeafComponent` (10, +2) · `RendererName` (4, +1).

---

## 10. Manifest sketch

```jsonc
{
  "id": "target-prep",
  "toolKey": "deeporigin.target-preparation",
  "toolVersion": "1",
  "identityHue": 210,
  "url": "/target-prep",
  "name": "Target Preparation",

  "steps": [
    { "id": "structure-report",
      "toolKey": "deeporigin.structure-report", "toolVersion": "latest",
      "visibility": "hidden",
      "inputSchema": { "type": "object", "properties": {
        "protein": { "type": "object", "x-data-type": "Protein", "properties": {
          "id":        { "type": "string", "x-data-type": "Protein.id" },
          "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } },
        "pdb_id":  { "type": "string", "x-data-type": "Protein.pdb_id" } } },
      "bind": { "structure_reports": "structureReport" } },

    { "id": "recommend",
      "toolKey": "deeporigin.protein-prep", "toolVersion": "2",
      "sync": true,                        // engine gap #1
      "visibility": "hidden",
      "inputSchema": { "type": "object", "properties": {
        "protein": { "type": "object", "x-data-type": "Protein", "properties": {
          "id":        { "type": "string", "x-data-type": "Protein.id" },
          "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } } } },
      "bind": { "recommendation": "recommendation" } }
      // `action: "recommend"` supplied via the ActionField's `inputs` override
  ],

  "layout": { "tiles": [
    { "id": "left", "component": "VerticalMosaic", "splitPercentage": 35, "tiles": [
      { "id": "protein-table",    "component": "Table" },
      { "id": "structure-filter", "component": "StructureFiltering" } ] },
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
                  "reportKey": "structureReport", "stepId": "recommend" } },
    "protein-viewer": { "component": "ProteinViewer", "title": "Mol* Viewer", "config": {
      "molstar": { "edit":    { "renderer": "renderStructureComponents" },
                   "results": { "renderer": "renderStructureAndPockets" } } } }
  },

  "parameters": {
    "resultLabel": { "singular": "prepared structure", "plural": "prepared structures" },
    "requiredEntities": [ { "type": "Protein",
      "message": "Select a protein before running Target Preparation" } ],
    "validations": [
      { "type": "state-resolved", "stateKey": "selection", "path": "decisions",
        "disallow": "review",
        "message": "Resolve every reviewed component before running" }
      // plus: exhaustive vs inventory, >=1 protein chain kept,
      //       no nucleic/mixed/other chain kept  (engine gap #3)
    ],
    "sections": [
      { "id": "structure-components", "label": "Structure Components",
        "fields": [ { "id": "components_summary", "type": "component-summary",
                      "optionsFrom": "selection", "labelsFrom": "recommendation" } ] },
      { "id": "preparation", "label": "Preparation", "fields": [
          { "id": "model_missing_loops", "type": "boolean", "default": true,
            "label": "Add missing loops", "helpText": "Models unresolved regions. Requires a PDB ID." },
          { "id": "pdb_id", "type": "string", "label": "PDB ID",
            "showIf": { "field": "model_missing_loops", "equals": true },
            "required": true, "pattern": "^[A-Za-z0-9]{4}$" } ] },
      { "id": "pocket", "label": "Pocket Detection", "collapsible": true, "fields": [
          { "id": "pocket_mode", "type": "enum", "default": "auto-find", "options": [
              { "label": "Auto-find",           "value": "auto-find" },
              { "label": "Define by selection", "value": "define-by-selection" },
              { "label": "From crystal ligand", "value": "from-crystal-ligand" } ] },
          { "id": "pocket_count", "type": "number", "default": 5,
            "showIf": { "field": "pocket_mode", "equals": "auto-find" } },
          { "id": "pocket_min_size", "type": "number", "default": 30,
            "showIf": { "field": "pocket_mode", "equals": "auto-find" } } ] },
      { "id": "run-details", "label": "Run Details", "fields": [
          { "id": "run_name", "type": "string", "label": "Run Name", "required": true,
            "x-body-key": "name", "x-run-name-default": "TargetPrep" },
          { "id": "action", "type": "enum", "default": "prepare", "hidden": true,
            "options": [ { "label": "Prepare", "value": "prepare" } ] } ] }
    ]
  },

  "inputSchema": { "type": "object",
    "required": ["action", "protein", "selection", "model_missing_loops", "pocket"],
    "properties": {
      "action":  { "type": "string",  "x-user-input": true },
      "protein": { "type": "object", "x-data-type": "Protein", "properties": {
        "id":        { "type": "string", "x-data-type": "Protein.id" },
        "file_path": { "type": "string", "x-data-type": "Protein.file_path" } } },
      "selection":           { "type": "object",  "x-from-state": "selection" },
      "model_missing_loops": { "type": "boolean", "x-user-input": true },
      "pdb_id":              { "type": "string",  "x-user-input": true,
                               "x-show-if": { "field": "model_missing_loops", "equals": true } },
      "pocket": { "type": "object", "properties": {
        "mode":            { "type": "string",  "x-from-form": { "field": "pocket_mode" } },
        "pocket_count":    { "type": "integer", "x-user-input": true,
                             "x-show-if": { "field": "pocket_mode", "equals": "auto-find" } },
        "pocket_min_size": { "type": "number",  "x-user-input": true,
                             "x-show-if": { "field": "pocket_mode", "equals": "auto-find" } } } } } }
}
```

Two things to watch when writing this for real:

- Root `additionalProperties: false` means `buildToolPayload` must emit **exactly** these
  keys. `x-show-if` (which drops a property entirely) is the right tool for `pdb_id` and the
  auto-find-only pocket fields — a `null` or `""` would be rejected.
- The app's `id` / `url` / filename (`target-prep`) are independent of the tool key.
  `/target-prep` vs `/target-preparation` is a product call; settle it before the route
  ships.

---

## 11. Test surface

- **Unit** — `buildToolPayload` emitting exactly the declared keys under each pocket mode;
  the decision reducer (exhaustiveness, `review` never auto-resolved, illegal keeps
  rejected); hex→Mol\*-colour mapping; `chain_id_mapping` application.
- **Storybook** — `StructureReportCard` across A–D and every `field_status` state;
  `StructureFiltering` with a realistic inventory (one chain, one ligand, two cofactors,
  **300 waters**) to prove the bulk/virtualized path.
- **E2E** — select a protein, assert both steps fire and both tiles populate; toggle a
  component and assert the viewer recolours; assert Run stays disabled while any `review`
  remains; assert a FAILED run with retained artifacts still renders them.
- **Contract** — assert the manifest's `inputSchema` keys are a subset of the published
  tool definition's `properties`, so a backend schema change fails CI rather than
  production.

---

## 12. Phasing

**Phase 0 — unblock (§8.1, §8.2).** Confirm the tool is registered in a usable environment,
and verify end to end that a `protein-prep`-produced Selection is accepted by
`target-preparation` `prepare`. Both are cheap; both invalidate the architecture if they
fail. Nothing below should start before §8.2 is answered.

**Phase 1 — the app.** Manifest + registration, proteins table, both pre-submit steps, both
new tiles, the Mol\* renderer, all four engine gaps, and the client-side validations in §7.

Parallelizable order:

1. Engine gaps #1 and #2 — everything depends on them, both small.
2. `StructureReportCard` and the Mol\* renderer — independent of each other and of the
   filtering tile; both Storybook-testable against fixtures with no backend.
3. `StructureFiltering` + engine gap #3. **The largest single piece** — §2.3, §2.4 and §2.5
   all live here.
4. Manifest, sidebar (gap #4), registration, results view (including the FAILED-with-
   artifacts path).

**Phase 2 — reports as table columns.** `results.detailColumns` / `aggregates` against
`structurereports`. No backend work; needs the `result_type` string from §8.3.

**Not scheduled — backend prerequisites.** §5a (prepared protein → `proteins` row) gates
Target Prep being useful to Docking and ABFE at all, and §2.1 / §2.3 are PRD decisions that
change the sidebar and the palette. None of them block phase 1 from starting; all of them
change what "done" means.
