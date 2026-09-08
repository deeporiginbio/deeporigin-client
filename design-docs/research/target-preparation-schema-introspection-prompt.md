# Agent prompt — introspect `deeporigin.target-preparation`

> **This prompt has been run.** Its output is
> [`target-preparation-backend-report.md`](./target-preparation-backend-report.md), and
> [`target-prep-do-studio-implementation-plan.md`](./target-prep-do-studio-implementation-plan.md)
> has been rewritten against it. Kept here as the reusable instrument: re-run it against a
> later tool version, or adapt it for the next tool that needs the same treatment.

Hand the block below to an agent that has **read access to `deeporiginbio/platform-toolbox`**
(the session that wrote it did not — `add_repo` requires push access there).

---

```
You are researching a backend tool so a frontend team can build a UI against it. Report
only what you can verify by reading files — never infer a field name from a PRD, a
directory name, or another tool's schema. Where something is genuinely absent, say
"not present" rather than guessing.

REPO: deeporiginbio/platform-toolbox
TOOL: deeporigin.target-preparation
ROOT: tools/target-preparation/

Read at least these, and follow anything they import:
  tools/target-preparation/workflow/tool-definition.json
  tools/target-preparation/workflow/workflow.yaml
  tools/target-preparation/workflow/preflight-service.yaml
  tools/target-preparation/README.md
  tools/target-preparation/RELEASE-NOTES.MD
  images/preflight/src/preflight_service/routes/target_preparation.py
  images/preflight/src/preflight_service/app.py
  tests/test_target_preparation.py
  tests/test_target_preparation_schema.py
  tests/test_target_preparation_composite.py
  images/preflight/tests/test_target_preparation_route.py
  CONTEXT.md   (whatever it says about this tool)

Answer these, in this order. Quote the file and line you got each answer from.

## 1. Identity
- The exact registered tool key. Confirm it is `deeporigin.target-preparation` and not
  something else — the directory name is not authoritative.
- Current version, and the MAJOR version a UI should pin (the platform resolves
  minor/patch at execution time, so the UI pins the major only).
- Is it registered/enabled in dev, staging, prod? Where is that declared?

## 2. Input schema — the critical part
Reproduce `inputs` from the tool definition VERBATIM as JSON. Then, as a table, one row
per top-level property:
  name | JSON type | required? | default | enum values (if any) | one-line meaning

Then answer specifically:
- What is the protein input called, and what is its shape? (`{id, file_path}`? a bare
  string? something else?)
- Is there a component-selection input — the keep/skip decisions per chain / ligand /
  cofactor / water? What is it called, and what is its exact shape? Does it carry a
  digest field (`source_sha256` or similar) and an analyzer version?
- Do these four toggles exist as inputs, and what are they actually called?
    * add missing atoms & residues
    * add missing loops  (compare `model_missing_loops` on deeporigin.protein-prep)
    * protonate
    * find pockets
  For each: present or absent, exact name, type, default.
- Is `pdb_id` an input? Is it conditionally required (e.g. only when loop modelling is
  on)? Is that enforced in the schema, in the workflow, or in the preflight route?
- Is there an output-naming input (the UI mockup shows an "Output Property Name" field)?
  What is it called and what does it actually name?
- Does the schema set `additionalProperties: false` anywhere? (If so, sending an extra
  key is a hard failure, not a warning — the UI must send exactly the declared keys.)

## 3. Output schema
Reproduce `outputs` VERBATIM. Then:
- Exact key names under `jobOutputs` for: the prepared structure, the pockets, the
  structure report.
- Does it write data-platform records? Specifically:
    * Does it CREATE A NEW ROW IN THE `proteins` ENTITY TABLE for the prepared structure?
      This is the single most important question in this report — if it only emits
      `preparedprotein` result rows, the prepared protein is invisible to Docking and
      ABFE, and the feature does not achieve its purpose. Show the code that writes it,
      or state clearly that nothing does.
    * What `result_type` values does it emit (e.g. `pocket`, `preparedprotein`)?
    * Does anything emit an indexed structure-report row (a `structurereport` result
      type)? Or does the report exist only in `jobOutputs`?
- Does the output PDB/CIF carry the prepared-protein stamp (`REMARK  99 DO_PREPARED`
  for PDB, `_deeporigin.prepared` for mmCIF)? Downstream tools skip cleanup based on it.

## 4. The preflight route — decides a UI architecture choice
Read images/preflight/src/preflight_service/routes/target_preparation.py closely.
- What is the HTTP method, path, request shape and response shape?
- DOES IT RETURN A COMPONENT INVENTORY — a list of chains / ligands / cofactors /
  waters with a keep/review/skip recommendation each? If yes, reproduce one example
  response, because the UI would call this instead of running a
  `deeporigin.protein-prep` `action: "recommend"` execution per row click.
- What does it reject, and with what error shape? Each rejection is a validation the UI
  should enforce client-side first, so the user sees a disabled button instead of a
  server error.
- Is it authenticated the same way as the tools API? Is it reachable from the browser,
  or only from inside the cluster? (If it is cluster-internal, the UI cannot call it and
  this whole option is void — say so explicitly.)

## 5. Workflow behaviour
From workflow.yaml and the composite tests:
- The ordered steps, and which tool/image each runs.
- Does it run its own component analysis internally? If so, does a caller-supplied
  selection still bind, or would an internal re-analysis override the user's choices?
  (If it overrides them, the UI's entire filtering panel is decorative — flag this
  loudly.)
- Are steps conditional on the input toggles, or do they always run?
- Partial failure: if pocket-finding or the report fails AFTER the protein was
  successfully prepared, what is the final execution status, and is the prepared
  protein still written?
- Is it synchronous or asynchronous? Does it accept a top-level `sync` body key?
- Is it billable, and does it quote as one line item or several?

## 6. Component identity — needed for 3D rendering
The UI colours each component in a Mol* viewer, which needs PDB addressing
(auth chain id, residue name, residue sequence number) per component.
- What identifies a component in the selection input — an opaque id string, a structured
  object, or both?
- If ids are strings, give the exact grammar per kind. In deeporigin.protein-prep the
  orders differ between kinds (`ligand:LIG:A:100` is resname:chain:resseq, but
  `water:A:HOH:310:` is chain:resname:resseq: with a trailing colon) — confirm whether
  target-preparation uses the same grammar or its own.
- WHAT DOES A COFACTOR COMPONENT ID LOOK LIKE? No fixture in do-dd-client contains one,
  and the UI mockup shows two (Mn2+, Zn2+). Find a real example.
- Is there an `author`-style block (chain_id / resname / resseq) alongside each
  component? That is what the renderer needs; without it the UI must parse ids.

## 7. Enum catalogues — closed or open?
For each of these, say whether the value set is CLOSED (fixed and safe for a frontend to
switch on exhaustively) or OPEN (tool-owned, may grow, so the frontend must accept any
string and render a fallback). Cite the definition, not an example.
- component `kind`
- component `subtype`   (do-dd-client accepts ANY string here — confirm whether that is
  deliberate)
- component `recommendation` / decision values
- component `reason_code`
- structure report `grade`, `metadata_source`, `field_status` values
- structure report `method_class`, `organism_class`
- anything else with a fixed value set

## 8. Anything that contradicts this summary of the intended UI
The UI is: user picks a protein from a table; a structure report and a component
inventory appear automatically; the user toggles Keep/Skip per component; a 3D viewer
colours kept vs excluded components; the user sets four boolean options and a run name;
one button submits ONE deeporigin.target-preparation execution. Call out anything in the
tool that makes this wrong, impossible, or more expensive than it looks.

## Output format
Markdown. Lead with a "Answers that change the UI plan" section of at most 6 bullets —
the things a frontend engineer would get wrong without this report. Then the numbered
sections above. Verbatim schema JSON in fenced blocks. Every claim carries a
file:line citation. End with an explicit list of anything you could NOT determine and
what would be needed to determine it.
```
