# `deeporigin.target-preparation` — backend research report for the frontend team

Repo `deeporiginbio/platform-toolbox`, branch `claude/upbeat-dijkstra-fdu502`
(tool content identical to `origin/main`). Every claim below is cited to a file
and line that I read. Nothing is inferred from a PRD, a directory name, or
another tool's schema.

---

## Answers that change the UI plan

1. **The mockup's four toggles do not exist. There is exactly one boolean input:
   `model_missing_loops`.** No "add missing atoms/residues", no "protonate", no
   "find pockets" input exists anywhere in the schema
   (`tools/target-preparation/workflow/tool-definition.json:207-432` — the full
   `properties` list is `action`, `model_missing_loops`, `pdb_id`, `pocket`,
   `protein`, `selection`). Protonation always runs
   (`CONTEXT.md:196`) and Pocket Finder always runs on `prepare`
   (`tools/target-preparation/workflow/workflow.yaml:191-199`, unconditional).

2. **It is two executions, not one.** `action: "recommend"` produces the
   component inventory; `action: "prepare"` consumes a digest-bound selection
   built from it. `prepare` requires `selection.source_sha256` +
   `selection.analyzer_version` + a *complete* `decisions` map
   (`tool-definition.json:399-431`), and the digest must match the analyzed
   file byte-for-byte
   (`packages/toolbox-core/src/toolbox_core/protein_component.py:821-826`).
   The single "Run" button cannot exist without a prior recommend execution.

3. **Nothing writes a row into the `proteins` entity table.** The prepared
   structure is emitted with `x-result-group: "preparedproteins"`
   (`tool-definition.json:723-724`), which is a dynamic `results__*` table, not
   the entity table. The only `create_protein` call site in the whole repo is
   `images/structure-report/src/pdb_import.py:259` (the `deeporigin.pdb-import`
   tool). Docking / ABFE will not see the prepared protein as a selectable
   Protein row. **Full detail in §3.**

4. **The component inventory has no indexed destination, and this tool always
   runs in workflow mode.** `recommendation` carries no `x-data-type` and no
   `x-result-group` (`tool-definition.json:726-832`), and
   `quote_target_preparation_method` hardcodes `"workflow"`
   (`images/preflight/src/preflight_service/validate_target_preparation.py:10-13`).
   The executions API "omits `jobOutputs` for workflow (Argo) mode by design"
   (`images/test-tool/src/test_tool/main.py:1065-1066`). By contrast
   `deeporigin.protein-prep` `recommend` quotes `"direct"` and therefore *does*
   return `jobOutputs`
   (`images/preflight/src/preflight_service/validate_protein_prep.py:141-143`).
   **The UI probably has to call `deeporigin.protein-prep` `recommend` to get an
   inventory it can actually read.**

5. **`keep` on a *ligand* means "extract to a separate file", not "retain in the
   protein".** `apply_selection` routes a kept `LIGAND` into
   `extract_ligand_keys`, and every other kept non-chain component into
   `keep_residue_keys`
   (`protein_component.py:1056-1063`); the reason string literally reads
   `"Recognized ligand {resname}; keep means extract"`
   (`protein_component.py:776-780`). A green "kept" swatch on a ligand in the
   Mol* viewer is wrong — the ligand is gone from the prepared structure.

6. **The preflight route is not a component-inventory endpoint.** It is a
   create-time billing/validation gate that returns `{method, counts,
   validations?}` and nothing else
   (`images/preflight/src/preflight_service/routes/target_preparation.py:32-72`),
   and it is a cluster-internal Knative service invoked by the platform, not by
   the browser (`tools/target-preparation/workflow/preflight-service.yaml:1-8`,
   `CONTEXT.md:1180-1192`). **That UI option is void — see §4.**

---

## 1. Identity

### Registered key

`"key": "deeporigin.target-preparation"`
— `tools/target-preparation/workflow/tool-definition.json:440`

Confirmed by test: `assert data["key"] == "deeporigin.target-preparation"`
— `tests/test_target_preparation.py:42`.

The directory name is `target-preparation`; the internal Argo/manifest name is
`"name": "deeporigin-target-preparation"` (`tool-definition.json:455`,
`workflow.yaml:4`). Neither is the registered key. The key is the dotted form.

### Version

| Field | Value | Citation |
|---|---|---|
| `version` | `1.0.0` | `tool-definition.json:1061` |
| `toolManifestVersion` | `6` | `tool-definition.json:1060` |
| `billingCode` | `DO_POCKET_FINDER` | `tool-definition.json:2` |
| `mcpFeatured` | `true` | `tool-definition.json:450` |
| `metadata` | `{"PocketFinder": "0.1.6", "SystemPrep": "0.3.31"}` | `tool-definition.json:451-454` |

**A UI should pin major `1`.** Per ADR-0021 the major only moves when the schema
contract breaks (`docs/adr/0021-tool-version-majors-track-schema-contract.md:30`
— `x-data-type`, `x-key`, and `x-result-group` changes count as breaking).

Release notes confirm this is the first release:
`## [1.0.0] - 2026-09-08 ... Initial deeporigin.target-preparation orchestration
Tool (DDOS-7091)` — `tools/target-preparation/RELEASE-NOTES.MD:5-10`.

### dev / staging / prod registration

**Not declared anywhere in the repo.** There is no per-environment enable list,
manifest, or flag. Registration is driven entirely by CI:

- `.github/workflows/deploy-toolbox.yml:3-8` — a push to `main` touching
  `tools/**` triggers the deploy.
- `.github/workflows/deploy-toolbox.yml:69-73` — a push always resolves
  `DEPLOY_ENV=dev`.
- `.github/workflows/deploy-toolbox.yml:34-36` — "Dev builds workflow bundles
  from git, registers on dev, and publishes to CodeArtifact. Staging/prod only
  download an existing CodeArtifact version and register." Staging/prod are
  `workflow_dispatch`-only (`deploy-toolbox.yml:31-41`, env choices
  `dev` / `stgn` / `prodn`).

The tool directory *is* present on `origin/main`
(`git ls-tree origin/main -- tools/target-preparation` returns all five files),
so a dev registration should have been attempted. But the repo's own integration
test still says otherwise:

> `pytest.skip("remote E2E waits until deeporigin.target-preparation is published")`
> — `tests/test_target_preparation.py:231`

and

> `"target-preparation E2E needs make start for preflight + child servings and a
> seeded tool registration"` — `tests/test_target_preparation.py:227-230`

**Verdict: cannot be determined from files.** Ask platform to run
`.github/workflows/dump-enabled-tools.yml` against dev/stgn/prodn — that
workflow exists precisely to dump the enabled-tool list per environment
(`.github/workflows/dump-enabled-tools.yml:37-52`).

---

## 2. Input schema

### Verbatim `inputs`

Source: `tools/target-preparation/workflow/tool-definition.json:3-439`
(re-serialized with 2-space indent and sorted keys, which is how it is stored on
disk — see `git log` commit `6e39b9b "fix: sort tool-definition.json keys"`).

```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "additionalProperties": false,
  "allOf": [
    {
      "if": {
        "properties": {
          "action": {
            "const": "recommend"
          }
        },
        "required": [
          "action"
        ]
      },
      "then": {
        "not": {
          "anyOf": [
            {
              "required": [
                "selection"
              ]
            },
            {
              "required": [
                "model_missing_loops"
              ]
            },
            {
              "required": [
                "pocket"
              ]
            }
          ]
        },
        "required": [
          "protein"
        ]
      }
    },
    {
      "if": {
        "properties": {
          "action": {
            "const": "prepare"
          }
        },
        "required": [
          "action"
        ]
      },
      "then": {
        "allOf": [
          {
            "required": [
              "protein",
              "selection",
              "model_missing_loops",
              "pocket"
            ]
          },
          {
            "else": {
              "required": [
                "pdb_id"
              ]
            },
            "if": {
              "properties": {
                "model_missing_loops": {
                  "const": false
                }
              },
              "required": [
                "model_missing_loops"
              ]
            },
            "then": {}
          },
          {
            "if": {
              "properties": {
                "pocket": {
                  "properties": {
                    "mode": {
                      "const": "auto-find"
                    }
                  }
                }
              }
            },
            "then": {
              "properties": {
                "pocket": {
                  "allOf": [
                    {
                      "required": [
                        "pocket_count",
                        "pocket_min_size"
                      ]
                    },
                    {
                      "not": {
                        "required": [
                          "crystal_ligand"
                        ]
                      }
                    }
                  ]
                }
              }
            }
          },
          {
            "if": {
              "properties": {
                "pocket": {
                  "properties": {
                    "mode": {
                      "const": "define-by-selection"
                    }
                  },
                  "required": [
                    "mode"
                  ]
                }
              }
            },
            "then": {
              "properties": {
                "pocket": {
                  "allOf": [
                    {
                      "required": [
                        "selections"
                      ]
                    },
                    {
                      "not": {
                        "required": [
                          "crystal_ligand"
                        ]
                      }
                    }
                  ]
                }
              }
            }
          },
          {
            "if": {
              "properties": {
                "pocket": {
                  "properties": {
                    "mode": {
                      "const": "from-crystal-ligand"
                    }
                  },
                  "required": [
                    "mode"
                  ]
                }
              }
            },
            "then": {
              "properties": {
                "pocket": {
                  "allOf": [
                    {
                      "required": [
                        "crystal_ligand"
                      ]
                    },
                    {
                      "not": {
                        "required": [
                          "selections"
                        ]
                      }
                    },
                    {
                      "not": {
                        "required": [
                          "pocket_count"
                        ]
                      }
                    },
                    {
                      "not": {
                        "required": [
                          "pocket_min_size"
                        ]
                      }
                    }
                  ]
                }
              }
            }
          }
        ]
      }
    }
  ],
  "description": "Target Preparation: recommend (source Structure Report + Protein Prep recommend) or prepare (Protein Prep prepare + prepared Structure Report + Pocket Finder).",
  "properties": {
    "action": {
      "description": "recommend: source report + component recommendations. prepare: apply Selection, prepare protein, final report, pockets.",
      "enum": [
        "recommend",
        "prepare"
      ],
      "type": "string",
      "x-display-name": "Action",
      "x-user-input": true
    },
    "model_missing_loops": {
      "description": "prepare only. Required boolean: when true, SystemPrep models missing loops and pdb_id is required.",
      "type": "boolean",
      "x-display-name": "Model missing loops",
      "x-user-input": true
    },
    "pdb_id": {
      "description": "4-character PDB ID for Structure Report metadata and loop templates. Required on prepare when model_missing_loops is true.",
      "pattern": "^[A-Za-z0-9]{4}$",
      "type": "string",
      "x-display-name": "PDB ID",
      "x-user-input": true
    },
    "pocket": {
      "additionalProperties": false,
      "description": "Nested Pocket Finder configuration (prepare only)",
      "properties": {
        "align_to_pocket": {
          "default": false,
          "description": "When mode is define-by-selection, run PCA on selection atoms for box.rotation_deg (oriented cube). Default false.",
          "type": "boolean",
          "x-display-name": "Align to Pocket",
          "x-user-input": true
        },
        "box_geometry": {
          "description": "from-crystal-ligand only. ligand-extents sizes the PCA-aligned box from ligand heavy-atom extents plus box_padding; fixed-radius uses a cube of edge 2 * pocket_radius. Runtime default is ligand-extents.",
          "enum": [
            "ligand-extents",
            "fixed-radius"
          ],
          "type": "string",
          "x-display-name": "Box Geometry",
          "x-user-input": true
        },
        "box_padding": {
          "description": "from-crystal-ligand only. Padding in angstroms added on each side of ligand PCA extents when box_geometry is ligand-extents. Runtime default is 4.0.",
          "minimum": 0,
          "type": "number",
          "x-display-name": "Box Padding",
          "x-user-input": true
        },
        "crystal_ligand": {
          "description": "from-crystal-ligand: identify a Protein Prep extracted ligand by component_id (resolved after prepare) or supply file_path/ligand_id",
          "properties": {
            "component_id": {
              "description": "Protein Prep extracted_ligands component_id",
              "type": "string",
              "x-display-name": "Component ID",
              "x-user-input": true
            },
            "file_path": {
              "description": "UFA path to the extracted ligand file (PDB, mmCIF, or SDF); omit when ligand_id is set",
              "format": "file",
              "type": "string",
              "x-data-type": "ExtractedLigand.file_path",
              "x-deeporigin-path": "/deeporigin/inputs/crystal_ligand.pdb"
            },
            "ligand_id": {
              "description": "Bare ligand/residue code (e.g. \"IBP\") to locate directly within the supplied protein file; omit when file_path is set. Errors if it matches zero or more than one residue instance.",
              "type": "string",
              "x-display-name": "Ligand ID",
              "x-user-input": true
            }
          },
          "type": "object",
          "x-display-name": "Crystal Ligand",
          "x-user-input": true
        },
        "mode": {
          "default": "auto-find",
          "description": "auto-find runs the PocketFinder classifier; define-by-selection builds one pocket from residue/ligand/cofactor selections plus radius; from-crystal-ligand builds one PCA-aligned pocket from an extracted ligand PDB.",
          "enum": [
            "auto-find",
            "define-by-selection",
            "from-crystal-ligand"
          ],
          "type": "string",
          "x-display-name": "Mode",
          "x-user-input": true
        },
        "pocket_count": {
          "description": "Number of pockets to find (auto-find). Required when self_test is false and mode is auto-find; when self_test is true, optional (defaults to 1).",
          "minimum": 1,
          "type": "integer",
          "x-display-name": "Pocket Count",
          "x-user-input": true
        },
        "pocket_min_size": {
          "description": "Minimum size of pockets to consider (auto-find). Required when self_test is false and mode is auto-find; when self_test is true, optional (defaults to 30).",
          "minimum": 1,
          "type": "number",
          "x-display-name": "Minimum Pocket Size",
          "x-user-input": true
        },
        "pocket_radius": {
          "default": 10.0,
          "description": "Half-edge of the docking cube in angstroms for define-by-selection, or for from-crystal-ligand when box_geometry is fixed-radius. Lab and box sizes are 2 * pocket_radius in those cases.",
          "exclusiveMinimum": 0,
          "type": "number",
          "x-display-name": "Pocket Radius",
          "x-user-input": true
        },
        "selections": {
          "description": "Residue, ligand, or cofactor selectors (define-by-selection). Author fields align with Protein Prep Component Recommendation.",
          "items": {
            "additionalProperties": false,
            "properties": {
              "author": {
                "additionalProperties": false,
                "description": "PDB/mmCIF author identity for the selected component",
                "properties": {
                  "chain_id": {
                    "description": "Author chain ID",
                    "type": "string"
                  },
                  "icode": {
                    "description": "Insertion code (empty when absent)",
                    "type": "string"
                  },
                  "resname": {
                    "description": "Residue or HETATM name",
                    "type": "string"
                  },
                  "resseq": {
                    "description": "Residue sequence number",
                    "type": "integer"
                  }
                },
                "required": [
                  "chain_id"
                ],
                "type": "object"
              },
              "kind": {
                "description": "Component kind being selected",
                "enum": [
                  "residue",
                  "ligand",
                  "cofactor"
                ],
                "type": "string"
              }
            },
            "required": [
              "kind",
              "author"
            ],
            "type": "object"
          },
          "minItems": 1,
          "type": "array",
          "x-display-name": "Selections",
          "x-user-input": true
        }
      },
      "type": "object",
      "x-display-name": "Pocket",
      "x-user-input": true
    },
    "protein": {
      "additionalProperties": false,
      "description": "Registered protein structure (id and file_path required)",
      "properties": {
        "file_path": {
          "description": "UFA path to a PDB or mmCIF structure",
          "format": "file",
          "type": "string"
        },
        "id": {
          "description": "Platform protein entity id (required for pocket indexing)",
          "type": "string",
          "x-data-type": "Protein.id"
        }
      },
      "required": [
        "file_path",
        "id"
      ],
      "type": "object",
      "x-data-type": "Protein"
    },
    "selection": {
      "additionalProperties": false,
      "description": "prepare only. Digest-bound complete keep/skip map for every analyzed component",
      "properties": {
        "analyzer_version": {
          "description": "Analyzer version that produced the recommendation",
          "type": "string"
        },
        "decisions": {
          "additionalProperties": {
            "enum": [
              "keep",
              "skip"
            ],
            "type": "string"
          },
          "description": "Map of component_id to keep or skip (no review)",
          "type": "object"
        },
        "source_sha256": {
          "description": "SHA-256 of the uploaded structure bytes",
          "type": "string"
        }
      },
      "required": [
        "source_sha256",
        "analyzer_version",
        "decisions"
      ],
      "type": "object",
      "x-display-name": "Selection",
      "x-user-input": true
    }
  },
  "required": [
    "action",
    "protein"
  ],
  "title": "Target Preparation input schema",
  "type": "object"
}```

### Top-level properties

Root `required`: `["action", "protein"]` — `tool-definition.json:433-436`.

| name | JSON type | required? | default | enum | meaning |
|---|---|---|---|---|---|
| `action` | `string` | **yes** (root) | none | `recommend`, `prepare` | Which branch of the composite runs. `tool-definition.json:208-217` |
| `protein` | `object` | **yes** (root) | none | — | The input structure; `{id, file_path}`, both required. `tool-definition.json:377-398` |
| `selection` | `object` | prepare only, **required** | none | — | Digest-bound complete keep/skip map. Forbidden on `recommend`. `tool-definition.json:399-431`, `:18-37`, `:56-63` |
| `model_missing_loops` | `boolean` | prepare only, **required** | **none — no default** | — | When true, SystemPrep models missing loops *and* `pdb_id` becomes required. Forbidden on `recommend`. `tool-definition.json:218-223`, `:26-30`, `:56-63` |
| `pocket` | `object` | prepare only, **required** | none | — | Nested Pocket Finder config. Forbidden on `recommend`. `tool-definition.json:231-376`, `:31-35`, `:56-63` |
| `pdb_id` | `string` | conditional (see below) | none | — | 4-char PDB ID, `pattern: ^[A-Za-z0-9]{4}$`. `tool-definition.json:224-230` |

Note there is **no schema-level `default`** on `action`, `model_missing_loops`,
or `protein`. The only `default`s in the whole input schema are nested inside
`pocket`: `pocket.align_to_pocket` → `false` (`:236-237`), `pocket.mode` →
`"auto-find"` (`:287-288`), `pocket.pocket_radius` → `10.0` (`:313-314`).
The Argo workflow has its own parameter defaults
(`workflow.yaml:9-21`: `action: "recommend"`, `model_missing_loops: "true"`),
but those are workflow-internal and are not the API contract.

### Nested `pocket` properties (prepare only)

| name | type | required | default | enum | meaning |
|---|---|---|---|---|---|
| `mode` | string | no (runtime default) | `auto-find` | `auto-find`, `define-by-selection`, `from-crystal-ligand` | Which pocket strategy. `:286-297` |
| `pocket_count` | integer ≥1 | **yes when `mode=auto-find`** | none | — | Number of pockets. `:298-304`, required at `:83-104` |
| `pocket_min_size` | number ≥1 | **yes when `mode=auto-find`** | none | — | Minimum pocket size. `:305-311`, required at `:83-104` |
| `selections` | array (minItems 1) | **yes when `mode=define-by-selection`** | none | — | Residue/ligand/cofactor selectors. `:320-371`, required at `:117-140` |
| `crystal_ligand` | object | **yes when `mode=from-crystal-ligand`** | none | — | `{component_id?, file_path?, ligand_id?}`. `:259-285`, required at `:153-176` |
| `pocket_radius` | number >0 | no | `10.0` | — | Half-edge of the docking cube. `:312-319` |
| `align_to_pocket` | boolean | no | `false` | — | PCA orientation for `define-by-selection`. `:235-241` |
| `box_geometry` | string | no | runtime `ligand-extents` | `ligand-extents`, `fixed-radius` | `from-crystal-ligand` only. `:242-251` |
| `box_padding` | number ≥0 | no | runtime `4.0` | — | `from-crystal-ligand` + `ligand-extents` only. `:252-258` |

Mode exclusivity is enforced by `not/required` clauses: `auto-find` forbids
`crystal_ligand` (`:105-110`); `define-by-selection` forbids `crystal_ligand`
(`:141-146`); `from-crystal-ligand` forbids `selections`, `pocket_count`, and
`pocket_min_size` (`:177-196`). Covered by
`tests/test_target_preparation_schema.py:109-143`.

### The protein input

**Name:** `protein`. **Shape:** an object with **both** `id` and `file_path`
required — it is *not* a bare string.

```json
"protein": {
  "additionalProperties": false,
  "description": "Registered protein structure (id and file_path required)",
  "properties": {
    "file_path": { "description": "UFA path to a PDB or mmCIF structure",
                   "format": "file", "type": "string" },
    "id":        { "description": "Platform protein entity id (required for pocket indexing)",
                   "type": "string", "x-data-type": "Protein.id" }
  },
  "required": ["file_path", "id"],
  "type": "object",
  "x-data-type": "Protein"
}
```
— `tool-definition.json:377-398`

Both are re-checked at create time in preflight with distinct check names
`target_preparation.protein.id` and `target_preparation.protein.file_path`
— `images/preflight/src/preflight_service/validate_target_preparation.py:36-51`.

### The component-selection input

**Name:** `selection`. **Prepare only** — supplying it on `recommend` is a hard
schema failure (`tool-definition.json:18-37`, test at
`tests/test_target_preparation_schema.py:65-70`).

```json
"selection": {
  "additionalProperties": false,
  "description": "prepare only. Digest-bound complete keep/skip map for every analyzed component",
  "properties": {
    "analyzer_version": { "description": "Analyzer version that produced the recommendation",
                          "type": "string" },
    "decisions": {
      "additionalProperties": { "enum": ["keep", "skip"], "type": "string" },
      "description": "Map of component_id to keep or skip (no review)",
      "type": "object"
    },
    "source_sha256": { "description": "SHA-256 of the uploaded structure bytes",
                       "type": "string" }
  },
  "required": ["source_sha256", "analyzer_version", "decisions"],
  "type": "object",
  "x-display-name": "Selection",
  "x-user-input": true
}
```
— `tool-definition.json:399-431`

- **Digest field: yes — `source_sha256`**, required. At runtime the prepare
  service recomputes the SHA-256 of the downloaded file and rejects a mismatch:
  `"selection.source_sha256 does not match the uploaded structure digest"`
  — `packages/toolbox-core/src/toolbox_core/protein_component.py:821-826`,
  called from `images/system-prep/src/prepare_protein_service.py:300-324`.
- **Analyzer version: yes — `analyzer_version`**, required, and checked against
  the analyzer that ran: `"selection.analyzer_version does not match the
  recommendation analyzer"` — `protein_component.py:828-838`. The current
  analyzer is `ANALYZER_VERSION = "1.1.0"` — `protein_component.py:37`.
- **`decisions` must be exhaustive and two-valued.** `review` is rejected:
  `"must be keep or skip (got ...); unresolved review is rejected"`
  (`protein_component.py:859-867`). Any missing component →
  `"selection.decisions is missing components: ..."`; any unknown id →
  `"selection.decisions has unknown components: ..."`
  (`protein_component.py:844-857`).
- **`decisions` is a flat `{component_id: "keep"|"skip"}` map.** There is no
  per-component object, no author block, no digest per row.

### The four toggles

| Mockup toggle | Present? | Exact name | Type | Default |
|---|---|---|---|---|
| add missing atoms & residues | **absent** | — | — | — |
| add missing loops | **present** | `model_missing_loops` | `boolean` | **no default; required on `prepare`** |
| protonate | **absent** | — | — | — |
| find pockets | **absent** | — | — | — |

- `model_missing_loops` — `tool-definition.json:218-223`; required on `prepare`
  at `:56-63`; forbidden on `recommend` at `:26-30`. Same name as on
  `deeporigin.protein-prep` (`tools/protein-prep/workflow/tool-definition.json`,
  and `CONTEXT.md:196-197`). Note the *parent* has no default, whereas
  `deeporigin.protein-prep` treats an omitted flag as on
  (`images/preflight/src/preflight_service/validate_protein_prep.py:128-132`).
- **Protonation is unconditional.** `CONTEXT.md:196`: "Protonation always runs
  on prepare." The `protonate_protein` field exists only as an **output**
  (`tool-definition.json:711-716`).
- **Pocket finding is unconditional on `prepare`.** The `pocket-finder` DAG task
  has no `when:` guard — `workflow.yaml:191-199`. The only choice is *which
  mode*, not *whether*.
- **The prepared Structure Report is also unconditional** —
  `workflow.yaml:178-190`, no `when:`.

### `pdb_id`

**Yes, it is an input** (`tool-definition.json:224-230`), typed `string` with
`pattern: "^[A-Za-z0-9]{4}$"`.

**Conditionally required: yes.** On `prepare`, `pdb_id` is required unless
`model_missing_loops` is exactly `false`:

```json
{ "if":   { "properties": { "model_missing_loops": { "const": false } },
            "required": ["model_missing_loops"] },
  "then": {},
  "else": { "required": ["pdb_id"] } }
```
— `tool-definition.json:64-81`

Enforced in **two** places:
1. **The JSON Schema** — `tool-definition.json:64-81`. Test:
   `tests/test_target_preparation_schema.py:82-86`.
2. **The preflight route** — `validate_target_preparation.py:78-87`, emitting
   check name `target_preparation.pdb_id` with detail
   `"pdb_id is required when model_missing_loops is true"`.

It is **not** enforced in the workflow: `workflow.yaml:14-15` defaults it to
`""` and `workflow.yaml:410-413` simply omits it from the child body when blank.

On `recommend`, `pdb_id` is optional and is passed to the source Structure
Report for RCSB metadata (`workflow.yaml:62-63`, `README.md:53-54`).

### Output-naming input

**Not present.** There is no run-name, output-property-name, result-name, or
label input. The complete property list is the six above
(`tool-definition.json:207-432`). The only "naming" fields anywhere are
`x-display-name` annotations, which name the *form controls*, not the results
(e.g. `"x-display-name": "PDB ID"` at `:228`).

### `additionalProperties: false`

**Yes, in six places** — an extra key is a hard schema rejection, not a warning:

| Location | Line |
|---|---|
| root of `inputs` | `tool-definition.json:5` |
| `pocket` | `:232` |
| `pocket.selections.items` | `:323` |
| `pocket.selections.items.author` | `:326` |
| `protein` | `:378` |
| `selection` | `:400` |

`selection.decisions` uses `additionalProperties` as a *value schema*
(`{"enum": ["keep","skip"]}`, `:407-417`) — arbitrary keys are allowed there,
but every value must be `keep` or `skip`.

**The one gap:** `pocket.crystal_ligand` (`:259-285`) does **not** set
`additionalProperties: false`, so extra keys pass schema there. Do not rely on
that.

Practical consequence for the UI: **send exactly the declared keys.** In
particular, do not send `sync`, `name`, `output_name`, `project_id`, or any
`recommend`-time `pocket`/`selection`/`model_missing_loops` scaffolding.

---

## 3. Output schema

### Verbatim `outputs`

Source: `tools/target-preparation/workflow/tool-definition.json:456-1058`.

```json
{
  "$schema": "http://json-schema.org/draft-07/schema",
  "description": "Union of Protein Prep, Structure Report, and Pocket Finder results published by child implementations under this Tool identity",
  "properties": {
    "audit_file_path": {
      "description": "UFA path to the transformation audit JSON (prepare)",
      "format": "file",
      "type": "string"
    },
    "extracted_ligands": {
      "description": "Kept ligands extracted as separate PDB files (prepare)",
      "items": {
        "additionalProperties": false,
        "properties": {
          "component_id": {
            "description": "Component id from the Selection",
            "type": "string"
          },
          "file_path": {
            "description": "UFA path to the extracted ligand PDB",
            "format": "file",
            "type": "string"
          }
        },
        "required": [
          "component_id",
          "file_path"
        ],
        "type": "object"
      },
      "type": "array",
      "x-data-type": "ExtractedLigand",
      "x-result-group": "extractedligands"
    },
    "pockets": {
      "description": "Detected binding pockets on the target protein",
      "items": {
        "additionalProperties": false,
        "properties": {
          "apolar_SASA": {
            "description": "Apolar Solvent Accessible Surface Area in square angstroms (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.apolar_SASA",
            "x-indexable": true
          },
          "box": {
            "additionalProperties": false,
            "description": "PCA-aligned docking box (OBB sizes + rotation for deeporigin.docking)",
            "properties": {
              "box_size_x": {
                "description": "PCA-aligned OBB X extent in angstroms",
                "type": "number"
              },
              "box_size_y": {
                "description": "PCA-aligned OBB Y extent in angstroms",
                "type": "number"
              },
              "box_size_z": {
                "description": "PCA-aligned OBB Z extent in angstroms",
                "type": "number"
              },
              "rotation_deg": {
                "description": "Euler [rx, ry, rz] degrees (Rz\u00b7Ry\u00b7Rx about pocket_center)",
                "items": {
                  "type": "number"
                },
                "maxItems": 3,
                "minItems": 3,
                "type": "array"
              }
            },
            "required": [
              "box_size_x",
              "box_size_y",
              "box_size_z",
              "rotation_deg"
            ],
            "type": "object",
            "x-data-type": "Pocket.box"
          },
          "box_size_x": {
            "description": "Lab-frame AABB docking box X dimension in angstroms (deprecated; prefer box)",
            "type": "number",
            "x-data-type": "Pocket.box_size_x"
          },
          "box_size_y": {
            "description": "Lab-frame AABB docking box Y dimension in angstroms (deprecated; prefer box)",
            "type": "number",
            "x-data-type": "Pocket.box_size_y"
          },
          "box_size_z": {
            "description": "Lab-frame AABB docking box Z dimension in angstroms (deprecated; prefer box)",
            "type": "number",
            "x-data-type": "Pocket.box_size_z"
          },
          "drugability_score": {
            "description": "Predicted drugability score indicating suitability as a drug target (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.drugability_score",
            "x-indexable": true
          },
          "file_path": {
            "description": "Path to the pocket PDB file on UFA",
            "type": "string",
            "x-data-type": "Pocket.file_path"
          },
          "hydrophobicity": {
            "description": "Hydrophobicity score of the pocket surface (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.hydrophobicity",
            "x-indexable": true
          },
          "pocket_center": {
            "description": "Center coordinates [x, y, z] of the pocket in angstroms",
            "items": {
              "type": "number"
            },
            "maxItems": 3,
            "minItems": 3,
            "type": "array",
            "x-data-type": "Pocket.center"
          },
          "pocket_count": {
            "description": "Number of pockets requested (1 for define-by-selection and from-crystal-ligand)",
            "type": "integer",
            "x-data-type": "Pocket.count"
          },
          "pocket_min_size": {
            "description": "Minimum pocket size (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.min_size"
          },
          "polar_SASA": {
            "description": "Polar Solvent Accessible Surface Area in square angstroms (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.polar_SASA",
            "x-indexable": true
          },
          "polar_apolar_SASA_ratio": {
            "description": "Ratio of polar to apolar Solvent Accessible Surface Area (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.polar_apolar_SASA_ratio",
            "x-indexable": true
          },
          "polarity": {
            "description": "Overall polarity measure of the pocket surface (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.polarity",
            "x-indexable": true
          },
          "protein_id": {
            "description": "ID of the parent protein this pocket belongs to",
            "type": [
              "string",
              "null"
            ],
            "x-data-type": "Protein.id",
            "x-key": true
          },
          "total_SASA": {
            "description": "Total Solvent Accessible Surface Area in square angstroms (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.total_SASA",
            "x-indexable": true
          },
          "volume": {
            "description": "Volume of the pocket in cubic angstroms (null for define-by-selection and from-crystal-ligand)",
            "type": [
              "number",
              "null"
            ],
            "x-data-type": "Pocket.volume",
            "x-indexable": true
          }
        },
        "required": [
          "protein_id",
          "file_path",
          "volume",
          "total_SASA",
          "polar_SASA",
          "apolar_SASA",
          "polar_apolar_SASA_ratio",
          "hydrophobicity",
          "drugability_score",
          "polarity",
          "pocket_center",
          "box",
          "box_size_x",
          "box_size_y",
          "box_size_z",
          "pocket_count",
          "pocket_min_size"
        ],
        "type": "object"
      },
      "type": "array",
      "x-data-type": "Pocket",
      "x-result-group": "pockets"
    },
    "protein": {
      "description": "Prepared protein artifacts (prepare)",
      "properties": {
        "force_field": {
          "description": "Protein force field used",
          "type": "string",
          "x-data-type": "PreparedProtein.force_field",
          "x-display-name": "Protein Force Field"
        },
        "pdb_id": {
          "description": "PDB ID supplied for loop modelling",
          "type": "string",
          "x-display-name": "PDB ID"
        },
        "ph": {
          "description": "pH used for protonation",
          "type": "number",
          "x-data-type": "PreparedProtein.ph",
          "x-display-name": "pH"
        },
        "protein_id": {
          "description": "ID of the protein (if provided on input)",
          "type": "string",
          "x-data-type": "Protein.id",
          "x-key": true
        },
        "protein_pdb_file_path": {
          "description": "Path to the cleaned, force-field-ready protein PDB file",
          "format": "file",
          "type": "string",
          "x-data-type": "PreparedProtein.protein_pdb_file_path"
        },
        "protonate_protein": {
          "description": "Whether the protein was protonated",
          "type": "boolean",
          "x-data-type": "PreparedProtein.protonate_protein",
          "x-display-name": "Protonate Protein"
        }
      },
      "required": [
        "protein_pdb_file_path"
      ],
      "type": "object",
      "x-data-type": "PreparedProtein",
      "x-result-group": "preparedproteins"
    },
    "recommendation": {
      "additionalProperties": false,
      "description": "Component inventory and recommendations (recommend)",
      "properties": {
        "analyzer_version": {
          "description": "Analyzer version that produced this recommendation",
          "type": "string"
        },
        "chain_id_mapping": {
          "additionalProperties": {
            "type": "string"
          },
          "description": "Author chain IDs remapped to single-character PDB chain IDs (mmCIF)",
          "type": "object"
        },
        "components": {
          "description": "Inventoried components with tri-state recommendations",
          "items": {
            "additionalProperties": false,
            "properties": {
              "author": {
                "additionalProperties": false,
                "properties": {
                  "chain_id": {
                    "type": "string"
                  },
                  "icode": {
                    "type": "string"
                  },
                  "label_asym_id": {
                    "type": "string"
                  },
                  "resname": {
                    "type": "string"
                  },
                  "resseq": {
                    "type": "integer"
                  }
                },
                "required": [
                  "chain_id"
                ],
                "type": "object"
              },
              "evidence": {
                "type": "object"
              },
              "id": {
                "type": "string"
              },
              "kind": {
                "enum": [
                  "chain",
                  "ligand",
                  "cofactor",
                  "water"
                ],
                "type": "string"
              },
              "label": {
                "type": "string"
              },
              "reason": {
                "type": "string"
              },
              "reason_code": {
                "type": "string"
              },
              "recommendation": {
                "enum": [
                  "keep",
                  "skip",
                  "review"
                ],
                "type": "string"
              },
              "subtype": {
                "type": "string"
              }
            },
            "required": [
              "id",
              "kind",
              "subtype",
              "label",
              "author",
              "recommendation",
              "reason_code",
              "reason"
            ],
            "type": "object"
          },
          "type": "array"
        },
        "source_sha256": {
          "description": "SHA-256 of the uploaded structure bytes",
          "type": "string"
        }
      },
      "required": [
        "source_sha256",
        "analyzer_version",
        "components",
        "chain_id_mapping"
      ],
      "type": "object"
    },
    "selection_file_path": {
      "description": "UFA path to the exact frozen Selection JSON (prepare)",
      "format": "file",
      "type": "string"
    },
    "structure_reports": {
      "description": "Structure Report rows",
      "items": {
        "additionalProperties": false,
        "properties": {
          "coverage": {
            "description": "Polymer residue coverage fraction (0-1)",
            "type": [
              "number",
              "null"
            ]
          },
          "coverage_score": {
            "description": "Coverage component score",
            "type": "number"
          },
          "field_status": {
            "additionalProperties": false,
            "description": "Per-field status: value, not_applicable, or unknown",
            "properties": {
              "coverage": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              },
              "inhibitor": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              },
              "method": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              },
              "organism": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              },
              "resolution": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              },
              "rfree": {
                "enum": [
                  "value",
                  "not_applicable",
                  "unknown"
                ],
                "type": "string"
              }
            },
            "required": [
              "organism",
              "method",
              "resolution",
              "coverage",
              "rfree",
              "inhibitor"
            ],
            "type": "object"
          },
          "grade": {
            "description": "Letter grade A-D",
            "enum": [
              "A",
              "B",
              "C",
              "D"
            ],
            "type": "string"
          },
          "has_ligand": {
            "description": "Whether the deposition has a non-water, non-ion ligand",
            "type": [
              "boolean",
              "null"
            ]
          },
          "inhibitor_score": {
            "description": "Inhibitor/holo component score",
            "type": "number"
          },
          "metadata_source": {
            "description": "Provenance of experimental metadata used for scoring",
            "enum": [
              "rcsb",
              "file_header",
              "file_header+rcsb"
            ],
            "type": "string"
          },
          "method": {
            "description": "Experimental method string",
            "type": [
              "string",
              "null"
            ]
          },
          "method_class": {
            "description": "Normalized method class",
            "type": [
              "string",
              "null"
            ]
          },
          "method_score": {
            "description": "Method component score",
            "type": "number"
          },
          "organism": {
            "description": "Source organism scientific name",
            "type": [
              "string",
              "null"
            ]
          },
          "organism_class": {
            "description": "Normalized organism class",
            "type": [
              "string",
              "null"
            ]
          },
          "organism_score": {
            "description": "Organism component score",
            "type": "number"
          },
          "pdb_id": {
            "description": "PDB ID used for experimental metadata",
            "type": [
              "string",
              "null"
            ]
          },
          "protein_id": {
            "description": "Platform protein entity id when provided",
            "type": [
              "string",
              "null"
            ],
            "x-data-type": "Protein.id",
            "x-key": true
          },
          "report_role": {
            "description": "Optional Target Preparation phase label (source or prepared)",
            "enum": [
              "source",
              "prepared"
            ],
            "type": "string"
          },
          "resolution": {
            "description": "Structure resolution in Angstroms",
            "type": [
              "number",
              "null"
            ]
          },
          "resolution_score": {
            "description": "Resolution component score",
            "type": "number"
          },
          "rfree": {
            "description": "Rfree (X-ray)",
            "type": [
              "number",
              "null"
            ]
          },
          "rfree_score": {
            "description": "Rfree component score",
            "type": "number"
          },
          "source_sha256": {
            "description": "SHA-256 of the uploaded structure bytes (omitted in pdb_id-only remote mode)",
            "type": "string"
          },
          "weighted_score": {
            "description": "Weighted Structure Report score",
            "type": "number"
          }
        },
        "required": [
          "metadata_source",
          "field_status",
          "resolution_score",
          "coverage_score",
          "rfree_score",
          "inhibitor_score",
          "method_score",
          "organism_score",
          "weighted_score",
          "grade"
        ],
        "type": "object"
      },
      "type": "array",
      "x-data-type": "StructureReport",
      "x-result-group": "structurereports"
    }
  },
  "type": "object"
}```

### There are no `jobOutputs` — read this first

The question assumes `jobOutputs`. For this tool that channel does not exist:

- `quote_target_preparation_method()` returns the literal `"workflow"` for every
  payload — `images/preflight/src/preflight_service/validate_target_preparation.py:10-13`
  (docstring: *"Target Preparation always runs as an Argo workflow."*).
  Confirmed by route tests asserting `payload["method"] == "workflow"` for both
  `recommend` and `prepare` —
  `images/preflight/tests/test_target_preparation_route.py:68` and `:100`.
- *"The executions API omits `jobOutputs` for workflow (Argo) mode by design"* —
  `images/test-tool/src/test_tool/main.py:1065-1066`, repeated at
  `images/test-tool/tests/test_user_logs.py:103-104`.

So the keys below are **output-schema keys / MQ result-group keys**, not
`jobOutputs` keys. Read them from the result tables named by `x-result-group`.

### Key names

| What | Output key | `x-data-type` | `x-result-group` | Line |
|---|---|---|---|---|
| Prepared structure | `protein` → `protein.protein_pdb_file_path` | `PreparedProtein` | `preparedproteins` | `:680-725`, path field `:706-711` |
| Pockets | `pockets` (array) → `pockets[].file_path` | `Pocket` | `pockets` | `:490-679` |
| Structure report | `structure_reports` (array) | `StructureReport` | `structurereports` | `:838-1055` |
| Extracted ligands | `extracted_ligands` (array) | `ExtractedLigand` | `extractedligands` | `:465-489` |
| Component inventory | `recommendation` (object) | **none** | **none** | `:726-832` |
| Frozen selection artifact | `selection_file_path` (string, UFA) | none | none | `:833-837` |
| Transformation audit | `audit_file_path` (string, UFA) | none | none | `:460-464` |

These match the standalone owners exactly, as composite tools are required to:
`tools/protein-prep/workflow/tool-definition.json:208` (`extractedligands`),
`:254` (`preparedproteins`); `tools/pocket-finder/workflow/tool-definition.json:562`
(`pockets`); `tools/structure-report/workflow/tool-definition.json:332`
(`structurereports`). Requirement: `.github/copilot-instructions.md:535`.

### Does it create a row in the `proteins` entity table?

**No. Nothing does.** This is the answer to your most important question, and it
is a clean negative:

1. **The prepared protein is routed to a dynamic results table, not an entity
   table.** The `protein` output declares
   `"x-data-type": "PreparedProtein"` / `"x-result-group": "preparedproteins"`
   — `tool-definition.json:723-724`. The repo's own routing table maps
   `x-result-group: proteins` to `x-data-type: Protein` (import) and
   `preparedproteins` to "Prepared ligand/protein (system-prep)" — two different
   rows — `.github/copilot-instructions.md:556-557`. The data-platform consumer
   "route[s] fixed-entity outputs (proteins/ligands) to their entity tables and
   everything else to the appropriate dynamic `results__*` table"
   — `images/import-dataset/src/import_dataset/main.py:571-573`.

2. **The repo forbids typing tool rows as `Protein`.** *"Array typed as `Protein`
   for tool rows → COPY to `proteins` fails"* —
   `.github/copilot-instructions.md:586` (regression table), and the authoring
   rule *"Do not reuse `Protein`, `Ligand`, or `Pose` unless rows are literally
   entity-table rows"* — `.github/copilot-instructions.md:510-513`.

3. **The only entity-creating call site in the repository:**

   ```python
   created = sdk.client.entities.create_protein(**create_kwargs)
   ```
   — `images/structure-report/src/pdb_import.py:259`
   (docstring at `:218`: *"Download or reuse cached coordinates, register a
   Protein, return jobOutputs"*). That is the `deeporigin.pdb-import` tool.
   A repo-wide grep for `create_protein` / `entities.create` returns only that
   line plus its tests (`images/structure-report/tests/test_pdb_import_server.py`)
   and the local test gateway stub (`tests/gateway/app.py:1020`).

4. **What target-preparation's children actually do:** they call
   `sdk.reporter.send_result(...)` and nothing else —
   `images/system-prep/src/prepare_protein_service.py:260` and `:468`,
   `images/structure-report/src/server.py:173`,
   `images/pocket-finder/src/core.py:495`, `:660`, `:817`. No entity API call
   appears in any of them.

5. The README states the design intent plainly: *"Child implementations remain
   the owners of their MQ result rows. This Tool declares the union of their
   output schemas and does not republish them."* —
   `tools/target-preparation/README.md:10-11`.

**Consequence:** after a successful `prepare`, the prepared structure exists as
(a) a UFA file at `protein.protein_pdb_file_path` and (b) a row in
`results__preparedproteins`. It does **not** appear in the `proteins` entity
table, so a protein-picker backed by that table will not offer it to Docking or
ABFE. If the feature's purpose is a prepared receptor that downstream tools can
select, **that wiring is missing and needs a backend change** (something has to
call `entities.create_protein`, as `pdb_import.py:259` does).

### `result_type` values

The declared, authoritative-in-repo values are the **`x-result-group`** strings
listed in the table above: `preparedproteins`, `pockets`, `structurereports`,
`extractedligands`.

A literal `result_type` field appears only in the **local test gateway
fixtures**, in singular form: `"result_type": "pocket"`
(`tests/gateway/state.py:139`), `"pose"` (`:178`), `"admetproperty"` (`:201`),
and a filter `{"column": "result_type", "op": "eq", "value": "pose"}`
(`images/system-prep/src/server.py:425`). **There is no `preparedprotein` or
`structurereport` literal anywhere in the repo.** I cannot tell you from these
files how the platform derives a singular `result_type` from a plural
`x-result-group`; that mapping lives outside this repository.

### Indexed structure-report rows

**Yes — the report is indexed, not only a job output.** `structure_reports`
carries `x-data-type: "StructureReport"` and `x-result-group: "structurereports"`
(`tool-definition.json:1053-1054`), and the producer publishes
`{"structure_reports": [row]}` via `send_result`
(`images/structure-report/src/server.py:171-175`).

Two rows are produced per `prepare` and one per `recommend`, distinguished by
`report_role`:

```json
"report_role": {
  "description": "Optional Target Preparation phase label (source or prepared)",
  "enum": ["source", "prepared"],
  "type": "string"
}
```
— `tool-definition.json:999-1006`

`recommend` → one `source` row (`workflow.yaml:58-59`).
`prepare` → one `prepared` row (`workflow.yaml:183-184`). Note `prepare` does
**not** re-emit a `source` row; only `recommend` does.
Tests: `tests/test_target_preparation.py:57`, `:74`;
`tests/test_target_preparation_composite.py:151-161`.

### The prepared-protein stamp

**Yes, PDB only.** SystemPrep writes it during `prepare`:

```python
stamp_prepared_protein_pdb(
    local_out,
    provenance=provenance_from_selection_and_applied(...),
)
```
— `images/system-prep/src/prepare_protein_service.py:408-421`, applied to
`protein_prepped.pdb` (`:401-407`).

The token constants:
`PREPARED_PROTEIN_STAMP_LINE = "REMARK  99 DO_PREPARED"` and
`PREPARED_PROTEIN_CIF_STAMP_KEY = "_deeporigin.prepared"` /
`PREPARED_PROTEIN_CIF_STAMP_VALUE = "DO_PREPARED"`
— `packages/toolbox-core/src/toolbox_core/prepared_protein_stamp.py:19-22`.

**Caveat you should surface to the backend team:** the module docstring says
*"Writers live outside this module for CIF (the Python client); toolbox detects
both forms and stamps PDB only"*
(`prepared_protein_stamp.py:4-6`). The prepare path always writes
`protein_prepped.pdb` (`prepare_protein_service.py:401-404`), so in practice the
target-preparation output is a stamped **PDB**. Downstream AUTO preprocess skips
cleanup on a stamped file — `CONTEXT.md:226-227`,
`docs/adr/0032-prepared-protein-stamp.md`.

---

## 4. The preflight route

File: `images/preflight/src/preflight_service/routes/target_preparation.py`
(80 lines, read in full).

### Method, path, request, response

- **Method / path:** `POST /target-preparation` —
  `images/preflight/src/preflight_service/app.py:230-236`. The path is selected
  by the annotation `deeporigin.io/serving-path: target-preparation` —
  `tools/target-preparation/workflow/preflight-service.yaml:7`.
- **Request:** the tool's own input payload, optionally wrapped —
  `payload = unwrap_nested_inputs_dict(body)`
  (`routes/target_preparation.py:54`). The route tests post the bare input
  object (`images/preflight/tests/test_target_preparation_route.py:59-65`).
- **Response:** always HTTP 200 with a billing-quote body.
  - `recommend` → `skip_billing_response(method="workflow")`, i.e.
    `{"method": "workflow", "counts": {"__SKIP": {"total": 1}}}` —
    `routes/target_preparation.py:46`; asserted at
    `test_target_preparation_route.py:68-69`.
  - `prepare` → `unit_billing_response(billing_code="DO_POCKET_FINDER",
    method="workflow", description="Target Preparation (Pocket Finder)")`, i.e.
    `{"method": "workflow", "counts": {"DO_POCKET_FINDER": {"total": 1}}}` —
    `routes/target_preparation.py:39-45`; asserted at
    `test_target_preparation_route.py:100-101`.
  - On validation failure the same 200 body gains
    `"validations": [<RFC 9457 problem>, ...]` —
    `routes/target_preparation.py:64-71`.

### Does it return a component inventory?

**No.** There is no inventory, no chain/ligand/cofactor/water list, and no
keep/review/skip recommendation anywhere in the route or its validator. The
entire response surface is `{method, counts, validations?}`. The route does not
download the structure at all — the validator's docstring says so explicitly:
*"Validate Target Preparation action branches **without UFA downloads**"* —
`images/preflight/src/preflight_service/validate_target_preparation.py:19`.

**The "call preflight instead of running a recommend execution" option is void.**
You cannot get an inventory from this endpoint.

### What it rejects, and the error shape

Every rejection is a `ValidationIssue(severity="error", check_name=..., detail=...)`
serialized to an RFC 9457 problem object via `validation_issue_to_problem`
(`routes/target_preparation.py:27-29`, `:69`), and also emitted as an ERROR
user log (`:65`). All from `validate_target_preparation.py`:

| `check_name` | Condition | `detail` | Line |
|---|---|---|---|
| `target_preparation.protein` | `protein` is not an object | `protein object with id and file_path is required` | `:27-34` |
| `target_preparation.protein.id` | `protein.id` missing/blank | `protein.id is required for Target Preparation` | `:36-43` |
| `target_preparation.protein.file_path` | `protein.file_path` missing/blank | `protein.file_path is required for Target Preparation` | `:44-51` |
| `target_preparation.model_missing_loops` | `prepare` and key absent | `model_missing_loops is required for action=prepare` | `:53-61` |
| `target_preparation.selection` | `prepare` and `selection` not an object | `selection is required for action=prepare` | `:62-69` |
| `target_preparation.pocket` | `prepare` and `pocket` not an object | `pocket configuration is required for action=prepare` | `:70-77` |
| `target_preparation.pdb_id` | `prepare` and `model_missing_loops is True` and `pdb_id` blank | `pdb_id is required when model_missing_loops is true` | `:78-87` |

All seven are trivially enforceable client-side — disable the submit button
rather than round-trip. Note the route deliberately does **not** duplicate JSON
Schema shape checks (`validate_target_preparation.py:21-23`), so the *schema*
rejections in §2 (unknown key, bad `pdb_id` pattern, mode-exclusivity, `review`
in decisions) surface from the tools-service instead, not from here.

Also note: **whitespace-only strings pass schema but fail preflight.**
`str(protein.get("id") or "").strip()` (`:36`) rejects `"  "`, while the schema's
bare `"type": "string"` accepts it. Trim client-side.

### Authentication and reachability

- Auth is by forwarded platform headers, read off the incoming request:
  `ToolsSDK.from_headers(dict(request.headers.items()))` —
  `routes/target_preparation.py:55`. Those are the `X-Do-*` headers
  (`X-Do-Auth-Token`, `X-Do-Org-Key`, `X-Do-Execution-Id`, …) enumerated at
  `images/served-tool-requester/call_served_tool.py:27-45`. This is **not** the
  tools API's user-facing auth.
- **It is cluster-internal.** The manifest is a Knative Service in namespace
  `deeporigin-executions` with `deeporigin.io/serving-role: preflight` —
  `tools/target-preparation/workflow/preflight-service.yaml:2-8`. `CONTEXT.md:1180-1192`
  describes it as the platform's create-time gate: *"Runs Validator checks and
  returns Estimator billing shape **before Argo scheduling**"* and *"Blocking
  failures … fail the execution at create time."* It is invoked by the platform
  harness during `executions.create`, not by a client.
- The FastAPI app does set permissive CORS (`allow_origins=["*"]`,
  `allow_credentials=False`) — `images/preflight/src/preflight_service/app.py:53-59`
  — but that is a property of the shared image, not evidence of an ingress. There
  is no public route, no gateway mapping, and no browser-facing hostname for it
  anywhere in this repo.

**Explicit verdict: the browser cannot call this endpoint. Design as if it does
not exist.** You will see its `validations[]` only as an execution-create error.

---

## 5. Workflow behaviour

Source: `tools/target-preparation/workflow/workflow.yaml` (772 lines, read in
full) and `tests/test_target_preparation.py`.

### Ordered steps

Entry DAG (`workflow.yaml:32-40`) branches on `action`.

**`action: recommend`** — `workflow.yaml:42-92`

| # | Task | Template | Image / tool | Line |
|---|---|---|---|---|
| 1 | `log-start` | `user-log` | `served-tool-requester:0.1.6` | `:45-52` |
| 2 | `source-report` | `call-structure-report` | `structure-report` serving `3.0.4`, `POST /structure-report`, `report_role: source`, progress suppressed | `:53-65`, template `:262-355` |
| 3 | `protein-recommend` | `protein-prep-served` | `sysprep` serving `2.0.10`, `POST /prepare-protein` with `action: recommend`, progress **not** suppressed (owns the 100% `__SKIP`) | `:66-74`, template `:357-478` |
| 4 | `fail-recommend` / `log-done` | `fail-child` / `user-log` | conditional on outcome | `:75-92` |

**`action: prepare`** — `workflow.yaml:94-208`

| # | Task | Template | Image / tool | Line |
|---|---|---|---|---|
| 1 | `log-start` | `user-log` | requester | `:97-104` |
| 2 | `try-prepare` | `protein-prep-served` | `sysprep` serving `2.0.10`, `action: prepare`, suppressed | `:105-113` |
| 3 | `heavy-prepare` | `protein-prep-heavy` | `functions/sysprep:2.0.10` pod, 8 CPU / 32 GiB, `python3 -m src.prepare_protein_main` — **only** when `try-prepare` timed out *and* loops are on | `:114-119`, template `:480-528` |
| 4 | `normalize-knative` / `normalize-heavy` → `resolve-prepared` | `normalize-protein-prep-outputs`, `resolve-prepare-artifacts` | requester (glue) | `:138-177`, templates `:530-653` |
| 5 | `prepared-report` | `call-structure-report` | `structure-report` `3.0.4`, `report_role: prepared`, suppressed | `:178-190` |
| 6 | `pocket-finder` | `call-pocket-finder` | `pocket-finder` serving `3.1.2`, `POST /find_pockets` (owns the 100% `DO_POCKET_FINDER`) | `:191-199`, template `:655-772` |
| 7 | `log-done` | `user-log` | requester | `:200-208` |

Image pins: `workflow.yaml:22-29`. Ordering asserted by
`tests/test_target_preparation.py:48-77`.

### Does it re-analyse components, and does the caller's selection survive?

**It re-analyses — but the caller's decisions are authoritative. Your filtering
panel is not decorative.**

The prepare service re-runs the analyzer on the freshly downloaded file:

```python
digest = sha256_file(Path(raw_local))
recommendation = recommend_protein_components(input_path=Path(raw_local),
                                              source_sha256=digest)
...
applied = apply_selection(input_path=pre_repaired, selection=selection,
                          recommendation=recommendation, ...,
                          source_sha256=digest)
```
— `images/system-prep/src/prepare_protein_service.py:300-324`

`apply_selection` uses the re-analysis **only** to (a) validate the digest and
analyzer version, (b) check the decisions map is exhaustive and legal, and
(c) resolve each `component_id` to its author identity. The keep/skip verdict
comes verbatim from `selection["decisions"]`:

```python
decisions: dict[str, str] = selection["decisions"]
by_id = {component.id: component for component in recommendation.components}
for component_id, decision in decisions.items():
    component = by_id[component_id]
    ...
    if decision == RecommendationDecision.SKIP.value:
        audit_skipped.append(entry); continue
```
— `packages/toolbox-core/src/toolbox_core/protein_component.py:1031-1063`

The recomputed `recommendation` field is **never** consulted in the branch that
decides keep vs skip. `CONTEXT.md:199-201` reinforces this as a design rule for
Protein Prep: *"Avoid … recomputing AUTO preprocess choices during prepare."*

**But three re-analysis-driven hard failures will hit your UI:**
- Digest drift (the file changed between recommend and prepare) →
  `protein_component.py:821-826`.
- Analyzer version drift (the serving was upgraded between your two calls) →
  `protein_component.py:828-838`. `ANALYZER_VERSION` is `"1.1.0"`
  (`protein_component.py:37`) — a minor bump of the sysprep image invalidates
  every in-flight selection.
- A decisions map that is not exactly the analyzed component set →
  `protein_component.py:844-857`.

All three surface as `ProteinComponentError` → HTTP 400 from the child
(`prepare_protein_service.py:325-326`), which the workflow classifies as
`error` → `fail-child` → **FAILED**.

Two more selection rules that must be enforced in the panel:
- Keeping a chain whose subtype is `nucleic`, `mixed`, or `other` →
  `"Cannot keep unsupported polymer chain {id} (subtype=...)"` —
  `protein_component.py:878-886`.
- Keeping zero protein chains →
  `"selection must keep at least one supported protein chain"` —
  `protein_component.py:892-893`.

### Are steps conditional on the toggles?

Only one step is conditional on an input, and it is not a toggle in the UI
sense:

- `heavy-prepare` runs only when the served prepare returned a 504 timeout
  **and** `model_missing_loops != 'false'` — `workflow.yaml:117-119`. Test:
  `tests/test_target_preparation.py:79-88`.
- `fail-timeout-no-heavy` is the mirror: a timeout with loops off is a hard
  failure with no fallback — `workflow.yaml:128-137`.
- `prepared-report` (`:178-190`) and `pocket-finder` (`:191-199`) have **no
  `when:`** — they always run on `prepare`.
- The pocket *mode* branch is not a DAG branch; it is `if/elif/else` inside the
  single `call-pocket-finder` script — `workflow.yaml:698-737`.

### Partial failure after the protein was prepared

**Execution status: FAILED. The prepared protein rows are still written.**

- The pocket step ends with
  `send-user-log --level error --message "TARGET_PREP_CHILD_FAILED: Pocket Finder failed"; exit 1`
  — `workflow.yaml:760-763`. The `prepared-report` step does the same
  (`:343-346`). A non-zero exit fails the task, which fails the `prepare-flow`
  DAG and the workflow.
- The prepared-protein result row was already published by the sysprep child via
  `sdk.reporter.send_result(result_data)`
  (`images/system-prep/src/prepare_protein_service.py:468`) *before* the
  workflow moved on, and the stamped PDB was already uploaded
  (`prepare_protein_service.py:422`). Nothing rolls it back.
- The README states the intended contract: *"ChildError[Child application
  failure] --> Failed[FAILED / ERROR user log / **completed artifacts
  retained**]"* — `tools/target-preparation/README.md:40`.
- **Zero pockets is *not* a failure.** README: *"Empty -->|no| NotReady[COMPLETED
  / ERROR user log / prepared artifacts retained]"* — `README.md:38`. This is
  consistent with `call_served_tool`'s success rule (2xx + parseable JSON +
  no top-level `error`/`detail` —
  `images/served-tool-requester/call_served_tool.py:76-86`); an empty `pockets`
  array is a success. I found **no explicit empty-pocket branch in
  `workflow.yaml`** — the ERROR user log the README promises would have to come
  from the pocket-finder image, and I could not locate it there.

Retries: every child template uses `retryPolicy: OnError` (transient pod/node
loss only, never an application failure) — `workflow.yaml:232-237`, `:350-355`,
`:473-478`, `:523-528`, `:767-772`; asserted by
`tests/test_target_preparation.py:114-126`.
Wall clock: `activeDeadlineSeconds: 18000` (5 h) — `workflow.yaml:7`; each child
HTTP call has `--timeout 660` (11 min) — `:339`, `:438`, `:757`.

### Synchronous or asynchronous? Is there a top-level `sync`?

**Asynchronous, always.** `method: "workflow"` for every payload —
`validate_target_preparation.py:10-13`.

**No, it does not accept a top-level `sync` key.** There is no `sync` property
in the input schema, and the root sets `additionalProperties: false`
(`tool-definition.json:5`), so sending one is a hard rejection. Contrast
`deeporigin.pocket-finder`, which *does* declare `sync`
(`tools/pocket-finder/workflow/tool-definition.json:344-347`).

The string `"sync": True` you will see at `workflow.yaml:696` is the *internal*
body the workflow sends to the pocket-finder serving — not your request.

### Billing

**One line item, always.** Both branches quote a single code with `total: 1`:

- `recommend` → `__SKIP`, total 1 —
  `routes/target_preparation.py:46`;
  `test_target_preparation_route.py:69` asserts
  `counts == {BILLING_SKIP_CODE: {"total": 1}}`.
- `prepare` → `DO_POCKET_FINDER`, total 1 —
  `routes/target_preparation.py:39-45`;
  `test_target_preparation_route.py:101` asserts
  `counts["DO_POCKET_FINDER"]["total"] == 1`.

Tool-level `billingCode` is `DO_POCKET_FINDER` (`tool-definition.json:2`).
The children never double-bill: their 100% progress is suppressed via
`X-Do-Suppress-Final-Progress: true`
(`workflow.yaml:64-65`, `:112-113`, `:189-190`;
`images/served-tool-requester/call_served_tool.py:71-72`), so exactly one child
publishes the final 100% — Protein Prep recommend (`__SKIP`) on the recommend
branch, Pocket Finder (`DO_POCKET_FINDER`) on prepare
(`tools/target-preparation/README.md:6-8`). Verified by
`tests/test_target_preparation_composite.py:73-108`.

**So a full user journey costs two executions but only one billed unit:** the
recommend execution is `__SKIP`.

---

## 6. Component identity for 3D rendering

### What identifies a component

**Both — but in different places, and that asymmetry matters.**

- In the **output** `recommendation.components[]`, every component carries an
  opaque `id` **and** a structured `author` block — `tool-definition.json:742-817`
  (both `id` and `author` are in `required`, `:806-815`).
- In the **input** `selection.decisions`, a component is identified by the
  **opaque id string only** — it is a flat `{id: "keep"|"skip"}` map with no
  author, no kind — `tool-definition.json:407-417`.

So the UI must hold the recommend response in memory (or re-fetch it) to render;
it sends back only strings.

### Exact grammar

Two builders, both in
`packages/toolbox-core/src/toolbox_core/protein_component.py`:

```python
def component_id_for_chain(chain_id: str) -> str:
    """Return the opaque transport id for a polymer/non-polymer chain."""
    return f"chain:{chain_id}"                                     # line 225-227

def component_id_for_residue(*, kind, chain_id, resname, resseq, icode="") -> str:
    """Return the opaque transport id for a residue instance."""
    return f"{kind.value}:{chain_id}:{resname}:{resseq}:{icode}"   # line 230-239
```

| kind | grammar | example |
|---|---|---|
| `chain` | `chain:<auth_chain_id>` | `chain:A` |
| `ligand` | `ligand:<chain>:<resname>:<resseq>:<icode>` | `ligand:A:IBP:100:` |
| `cofactor` | `cofactor:<chain>:<resname>:<resseq>:<icode>` | `cofactor:A:MN:501:` |
| `water` | `water:<chain>:<resname>:<resseq>:<icode>` | `water:A:HOH:310:` |

`icode` is `""` for the overwhelmingly common no-insertion-code case, which is
why non-chain ids end in a trailing colon. Call sites confirming each kind:
water `:622-628`, cofactor/metal `:662-668`, cofactor/organic `:718-724`,
cofactor/artifact `:740-746`, ligand `:762-768`, chain `:424`, `:449`, `:537`,
`:555`.

### Correcting the premise in the brief

**`deeporigin.protein-prep` does *not* use different field orders per kind.**
There is one residue builder shared by ligand, cofactor, and water
(`protein_component.py:230-239`), so the order is always
`kind : chain : resname : resseq : icode`.

- Your `water:A:HOH:310:` example **is** correct — and it is
  `chain:resname:resseq:icode`, not a special case.
- Your `ligand:LIG:A:100` example (resname:chain:resseq) is **not** a shape this
  code can emit. The correct form is `ligand:A:LIG:100:`.

**`deeporigin.target-preparation` uses exactly the same grammar** — its
`selection.decisions` keys are validated against
`{component.id for component in recommendation.components}`
(`protein_component.py:844-848`), i.e. the ids this builder produced. There is
no target-preparation-specific grammar.

⚠️ **The repo's own test fixtures use a wrong shape.**
`"ligand:A:120:IBP"` appears at `tests/test_target_preparation_schema.py:103`,
`:112`, `:123`, `:134` and `tests/test_target_preparation.py:163`, `:190-191`.
That is `kind:chain:resseq:resname` — resseq and resname swapped. Those tests
pass because the strings are never round-tripped through the analyzer (they are
compared as opaque keys). **Do not copy component ids out of these tests.**
Build them from the analyzer, or better, echo back exactly the `id` strings the
recommend response gave you.

### What a cofactor component id looks like

I searched the whole repository for a literal `cofactor:` id and found **none** —
no fixture, no test, no doc contains one. So here is the derivation from the
producing code, which is authoritative.

A metal cofactor is created at `protein_component.py:660-712`:

```python
cid = component_id_for_residue(
    kind=ComponentKindV2.COFACTOR,
    chain_id=chain.id, resname=resname, resseq=resseq, icode=icode,
)
```

`resname` is the uppercased, stripped PDB residue name
(`protein_component.py:600`). For your mockup's two metals — Mn²⁺ and Zn²⁺ — the
PDB residue names are `MN` and `ZN`, both members of `METAL_ELEMENTS`
(`packages/toolbox-core/src/toolbox_core/protein_preprocess.py:87`, `:92`).

```
cofactor:A:MN:501:      # Mn2+, auth chain A, resseq 501, no insertion code
cofactor:A:ZN:302:      # Zn2+, auth chain A, resseq 302, no insertion code
cofactor:B:HEM:401:     # heme (ORGANIC_COFACTORS, protein_component.py:42)
cofactor:A:FAD:600:     # FAD  (ORGANIC_COFACTORS, protein_component.py:49)
```

Which HETATM records become `cofactor` rather than `ligand`:
- any residue whose element is in `METAL_ELEMENTS`
  (`protein_preprocess.py:77-104`, includes `MN`, `ZN`, `MG`, `CA`, `FE`, `CU`,
  `NI`, `CO`, `CD`, `HG`, …) — `protein_component.py:660-712`;
- any resname in `ORGANIC_COFACTORS` = `HEM HEA HEB HEC NAD NAP NDP FAD FMN PLP
  COA SAM SAH TPP THF FES SF4 F3S` — `protein_component.py:40-61`, handled at
  `:717-735`;
- any resname in `ARTIFACT_RESNAMES` (`SO4 PO4 GOL EDO PEG … CL BR NA K …`) —
  `protein_component.py:64-96`, handled at `:737-758`, always `skip`.

### Is there an `author` block alongside each component?

**Yes — and it is exactly what a renderer needs.**

```json
"author": {
  "additionalProperties": false,
  "properties": {
    "chain_id":      { "type": "string"  },
    "icode":         { "type": "string"  },
    "label_asym_id": { "type": "string"  },
    "resname":       { "type": "string"  },
    "resseq":        { "type": "integer" }
  },
  "required": ["chain_id"],
  "type": "object"
}
```
— `tool-definition.json:746-769`

Only `chain_id` is required by the schema; the serializer omits `resname`,
`resseq`, `icode`, and `label_asym_id` when they are `None`/empty
(`protein_component.py:243-255`). In practice a chain component has only
`chain_id`, and every residue component has `chain_id` + `resname` + `resseq`
(the analyzer always sets all three — `protein_component.py:603-608`), with
`icode` present only when non-blank.

**Do not parse the ids.** Use `author` for Mol* addressing and treat `id` as an
opaque key you echo back in `decisions`.

Two more renderer-relevant fields on each component:
`label` (a human string like `"MN A 501"` / `"chain A"`,
`protein_component.py:674`, `:539`) and `evidence` (a free-form object, e.g.
`{"coordinating_residues": 4}` — `tool-definition.json:770-772`,
`protein_component.py:681`).

And one at the top of `recommendation`: `chain_id_mapping`, *"Author chain IDs
remapped to single-character PDB chain IDs (mmCIF)"* —
`tool-definition.json:734-740`, required (`:825-830`). **For an mmCIF input the
`author.chain_id` in the inventory may not be the chain id in the original
file.** Apply this map before addressing the viewer.

---

## 7. Enum catalogues — closed or open?

"CLOSED" below means the *tool definition* pins the value set with a JSON Schema
`enum`. "OPEN" means the schema declares a bare `string` and the value set lives
only in producer code, where it can grow without a schema change.

| Value set | Verdict | Definition (not an example) |
|---|---|---|
| component `kind` | **CLOSED** — `chain`, `ligand`, `cofactor`, `water` | `tool-definition.json:776-784`; matches `ComponentKindV2` (`protein_component.py:146-152`) |
| component `subtype` | **OPEN** — `{"type": "string"}`, no enum | `tool-definition.json:802-804` |
| component `recommendation` | **CLOSED** — `keep`, `skip`, `review` | `tool-definition.json:794-801`; matches `RecommendationDecision` (`protein_component.py:155-160`) |
| `selection.decisions` values | **CLOSED** — `keep`, `skip` only (no `review`) | `tool-definition.json:407-417`; runtime rejection at `protein_component.py:859-867` |
| component `reason_code` | **OPEN** — `{"type": "string"}`, no enum | `tool-definition.json:791-793` |
| report `grade` | **CLOSED** — `A`, `B`, `C`, `D` | `tool-definition.json:917-926` |
| report `metadata_source` | **CLOSED** — `rcsb`, `file_header`, `file_header+rcsb` | `tool-definition.json:938-946` |
| report `field_status.*` | **CLOSED** — `value`, `not_applicable`, `unknown`, on all six of `organism`/`method`/`resolution`/`coverage`/`rfree`/`inhibitor` | `tool-definition.json:854-915` |
| report `method_class` | **OPEN in the schema** — `["string","null"]` | `tool-definition.json:954-960` |
| report `organism_class` | **OPEN in the schema** — `["string","null"]` | `tool-definition.json:972-978` |
| report `report_role` | **CLOSED** — `source`, `prepared` | `tool-definition.json:999-1006` |
| input `action` | **CLOSED** — `recommend`, `prepare` | `tool-definition.json:208-217` |
| input `pocket.mode` | **CLOSED** — `auto-find`, `define-by-selection`, `from-crystal-ligand` | `tool-definition.json:286-297` |
| input `pocket.box_geometry` | **CLOSED** — `ligand-extents`, `fixed-radius` | `tool-definition.json:242-251` |
| input `pocket.selections[].kind` | **CLOSED** — `residue`, `ligand`, `cofactor` | `tool-definition.json:351-360` |

### On `subtype`

**do-dd-client accepting any string is correct, and it matches the contract.**
The schema declares a bare `string` (`tool-definition.json:802-804`) while the
producer emits ten values from inline literals with no enum type behind them:

| kind | subtype values | source |
|---|---|---|
| `chain` | `protein`, `nucleic`, `mixed`, `other` | `_chain_polymer_subtype` — `protein_component.py:380-390`; used at `:506`, `:539` |
| `chain` | `short_peptide` | `protein_component.py:557` |
| `water` | `coordinating`, `crystal` | `protein_component.py:634`, `:647` |
| `cofactor` | `metal`, `ion`, `organic`, `other` | `protein_component.py:674`, `:690`, `:703`, `:729`, `:751` |
| `ligand` | `organic` | `protein_component.py:774`, `:790` |

I cannot tell from these files whether "any string" was a deliberate client
decision or a coincidence, but it is the right behaviour: render a fallback for
unknown subtypes. Note that `subtype` is **semantically load-bearing** despite
being open — keeping a chain with subtype `nucleic`/`mixed`/`other` is a hard
error (`protein_component.py:878-886`), so the UI must special-case at least
those three by string comparison.

### On `reason_code`

Also OPEN in the schema (`tool-definition.json:791-793`). Twelve values exist in
the producer today, all inline literals with no enum:

`keep_chain` (`:430`), `duplicate_chain` (`:455`), `non_protein_chain` (`:543`),
`short_peptide` (`:561`), `coordinating_water` (`:638`), `water_review` (`:651`),
`coordinated_metal` (`:678`), `crystallization_artifact` (`:694`, `:755`),
`under_coordinated_metal` (`:707`), `organic_cofactor` (`:733`),
`recognized_ligand` (`:778`), `ambiguous_ligand` (`:794`) — all in
`packages/toolbox-core/src/toolbox_core/protein_component.py`.

Switch on these for iconography if you like, but always fall back to the
free-text `reason` field (required, `tool-definition.json:788-790`), which the
producer always populates with a human sentence.

### On `method_class` / `organism_class`

The schema leaves them open, but the producer *is* closed today:

```python
FieldStatus   = Literal["value", "not_applicable", "unknown"]
MethodClass   = Literal["cryo-em", "x-ray", "nmr", "other", "unknown"]
OrganismClass = Literal["human", "mammal", "vertebrate", "other", "unknown"]
```
— `images/structure-report/src/scoring.py:7-9`

with one twist: the row writer converts `"unknown"` to `null` before emitting —
`"organism_class": organism if organism != "unknown" else None` and the same for
method (`images/structure-report/src/scoring.py:290-292`). **So over the wire the
value set is `human|mammal|vertebrate|other|null` and
`cryo-em|x-ray|nmr|other|null`.** Because the schema does not pin them, treat as
OPEN: accept any string, and handle `null` explicitly.

`metadata_source` is `MetadataSource = str  # rcsb | file_header | file_header+rcsb`
in code (`images/structure-report/src/server.py:46`) — a comment, not a type —
but the *schema* pins it to those three (`tool-definition.json:938-946`), so it
is CLOSED for your purposes.

### Other fixed sets worth knowing

- `pdb_id` is not an enum but is pattern-constrained: `^[A-Za-z0-9]{4}$`
  (`tool-definition.json:226`).
- `pockets[].box.rotation_deg` is a fixed-length array: `minItems: 3`,
  `maxItems: 3` (`tool-definition.json:520-529`); same for
  `pockets[].pocket_center` (`:583-592`).
- Analyzer version is a single pinned constant, not an enum:
  `ANALYZER_VERSION = "1.1.0"` (`protein_component.py:37`).

---

## 8. Contradictions with the intended UI

Restating your plan: *user picks a protein from a table → a structure report and
a component inventory appear automatically → user toggles Keep/Skip per
component → a 3D viewer colours kept vs excluded → user sets four booleans and a
run name → one button submits ONE `deeporigin.target-preparation` execution.*

Here is everything in the tool that breaks that.

**1. "One button, one execution" is impossible. It is two executions,
sequentially dependent.**
`action: "recommend"` and `action: "prepare"` are mutually exclusive branches
(`tool-definition.json:6-42`, `:43-204`; `workflow.yaml:32-40`), and `prepare`
cannot be constructed without three values that only `recommend` produces:
`selection.source_sha256`, `selection.analyzer_version`, and the exact component
id set. Budget for a two-phase flow with a visible wait between them.

**2. The inventory may be unreadable from a target-preparation execution.**
Covered in §3 and in the lead bullets: `recommendation` has no `x-result-group`
(`tool-definition.json:726-832`), and workflow-mode executions have no
`jobOutputs` (`images/test-tool/src/test_tool/main.py:1065-1066`). Meanwhile
`deeporigin.protein-prep` `recommend` quotes `"direct"`
(`validate_protein_prep.py:141-143`) and therefore returns a synchronous body.
**Recommendation: get the inventory from `deeporigin.protein-prep`
`action: "recommend"` (fast, synchronous, `__SKIP`-billed), and use
`deeporigin.target-preparation` only for `action: "prepare"`.** Confirm with
platform before committing — this is the one claim that depends on runtime
behaviour I cannot read from these files.

**3. Three of the four booleans do not exist; the fourth drags a mandatory text
field with it.** Only `model_missing_loops` exists (§2). And because it defaults
to *on* in every other surface, your form's default state makes `pdb_id`
**required** (`tool-definition.json:64-81`,
`validate_target_preparation.py:78-87`), constrained to `^[A-Za-z0-9]{4}$`. A
user-uploaded structure with no PDB entry therefore **cannot** run with loop
modelling — they must switch it off. Make that trade explicit in the UI.

**4. There is no run name, and inventing one is a hard failure.** No
output-naming input exists, and `additionalProperties: false` at the root
(`tool-definition.json:5`) turns an extra `name` key into a schema rejection at
execution-create, not a warning.

**5. The Keep/Skip panel must be exhaustive, including every water.** A partial
decisions map is rejected outright (`protein_component.py:844-857`), and every
water in the file is inventoried as its own component
(`protein_component.py:622-656`) with `recommendation: "review"` by default for
non-coordinating waters (`:647-655`). For a typical crystal structure that is
hundreds of rows the user must resolve to `keep` or `skip` before the button can
enable — `review` is rejected (`protein_component.py:859-867`). You need
bulk "skip all crystal waters" affordances, not a per-row toggle list.

**6. "Kept" does not mean "in the prepared structure" for ligands.** A kept
ligand is *extracted* into a separate PDB and removed from the receptor
(`protein_component.py:1056-1063`, reason text at `:776-780`), and surfaces as
an `extracted_ligands[]` row (`tool-definition.json:465-489`). Colour ligands as
a third state ("extract"), not as "keep".

**7. Some Keep choices are illegal and must be blocked client-side.** Keeping a
`nucleic`/`mixed`/`other` chain, or keeping zero protein chains, are 400s
(`protein_component.py:878-893`). The inventory already tells you: those chains
arrive with `recommendation: "review"` and `reason_code: "non_protein_chain"`
(`protein_component.py:534-550`).

**8. `prepare` needs a whole pocket configuration the mockup has no room for.**
`pocket` is required (`tool-definition.json:56-62`) with three modes and
mutually exclusive required fields per mode (§2). The cheapest viable default is
`{"mode": "auto-find", "pocket_count": N, "pocket_min_size": M}` — but
`pocket_count` and `pocket_min_size` have **no schema defaults**
(`tool-definition.json:298-311`), so the UI must supply numbers.

**9. Your selection can go stale between the two executions.** A sysprep serving
upgrade that bumps `ANALYZER_VERSION` (`protein_component.py:37`) invalidates
every open selection with
`"selection.analyzer_version does not match the recommendation analyzer"`
(`protein_component.py:828-838`). Handle that error by re-running recommend and
diffing, not by retrying.

**10. The prepared protein does not land in the `proteins` table** — §3. If the
user's mental model is "prepare it, then dock it", that link does not exist yet.

**11. `prepare` fails hard if pockets fail, after the protein is already
prepared** — §5. The execution reads FAILED while a perfectly good prepared
structure sits in `results__preparedproteins`. Do not treat FAILED as "nothing
happened"; surface the retained artifacts.

**12. The preflight route is not a browser-callable inventory service** — §4.

**13. `mcpFeatured: true`** (`tool-definition.json:450`) means this tool is also
exposed through the MCP surface. Any change to the two-step contract has a second
consumer.

### What is fine

- One protein picked from a table maps cleanly to
  `protein: {id, file_path}` — both required, both available from an entity row
  (`tool-definition.json:377-398`).
- The structure report really does appear automatically: `recommend` runs it
  first, unprompted, with `report_role: "source"`
  (`workflow.yaml:53-65`), and `prepare` adds a `"prepared"` row
  (`workflow.yaml:178-190`). No extra call needed. Note these are *indexed*
  rows, so they are readable even in workflow mode.
- The `author` block gives the viewer everything it needs — chain id, resname,
  resseq, icode (`tool-definition.json:746-769`) — with no id parsing, provided
  you apply `chain_id_mapping` for mmCIF (`:734-740`).
- Billing is one line item per execution, and the inventory step is free (§5).

---

## What I could not determine

| Question | Why not | What would settle it |
|---|---|---|
| Whether the tool is actually registered/enabled on dev, staging, prod | No per-environment declaration exists in the repo; `tests/test_target_preparation.py:231` says publication was still pending | Run `.github/workflows/dump-enabled-tools.yml` for each env, or ask platform |
| The literal `result_type` string the platform derives from each `x-result-group` | Only the local test gateway has literals (`tests/gateway/state.py:139`, `:178`, `:201`), and none of them cover prepared proteins or structure reports | Query a real completed execution's results, or read the data-platform MQ consumer (outside this repo) |
| Whether a workflow-mode `send_result` payload with **no** `x-result-group` (i.e. `recommendation`) is stored anywhere at all | The repo documents the failure mode (*"Missing `x-result-group` → No dynamic result table"*, `.github/copilot-instructions.md:589`) but not the MQ consumer's actual behaviour | Run one `recommend` execution and inspect what the API returns |
| Whether the mmCIF stamp `_deeporigin.prepared` is ever written on this path | `prepared_protein_stamp.py:4-6` says toolbox "stamps PDB only" and the prepare path writes `protein_prepped.pdb` (`prepare_protein_service.py:401-404`); the CIF writer is "the Python client", which is not in this repo | Check `do-dd-client` for the CIF stamp writer |
| Whether the README's "zero pockets → COMPLETED + ERROR user log" branch is implemented | `README.md:36-38` asserts it; I found no such branch in `workflow.yaml` and could not locate it in `images/pocket-finder/src/core.py` | Ask the pocket-finder owner, or run an auto-find that finds nothing |
| A real, analyzer-produced cofactor component id from a live structure | No fixture in this repo contains one (grep for `cofactor:` returns zero hits) | Run `deeporigin.protein-prep` `recommend` on a metalloprotein (e.g. a Zn/Mn structure) and read `recommendation.components[]`. The derivation in §6 is from the id builder itself and should hold. |
| Whether do-dd-client's "accept any `subtype` string" is a deliberate decision | That client is not in this repository | Ask its owner — though it is the correct behaviour either way (§7) |
| How the UI is expected to obtain `protein.id` for an uploaded (non-PDB) file | Only `deeporigin.pdb-import` creates Protein entities (`images/structure-report/src/pdb_import.py:259`); no upload path appears here | Ask platform for the upload/registration flow |
