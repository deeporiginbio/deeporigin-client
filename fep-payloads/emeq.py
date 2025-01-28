payload = json.dumps(
    {
        "toolId": "deeporigin/md-suite-emeq",
        "inputs": {
            "input": {
                "columnId": "_column:5rAPtk2oPl9psZZU3IeS0",
                "rowId": "fma-1",
                "databaseId": "fep-mason-all",
            },
            "force": 1,
            "test_run": 0,
            "system": "complex",
            "run_name": "test-run",
            "threads": 0,
            "em_solvent": True,
            "em_all": True,
            "nvt_heating_ns": 0.1,
            "npt_reduce_restraints_ns": 0.2,
            "from_run": "__USE_SYSTEM",
            "emeq_md_options": {
                "Δt": 0.004,
                "T": 298.15,
                "cutoff": 0.9,
                "fourier_spacing": 0.12,
                "hydrogen_mass": 2,
            },
        },
        "outputs": {
            "output_file": {
                "columnId": "_column:9YNlqWoctV1i3I04o3HMM",
                "rowId": "fma-1",
                "databaseId": "fep-mason-all",
            }
        },
        "clusterId": "d638d3cd-833f-476b-87fe-ca6d76e66049",
    }
)
