import warnings

import cement

from ..feature_flags import (
    FeatureNotAvailableWarning,
)
from ..feature_flags import (
    get_value as get_feature_flags,
)
from ..variables import (
    VariableStatus,
    VariableType,
    disable_variable_auto_updating,
    enable_variable_auto_updating,
    get_variable_types_by_values,
    install_variables,
    uninstall_variables,
)
from ..variables.types import EnvironmentVariable

__all__ = [
    "CONTROLLERS",
    "VariablesController",
    "InstallVariablesController",
    "AutoInstallVariablesController",
    "UninstallVariablesController",
]


class VariablesController(cement.Controller):
    """controller to retrieve and install variables and secrets from the Deep Origin platform"""

    class Meta:
        label = "variables"
        stacked_on = "base"
        stacked_type = "nested"
        help = (
            "Retrieve and install variables and secrets from the Deep Origin platform"
        )
        description = (
            "Retrieve the variables and secrets for your bench from the Deep Origin platform and "
            "install them into the bench. Includes your variables and secrets, as well as those of "
            "the parent organization of the bench."
        )


class InstallVariablesController(cement.Controller):
    """controller to retrieve and install variables and secrets from the Deep Origin platform"""

    class Meta:
        label = "variables-install"
        aliases = ["install"]
        stacked_on = "variables"
        stacked_type = "nested"
        help = "Retrieve variables and secrets from the Deep Origin platform and install them"
        description = (
            "Retrieve the variables and secrets for your bench from the Deep Origin platform and "
            "install them into the bench. Includes your variables and secrets, as well as those of "
            "the parent organization of the bench.\n"
            "\n"
            "By default, this will not overwrite direct changes to variables. To overwrite direct "
            "changes, include the `--overwrite` argument."
        )
        arguments = [
            # TODO: Enable updating variables only for users or only for organizations
            # (
            #     ["--no-user"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install your user variables and secrets [default: False]",
            #     },
            # ),
            # (
            #     ["--no-org"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install the variables and secrets for the parent organization of your bench [default: False]",
            #     },
            # ),
            (
                ["--type"],
                {
                    "dest": "types",
                    "type": str,
                    "nargs": "+",
                    "default": list(VariableType.__members__.keys()),
                    "help": f"Type of variable/secret to install [default: {', '.join(sorted(VariableType.__members__.keys()))}]",
                },
            ),
            (
                ["--overwrite"],
                {
                    "action": "store_true",
                    "help": "Whether to overwrite direct changes to variables and secrets.",
                },
            ),
        ]

    @cement.ex(hide=True)
    def _default(self):
        """default action when no sub-command is passed"""

        feature_flags = get_feature_flags()
        if not feature_flags.variables:
            msg = "Updating variables is not yet available. For beta access, please contact support at support@deeporigin.com."
            warnings.warn(msg, FeatureNotAvailableWarning)
            return

        args = self.app.pargs

        # TODO: Enable updating variables only for users or only for organizations
        # user = not args.no_user
        # org = not args.no_org
        user = True
        org = True

        types = get_variable_types_by_values(args.types)
        overwrite = args.overwrite

        variable_modifications = install_variables(
            user=user, org=org, types=types, overwrite=overwrite
        )

        modified_variable_names = []
        added_variable_names = []
        deleted_variable_names = []
        unmodified_variable_names = []
        for variable_modification in variable_modifications.values():
            variable_type = variable_modification["type"]
            variable_name = variable_modification["name"] or ""
            variable_status = variable_modification["status"]

            if variable_status == VariableStatus.modified:
                modified_variable_names.append(f"{variable_type}: {variable_name}")
            elif variable_status == VariableStatus.added:
                added_variable_names.append(f"{variable_type}: {variable_name}")
            elif variable_status == VariableStatus.deleted:
                deleted_variable_names.append(f"{variable_type}: {variable_name}")
            else:
                unmodified_variable_names.append(f"{variable_type}: {variable_name}")

        modified_variable_names.sort()
        added_variable_names.sort()
        deleted_variable_names.sort()
        unmodified_variable_names.sort()

        num_modified_variable_names = len(modified_variable_names)
        num_added_variable_names = len(added_variable_names)
        num_deleted_variable_names = len(deleted_variable_names)
        num_unmodified_variable_names = len(unmodified_variable_names)

        modified_variable_names = "\n  ".join(modified_variable_names)
        added_variable_names = "\n  ".join(added_variable_names)
        deleted_variable_names = "\n  ".join(deleted_variable_names)
        unmodified_variable_names = "\n  ".join(unmodified_variable_names)

        if num_modified_variable_names:
            print(
                f"{num_modified_variable_names} variables were modified:\n  {modified_variable_names}"
            )
        else:
            print("No variables were modified")

        print("")

        if num_added_variable_names:
            print(
                f"{num_added_variable_names} variables were added:\n  {added_variable_names}"
            )
        else:
            print("No variables were added")

        if num_deleted_variable_names:
            print(
                f"{num_deleted_variable_names} variables were deleted:\n  {deleted_variable_names}"
            )
        else:
            print("No variables were deleted")

        if num_unmodified_variable_names:
            print(
                f"{num_unmodified_variable_names} variables were unmodified:\n  {unmodified_variable_names}"
            )
        else:
            print("No variables were unmodified")

        n_env_changes = 0
        for variable_modification in variable_modifications.values():
            if (
                variable_modification["type"]
                in [
                    VariableType.EnvironmentVariable.name,
                    VariableType.SecretEnvironmentVariable.name,
                ]
                and variable_modification["status"] != VariableStatus.unmodified
            ):
                n_env_changes += 1

        if n_env_changes:
            env_filename = EnvironmentVariable.get_env_filename()
            print("")
            print(
                (
                    f"{n_env_changes} environment variables have changed. "
                    "To install their new values into this shell, run the following command:"
                )
            )
            print(f"  set -o allexport && source {env_filename} && set +o allexport")


class AutoInstallVariablesController(cement.Controller):
    """controller to automatically install variables and secrets from the Deep Origin platform"""

    class Meta:
        label = "variables-auto-install"
        aliases = ["auto-install"]
        stacked_on = "variables"
        stacked_type = "nested"
        help = "Enable variables and secrets to be automatically retrieved from the Deep Origin platform and installed"
        description = (
            "Enable or disable the variables and secrets for your bench from the Deep Origin "
            "platform to be automatically installed once added or modified. Includes your "
            "variables and secrets, as well as those of the parent organization of the bench.\n"
            "\n"
            "By default, this will not overwrite direct changes to variables. To overwrite direct "
            "changes, include the `--overwrite` argument."
        )
        arguments = [
            # TODO: Enable updating variables only for users or only for organizations
            # (
            #     ["--no-user"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install your user variables and secrets [default: False]",
            #     },
            # ),
            # (
            #     ["--no-org"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install the variables and secrets for the parent organization of your bench [default: False]",
            #     },
            # ),
            (
                ["--type"],
                {
                    "dest": "types",
                    "type": str,
                    "nargs": "+",
                    "default": list(VariableType.__members__.keys()),
                    "help": f"Type of variable/secret to install [default: {', '.join(sorted(VariableType.__members__.keys()))}]",
                },
            ),
            (
                ["--overwrite"],
                {
                    "action": "store_true",
                    "help": "Whether to overwrite direct changes to variables and secrets.",
                },
            ),
            (
                ["--time"],
                {
                    "dest": "time_period_min",
                    "type": int,
                    "default": 30,
                    "help": "Time period for updating variables from the Deep Origin platform [default: 30 min, units: min]",
                },
            ),
            (
                ["--disable"],
                {
                    "action": "store_true",
                    "help": "Whether to disable auto updating [default: False]",
                },
            ),
        ]

    @cement.ex(hide=True)
    def _default(self):
        """default action when no sub-command is passed"""

        feature_flags = get_feature_flags()
        if not feature_flags.variables:
            msg = "Updating variables is not yet available. For beta access, please contact support at support@deeporigin.com."
            warnings.warn(msg, FeatureNotAvailableWarning)
            return

        args = self.app.pargs

        enable = not args.disable

        if enable:
            # TODO: Enable updating variables only for users or only for organizations
            # user = not args.no_user
            # org = not args.no_org
            user = True
            org = True

            types = get_variable_types_by_values(args.types)
            overwrite = args.overwrite
            time_period_min = args.time_period_min

            enable_variable_auto_updating(
                user=user,
                org=org,
                types=types,
                time_period_min=time_period_min,
                overwrite=overwrite,
            )
        else:
            disable_variable_auto_updating()


class UninstallVariablesController(cement.Controller):
    """controller to uninstall variables and secrets from the Deep Origin platform"""

    class Meta:
        label = "variables-uninstall"
        aliases = ["uninstall"]
        stacked_on = "variables"
        stacked_type = "nested"
        help = "Uninstall variables and secrets retrieved from the Deep Origin platform"
        description = (
            "Remove variables and secrets retrieved from the Deep Origin platform, including "
            "uninstalling variables and secrets.\n"
            "\n"
            "By default, this will not overwrite direct changes to variables. To overwrite direct "
            "changes, include the `--overwrite` argument."
        )
        arguments = [
            # TODO: Enable updating variables only for users or only for organizations
            # (
            #     ["--no-user"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install your user variables and secrets [default: False]",
            #     },
            # ),
            # (
            #     ["--no-org"],
            #     {
            #         "action": "store_true",
            #         "help": "Whether to not pull and install the variables and secrets for the parent organization of your bench [default: False]",
            #     },
            # ),
            (
                ["--type"],
                {
                    "dest": "types",
                    "type": str,
                    "nargs": "+",
                    "default": list(VariableType.__members__.keys()),
                    "help": f"Type of variable/secret to install [default: {', '.join(sorted(VariableType.__members__.keys()))}]",
                },
            ),
            (
                ["--overwrite"],
                {
                    "action": "store_true",
                    "help": "Whether to overwrite direct changes to variables and secrets.",
                },
            ),
        ]

    @cement.ex(hide=True)
    def _default(self):
        """default action when no sub-command is passed"""

        feature_flags = get_feature_flags()
        if not feature_flags.variables:
            msg = "Updating variables is not yet available. For beta access, please contact support at support@deeporigin.com."
            warnings.warn(msg, FeatureNotAvailableWarning)
            return

        args = self.app.pargs

        # TODO: Enable updating variables only for users or only for organizations
        # user = not args.no_user
        # org = not args.no_org
        user = True
        org = True

        types = get_variable_types_by_values(args.types)
        overwrite = args.overwrite

        uninstall_variables(
            user=user,
            org=org,
            types=types,
            overwrite=overwrite,
        )


CONTROLLERS = [
    UninstallVariablesController,
    AutoInstallVariablesController,
    InstallVariablesController,
    VariablesController,
]
