#!/bin/bash

# This generates file service client file and copies the the resulting files to
# deeporigin-client/src/files.

# Get script directory and save current directory
SCRIPT_DIR="$(dirname "$0")"
CURRENT_DIR="$(pwd)"

# Change to script directory for generation
cd "$SCRIPT_DIR"

OPENAPI_YAML_PATH="./gen_file_openapi.yaml"

# Check if "platform" argument is passed. Platfrom auto-generated OpenAPI spec
# had bugs so we keep local modified copy for now.
if [[ "$1" == "platform" ]]; then    
    OPENAPI_YAML_PATH="../../platform/docs/file-service/openapi.yaml"
    echo "Using OpenAPI spec from platform path: $OPENAPI_YAML_PATH"
else
    echo "Using local OpenAPI spec: $OPENAPI_YAML_PATH"
fi

# Generate using openapi-python-client
openapi-python-client generate \
--path "$OPENAPI_YAML_PATH" \
--config .openapi-python-client.yml
#--custom-template-path "$SCRIPT_DIR"/templates

# Create target directory if it doesn't exist
mkdir -p ../src/files/file_service

# Copy generated files to target location
cp -r file-service/file_service/* ../src/files/file_service/
cp file-service/README.md ../src/files/file_service/

# Clean up generated files,
rm -rf file-service
# Return to original directory
cd "$CURRENT_DIR"


<<'COMMENT'
# THis was used to generate the file service client before openapi-python-client was used.
# Keeping it for now in case we need it again.

openapi-generator-cli generate \
-i "$SCRIPT_DIR"/../../platform/docs/file-service/openapi.yaml \
-g python \
-o "$SCRIPT_DIR"/../src/files \
--package-name file_service \
--global-property skipDocs=true \
--global-property skipTests=true \
--global-property generateAliasAsModel=false \
--additional-properties generateSourceCodeOnly=true

# Doesn't seem to work:
# --additional-properties singleSourceFile=true

# Fix imports to be relative to 'files' package
find "$SCRIPT_DIR"/../src/files/file_service -type f -name "*.py" -exec sed -i 's/from file_service\./from files.file_service./g' {} \;
find "$SCRIPT_DIR"/../src/files/file_service -type f -name "*.py" -exec sed -i 's/import file_service\./import files.file_service./g' {} \;

# Clean up individual files and test directories
rm -rf "$SCRIPT_DIR"/../src/files/file_service/models/
rm -rf "$SCRIPT_DIR"/../src/files/file_service/test
COMMENT
