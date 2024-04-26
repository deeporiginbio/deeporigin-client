# The purpose of this script is to catch any warnings
# during the docs built process and surface them as errors
# so that there are no missing links or other issues



if [ "$CI" = "true" ]; then
  echo "Running in GitHub Actions runner, installing repo."
  echo "🚧 Installing repo using pip..."
  pip install -qqq --upgrade pip
  pip install -qqq -e .[docs]
  echo "Installed using pip."
  if ! command -v mkdocs &> /dev/null
  then
      echo "❌ mkdocs could not be found. Fatal"
      exit 2;
  fi

  MKDOCS_OUT="$(mkdocs build -s 2>&1)"

else
  echo "Running Locally, will not install."
fi


if [ "$?" -gt 0 ]; then
  echo "Something went wrong building docs. The error is:";
  echo $MKDOCS_OUT
  exit 3;
fi
warnings=$(echo $MKDOCS_OUT | grep "WARNING" | wc -l)
if [ "$warnings" -gt 0 ]; then
  echo "WARNINGS were found when making docs; aborting.";
  exit 4;
fi

echo "Built docs successfully"





