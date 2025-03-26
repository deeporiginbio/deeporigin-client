# install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh

# If uv tells you to add something to PATH,
# you may need to do that first
uv init

# we use the latest stable version of python
uv python install 3.13

# install deeporigin with the tools extra
uv add --upgrade deeporigin --extra tools

# also install as a uv tool to run from the command line
uv tool install deeporigin

