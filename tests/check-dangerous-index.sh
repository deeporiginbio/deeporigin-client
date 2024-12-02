#!/bin/bash

# Search for the pattern and count the matches
if [ "$(grep -rE 'if [a-zA-Z_][a-zA-Z0-9_]*\.[a-zA-Z_][a-zA-Z0-9_]* is None:' src | wc -l)" -gt 0 ]; then
    echo "foo.bar is None pattern found! Replace with \`if getattr(foo, 'bar', None) is None:\`"
    grep -rE 'if [a-zA-Z_][a-zA-Z0-9_]*\.[a-zA-Z_][a-zA-Z0-9_]* is None:' src
    exit 1
else
    echo "No matches found."
    exit 0
fi