#!/bin/bash

python3 -m build
twine upload dist/*$1*.whl dist/*$1*.tar.gz