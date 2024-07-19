#!/bin/bash

python3 setup.py sdist bdist_wheel
twine upload dist/*$1*.whl dist/*$1*.tar.gz