"""Pytest configuration for the MIT1002 model repo.

This file exists so that tests can import repo modules, e.g.
``from tools.deprecate import read_records``. Pytest prepends the directory
containing the rootdir ``conftest.py`` to ``sys.path``, which makes ``scripts``
and ``test`` importable as packages regardless of where pytest is invoked from.

The tests also read ``model.xml`` by relative path, so they assume the repo root
is the working directory. That is what CI does and what ``pytest`` from the repo
root does.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
