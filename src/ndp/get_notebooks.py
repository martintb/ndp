#!/usr/bin/env python
from importlib.resources import files
import shutil
import os
import pathlib


if __name__=="__main__":
    current_directory = pathlib.Path()

    for filepath in files('ndp.jupyter').joinpath('').glob('*ipynb'):
        shutil.copy(filepath,current_directory)
