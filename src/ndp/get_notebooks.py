#!/usr/bin/env python
''' Module function for easy access to repository notebooks
'''
from importlib.resources import files
import shutil
import os
import pathlib


if __name__=="__main__":
    current_directory = pathlib.Path()

    for filepath in files('ndp.jupyter').iterdir():
        if filepath.is_file() and ('.ipynb' in filepath.parts[-1]):
            shutil.copy(filepath,current_directory)
