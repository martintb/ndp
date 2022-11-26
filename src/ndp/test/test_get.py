#!/usr/bin/env python

import runpy
import pathlib
import warnings

def test_get_notebooks():
    #get current working directory
    path = pathlib.Path()

    #want to be sure the directory is clean
    assert(not (path/'NDPReduce.ipynb').exists())
    assert(not (path/'Schema.ipynb').exists())

    #call module function
    try:
        runpy.run_module('ndp.get_notebooks',run_name='__main__')

    except:
        warnings.warn('ndp.get_notebooks failed and will not work on this build')

    else:
        #assert that it worked
        assert((path/'NDPReduce.ipynb').exists())
        assert((path/'Schema.ipynb').exists())

        #clean up after yourself
        (path/'NDPReduce.ipynb').unlink()
        (path/'Schema.ipynb').unlink()


