#!/usr/bin/env python

import runpy
import pathlib

def test_get_notebooks():
    #get current working directory
    path = pathlib.Path()

    #want to be sure the directory is clean
    assert(not (path/'NDPReduce.ipynb').exists())
    assert(not (path/'Schema.ipynb').exists())

    #call module function
    runpy.run_module('ndp.get_notebooks',run_name='__main__')

    #assert that it worked
    assert((path/'NDPReduce.ipynb').exists())
    assert((path/'Schema.ipynb').exists())

    #clean up after yourself
    (path/'NDPReduce.ipynb').unlink()
    (path/'Schema.ipynb').unlink()


