#!/usr/bin/env python

# read version from installed package
from importlib.metadata import version
__version__ = version(__name__)


from ndp.reduce import *
from ndp.schema import *
