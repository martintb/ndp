#!/usr/bin/env python
import ndp
import json
import os
from pprint import pprint
from ndp import schema

def test_schema():
    '''Just try to instantiate and load some values into a schema

    A more in-depth test that actually loads files is needed elsewhere
    '''
    heat0x = ndp.schema()

    heat0x.schema['Operations'] = ['Eval', 'Load', 'Norm', 'Correct', 'Absolute', 'Bin', 'Save']
    heat0x.schema['Eval'] = ['TRIM']
    heat0x.schema['Load'] = ['Sam Dat', 'Sam Mon', 'Bgd Dat', 'Bgd Mon', 'Ref Dat', 'Ref Mon']
    heat0x.schema['Norm'] = ['Sam', 'Ref', 'Bgd']
    heat0x.schema['Corr'] = ['Sam', 'Ref']
    heat0x.schema['Absolute']['Atom'] = 'B'
    heat0x.schema['Absolute']['Cross Sec'] = 3600.48
    heat0x.schema['Absolute']['Abundance'] = 0.196
    heat0x.schema['Absolute']['Ref Conc'] = 5.22e15
    heat0x.schema['Absolute']['Ref Conc Uncert'] = 3e13
    heat0x.schema['Absolute']['Branch Frac'] = 0.94
    heat0x.schema['Bin'] = 21
    heat0x.schema['Save']['Columns'] = ['Channels', 'Energy', 'Depth', 'Counts', 'Atoms/cm2', 'Atoms/cm2 Uncert', 'Atoms/cm3', 'Atoms/cm3 Uncert']


