# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:48:13 2022

@author: rljones
"""

import json
import os

class schema():
    """
    Class that defines a data object for neutron depth profiling (NDP) data analysis.
    Data object describes values collected for a single sample, including
    sample, background, and reference with associated monitors and
    relevant parameters required for data processing.
    
    ndp.data is intended to provide enough information to recreate the entire data reduction
    """
    
    def __init__(self):
       self.schema = {
        'Operations' : ['Load', 'Bin', 'Save'],
        'Eval' : {
            'Type': 'TRIM',
            },            
        'Load' : ['Sam Dat'],
        'Norm' : [],
        'Corr' : [],
        'Absolute' : {
            'Ref Atom' : 'B',
            'Ref Cross Sec' : 3600.48,
            'Ref Abundance' : 0.196,
            'Ref Conc' : 5.22e15,
            'Ref Conc Uncert' : 3e13,
            'Atom' : '',
            'Cross Sec' : 0.0,
            'Abundance' : 0.0,
            'Branch Frac' : 0.0,
            },
        'Bin' : 1,
        'Save' : {
            'Columns' : ['Channels', 'Counts'],
            'Path' : '',
            'Filename' : 'ndp_default.csv'
        },
        'TRIM' : [],
        'Sam Dat' : {
            'Path' : '',
            'Files' : ''
        },
        'Sam Mon' : {
            'Path' : '',        
            'Files' : ''
        },
        'Bgd Dat' : {
            'Path' : '',
            'Files' : ''
        },
        'Bgd Mon' : {
            'Path' : '',
            'Files' : ''
        },
        'Ref Dat' : {
            'Path': '',
            'Files': ''
        },
        'Ref Mon' : {
            'Path' : '',
            'Files': ''
        },
    }

 

    def load_schema(self, schemafile):
        """
        Load a previously written schema file
        
        Returns
        A schema data object
        """
        
        with open(schemafile, 'r') as f:
            schema = json.load(f)
            
        return schema
    

    def save_schema(self, schemafile):
        """
        Save a schema data object to a JSON file
        """

        with open(schemafile, 'w') as f:
            json.dump(self.schema, f, indent=4)
            
    def get_filelist(self, path, tag):
        """
        Generate a list of filenames in path that contain (tag) in the name

        Parameters
        ----------
        path : string
            Path to the directory containing the files
        tag : string
            Text common to all filenames of interest. Use a null string if all 
            files in a directory are desired.

        Returns
        -------
        List of strings with filenames

        """
        
        dirlist = os.listdir(path)
        return [x for x in dirlist if tag in x]

    
    def add_TRIMlayer(self, path, tag):
        """
        Add a material layer to the TRIM schema. Function  modifies a schema data object.
        
        Parameters
        ----------
        path : string
            Path to the TRIM files defining the energy to depth relationship for a given material layer
        tag : string
            Text common to all filenames of interest. Use a null string if all 
            files in a directory are desired.

        Returns
        -------
        None

        """
        
        new_layer = {
            'Path' : path,
            'Files' : self.get_filelist(path,tag)
        }
        self.schema['TRIM'].append(new_layer)
        

        
