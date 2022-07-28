# -*- coding: utf-8 -*-
"""
Neutron Depth Profiling Data Reduction
"""

import numpy as np
import math
import re
import os
import csv
import json
from datetime import datetime


class NdpData():
    """
    Class that defines the data object for neutron depth profiling (NDP).
    Data object describes values collected for a single sample, including
    sample, background, and reference with associated monitors and
    relevant parameters required for data processing.
    
    ndp.data is intended to provide enough information to recreate the entire data reduction
    """
    
    def __init__(self):
        
        
        self.instrument = {
            "Configuration" : "Default",
            "Beam Energy": 1472.35,
            "Num Channels": 4096,
            "Zero Channel": 2077,
            "Mon Peak Channels": [
                1900,
                2901
            ],
            "Calib Coeffs": [
                0.7144,
                -12.45
            ]
        }
        
        self.detector = {
            "Name" : "Lynx",
            "Channels" : np.arange(0,4096),
            "Energy" : np.zeros(4096),
            "Depth" : np.zeros(4096),
            "Corr Depth" : np.zeros(4096),
            "dDepth" : np.zeros(4096),
            "dDepth Uncert" : np.zeros(4096)
        }
        

        self.data = {
            "Sam Dat" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "Sam Mon" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "Bgd Dat" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "Bgd Mon" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "Ref Dat" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "Ref Mon" : {
                "Files" : [],
                "Labels" : [],
                "Detector" : [],
                "Live Time" : 0.0,
                "Real Time" : 0.0,
                "Operations" : []
            },
            "TRIM" : {
                "Files" : []
            },
        }
        
        self.atom = {
            "He" : {
                "Cross Sec" : 5322.73,
                "Abundance" : 0.00000134
                },
            "Li" : {
                "Cross Sec" : 939.09,
                "Abundance" : 0.0759
                },
            "B" : {
                "Cross Sec" : 3600.48,
                "Abundance" : 0.196
                },
            "N" : {
                "Cross Sec" : 1.86,
                "Abundance" : 0.99636
                },
            }

    
    def readconfig(self, config_filename = "NDPInstrumParms.dat"):
        """
        Read NDPReduce configuration file
        """
        
        with open(config_filename) as f:
            self.instrument = json.load(f)
        return


    def runschema(self, schemafilename = "schema.txt"):
        """
        Run the operations listed in the schema file
        """

        bin_flag = True #Bin op must run, even if bin size is 1
        self.readschema(schemafilename)
        
        ops = self.schema['Operations']
        
        for op in ops:
            if 'Eval' in op:
                self.data['TRIM']["Path"] = self.schema['TRIM']['Path']
                filelist = os.listdir(self.data['TRIM']["Path"])
                self.data['TRIM']["Files"] = [x for x in filelist if self.schema['TRIM']['Tag'] in x] 
                self.evalTRIM(self.schema['TRIM']['Path'])
                self.chan2depth()
            if 'Load' in op:
                for dt in self.schema['Load']:
                    self.data[dt]["Path"] = self.schema[dt]['Path']
                    filelist = os.listdir(self.data[dt]["Path"])
                    self.data[dt]["Files"] = [x for x in filelist if self.schema[dt]['Tag'] in x] 
                    self.load_NISTNDP(dt)
            if 'Norm' in op:
                for dt in self.schema['Norm']:
                    self.normalize(dt)
            if 'Corr' in op:
                for dt in self.schema['Corr']:
                    self.correct(dt)
            if 'Absolute' in op:
                for dt in self.schema['Absolute']:
                    self.cross_sec("B", "Sam Dat")
                    self.cross_sec("B", "Ref Dat")
                    self.ref_integrate()
                    self.scale2ref("Sam Dat")                    
            if 'Bin' in op:
                bin_size = int(self.schema['Bin'])
                self.bin_channels(bin_size)
                bin_flag = False
            if 'Save' in op:
                filename = self.schema['Save']['Filename']
                path = self.schema['Save']['Path']
                columns = self.schema['Save']['Columns']
                self.saveAtoms(path, filename, columns)
        
        #If user did not bin, then run at the end to get data into binned arrays (binsize = 1)
        if bin_flag: 
            bin_size = 1
            self.bin_channels(bin_size)
            bin_flag = False
            
        return

        

