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


class ndpData():
    """
    Modules that process neutron depth profiling (NDP) data.
    schema = a JSON file with parameters describing input files and other parameters.
    instrument = data structure defining parameters related to the instrument
    data = object that stores the entire data reduction process

    output = user defined comma separated value text file with 2 to 6 columns    
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
        
        self.TRIM = {
            'Type' : 'TRIM',
            'Layers' : [],
            'Coeffs' : np.zeros(3),
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
        }
        
    
    def readconfig(self, config_filename="instrument.dat"):
        """
        Read NDPReduce configuration file
        """
        
        with open(config_filename) as f:
            self.instrument = json.load(f)
        return


    def runschema(self, schemafilename = "schema.json"):
        """
        Run the operations listed in the schema file
        """

        bin_flag = True #Bin op must run, even if bin size is 1
        self.readschema(schemafilename)
        
        ops = self.schema['Operations']
        
        for op in ops:
            if 'Eval' in op:
                kev, depth, mat = self.evalTRIM(self.schema['TRIM'], self.instrument['Beam Energy'])
                self.TRIM['Layers'] = self.schema['TRIM']
                self.TRIM['Coeffs'] = self.fit_TRIM(kev, depth)
                self.TRIM['Energy'] = kev
                self.TRIM['Depth'] = depth
                self.TRIM['Material'] = mat
                self.chan2depth()
            if 'Load' in op:
                for dt in self.schema['Load']:
                    self.data[dt]["Path"] = self.schema[dt]['Path']
                    self.data[dt]["Files"] = self.schema[dt]["Files"]
                    self.loadfiles(dt)
            if 'Norm' in op:
                for dt in self.schema['Norm']:
                    self.normalize(dt)
            if 'Corr' in op:
                for dt in self.schema['Corr']:
                    self.correct(dt)
            if 'Absolute' in op:
                self.data['Ref Dat']["Atom"] = self.schema['Absolute']['Ref Atom'] 
                self.data['Ref Dat']["Cross Sec"] = self.schema['Absolute']['Ref Cross Sec'] 
                self.data['Ref Dat']["Abundance"] = self.schema['Absolute']['Ref Abundance']
                self.data['Ref Dat']["Conc"] = self.schema['Absolute']['Ref Conc'] 
                self.data['Ref Dat']["Conc Uncert"] = self.schema['Absolute']['Ref Conc Uncert']
                self.data['Sam Dat']["Atom"] = self.schema['Absolute']['Atom'] 
                self.data['Sam Dat']["Cross Sec"] = self.schema['Absolute']['Cross Sec'] 
                self.data['Sam Dat']["Abundance"] = self.schema['Absolute']['Abundance']
                self.data['Sam Dat']["Branch Frac"] = self.schema['Absolute']['Branch Frac']                    
                self.ref_integrate()
                self.scale2ref()                    
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

    def loadfiles(self, dt):
        """
        Function to load a list of NDP data files of a given datatype (Sam Dat, Sam Mon, etc),
        load header info and sum of counts/channel into ndp.data
        
        """
    
        path = self.data[dt]["Path"]
        filelist = self.data[dt]["Files"]
        numfiles = len(filelist)
        numchannels = self.instrument["Num Channels"]
        
        if("Channel Sum" not in self.data[dt]["Operations"]):
                self.data[dt]["Operations"].append("Channel Sum")
        
        for filenum in range(numfiles):
            ndp_file = path + filelist[filenum]
            with open(ndp_file) as f:
                lines = f.readlines() #reads all of the file into a numbered list of strings
                self.data[dt]["Detector"].append(lines[0][12:-1])
                self.data[dt]["Labels"].append(lines[1][8:-1])
                self.data[dt]["Live Time"] += float(lines[3][12:-1])
                self.data[dt]["Real Time"] += float(lines[4][12:-1])
            
                if(filenum<1): 
                    self.data[dt]["Datetime"] = \
                        datetime.strptime(lines[2][12:-10],'%a %b %y %H:%M:%S')
                    self.data[dt]["Counts"] = np.zeros(numchannels)
                
                for channel in range(numchannels):
                    counts = lines[channel+8].split()
                    self.data[dt]["Counts"][channel] += float(counts[1])
    
        self.deadtime(dt)
        
        return
        
    def readschema(self, filename):
        """
        Read data schema file
        """
                   
        with open(filename) as f:
            self.schema = json.load(f)
        
        
        
        return
    
    
    def evalTRIM(self, trim_list, beam_energy):
        """
        Read in blocks of TRIM simulation data files to calculate polynomial
        coefficients that relate depth to energy of an ion escaping a sample
        of multiple material layers

        Parameters
        ----------
        trim_list : ordered list of data defining a material layer. Each list item
                    is a dictionary with keys <Path> and <Tag>. 
        
        Returns
        -------
        data object containing trim_list and an additional 

        """

        # Regex to extract the material label from TRIM header (line 10)
        p1 = re.compile(r'\>.*\)')    
        # Regex to extract energy from the third or fourth column of the TRIM file
        p2 = re.compile(r'\.[0-9]*E\+[0-9]*')
        # Regex to extract depth from the fourth or fifth column of the TRIM file
        p3 = re.compile(r'[0-9]*E-[0-9]*')
        

        # TRIM calculations are done for each material layer, with multiple
        # files within a layer. TRIM provides the relative energy and depth to 
        # the top interface. Offset determines the energy and depth of the
        # top interface of the current layer.
        ev_offset = 0.0
        depth_offset = 0.0        
        energy = np.array([beam_energy])
        thick = np.array([0.0])
        coeffs = np.zeros(3)
        material = []
        
        for layer in trim_list:
            for filename in layer['Files']:
                trim_file = layer['Path'] + filename
        
                with open(trim_file) as f:
                    lines = f.readlines()
                    num_lines = len(lines)-12
                    ev = np.zeros(num_lines)
                    depth = np.zeros(num_lines)
                    m = p1.search(lines[9])
                    material.append(m.group()[1:])
                    i=0
                    for line in lines[12:]:        
                        m = p2.search(line)
                        ev[i] = float(m.group())/1000.0 - ev_offset
                        m = p3.search(line)
                        depth[i] = float(m.group())/10 + depth_offset
                        i += 1
                    
                energy = np.append(energy, np.median(ev))
                thick = np.append(thick, np.average(depth))

            ev_offset = energy[0] - np.min(energy)
            depth_offset = thick[0] + np.max(thick)
                
        return energy, thick, material
    
    
    def fit_TRIM(self, energy, thick):
        """
        Fit a polynomial to thickness and energy data and return the three coefficients
        
        Parameters
        ----------
        median_KeV : Escape energy of each sub-layer
            DESCRIPTION.
        thick : 1-D numpy array
            Thickness of each sub-layer

        Returns
        -------
        1x3 numpy array of coeffiients for the 2nd order polynomial fit

        """
        return np.polyfit(energy, thick, deg=2)
        

    def chan2depth(self):

    
        # These values change infrequently and are provided by the instrument scientist
        m, b = self.instrument["Calib Coeffs"]
        self.detector["Energy"] = (m*self.detector["Channels"]) + b
    
        # These values are derived from SRIM/TRIM, freeware used to calculate energy of generated ions in matter
        # Depth is in nanometers
        a, b, c = self.TRIM["Coeffs"]
        self.detector["Depth"] = a*np.power(self.detector["Energy"],2) \
            + b*self.detector["Energy"] + c
        
        # Zero channel is defined through the experimental setup
        zerochan = self.instrument["Zero Channel"]
        self.detector["Corr Depth"] = self.detector["Depth"] \
            - self.detector["Depth"][zerochan]
        
        # Del Depth in centimeters
        self.detector["Del Depth"] = np.zeros(len(self.detector["Corr Depth"]))
        for x in range(len(self.detector["Corr Depth"])-1):
            self.detector["Del Depth"][x] = 1e-7*(self.detector["Corr Depth"][x-1] - self.detector["Corr Depth"][x])
    
        self.detector["Del Depth Uncert"] = 0.05 * self.detector["Del Depth"]
        
        return
    
           

    def deadtime(self, dt):
        """
        Returns a copy of ndp.data with deadtime corrected counts for each of the
        sample types. Optional argument to specify the datatypes.
        
        """
    
        old_settings = np.seterr(all='ignore')  #seterr to known value
        np.seterr(all='ignore')
    
    
        if("Deadtime Scaled" not in self.data[dt]["Operations"]):
            self.data[dt]["Operations"].append("Deadtime Scaled")
        livetime = self.data[dt]["Live Time"]
        realtime = self.data[dt]["Real Time"]
        self.data[dt]["Dt ratio"] = livetime/realtime
        self.data[dt]["Counts/Dt"] = self.data[dt]["Counts"]*realtime/livetime
        self.data[dt]["Counts/Dt Uncert"] = np.nan_to_num(np.divide(
            self.data[dt]["Counts/Dt"],np.sqrt(self.data[dt]["Counts"])))
    
        np.seterr(**old_settings)
    
        return
    
    
    def normalize(self, dt):
        """
        Calculate (data file counts)/(monitor file counts)
    
        returns ndp_norm
        """
    
        old_settings = np.seterr(all='ignore')  #seterr to known value
        np.seterr(all='ignore')
    
        dt_dat = dt + " Dat"
        dt_mon = dt + " Mon"
        
        if("Normalized" not in self.data[dt_dat]["Operations"]):
            self.data[dt_dat]["Operations"].append("Normalized")
    
        #Sum over range set to capture 10B alpha peaks (channels 1900-2900)
        lowchan, hichan = self.instrument["Mon Peak Channels"]
        mon_sum = np.sum(self.data[dt_mon]["Counts/Dt"][lowchan:hichan])        
        self.data[dt_dat]["Monitor"] = mon_sum
        self.data[dt_dat]["Monitor Uncert"] = math.sqrt(mon_sum)    
    
        self.data[dt_dat]["Norm Cts"] = self.data[dt_dat]["Counts/Dt"]/self.data[dt_dat]["Monitor"]
        x2 = np.nan_to_num(np.power(self.data[dt_dat]["Counts/Dt Uncert"]/self.data[dt_dat]["Counts/Dt"],2))
        y2 = math.pow(self.data[dt_dat]["Monitor Uncert"]/self.data[dt_dat]["Monitor"],2)
        self.data[dt_dat]["Norm Cts Uncert"] = self.data[dt_dat]["Norm Cts"]*np.sqrt(x2+y2)
    
        np.seterr(**old_settings)
    
        return
    
    
    def correct(self, dt):
        """
        Subtract background and return a corrected data file
        
        data_norm, bkgd_norm are normalized data objects from ndp_norm()
        """
        
        dt_dat = dt + " Dat"
        if("Corrected" not in self.data[dt_dat]["Operations"]):
            self.data[dt_dat]["Operations"].append("Corrected")
    
        self.data[dt_dat]["Corr Cts"] = self.data[dt_dat]["Norm Cts"]-self.data["Bgd Dat"]["Norm Cts"]
        x2 = np.power(self.data[dt_dat]["Norm Cts Uncert"],2)
        y2 = np.power(self.data["Bgd Dat"]["Norm Cts Uncert"],2)
        self.data[dt_dat]["Corr Cts Uncert"] = np.sqrt(x2+y2)
    
        return
    
            
    def ref_integrate(self):
        """
        Integrate the alpha peaks of the reference data set
        Also, set atomic concentration field here for now
        """
    
        dt = 'Ref Dat'
        if("Integrated Peaks" not in self.data[dt]["Operations"]):
            self.data[dt]["Operations"].append("Integrated Peaks")
    
        self.data[dt]["alpha*"] = np.sum(self.data[dt]["Corr Cts"][1791:2142])        
        self.data[dt]["alpha"] = np.sum(self.data[dt]["Corr Cts"][2291:2592])        
        cts_uncert2 = np.power(self.data[dt]["Corr Cts Uncert"], 2)
        self.data[dt]["alpha* Uncert"] = math.sqrt(np.sum(cts_uncert2[1791:2142]))
        self.data[dt]["alpha Uncert"] = math.sqrt(np.sum(cts_uncert2[2291:2592]))
            
        return
    
    
    def scale2ref(self):
        """
        Use reference sample data to convert counts to number of atoms
        """
    
        old_settings = np.seterr(all='ignore')  #seterr to known value
        np.seterr(all='ignore')
    
        if("Scaled to Reference" not in self.data['Sam Dat']["Operations"]):
            self.data['Sam Dat']["Operations"].append("Scaled to Reference")
            
        alpha_cts = self.data["Ref Dat"]["alpha*"]
        sample_cross = self.data['Sam Dat']["Cross Sec"]
        ref_cross = self.data["Ref Dat"]["Cross Sec"]
        ref_conc = self.data['Ref Dat']["Conc"]
        branch_frac = self.data['Sam Dat']["Branch Frac"]
        abundance = self.data['Sam Dat']['Abundance']
    
        scale_coeff = (ref_conc * ref_cross) / (alpha_cts * sample_cross)
        self.data['Sam Dat']["Atoms/cm2"] = scale_coeff * self.data['Sam Dat']["Corr Cts"]
        self.data['Sam Dat']["Atoms/cm2"] /= (branch_frac*abundance)
        
        ratio1 = math.pow(self.data["Ref Dat"]["alpha* Uncert"]/alpha_cts,2)
        ratio2 = math.pow(self.data["Ref Dat"]["Conc Uncert"]/ref_conc,2)
        ratio3 = np.nan_to_num(np.power(self.data['Sam Dat']["Corr Cts Uncert"]/self.data['Sam Dat']["Corr Cts"], 2))
        self.data['Sam Dat']["Atoms/cm2 Uncert"] = self.data['Sam Dat']["Atoms/cm2"]*np.sqrt(ratio1 + ratio2 + ratio3)
        
        self.data['Sam Dat']["Atoms/cm3"] = np.nan_to_num(self.data['Sam Dat']["Atoms/cm2"]/self.detector["Del Depth"])
    
        ratio1 = np.nan_to_num(np.power(self.data['Sam Dat']["Atoms/cm2 Uncert"]/self.data['Sam Dat']["Atoms/cm2"],2))
        ratio2 = np.nan_to_num(np.power(self.detector["Del Depth Uncert"]/self.detector["Del Depth"],2))
        self.data['Sam Dat']["Atoms/cm3 Uncert"] = self.data['Sam Dat']["Atoms/cm3"]*np.sqrt(ratio1 + ratio2)
    
        np.seterr(**old_settings)
    
        return
    
    
    def bin_channels(self, bin_size=21):
        """
        Bin channels from ndp.detector
        """

        num_channels = self.instrument["Num Channels"]
        num_bins = int(num_channels/bin_size)+1
                
        self.detector["Channels Binned"] = np.arange(num_bins)
        self.detector["Energy Binned"] = np.zeros(num_bins)
        self.detector["Depth Binned"] = np.zeros(num_bins)
        self.data["Sam Dat"]["Counts Binned"] = np.zeros(num_bins)
        self.data["Sam Dat"]["Atoms/cm2 Binned"] = np.zeros(num_bins)
        self.data["Sam Dat"]["Atoms/cm2 Binned Uncert"] = np.zeros(num_bins)
        self.data["Sam Dat"]["Atoms/cm3 Binned"] = np.zeros(num_bins)
        self.data["Sam Dat"]["Atoms/cm3 Binned Uncert"] = np.zeros(num_bins)
    
        uncert2_1 = np.power(self.data["Sam Dat"]["Atoms/cm2 Uncert"],2)
        uncert2_2 = np.power(self.data["Sam Dat"]["Atoms/cm3 Uncert"],2)
    
        for bin in range(num_bins-1):
            lo = bin*bin_size
            hi = (bin+1)*bin_size
            if hi > num_channels:
                hi = num_channels
                bin_size = hi - lo
            self.detector["Energy Binned"][bin] = np.median(self.detector["Energy"][lo:hi])
            self.detector["Depth Binned"][bin] = np.median(self.detector["Corr Depth"][lo:hi])
            self.data["Sam Dat"]["Counts Binned"][bin] = np.average(self.data["Sam Dat"]["Counts"][lo:hi])
            self.data["Sam Dat"]["Atoms/cm2 Binned"][bin] = \
                np.average(self.data["Sam Dat"]["Atoms/cm2"][lo:hi])
            self.data["Sam Dat"]["Atoms/cm2 Binned Uncert"][bin] = \
                math.sqrt(np.sum(uncert2_1[lo:hi]))/bin_size
            self.data["Sam Dat"]["Atoms/cm3 Binned"][bin] = \
                np.average(self.data["Sam Dat"]["Atoms/cm3"][lo:hi])
            self.data["Sam Dat"]["Atoms/cm3 Binned Uncert"][bin] = \
                math.sqrt(np.sum(uncert2_2[lo:hi]))/bin_size
            
        return
    
    

 
    def saveAtoms(self, path, filename, data_cols):
        """
        Write six column CSV for atoms/cm2 and atoms/cm3
        """
        
        num_cols = len(data_cols)
        len_cols = self.detector['Channels Binned'].size
    
        i=0
        columns = np.zeros((num_cols, len_cols))
        
        for dkey in data_cols:
            if dkey == 'Channels':
                columns[i][:] = self.detector['Channels Binned']
            if dkey == 'Energy':
                columns[i][:] = self.detector['Energy Binned']
            if dkey == 'Depth':
                columns[i][:] = self.detector['Depth Binned']
            if dkey == 'Counts':
                columns[i][:] = self.data['Sam Dat']['Counts Binned']
            if dkey == 'Atoms/cm2':
                columns[i][:] = self.data['Sam Dat']['Atoms/cm2 Binned']
            if dkey == 'Atoms/cm2 Uncert':
                columns[i][:] = self.data['Sam Dat']['Atoms/cm2 Binned Uncert']
            if dkey == 'Atoms/cm3':
                columns[i][:] = self.data['Sam Dat']['Atoms/cm3 Binned']
            if dkey == 'Atoms/cm3 Uncert':
                columns[i][:] = self.data['Sam Dat']['Atoms/cm3 Binned Uncert']
            i += 1
        
        
        header = [['NIST Neutron Depth Profiling Data File'],
                  ['Sample Data Files'],
                  [self.data['Sam Dat']["Files"]],
                  ['Sample Monitor Files'],
                  [self.data['Sam Mon']["Files"]],
                  ['Background Data Files'],
                  [self.data['Bgd Dat']["Files"]],
                  ['Background Monitor Files'],
                  [self.data['Bgd Mon']["Files"]],
                  ['Reference Data Files'],
                  [self.data['Ref Dat']["Files"]],
                  ['Reference Monitor Files'],
                  [self.data['Ref Mon']["Files"]],
                  ['Sample Data Operations'],
                  [self.data['Sam Dat']['Operations']],
                  [self.data['Sam Dat']["Datetime"]],
                  [' '],
                  data_cols
                  ]
        
        numlines = len(self.detector['Channels Binned'])
    
        with open((path+filename), 'w', newline='') as csvfile:
            #using excel comma separated value format
            writer = csv.writer(csvfile, dialect = 'excel')
            for x in range(len(header)):
                writer.writerow(header[x]) 
            writer.writerows(np.transpose(columns)) 
            
        return
