# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:45:44 2022

@author: rljones
"""

def readschema(self, filename):
    """
    Read data schema file
    """
               
    with open(filename) as f:
        self.schema = json.load(f)
    
    
    
    return


def evalTRIM(self, path):

    numfiles = len(self.data["TRIM"]["Files"])
    self.data["TRIM"]["Median KeV"] = np.zeros(numfiles+1)
    self.data["TRIM"]["Thick"] = np.zeros(numfiles+1)
    self.data["TRIM"]["Median KeV"][0] = self.instrument["Beam Energy"]
    self.data["TRIM"]["Thick"][0] = 0.0
    self.data["TRIM"]["Coeffs"] = np.zeros(3)
    
    # Regex to extract layer thickness in Angstroms from TRIM header (line 10)
    # p = re.compile(r'\(\s*[0-9]+\s*A\)')

    
    # Regex to extract energy from the third or fourth column of the TRIM file
    p = re.compile(r'\.[0-9]*E\+[0-9]*')
    
    # Regex to extract depth from the fourth or fifth column of the TRIM file
    q = re.compile(r'[0-9]*E-[0-9]*')
    
    
    for filenum in range(numfiles):
        trim_file = path + self.data["TRIM"]["Files"][filenum]

        with open(trim_file) as f:
            lines = f.readlines()
            num_lines = len(lines)-12
            ev = np.zeros(num_lines)
            depth = np.zeros(num_lines)

            i=0
            for line in lines[12:]:        
                m = p.search(line)
                ev[i] = float(m.group())
                m = q.search(line)
                depth[i] = float(m.group())
                i += 1
                
        self.data["TRIM"]["Median KeV"][filenum+1] = np.median(ev)/1000
        self.data["TRIM"]["Thick"][filenum+1] = np.average(depth)/10
        self.data["TRIM"]["Coeffs"] = np.polyfit(\
            self.data["TRIM"]["Median KeV"], self.data["TRIM"]["Thick"], 2)
    
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


def cross_sec(self, atom, dt):
    """
    Define the cross section of the sample
    """
    
    self.data[dt]["Atom"] = atom
    self.data[dt]["Atom Cross Sec"] = self.atom[atom]['Cross Sec']
    self.data[dt]["Atom Abundance"] = self.atom[atom]['Abundance']
    self.data[dt]["Atom Conc"] = 5.22e15
    self.data[dt]["Atom Conc Uncert"] = 3e13
    self.data[dt]["Atom Branch Frac"] = 0.94

        
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


def scale2ref(self, dt):
    """
    Use reference sample data to convert counts to number of atoms
    """

    old_settings = np.seterr(all='ignore')  #seterr to known value
    np.seterr(all='ignore')

    if("Scaled to Reference" not in self.data[dt]["Operations"]):
        self.data[dt]["Operations"].append("Scaled to Reference")
        
    alpha_cts = self.data["Ref Dat"]["alpha*"]
    sample_cross = self.data[dt]["Atom Cross Sec"]
    ref_cross = self.data["Ref Dat"]["Atom Cross Sec"]
    ref_conc = self.data[dt]["Atom Conc"]
    branch_frac = self.data[dt]["Atom Branch Frac"]
    abundance = self.data[dt]["Atom Abundance"]

    scale_coeff = (ref_conc * ref_cross) / (alpha_cts * sample_cross)
    self.data[dt]["Atoms/cm2"] = scale_coeff * self.data[dt]["Corr Cts"]
    self.data[dt]["Atoms/cm2"] /= (branch_frac*abundance)
    
    ratio1 = math.pow(self.data["Ref Dat"]["alpha* Uncert"]/alpha_cts,2)
    ratio2 = math.pow(self.data["Ref Dat"]["Atom Conc Uncert"]/ref_conc,2)
    ratio3 = np.nan_to_num(np.power(self.data[dt]["Corr Cts Uncert"]/self.data[dt]["Corr Cts"], 2))
    self.data[dt]["Atoms/cm2 Uncert"] = self.data[dt]["Atoms/cm2"]*np.sqrt(ratio1 + ratio2 + ratio3)
    
    self.data[dt]["Atoms/cm3"] = np.nan_to_num(self.data[dt]["Atoms/cm2"]/self.detector["Del Depth"])

    ratio1 = np.nan_to_num(np.power(self.data[dt]["Atoms/cm2 Uncert"]/self.data[dt]["Atoms/cm2"],2))
    ratio2 = np.nan_to_num(np.power(self.detector["Del Depth Uncert"]/self.detector["Del Depth"],2))
    self.data[dt]["Atoms/cm3 Uncert"] = self.data[dt]["Atoms/cm3"]*np.sqrt(ratio1 + ratio2)

    np.seterr(**old_settings)

    return


def bin_channels(self, bin_size=21):
    """
    Bin channels from ndp.detector
    """

    # Note that the last bin is not handled correctly as it will not necessarily have
    # the same number of channels as the other bins. To fix this would involve
    # some time and code testing. Assuming that the last bin is never interesting, but
    # perhaps we should eventually fix this just in case.
    #
    # Fix is in recalculating bin_size each time, or at least in the final bin.
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


def chan2depth(self):
    """
    Convert channels to energy and then to a relative depth with uncertainties
    Depth is based on the TRIM simulation but will have the incorrect origin
    """

    # These values change infrequently and are provided by the instrument scientist
    m, b = self.instrument["Calib Coeffs"]
    self.detector["Energy"] = (m*self.detector["Channels"]) + b

    # These values are derived from SRIM/TRIM, freeware used to calculate energy of generated ions in matter
    # Depth is in nanometers
    a, b, c = self.data["TRIM"]["Coeffs"]
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


def bin_channels(self, bin_size=21):
    """
    Bin channels from ndp.detector
    """

    # Note that the last bin is not handled correctly as it will not necessarily have
    # the same number of channels as the other bins. To fix this would involve
    # some time and code testing. Assuming that the last bin is never interesting, but
    # perhaps we should eventually fix this just in case.
    #
    # Fix is in recalculating bin_size each time, or at least in the final bin.
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
            print('Channels')
            columns[i][:] = self.detector['Channels Binned']
        if dkey == 'Energy':
            print('Energy')
            columns[i][:] = self.detector['Energy Binned']
        if dkey == 'Depth':
            print('Depth')
            columns[i][:] = self.detector['Depth Binned']
        if dkey == 'Counts':
            print('Counts')
            columns[i][:] = self.data['Sam Dat']['Counts Binned']
        if dkey == 'Atoms/cm2':
            print('A/cm2')
            columns[i][:] = self.data['Sam Dat']['Atoms/cm2 Binned']
        if dkey == 'Atoms/cm2 Uncert':
            print('A/cm2 un')
            columns[i][:] = self.data['Sam Dat']['Atoms/cm2 Binned Uncert']
        if dkey == 'Atoms/cm3':
            print('A/cm3')
            columns[i][:] = self.data['Sam Dat']['Atoms/cm3 Binned']
        if dkey == 'Atoms/cm3 Uncert':
            print('A/cm3 un')
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
