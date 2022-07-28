# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:50:16 2022

@author: rljones
"""

def load_NISTNDP(self, dt):
    """
    Function to load a list of NDP data files of a given datatype (Sam Dat, Sam Mon, etc),
    load header info and sum of counts/channel into ndp.data
    
    """

    filelist = self.data[dt]["Files"]
    numfiles = len(filelist)
    numchannels = self.instrument["Num Channels"]
    
    if("Channel Sum" not in self.data[dt]["Operations"]):
            self.data[dt]["Operations"].append("Channel Sum")
    
    for filenum in range(numfiles):
        ndp_file = filelist[filenum]
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

    return