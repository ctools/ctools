#!/usr/bin/env python

import gammalib
import sys
import os
import glob
import json

def print_help():
    print("")
    print("tsmerge.py expects two keyword argumens: use via command line as follows")
    print("tsmerge.py files=<list of files> outfile=<desired outfile name>")
    print("-------------------------------------------------------------------------")
    print("*** \"outfile\" should be a path to a fits file")
    print("*** \"files\" should be either be:")
    print("     - an ascii-file listing the files your want to merge")
    print("  OR - a wildcard string describing the files for merging")
    print("-------------------------------------------------------------------------")
    print("")
    print("Example: python tsmerge,py files=tsmap*.fits outfile=tsmap.fits")
    print("")
    sys.exit(-1)


# Maps class holding the different sky maps 
# created by cttsmap
class maps:
    
    # Constructor
    def __init__(self,fitsfile):
        
        self.filename = fitsfile
        if not os.path.isfile(fitsfile):
            self.error("Maps instance cannot be constructed: File \""+self.filename+"\" does not exist")
        
        # open FITS file
        fits = gammalib.GFits(fitsfile)
        
        # get TS map
        self.tsmap = gammalib.GSkymap()
        self.tsmap.read(fits[0])
        
        # Get status map
        if not fits.contains("STATUS MAP"):
            self.error("No extension \"STATUS MAP\" found in \""+self.filename+"\"")
        self.statusmap = gammalib.GSkymap()
        self.statusmap.read(fits["STATUS MAP"])
              
        # Get other maps 
        self.maps = []
        self.mapnames = []
        for hdu in fits:
            if hdu.extname() != "IMAGE" and hdu.extname() != "STATUS MAP":
                skymap = gammalib.GSkymap()
                skymap.read(hdu)
                self.maps.append(skymap)
                self.mapnames.append(hdu.extname())

    # len operator
    def __len__(self):
        return len(self.maps)

    # operator []
    def __getitem__(self,index):
        return self.maps[index]

    # Check the status map for completeness
    def check_status(self):
        check = True
        for pix in self.statusmap:
            if pix < 0.5:
                check = False
                break
        return check
    
    # Save maps to outfile
    def save(self,outfile):
        
        fits = gammalib.GFits()
        self.tsmap.write(fits)
        
        for i in range(len(self.maps)):
            self.maps[i].write(fits)
        
        for i in range(len(self.mapnames)):   
            fits[i+1].extname(self.mapnames[i])
        
        if not self.check_status():
            self.statusmap.write(fits)
            fits[fits.size()-1].extname("STATUS MAP")
        fits.saveto(outfile,True)
           
    # Add map instances to each other
    # in this way the TS maps get merged  
    def add(self,addmaps):
        if not len(addmaps) == len(self.maps):
            self.error("Maps cannot be added: Maps from \""+self.filename+"\" do not match from \""+addmaps.filename+"\"")
            
        # Loop over bins    
        for i in range(self.tsmap.npix()):
            if addmaps.statusmap[i] > 0.5:
                self.tsmap[i] = addmaps.tsmap[i]
                self.statusmap[i] = addmaps.statusmap[i]
                for j in range(len(self.maps)):
                    self.maps[j][i] = addmaps[j][i]
    # Unify error messages
    def error(self,message):  
        print("*** ERROR ***: "+message)
        sys.exit(-1)
        
        
def pars_from_argv(argv):

    if argv is None:
        argv = sys.argv

    args = []
    kwargs = {}
    for s in argv[1:]:
        if s.count('=') == 1:
            key, value = s.split('=', 1)
        else:
            key, value = None, s
        try:
            value = json.loads(value) 
        except ValueError:
            pass
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return args, kwargs


def get_kwargs(argv):
    kwargs = {}

    for s in argv:
        if s.count('=') == 1:
            key, value = s.split('=', 1) 
            kwargs[key] = value
            
        else:
            print_help()         
            
    return kwargs   

if __name__ == "__main__":
    # Initialise list of files

    kwargs = get_kwargs(sys.argv[1:])
    if "outfile" not in kwargs or "files" not in kwargs:
        print_help()
    outfile = kwargs["outfile"]
    files = kwargs["files"]
    allfiles = []
    if os.path.isfile(files):
        allfiles = open(files).read().splitlines()   
    else:
        allfiles = glob.glob(files)
        
    # Initialise files to be merged
    mergefiles = []
    
    # Test files for the entry status map
    # and add them to the file list
    for fitsfile in allfiles:    
        fits = gammalib.GFits(fitsfile)
        if fits.contains("STATUS MAP"):
            mergefiles.append(fitsfile)
        fits.close()
        
    # create a maps instance from the first file
    finalmaps = maps(mergefiles[0])
    
    # add other files to the maps instance
    for fitsfile in mergefiles[1:]:
        finalmaps.add(maps(fitsfile))
    
    # save the filled maps to file
    finalmaps.save(outfile=outfile)
    
