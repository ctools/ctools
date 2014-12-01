# ==========================================================================
# This script provides a number of functions that are useful for CTA data
# analysis.
#
# Copyright (C) 2014 David Sanchez
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import ctools as ct
import gammalib as gl
from math import sqrt,log10,floor
import sys,string

class base(object):
    '''This class is used to defined the colors of the print''' 
    def __init__(self):
        self.classname = self.__class__.__name__
        self.errorcolor = "\033[31m"#red
        self.infocolor = "\033[34m"#blue
        self.warningcolor = "\033[33m"#yellow
        self.successcolor = "\033[32m"#green
        self.endcolor = "\033[0m"#reset
        self.prependfunction = True
        self.basemembers = ["classname","errorcolor","infocolor","warningcolor","successcolor","endcolor","prependfunction"]
        
    def error(self,message, functionname = ""):
    '''Print the error message'''
        printstring = ""
        if functionname == "":
            printstring = "\n"+self.errorcolor+"*** Error ["+self.classname+"]: "+message+" ***\n"+self.endcolor
        else:
            printstring = "\n"+self.errorcolor+"*** Error ["+self.classname+"::"+functionname+"]: "+message+" ***\n"+self.endcolor
        sys.exit(printstring)
    
    def info(self,message,newline=True):
    '''Print the information message'''
        printstring = self.infocolor+"["+self.classname+"]: "+message+self.endcolor     
        if newline:
            print(printstring)
        else:
            print(self.infocolor+message+self.endcolor)
            sys.stdout.flush()  
        
    def warning(self,message,functionname = ""):*
    '''Print the warning message'''
        printstring = ""
        if functionname == "":
            printstring = self.warningcolor+"["+self.classname+"] Warning: "+message+self.endcolor
        else:
            printstring = self.warningcolor+"["+self.classname+"::"+functionname+"] Warning: "+message+self.endcolor
        print(printstring)
        
    def success(self,message):
    '''Print the success message'''
        printstring = self.successcolor+"["+self.classname+"]: "+message+self.endcolor
        print(printstring)
        
    def progress(self,message="."):
    '''Print the progess message'''
        string = self.infocolor+message+self.endcolor
        print(string)
        
class ResultsStorage(dict,base):
    """class to store the results in the good format"""
    def __init__(self,*arg,**kw):
        super(ResultsStorage, self).__init__(*arg, **kw)
        base.__init__(self)
        # Can store the differential flux, upper limit, time, energy and so one
        self["dnde"]      = {}
        self["ulim_dnde"] = []
        self["time"]      = {}
        self["ener"]      = {}
        self["TS"]        = []
        self["Iflux"]     = {}
        
        #store the time : Min, Max and units of the variable
        self["time"]["tmin_value"] = []
        self["time"]["tmax_value"] = []
        self["time"]["unit"]       = "sec"
        
        #store the energy : Value min,max (ed_* and eu_*) and units of the variable
        self["ener"]["value"]      = []
        self["ener"]["ed_value"]   = []
        self["ener"]["eu_value"]   = []
        self["ener"]["unit"]       = "MeV"
        
        #store the differential flux : Value min,max (ed_* and eu_*) and units of the variable
        self["dnde"]["value"]      = []
        self["dnde"]["ed_value"]   = []
        self["dnde"]["eu_value"]   = []
        self["dnde"]["unit"]       = "ph/cm2/s/MeV"
        
        #store the integral flux : Value min,max (ed_* and eu_*) and units of the variable
        self["Iflux"]["value"]     = []
        self["Iflux"]["ed_value"]  = []
        self["Iflux"]["eu_value"]  = []
        self["Iflux"]["unit"]      = "ph/cm2/s"
        
    def Print(self):
    """ Print the value stored in the class"""
        self.success("Results are :")
        for k in self.iterkeys():
            self.info("Key name : "+k)
            try :
                for ki in self[k].iterkeys():
                    print("\t",ki,"\t",self[k][ki])
            except:
                print(self[k])
        print()
        
class Analyser(base):
    '''Class to run the analysis of the data. Currently the parameters are hard coded but this will change once
        the method to feed the class will be defined
        The classe can work with GObservation as input (simulation) or fits files (simulation+real data)'''
    def __init__(self):
        super(Analyser,self).__init__()
        self.info("Creating "+self.classname+" object")
        #Observation and ctlike object
        self.m_obs = None
        self.like  = None
        
        #parameters
        self.m_ebinalg  = "LOG"
        self.m_eunit    = "TeV"
        self.m_enumbins = 10#nbins
        self.m_usepnt   = True # Use pointing for map centre
        self.m_coordsys = "CEL"#coord
        self.m_proj     = "TAN"#proj
        self.m_binned   = False
        self.m_edisp    = False
        self.m_stat     = "POISSON"#stat

        # booleans to ensure that input files and parameters have been defined
        # set to false at the begining and turn to true while a file/parameter is defined.
        self.m_source_set = False
        self.m_files_set  = False
        self.m_time_set   = False
        self.m_energy_set = False
        self.m_irfs_set   = False
        self.m_roi_set    = False
        self.m_xml_set    = False
        
    def set_xml(self,xml):
        '''Set the model XML file'''
        self.m_xml     = xml
        ind = xml.find(".xml")
        if not ind==-1:
            self.m_xml_out = xml[:xml.find(".xml")]+"_results"+xml[xml.find(".xml"):]
        self.m_xml_set = True
        
    def set_roi(self,roi):
        '''Set the ROI parameter and compute the number of pixel, bin size of the map'''
        self.m_roi      = float(roi)
        self.m_nxpix    = int(self.m_roi*100) #nxpix
        self.m_nypix    = int(self.m_roi*100) #nypix
        self.m_binsz    = 0.05                #binsz
        self.m_roi_set  = True

    def set_irfs(self,caldb,irf):
        '''Set the IRFs definition'''
        self.m_caldb    = caldb
        self.m_irf      = irf
        self.m_irfs_set = True

    def set_Fitsfiles(self,evtfile,tag="CTA_analysis"):
        '''Set the input and output fits files'''
        self.m_evtfile   = "event_selected_"+tag+".fits"    
        self.m_raw_evt   = evtfile
        self.m_cntfile   = "countmap_"+tag+".fits"
        self.m_files_set = True

    def set_source(self,name,ra,dec):
        '''Set the source name and the RA/DEC position of the source'''
        self.m_ra         = ra
        self.m_dec        = dec
        self.m_name       = name
        self.m_source_set = True

    def set_time_boundary(self,tmin,tmax):
        """Set the start and stop time of the observation"""
        self.m_tmin        = float(tmin)
        self.m_tmax        = float(tmax)
        self.m_time_set    = True

    def set_energy_boundary(self,emin,emax):
        """Set the energy bounds of the observation"""
        self.m_emin = emin
        self.m_emax = emax
        ebounds = gl.GEbounds(gl.GEnergy(emin, self.m_eunit), \
                                gl.GEnergy(emax, self.m_eunit))

        if self.m_obs:
            for obs in self.m_obs:
                obs.events().ebounds(ebounds)
        self.m_energy_set = True

    def validate(self):
        '''Check is all the parameters has been set'''
        if not self.m_time_set :
            self.error("Time range of the observations not set")

        if not self.m_files_set :
            self.error("FITS files of the observations not set")

        if not self.m_source_set :
            self.error("Source not set")

        if not self.m_energy_set :
            self.error("Energy range of the observations not set")
            
        if not self.m_irfs_set :
            self.error("IRFs to use not set")

        if not self.m_roi_set :
            self.error("ROI to use not set")

        if not self.m_xml_set :
            self.error("XML files for the sky model not set")

    def copy(self,analyse,tag="CTA_analysis_copy"):
        '''Make a copy of the class'''
        self.set_xml(analyse.m_xml)
        self.set_irfs(analyse.m_caldb,analyse.m_irf)
        self.set_roi(analyse.m_roi)
        self.set_time_boundary(analyse.m_tmin,analyse.m_tmax)
        self.set_energy_boundary(analyse.m_emin,analyse.m_emax)
        self.set_source(analyse.m_name,analyse.m_ra,analyse.m_dec)
        self.set_Fitsfiles(analyse.m_evtfile,tag)


    def set_obs(self,obs):
        '''set the GObservation container'''
        self.m_obs = obs 

    def ctselect(self):
        '''ctselect application and set parameters'''
        self.info("Running ctselect to cut on events")
        self.validate()# check the parameters
        
        if self.m_obs:
            filter = ct.ctselect(self.m_obs)
        else:
            filter = ct.ctselect()
        
        #set the parameters for running the tool
        filter["infile"] = self.m_raw_evt               
        filter["outfile"] = self.m_evtfile  
        filter["usepnt"].boolean(self.m_usepnt) # Use pointing for map centre
        filter["ra"]  = self.m_ra
        filter["dec"]  = self.m_dec
        filter["rad"]  = self.m_roi
        filter["tmin"] = self.m_tmin
        filter["tmax"] = self.m_tmax
        filter["emin"].real(self.m_emin)
        filter["emax"].real(self.m_emax)
        filter["expr"] = ""
        filter.logFileOpen()

        filter.run()
        # if no Gobservation have been defined, then save the fits file.
        if not(self.m_obs):
            filter.save()

        if self.m_obs:
            # Make a deep copy of the observation that will be returned
            # (the ctbin object will go out of scope one the function is
            # left)
            self.m_obs = filter.obs().copy()
        
    def ctbin(self,log=False,debug=False):
        '''ctbin application and set parameters'''
        self.info("Running ctbin to create count map")
        self.validate()# check the parameters

        if self.m_obs:
            bin = ct.ctbin(self.m_obs)
        else:
            bin = ct.ctbin()
            bin["evfile"] = self.m_evtfile

        #set the parameters for running the tool
        bin["outfile"] = self.m_cntfile
        bin["ebinalg"].string(self.m_ebinalg)
        bin["emin"].real(self.m_emin)
        bin["emax"].real(self.m_emax)
        Nbdecade = log10(self.m_emax)-log10(self.m_emin)#Compute the number of decade
        bin["enumbins"].integer(int(self.m_enumbins*Nbdecade))
        bin["usepnt"].boolean(self.m_usepnt) # Use pointing for map centre
        bin["nxpix"].integer(self.m_nxpix)
        bin["nypix"].integer(self.m_nypix)
        bin["binsz"].real(self.m_binsz)
        bin["coordsys"].string(self.m_coordsys)
        bin["proj"].string(self.m_proj)
        
        # Optionally open the log file
        if log:
            bin.logFileOpen()

        # Optionally switch-on debugging model
        if debug:
            bin["debug"].boolean(True)

        # Run ctbin application. This will loop over all observations in
        # the container and bin the events in counts maps
        bin.execute()
        if self.m_obs:
            # Make a deep copy of the observation that will be returned
            # (the ctbin object will go out of scope one the function is
            # left)
            self.m_obs = bin.obs().copy()
            
    def create_fit(self,log=False,debug=False):
        '''Create the like object from ctlike'''
        # create ctlike instance with given parameters
        self.info("Fitting Data using ctlike")
        self.validate()
        
        #set the parameters of the tools
        if self.m_obs:
            self.like = ct.ctlike(self.m_obs)
        else:
            self.like = ct.ctlike()
            if self.m_binned#binned or unbinned analysis
                self.like["infile"] = self.m_cntfile
            else:
                self.like["infile"] = self.m_evtfile
                
            self.like["stat"].string(self.m_stat)
            self.like["caldb"].string(self.m_caldb)
            self.like["irf"].string(self.m_irf)
        self.like["srcmdl"] = self.m_xml 
        self.like["edisp"].boolean(self.m_edisp)
        self.like["outmdl"] = self.m_xml_out

        # Optionally open the log file
        if log:
            self.like.logFileOpen()
            
        # Optionally switch-on debugging model
        if debug:
            self.like["debug"].boolean(True)
            
    def compute_butterfly(self,name,log=False,debug=False):
        '''create ctlike instance with given parameters'''
        self.info("Computing the data with ctbutterfly")
        self.validate()
        
        if self.m_obs:
            butterfly = ct.ctbutterfly(self.m_obs)
        else:#parameters of the tools
            butterfly = ct.ctbutterfly()
            butterfly["infile"] = self.m_evtfile
            butterfly['srcname'] = name
            butterfly["emin"].real(self.m_emin)
            butterfly["emax"].real(self.m_emax)
            butterfly["caldb"].string(self.m_caldb)
            butterfly["irf"].string(self.m_irf)
            butterfly["srcmdl"] = self.m_xml
            filename = name+"_butterfly.txt"
            butterfly["outfile"] = filename

        # Optionally open the log file
        if log:
            butterfly.logFileOpen()
        # Optionally switch-on debugging model
        if debug:
            butterfly["debug"].boolean(True)
        
        butterfly.run() #Compute the butterfly
        butterfly.save() # save the txt files

        #Store the results in the Storage class
        self.m_butterfly = ResultsStorage()
        data = open(filename,'r').readlines()
        for line in data:
            words = string.split(line)
            self.m_butterfly["ener"]["value"].append(float(words[0]))
            self.m_butterfly["dnde"]["value"].append(float(words[1]))
            self.m_butterfly["dnde"]["ed_value"].append(float(words[1]) - float(words[2]))
            self.m_butterfly["dnde"]["eu_value"].append(float(words[1]) + float(words[2]))
        
    def fit(self,log=False,debug=False):
        self.validate()
        if not(self.like):
            self.warning("ctlike object not created, creating now")
            self.create_fit(log,debug)

        # Run ctlike application.
        self.like.run()
        # Save the results in XML
        self.like.save()
        self.success("Fit performed")
        
    def PrintResults(self,srcname = ""):
        ''' Print results of the fit '''
        self.info("Results of the Fit")
        for m in self.like.obs().models():
            if srcname == m.name() or srcname=="":
                print("Model : "+m.name())
                print(m)
                
    def GetSrcResuls(self,Res,srcname,E0=None,factor = 1e6):
        """function to get and store the results of a fit. Store the results in the Storage class element"""
        results = self.GetSrcParams(srcname)

        #Store the time
        Res["time"]["tmin_value"].append(self.m_tmin)
        Res["time"]["tmax_value"].append(self.m_tmax)
        Res["time"]["unit"] = "sec"

        #Store the energy
        if (E0 is not None):#if a 'pivot' energy is given
            Res["ener"]["value"].append(E0*factor)
            Res["ener"]["ed_value"].append((E0-self.m_emin)*factor)
            Res["ener"]["eu_value"].append((self.m_emax-E0)*factor)
            Res["ener"]["unit"] = "MeV"
        else:#if not, the geometric mean energy is used
            Res["ener"]["value"].append(sqrt(self.m_emin*self.m_emax))
            Res["ener"]["ed_value"].append(-self.m_emin+sqrt(self.m_emin*self.m_emax))
            Res["ener"]["eu_value"].append(self.m_emax-sqrt(self.m_emin*self.m_emax))
            Res["ener"]["unit"] = self.m_eunit

        #save the flux values
        Res["Iflux"]["value"].append(results["flux"])
        Res["Iflux"]["eu_value"].append(results["eflux"])
        Res["Iflux"]["ed_value"].append(results["eflux"])
        
        #save the differential flux in a src model independent maner
        srcmodel = self.GetModel(srcname)
        Res["dnde"]["value"].append(srcmodel.spectral()[0].value())
        Res["dnde"]["eu_value"].append(srcmodel.spectral()[0].error())
        Res["dnde"]["ed_value"].append(srcmodel.spectral()[0].error())

        Res["TS"].append(100) ## TODO use the real TS value

        UL = UpperLimitComputer(self,srcname) #Compute the UL. To be updated
        UL.run()
        Res["ulim_dnde"].append(0)

        return Res

    def GetSrcParams(self,srcname):
        '''Retrive the source parameters from the like object and store it in a dictonnary'''
        results = {}
        found = False
        m = self.GetModel(srcname)
        
        #get the integral flux
        results["flux"] = m.spectral().flux(gl.GEnergy(self.m_emin,self.m_eunit),gl.GEnergy(self.m_emax,self.m_eunit))
        results["eflux"] = m.spectral().eflux(gl.GEnergy(self.m_emin,self.m_eunit),gl.GEnergy(self.m_emax,self.m_eunit))

        # Get free parameters
        for par in m.spectral():
            if par.is_free():
                results[par.name()] = par.value()
                results["e"+par.name()] = par.error()

        return results

    def GetModel(self,srcname):
        ''' Get the spectral parameters of the sources stored in the like object''' 
        found = False
        for m in self.like.obs().models():
            if m.name() == srcname:
                return m
        self.error("source "+srcname+" not in the source list")
    
def MakeEbin(analyse,nbins,srcname,binned = False):
    #Store the results
    Res = ResultsStorage()
    
    #check the source name, exit is the srcname is false
    analyse.GetModel(srcname)
    
    #compute the conversion factor to have energy in MeV
    Energy_list = {"MeV":1,"GeV":1e3,"TeV":1e6} #TODO
    factor = Energy_list[analyse.m_eunit]
    
    for i in xrange(nbins):
        obs_copy = analyse.like.obs().copy()
        emin = pow(10,log10(analyse.m_emin) + (log10(analyse.m_emax)-log10(analyse.m_emin))/(nbins)*i)
        emax = pow(10,log10(analyse.m_emin) + (log10(analyse.m_emax)-log10(analyse.m_emin))/(nbins)*(i+1))
        E0 = pow(10, (log10(emin) + log10(emax)) / 2)
        print("Making energy bin between %2.1e"%(emin)+analyse.m_eunit+" and %2.1e"%(emax)+analyse.m_eunit)

        # change the PivotEnergy to E0
        m = obs_copy.models()[srcname]
        for par in m.spectral():
            if par.name() == "PivotEnergy":
                par.value(E0*factor)
                            
        ana_bin = Analyser()

        #copy the obs object
        ana_bin.set_obs(obs_copy)

        #copy the analyser parameters
        ana_bin.copy(analyse)

        # Set energy boundaries
        ana_bin.set_energy_boundary(emin,emax)
        ana_bin.ctselect()

        # if binned:
            # ana_bin.ctbin()

        # Create the fit object
        ana_bin.create_fit(log=False,debug=False)
        
        # Run ctlike application.
        ana_bin.like.run()

        #Retrive the results in the good format
        Res = ana_bin.GetSrcResuls(Res,srcname,E0,factor)
        
        del ana_bin
    return Res
    
def MakeLC(analyse,nbins,srcname,binned = False):
    #Store the results
    Res = ResultsStorage()

    #check the source name, exit is the srcname is false
    analyse.GetModel(srcname)
    
    for i in xrange(nbins):
        obs_copy = analyse.like.obs().copy()
        tmin = analyse.m_tmin + (analyse.m_tmax-analyse.m_tmin)/(nbins)*i
        tmax = analyse.m_tmin + (analyse.m_tmax-analyse.m_tmin)/(nbins)*(i+1)
        print("Making time bin between %2.1e"%(tmin)+" and %2.1e"%(tmax))
                            
        ana_bin = Analyser()

        #copy the obs object
        ana_bin.set_obs(obs_copy)
        
        #copy the analyser parameters
        ana_bin.copy(analyse)

        # Set energy boundaries
        ana_bin.set_time_boundary(tmin,tmax)
        ana_bin.ctselect()

        # if binned:
            # ana_bin.ctbin()

        # Create the fit object
        ana_bin.create_fit(log=False,debug=False)
        
        # Run ctlike application.
        ana_bin.like.run()

        #Retrive the results in the good format
        Res = ana_bin.GetSrcResuls(Res,srcname,tmin=tmin,tmax=tmax)
        
        del ana_bin
    return Res
        
class UpperLimitComputer(base):
    '''This class compute the ul and will be removed'''
    def __init__(self,analyser,srcname,parname = "Prefactor"):
        super(UpperLimitComputer,self).__init__()
        self.info("Creating "+self.classname+" object")

        obs = analyser.like.obs().copy()
        self.like = ct.ctlike(obs)
        self.succes = False

        # get spectral model
        self.model = self.like.obs().models()[srcname]
        self.spec = self.model.spectral()
        self.parname = parname
        self.bestloglike = analyser.like.opt().value()
        self.dloglike = 3.84/2.0 #95% CL now

    def _bisec(self,a,b,tol = 1e-12):
        """iterative function to find the roots"""
        if (self._func(a)==0.):
            return a
        elif (self._func(b)==0):
            return b
        elif (self._func(b)*self._func(a)>0.):
            self.warning("Finding root for upper limit calculation failed")
            return 0
        mid = (a+b)/2
        err = abs(a-b)/2
        fm = self._func(mid)
        if (err<tol or self._func(mid)==0.):
            self.succes = True
            return mid
        elif (fm*self._func(a)<0):
            b = mid
        else :
            a = mid
        return self._bisec(a,b)
           
    def _func(self,value): 
        """function to be minimized"""
        # get parmeter value
        val = value*self.spec[self.parname].scale()
        
        # ensure value to be above minimum boundary
        if val <= self.spec[self.parname].min():
            return 0.0-self.self.bestloglike-self.dloglike
        
        # set value
        self.spec[self.parname].value(val)
        
        # Recalculate likelihood
        self.like.run()
        logl = self.like.opt().value()   
        
        # Return likelihood difference
        return (logl-self.bestloglike-self.dloglike)
    
    def run(self):
        
        # get spectral parameters and set start value
        value = self.spec[self.parname].value()
        error = self.spec[self.parname].error()
        scale = self.spec[self.parname].scale()

        #Fix all parameters
        for model in self.like.obs().models():
            for par in model:
                par.fix()
                
        
        # set start value to 2 sigma above fit value
        self.startpar = value / scale + 2* error/scale

        # set start (a) and (b) depend on the CL asked
        a= value / scale +  (self.dloglike/4.)*error/scale
        b=value / scale + (self.dloglike*2)* error/scale
        ul = self._bisec(a,b)
        
        self.ulimit = ul*self.spec[self.parname].scale()
        # Test for convergence
        if self.succes:
            self.success("Upper limit of parameter \""+self.parname +"\"  determined to "+str(self.ulimit))
