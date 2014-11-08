# ==========================================================================
# This script provides a number of functions that are useful for plotting.
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
from math import log10,pow,sqrt
import os

try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        has_matplotlib = True
except:
        has_matplotlib = False

try:
        import aplpy
        has_aplpy = True
except:
        has_aplpy = False

class options(object):
    def __init__(self):
        self.fmt         = 'o'
        self.color       = BLUE
        self.linecolor   = 'black'
        self.markercolor = self.color
        self.tslim       = 4.0
        self.elinewidth  = 2.0
        self.linewidth   = 2.0
        self.linestyle   = '--'
        self.capsize     = 5.0
        self.markersize  = 5.0
        self.limlength   = 0.4 
        self.label       = "_nolegend_"

import analysisutils 

class SpectralBase(object):
    def __init__(self,sed=True,sed_factor=1.):
        self.m_sed_factor = sed_factor
        self.m_sed        = sed
        
    def _convertSED(self,tab,energy):
        """convert the a table storing dN/dE in sed with right units"""

        if len(tab) != len(energy):
            self.error("input table and energy table does not match")
            return
        if self.m_sed:
            for i in xrange(len(tab)):
                tab[i] *= energy[i]**2*self.m_sed_factor
        return tab

BLUE="#337EF5"
GREY='grey'
GREEN = '#21851D'
RED = '#BD0B0B'
YELLOW = '#D4CF44'
ORANGE = '#D19534'


class residuals(analysisutils.base,options,SpectralBase):
    def __init__(self,spectrum,datapoint,sed,sed_factor,eunit="TeV"):
        super(residuals,self).__init__()
        options.__init__(self)
        SpectralBase.__init__(self,sed,sed_factor)

        # load the gammalib spectrum object
        self.spectrum  = spectrum
        self.SetDataPoint(datapoint)
        
        self.xerrors = True
        self.refUnit = eunit     # reference unit for conversion
        
    def SetDataPoint(self,datapoint):
        self.m_point_info = datapoint
        self.Npt = len(self.m_point_info["dnde"]["value"])

    def _plt_residuals(self):
        th_val = []
        index  = []
        ed_y   = []
        eu_y   = []
        y      = []
        
        Energy_list = {"MeV":1,"GeV":1e-3,"TeV":1e-6} #TODO
        for i in range(self.Npt):
            if self.m_point_info["TS"][i] >= self.tslim:
                th_val.append(self.spectrum.eval(gl.GEnergy(self.m_point_info["ener"]["value"][i],self.m_point_info["ener"]["unit"]),gl.GTime(0)))
                index.append(i)
            else:
                th_val.append(0)
        if self.m_sed:
            self._convertSED(th_val,self.m_point_info["ener"])
        for i in index:
                y.append((self.m_point_info["dnde"]["value"][i] - th_val[i]) / th_val[i])
                ed_y.append(self.m_point_info["dnde"]["ed_value"][i]/th_val[i])
                eu_y.append(self.m_point_info["dnde"]["eu_value"][i]/th_val[i])
                
        if self.xerrors:
            plt.errorbar(self.m_point_info["ener"]["value"], y, xerr=[self.m_point_info["ener"]["ed_value"],self.m_point_info["ener"]["ed_value"]],yerr=[ed_y,eu_y],fmt=self.fmt,color=self.markercolor,elinewidth=self.elinewidth) 
        else:
            plt.errorbar(self.m_point_info["ener"]["value"], y, yerr=[ed_y,eu_y],fmt=self.fmt,color=self.markercolor,elinewidth=self.elinewidth) 
        
        plt.axhline(0.0,color=self.linecolor,lw=self.linewidth,ls=self.linestyle)


class ulimgraph(analysisutils.base,options,SpectralBase):
    def __init__(self,datapoint,sed,sed_factor,eunit='TeV'):
        super(ulimgraph,self).__init__()
        #~ super(ulimgraph,self).__init__()
        options.__init__(self)
        SpectralBase.__init__(self,sed,sed_factor)
        
        self.xerrors=True  
        self.nlims = -1
        self.refUnit = eunit
        self.SetDataPoint(datapoint)
    
    def SetDataPoint(self,datapoint):
        self.m_point_info = datapoint
        self.Npt = len(self.m_point_info["dnde"]["value"])
       
    def _plt_points(self):
        results = analysisutils.ResultsSorage()
        
        Energy_list = {"MeV":1,"GeV":1e-3,"TeV":1e-6} #TODO
        for i in range(self.Npt):

            self.m_point_info["ener"]["value"][i]*=Energy_list[self.refUnit]
            self.m_point_info["ener"]["ed_value"][i]*=Energy_list[self.refUnit]
            self.m_point_info["ener"]["eu_value"][i]*=Energy_list[self.refUnit]
            self.m_point_info["ener"]["unit"] = self.refUnit

            if self.m_point_info["TS"][i] >= self.tslim:

                if self.m_point_info["dnde"]["ed_value"][i]>self.m_point_info["dnde"]["value"][i]:
                    self.warning("\tflux error: "+str(self.m_point_info["dnde"]["ed_value"][i])+" is larger than flux value "+str(self.m_point_info["dnde"]["ed_value"][i]))
                    self.warning("\tTS value is however TS="+str(self.m_point_info["TS"][i])+", weird!")
            else :
                self.m_point_info["dnde"]["ed_value"][i] = 0
                self.m_point_info["dnde"]["eu_value"][i] = 0
            
        if self.xerrors:
            plt.errorbar(self.m_point_info["ener"]["value"], self.m_point_info["dnde"]["value"], xerr=[self.m_point_info["ener"]["ed_value"], self.m_point_info["ener"]["eu_value"]],yerr=[self.m_point_info["dnde"]["ed_value"], self.m_point_info["dnde"]["eu_value"]],fmt=self.fmt,color=self.markercolor,elinewidth=self.elinewidth) 
        else:
            plt.errorbar(self.m_point_info["ener"]["value"], self.m_point_info["dnde"]["value"],yerr=[results["dnde"]["ed_value"],results["dnde"]["eu_value"]],fmt=self.fmt,color=self.markercolor,elinewidth=self.elinewidth) 
        
    def _plt_limits(self):
        ed_flux = []
        eu_flux = []
        flux    = []
        ener    = []
        
        drawn_lims = 0
        Energy_list = {"MeV":1,"GeV":1e-3,"TeV":1e-6} #TODO
        for i in range(self.Npt):
            if self.m_point_info["TS"][i] < self.tslim:
                self.info("Plot an upper limit for E = "+self.m_point_info["ener"]["value"][i]+" with a TS = "+self.m_point_info["TS"][i])
                if drawn_lims >= self.nlims and self.nlims!=-1:
                    break

                ed_flux.append(self.limlength*self.m_point_info["ulim_dnde"][i])
                eu_flux.append(0.0)
                flux.append(self.m_point_info["ulim_dnde"][i])
                drawn_lims+=1
                
        # if self.xerrors:
            # plt.errorbar(self.m_point_info["ener"]["value"], self.m_point_info["dnde"]["value"], xerr=[self.m_point_info["ener"]["ed_value"], self.m_point_info["ener"]["eu_value"]],yerr=[self.m_point_info["dnde"]["ed_value"],self.m_point_info["dnde"]["eu_value"]],fmt=self.fmt,markersize=0.0,elinewidth=self.elinewidth,lolims=True,capsize=self.capsize,label='_nolegend_',mfc = self.color,ecolor=self.color,ms =self.markersize)      
        # else:
        plt.errorbar(ener, flux, yerr=[ed_flux, eu_flux],fmt=self.fmt,markersize=0.0,elinewidth=self.elinewidth,lolims=True,capsize=self.capsize,label='_nolegend_',mfc = self.color,ecolor=self.color,ms =self.markersize)      
            

class LightCurvePlotter(analysisutils.base,options):
    """ Class to plot the Light curve"""
    def __init__(self,srcname,datapoint):
        super(LightCurvePlotter,self).__init__()
        options.__init__(self)
        self.m_name = srcname
        self.m_data = datapoint #data point in a Results Storage class

    def draw(self):
        if not(has_matplotlib):
            self.warning("matplotlib module not found, can draw")
            return

        time     = []
        dtime    = []
        flux     = []
        dflux_ed = []
        dflux_eu = []
        for i in range(len(self.m_data["time"]["tmin_value"])):
            time.append((self.m_data["time"]["tmin_value"][i]+self.m_data["time"]["tmax_value"][i])/2)
            dtime.append((self.m_data["time"]["tmax_value"][i]-self.m_data["time"]["tmin_value"][i])/2)
            flux.append(self.m_data["Iflux"]["value"][i])
            dflux_eu.append(self.m_data["Iflux"]["eu_value"][i])
            dflux_ed.append(self.m_data["Iflux"]["ed_value"][i])

        plt.figure("results",figsize=(12,7),edgecolor="w")
        unit = self.m_data["time"]["unit"]
        plt.xlabel("Time ["+ unit+"]",fontsize=15)
        unit = self.m_data["Iflux"]["unit"]
        plt.ylabel("E$^{2}$ dN/dE ["+unit+"]",fontsize=15)

        print(time)
        print(dtime)
        print(flux)
        print(dflux_ed)
        print(dflux_eu)
        plt.errorbar(time,flux,xerr=[dtime,dtime],yerr=[dflux_ed,dflux_eu],fmt=self.fmt,color=self.markercolor,elinewidth=self.elinewidth) 
        plt.show()

class SpectrumPlotter(analysisutils.base,SpectralBase):
    """ Class to compute and plot the spectra (i.e best fit model and butterfly).
        if data points are provided, only change applied  is energy unit convertion to match the butterfly one
        (i.e. no dN/dE to SED conversion).
        It is currently assumed that the unit is MeV for the energy of the data point. this will be change 
        once the I/O container format will be defined."""
    def __init__(self,srcname,like,datapoint,emin=0.1,emax=100,energy = "TeV",SED = True ,npt = 50):
        super(SpectrumPlotter,self).__init__()

        self.m_ebound   = [emin,emax]
        self.m_name     = srcname
        self.m_like     = like
        self.m_spectral = like.obs().models()[srcname].spectral()
        self.m_energy   = []
        self.m_flux     = []
        self.m_error    = []
        self.m_but      = []
        self.m_enebut   = []
        self.m_npt      = npt
        self.m_covar    = gl.GMatrix(1,1)
        self.covar()
        #~ this part take care of the style : SED or differential flux
        #~ in the case of SED, the y-unit is erg/cm2/s
        self.m_sed = SED
        self.eunit = energy
        self.m_sed_factor = 1.
        self._validateUnits()
        SpectralBase.__init__(self,self.m_sed,self.m_sed_factor)

        #Assume MeV for the energy #TODO
        self.points = ulimgraph(datapoint,self.m_sed,gl.MeV2erg,energy)
        self.residuals = residuals(self.m_spectral,datapoint,self.m_sed,gl.MeV2erg,energy)

    def _validateUnits(self):
        Energy_list = {"MeV":1,"GeV":1e3,"TeV":1e6}
        if not(Energy_list.has_key(self.eunit)):
            self.warning("Change energy unit to TeV")
            self.eunit = "TeV"
        
        if self.m_sed:
            self.m_sed_factor = gl.MeV2erg*Energy_list[self.eunit]**2

    def covar(self):
        # Get covariance matrix if fit has converged
        try:
            fullcovar = self.m_like.obs().function().curvature().invert()
        except:
            self.warning("Covariance matrix not determined successfully")
            return
        #~ get the matrix element of the source of interest
        #~ 1) get the right element indices
        idx = 0
        par_index_map = []
        for m in self.m_like.obs().models() :
            if m.name()!= self.m_name:
                idx += m.size()
                continue
            idx += m.spatial().size()
            for par in m.spectral():
                if par.is_free():
                    par_index_map.append(idx)
                    idx += 1

        #~ 2) get the elemets and store them in a matrix
        self.m_covar = gl.GMatrix(len(par_index_map),len(par_index_map))
        i = 0
        for xpar in par_index_map:
            j = 0
            for ypar in par_index_map:
                self.m_covar[i,j] = fullcovar[xpar,ypar]
                j += 1
            i += 1
        self.success("Covariance Matrix succesfuly computed")

    def _makespectrum(self):
        for i in xrange(self.m_npt):
            #Compute the energy
            lene = log10(self.m_ebound[0])+log10(self.m_ebound[1])*i/(self.m_npt-1)
            self.m_energy.append(pow(10,lene))
            #Compute the flux
            self.m_flux.append(self.m_spectral.eval(gl.GEnergy(pow(10,lene),self.eunit),gl.GTime(0)))

            #Compute the gradients
            if self.m_covar.size() != 1:
                self.m_spectral.eval_gradients(gl.GEnergy(pow(10,lene), self.eunit),gl.GTime(0))
                #store them into a GVector
                Derivative = gl.GVector(self.m_covar.size()/2)
                j = 0
                for par in self.m_spectral:
                    if par.is_free():
                        Derivative[j] = par.factor_gradient()
                    j += 1
                #~ computed the error
                self.m_error.append(sqrt(Derivative * (self.m_covar*Derivative)))

        #~ convert the flux in sed if asked
        self.m_flux = self._convertSED(self.m_flux,self.m_energy)
        if self.m_covar.size() != 1:
            self.m_error = self._convertSED(self.m_error,self.m_energy)
            self._makebutterfly() # make the butterfly
        self.success("Spectrum computed")

    def _makebutterfly(self):
        """make the butterfly by appending element in a table"""
        for i in xrange(self.m_npt):
            self.m_but.append(self.m_flux[i]+self.m_error[i])
            self.m_enebut.append(self.m_energy[i])
        for i in xrange(self.m_npt):
            idx = self.m_npt-i-1
            self.m_but.append(self.m_flux[idx]-self.m_error[idx])
            self.m_enebut.append(self.m_energy[idx])
        self.m_but.append(self.m_but[0])
        self.m_enebut.append(self.m_enebut[0])

    def write(self):
        self.warning("write function not yet implemented")

    def draw(self):
        """plot the results, provide a file named ptfile for the flux points"""
        gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
        self._makespectrum()
        self.write()

        if not(has_matplotlib):
            self.warning("matplotlib module not found, can draw")
            return

        plt.figure("results",figsize=(12,7),edgecolor="w")

        sb1 = plt.subplot(gs[0])
        plt.loglog()
        plt.xlabel("Energy ["+ self.eunit+"]",fontsize=15)
        
        if self.m_sed:
            plt.ylabel("E$^{2}$ dN/dE [erg cm$^{2}$ s$^{-1}$]",fontsize=15)
        else:
            plt.ylabel("dN/dE [cm$^{-2}$ s$^{-1}$ "+ self.eunit+"$^{-1}$]",fontsize=15)

        plt.plot(self.m_energy,self.m_flux)
        plt.plot(self.m_enebut,self.m_but)
        
        #~ plot the data points on top of the buterfly
        self.points._plt_points()
        self.points._plt_limits()

        #~ plot the residual points
        sb2 = plt.subplot(gs[1])     
        plt.xscale('log')     
        self.residuals._plt_residuals()

        plt.show()


# ============= #
# Plot skymaps  #
# ============= #
class MapPlotter(analysisutils.base):
    def __init__(self,modelmap, countmap):
        super(MapPlotter,self).__init__()
        self.modelmap = modelmap
        self.countmap = countmap
        
        if not(has_aplpy) :
            self.error("Module aplypy needed")
        if not(has_matplotlib) :
            self.error("Module matplotlib needed")
            
    def draw(self):
        """ Plot the model map.    """
        # Load model map
        fig = plt.figure()
        f1 = aplpy.FITSFigure(self.countmap, figure=fig, subplot=[0.1,0.1,0.35,0.8])
        f1.set_tick_labels_font(size='small')
        f1.set_axis_labels_font(size='small')
        f1.show_colorscale()

        f1.tick_labels.set_yformat('dd:mm:ss')
        f1.tick_labels.set_xformat('hh:mm:ss')
        f1.axis_labels.set_xtext('Right Ascension (J2000)')
        f1.axis_labels.set_ytext('Declination (J2000)')
        f1.ticks.set_length(10.5, minor_factor=0.5)

        f2 = aplpy.FITSFigure(self.modelmap, figure=fig, subplot=[0.5,0.1,0.35,0.8])
        f2.set_tick_labels_font(size='small')
        f2.set_axis_labels_font(size='small')
        f2.show_colorscale()

        f2.hide_yaxis_label()
        f2.hide_ytick_labels()

        fig.canvas.draw()
        fig.show()
