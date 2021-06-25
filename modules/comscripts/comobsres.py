#! /usr/bin/env python
# ==========================================================================
# Generate residuals of COMPTEL observations
#
# Copyright (C) 2019-2021 Juergen Knoedlseder
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
import sys
import glob
import os
import math
import gammalib
import ctools


# ======================= #
# Gaussian function class #
# ======================= #
class gaussian(gammalib.GPythonOptimizerFunction):

    # Constructor
    def __init__(self, x_vals, y_vals):

        # Call base class constructor
        gammalib.GPythonOptimizerFunction.__init__(self)

        # Set eval method
        self._set_eval(self.eval)

        # Set data
        self._x_vals = x_vals
        self._y_vals = y_vals

    # Methods
    def eval(self):
        """
        Evaluate function
        """
        # Recover parameters
        pars  = self._pars()
        norm  = pars[0].value()
        mean  = pars[1].value()
        sigma = pars[2].value()

        # Evaluate function values
        y = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in self._x_vals]

        # Compute weights (1/sqrt(y))
        weight = []
        for val in self._y_vals:
            if val > 0.0:
                weight.append(1.0/val)
            else:
                weight.append(0.0)

        # Compute Chi Square
        value = 0.0
        for i in range(len(self._x_vals)):
            arg    = self._y_vals[i] - y[i]
            value += arg * arg * weight[i]

        # Evaluate gradient and curvature
        sigma2 = sigma  * sigma
        sigma3 = sigma2 * sigma
        for i in range(len(self._x_vals)):

            # Evaluate function gradients
            dx     = self._x_vals[i] - mean
            dnorm  = y[i]         / norm   * pars[0].scale()
            dmean  = y[i] * dx    / sigma2 * pars[1].scale()
            dsigma = y[i] * dx**2 / sigma3 * pars[2].scale()

            # Setup gradient vector
            arg                 = (self._y_vals[i] - y[i]) * weight[i]
            self.gradient()[0] -= arg * dnorm
            self.gradient()[1] -= arg * dmean
            self.gradient()[2] -= arg * dsigma

            # Setup curvature matrix
            self.curvature()[0,0] +=  dnorm  * dnorm   * weight[i]
            self.curvature()[0,1] +=  dnorm  * dmean   * weight[i]
            self.curvature()[0,2] +=  dnorm  * dsigma  * weight[i]
            self.curvature()[1,0] +=  dmean  * dnorm   * weight[i]
            self.curvature()[1,1] +=  dmean  * dmean   * weight[i]
            self.curvature()[1,2] +=  dmean  * dsigma  * weight[i]
            self.curvature()[2,0] +=  dsigma * dnorm   * weight[i]
            self.curvature()[2,1] +=  dsigma * dmean   * weight[i]
            self.curvature()[2,2] +=  dsigma * dsigma  * weight[i]

        # Set value
        self._set_value(value)

        # Return
        return


# =============== #
# comobsres class #
# =============== #
class comobsres(ctools.csobservation):
    """
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._outmap  = gammalib.GSkyMap()
        self._ebounds = gammalib.GEbounds()

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self.obs().load(self['inobs'].filename())

        # Set models if we have none
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        # Query parameters
        self['algorithm'].string()
        self['armmin'].real()
        self['armmax'].real()
        self['margin'].real()

        # Query PNG parameters
        if self['png'].boolean():
            self['outfolder'].string()

        # Query DRI parameters
        if self['dri'].boolean():
            self['outfolder'].string()
            self['algorithm'].string()
            self['grouping'].integer()

        # Query ahead output skymap filename
        if self._read_ahead():
            self['outmap'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _get_map(self, inobs):
        """
        Get sky map with size determined from observations

        Parameters
        ----------
        inobs : `~gammalib.GObservations`
            Observation container

        Returns
        -------
        map : `~gammalib.GSkyMap`
            Sky map
        """
        # Initialise boundaries
        glon_min = None
        glon_max = None
        glat_min = None
        glat_max = None

        # Loop over all input observations
        for obs in inobs:

            # Get cube as sky map
            map = obs.events().dre().map()

            # Get step size in longitude
            cdelt1 = map.projection().cdelt(0)

            # Get sky map limits
            if cdelt1 > 0.0:
                l_min = map.pix2dir(gammalib.GSkyPixel(-0.5,         map.ny()/2.0)).l_deg()
                l_max = map.pix2dir(gammalib.GSkyPixel(map.nx()-0.5, map.ny()/2.0)).l_deg()
            else:
                l_min = map.pix2dir(gammalib.GSkyPixel(map.nx()-0.5, map.ny()/2.0)).l_deg()
                l_max = map.pix2dir(gammalib.GSkyPixel(-0.5,         map.ny()/2.0)).l_deg()
            b_min = map.pix2dir(gammalib.GSkyPixel(map.nx()/2.0,         -0.5)).b_deg()
            b_max = map.pix2dir(gammalib.GSkyPixel(map.nx()/2.0, map.ny()-0.5)).b_deg()
            if l_min > l_max:
                l_min -= 360.0

            # Set largest boundary
            if glon_min == None:
                glon_min = l_min
                glon_max = l_max
                glat_min = b_min
                glat_max = b_max
            else:
                if l_min < glon_min:
                    glon_min = l_min
                if l_max < glon_max:
                    glon_max = l_max
                if b_min < glat_min:
                    glat_min = b_min
                if b_max > glat_max:
                    glat_max = b_max

        # Compute map centre
        glon = 0.5 * (glon_min + glon_max)
        glat = 0.5 * (glat_min + glat_max)

        # Add margins
        glon_min -= self['margin'].real()
        glat_min -= self['margin'].real()
        glon_max += self['margin'].real()
        glat_max += self['margin'].real()

        # Restrict map to sky
        if glat_min < -90.0:
            glat_min = -90.0
        if glat_max > 90.0:
            glat_max = 90.0

        # Compute width and number of pixels
        glon_width = glon_max - glon_min
        glat_width = glat_max - glat_min
        nx         = int(glon_width+0.5)
        ny         = int(glat_width+0.5)
        if nx > 360:
            nx = 360
        if ny > 180:
            ny = 180

        # Allocate map
        map = gammalib.GSkyMap('TAN', 'GAL', glon, glat, -1.0, 1.0, nx, ny)

        # Return map
        return map

    def _get_ebounds(self, inobs):
        """
        Get energy boundaries of observations

        Parameters
        ----------
        inobs : `~gammalib.GObservations`
            Observation container

        Returns
        -------
        ebounds : `~gammalib.GEbounds`
            Energy boundaries
        """
        # Initialise list of energy boundaries
        emins = []
        emaxs = []

        # Loop over all input observations
        for obs in inobs:

            # Get energy boundaries for observation
            ebounds = obs.events().ebounds()
            emin    = ebounds.emin().MeV()
            emax    = ebounds.emax().MeV()

            # If energy is not yet in list then add it now
            if emin not in emins and emax not in emaxs:
                emins.append(emin)
                emaxs.append(emax)

        # Build energy boundaries
        ebounds = gammalib.GEbounds()

        # Append boundaries
        for i in range(len(emins)):
            ebounds.append(gammalib.GEnergy(emins[i], 'MeV'),
                           gammalib.GEnergy(emaxs[i], 'MeV'))

        # Return energy boundaries
        return ebounds

    def _fill_map(self, dre, drm, imap):
        """
        Fill residual map

        Parameters
        ----------
        dre : `~gammalib.GCOMDri`
            Event cube
        drm : `~gammalib.GCOMDri`
            Model cube
        imap : int
            Energy index of sky map

        Returns
        -------
        sig : list of float
            List of residual significances
        """
        # Get residual algorithm
        algorithm = self['algorithm'].string()

        # Get ARM range
        armmin = self['armmin'].real()
        armmax = self['armmax'].real()

        # Initialise list of residual significances
        sig = []

        # Loop over all sky map pixels
        for i in range(self._outmap.npix()):

            # Get sky direction of sky map pixel
            dir = self._outmap.inx2dir(i)

            # Compute events and model value
            events = dre.cone_content(dir, armmin, armmax)
            model  = drm.cone_content(dir, armmin, armmax)

            # Compute residual
            if algorithm == 'SUB':
                residual  = events - model
            elif algorithm == 'SUBDIV':
                if model > 0.0:
                    residual = (events - model) / model
                else:
                    residual = 0.0
            elif algorithm == 'SUBDIVSQRT':
                if model > 0.0:
                    residual = (events - model) / math.sqrt(model)
                else:
                    residual = 0.0
            elif algorithm == 'SIGNIFICANCE': # Significance from Li&Ma
                if events >= model:
                    sign = +1.0
                else:
                    sign = -1.0
                if model > 0.0:
                    if events > 0.0:
                        log_val  = math.log(events / model)
                        residual = (events * log_val) + model - events
                        if residual < 0.0:   # Cope with rounding errors
                            residual = 0.0
                    else:
                        residual = model
                else:
                    residual = 0.0
                residual *= 2.0
                residual  = math.sqrt(residual)
                residual *= sign
            else:
                raise TypeError('Algorithm "' + algorithm + '" not known')

            # Store residual
            self._outmap[i, imap] = residual

            # Append residual to list of significances
            if model > 0.0:
                sig.append(residual)

        # Write residual range
        if len(sig) > 0:
            self._log_value(gammalib.NORMAL, 'Minimum map residual', min(sig))
            self._log_value(gammalib.NORMAL, 'Maximum map residual', max(sig))
        else:
            self._log_value(gammalib.NORMAL, 'Minimum map residual', 'empty map')
            self._log_value(gammalib.NORMAL, 'Maximum map residual', 'empty map')

        # Return
        return sig

    def _compute_residual_dri(self, dre, drm, imap, energy_bin, group):
        """
        Compute residual RDI

        Parameters
        ----------
        dre : `~gammalib.GCOMDri`
            Event cube
        drm : `~gammalib.GCOMDri`
            Model cube
        imap : int
            Energy index of sky map
        energy_bin : str
            Energy bin string

        Returns
        -------
        sig : list of float
            List of residual significances
        """
        # Get residual algorithm
        algorithm = self['algorithm'].string()

        # Set filename of residual FITS file
        filename = 'residual_%s.fits' % (energy_bin)

        # Initialise result
        drs = dre.copy()
 
        # Initialise list of residual significances
        sig = []

        # Loop over data space
        for iphibar in range(drs.nphibar()):
            for iy in range(0, drs.npsi(), group):
                for ix in range(0, drs.nchi(), group):

                    # Compute number of events and model in sub bin
                    events = 0.0
                    model  = 0.0
                    for jy in range(group):
                        y = iy + jy
                        if y < drs.npsi():
                            for jx in range(group):
                                x = ix + jx
                                if x < drs.nchi():

                                    # Compute index
                                    inx = x + y * drs.nchi() + iphibar * drs.nchi() * drs.npsi()

                                    # Sum events and model value
                                    events += dre[inx]
                                    model  += drm[inx]

                    # Compute residual
                    if algorithm == 'SUB':
                         residual  = events - model
                    elif algorithm == 'SUBDIV':
                         if model > 0.0:
                              residual = (events - model) / model
                         else:
                              residual = 0.0
                    elif algorithm == 'SUBDIVSQRT':
                         if model > 0.0:
                              residual = (events - model) / math.sqrt(model)
                         else:
                              residual = 0.0
                    elif algorithm == 'SIGNIFICANCE': # Significance from Li&Ma
                         if events >= model:
                              sign = +1.0
                         else:
                              sign = -1.0
                         if model > 0.0:
                              if events > 0.0:
                                   log_val  = math.log(events / model)
                                   residual = (events * log_val) + model - events
                              else:
                                   residual = model
                         else:
                              residual = 0.0
                         residual *= 2.0
                         residual  = math.sqrt(residual)
                         residual *= sign
                    else:
                         raise TypeError('Algorithm "' + algorithm + '" not known')

                    # Store residuals
                    for jy in range(group):
                        y = iy + jy
                        if y < drs.npsi():
                            for jx in range(group):
                                x = ix + jx
                                if x < drs.nchi():
                                    inx = x + y * drs.nchi() + iphibar * drs.nchi() * drs.npsi()
                                    drs[inx] = residual

                    # Append residual to list of significances
                    if model > 0.0:
                        sig.append(residual)

        # Save residual DRI
        drs.save(filename, True)

        # Write residual range
        self._log_value(gammalib.NORMAL, 'Minimum DRI residual', min(sig))
        self._log_value(gammalib.NORMAL, 'Maximum DRI residual', max(sig))

        # Return
        return sig

    def _plot_sighist(self, sighist, nbins=80, sigma_min=-6.0, sigma_max=6.0, fontsize=11, title=None):
        """
        Plot the significance distribution

        Parameters
        ----------
        sighist : list of floats
            Significance values
        nbins : int, optional
            Number of histogram bins
        sigma_min : float, optional
            Lower limit of the x axis in the plot
        sigma_max : float, optional
            Upper limit of the x axis in the plot
        fontsize : int, optional
              Font size for axes
        title : str, optional
            Title string
        """
        # Plot only if matplotlib is available
        try:

            # Import matplotlib
            import matplotlib.pyplot as plt

            # Plot histogram
            self._plot_histogram(sighist, nbins=nbins, sigma_min=sigma_min,
                                 sigma_max=sigma_max, fontsize=fontsize,
                                 title=title)

        except (ImportError, RuntimeError):

            # Inform about impossibility to plot
            self._log_string(gammalib.NORMAL, 'matplotlib not available, '+
                             'no significance histogram plotted.')

        # Return
        return

    def _plot_histogram(self, sighist, nbins=80, sigma_min=-6.0, sigma_max=6.0, fontsize=11, title=None):
        """
        Plot the significance distribution

        Parameters
        ----------
        sighist : list of floats
            Significance values
        nbins : int, optional
            Number of histogram bins
        sigma_min : float, optional
            Lower limit of the x axis in the plot
        sigma_max : float, optional
            Upper limit of the x axis in the plot
        fontsize : int, optional
            Font size for axes
        title : str, optional
            Title string
        """
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Compute bin edges, hence use nbins+1
        binwidth  = (sigma_max - sigma_min) / float(nbins)
        bin_edges = [sigma_min + binwidth*i for i in range(nbins+1)]

        # Draw significance histogram
        f, ax   = plt.subplots()
        y, _, _ = plt.hist(sighist, bins=bin_edges, histtype='step', color='k')

        # Set initial Gaussian parameters
        y_max = float(y.max())
        par1  = gammalib.GOptimizerPar('Norm',  y_max)
        par2  = gammalib.GOptimizerPar('Mean',  0.0)
        par3  = gammalib.GOptimizerPar('Sigma', 1.0)
        pars  = gammalib.GOptimizerPars()
        pars.append(par1)
        pars.append(par2)
        pars.append(par3)

        # Set fit function
        x = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(nbins)]
        fct = gaussian(x, y)

        # Optimize function and compute errors
        opt = gammalib.GOptimizerLM()
        opt.optimize(fct, pars)
        opt.errors(fct, pars)

        # Recover parameters and errors
        norm    = pars[0].value()
        e_norm  = pars[0].error()
        mean    = pars[1].value()
        e_mean  = pars[1].error()
        sigma   = pars[2].value()
        e_sigma = pars[2].error()

        # Draw text box
        msg = 'mean: %.3f +/- %.3f \n width: %.3f +/- %.3f' % \
              (mean, e_mean, sigma, e_sigma)
        ax.text(0.97, 0.95, msg, ha='right', va='top',
                bbox=dict(edgecolor='red', facecolor='white'),
                transform=ax.transAxes, fontsize=8)

        # Plot the normal
        yvals = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in bin_edges]
        ax.plot(bin_edges, yvals, 'r-')

        # Configure the plot
        ax.set_ylim(0.5, ax.get_ylim()[1]*2.5)
        ax.set_xlim(sigma_min, sigma_max)
        ax.set_xlabel('Significance', fontsize=fontsize)
        ax.set_ylabel('Entries', fontsize=fontsize)
        ax.grid()
        ax.set_yscale('log')

        # Add title
        if title != None:
            ax.set_title(title)

        # Save figure
        plt.savefig(title+'.png')

        # Return
        return


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Setup residual map')

        # Get map
        self._outmap = self._get_map(self.obs())

        # Get energy boundaries
        self._ebounds = self._get_ebounds(self.obs())

        # Set number of maps to number of energy bins
        self._outmap.nmaps(self._ebounds.size())

        # Log map and energy boundaries
        self._log_string(gammalib.NORMAL, str(self._outmap))
        self._log_string(gammalib.NORMAL, str(self._ebounds))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Compute residual map')

        # Get grouping for DRI residuals
        if self['dri'].boolean():
            group = self['grouping'].integer()
            self._log_value(gammalib.NORMAL, 'Chi/Psi grouping', group)

        # Get output folder and create folder
        if self['png'].boolean() or self['dri'].boolean():
            outfolder = self['outfolder'].string()
            try:
                os.makedirs(gammalib.expand_env(outfolder))
            except OSError:
                pass
        else:
            outfolder = '.'

        # Loop over all input observations
        for obs in self.obs():

            # Write header
            self._log_header3(gammalib.NORMAL, self._get_obs_header(obs))

            # Get sky map energy index
            imap = self._ebounds.index(obs.events().ebounds().elogmean(0))

            # Throw an exception if index is invalid
            if imap < 0 or imap >= self._outmap.nmaps():
                msg = 'comobsres.run: Invalid sky map index %d encountered.' % imap
                raise RuntimeError(msg)

            # Get DRE
            dre = obs.events().dre()

            # Compute DRM
            drm = obs.drm(self.obs().models())

            # Recover energy bins
            emin = str(obs.events().ebounds().emin().MeV())
            emax = str(obs.events().ebounds().emax().MeV())

            # Build filename without whitespace
            fname = '%s/%s_%s-%sMeV' % (outfolder, obs.id(), emin, emax)
            fname = fname.replace(' ', '_')

            # Fill residual map
            sig = self._fill_map(dre, drm, imap)

            # Generate significance histogram
            if self['png'].boolean():
                self._plot_sighist(sig, title='%s_map' % (fname))

            # Build DRI residuals
            if self['dri'].boolean():

                # Append grouping to filename
                fname = '%s_g%d' % (fname, group)

                # Generate residual data space
                sig = self._compute_residual_dri(dre, drm, imap, '%s' % fname, group)

                # Generate significance histogram
                if self['png'].boolean():
                    self._plot_sighist(sig, title='%s_dri' % fname)

        # Return
        return

    def save(self):
        """ 
        Save residual map file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save residual map')

        # Get output filename
        outmap = self['outmap'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outmap.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmap.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log filename
            self._log_value(gammalib.NORMAL, 'Residual map file',
                                             outmap.url())

            # Save observations
            self._outmap.save(outmap, self['clobber'].boolean())

            # Save energy boundaries
            self._ebounds.save(outmap, self['clobber'].boolean())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsres(sys.argv)

    # Execute application
    app.execute()
