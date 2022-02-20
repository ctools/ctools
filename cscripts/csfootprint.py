#! /usr/bin/env python
# ==========================================================================
# Carbon footprint report script
#
# Copyright (C) 2022 Juergen Knoedlseder
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
import gammalib
import ctools


# ================= #
# csfootprint class #
# ================= #
class csfootprint(ctools.cscript):
    """
    Carbon footprint report script
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor

        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._statistics = {}

        # Return
        return

    # State methods for pickling
    def __getstate__(self):
        """
        Extend ctools.cscript __getstate__ method

        Returns
        -------
        state : dict
            Pickled instance
        """
        # Set pickled dictionary
        state = {'base'       : ctools.cscript.__getstate__(self),
                 'statistics' : self._statistics}

        # Return pickled dictionary
        return state

    def __setstate__(self, state):
        """
        Extend ctools.cscript __setstate__ method

        Parameters
        ----------
        state : dict
            Pickled instance
        """
        # Set state
        ctools.cscript.__setstate__(self, state['base'])
        self._statistics = state['statistics']

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query input parameters
        self['infile'].filename()
        if self['tmin'].is_valid():
            self['tmin'].time()
        if self['tmax'].is_valid():
            self['tmax'].time()

        # Query ahead output model filename
        if self._read_ahead():
            if self['outfile'].is_valid():
                self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    # Load statistics
    def _load_statistics(self):
        """
        Load statistics
        """
        # Initialise statistics
        statistics = {}

        # Get filename
        infile = self['infile'].filename()

        # Load XML file
        xml = gammalib.GXml()
        xml.load(infile)

        # Get useful nodes
        header    = xml.element('statistics > header')
        data      = xml.element('statistics > data')
        dates     = header.element('dates')
        countries = data.element('countries')
        versions  = data.element('versions')
        daily     = data.element('daily')

        # Extract header dates information
        statistics['creation'] = dates.element('creation').string()
        statistics['modified'] = dates.element('modified').string()
        statistics['start']    = dates.element('start').string()
        statistics['stop']     = dates.element('stop').string()

        # Set used time interval [start,stop]
        if self['tmin'].is_valid():
            start = self['tmin'].time().utc()
            if start < statistics['start']:
                start = statistics['start']
        else:
            start = statistics['start']
        if self['tmax'].is_valid():
            stop = self['tmax'].time().utc()
            if stop > statistics['stop']:
                stop = statistics['stop']
        else:
            stop = statistics['stop']

        # Set used time interval
        statistics['use_start'] = start
        statistics['use_stop']  = stop

        # Extract list of countries
        statistics['countries'] = []
        num = countries.elements()
        for i in range(num):
            element = countries.element(i)
            country = element.name()
            calls   = element.attribute('calls')
            wall    = element.attribute('wall')
            cpu     = element.attribute('cpu')
            gCO2e   = element.attribute('gCO2e')
            entry   = {'country': country, 'calls': calls, 'wall': wall,
                       'cpu': cpu, 'gCO2e': gCO2e}
            statistics['countries'].append(entry)

        # Extract list of versions
        statistics['versions'] = []
        num = versions.elements()
        for i in range(num):
            element = versions.element(i)
            version = element.name()
            calls   = element.attribute('calls')
            wall    = element.attribute('wall')
            cpu     = element.attribute('cpu')
            gCO2e   = element.attribute('gCO2e')
            entry   = {'version': version, 'calls': calls, 'wall': wall,
                       'cpu': cpu, 'gCO2e': gCO2e}
            statistics['versions'].append(entry)

        # Extract daily and tools statistics
        statistics['daily'] = []
        statistics['tools'] = []
        num                 = daily.elements()
        total_calls         = 0
        total_wall          = 0.0
        total_cpu           = 0.0
        total_gCO2e         = 0.0
        for i in range(num):

            # Get date
            element = daily.element(i)
            date    = element.element('value',0).string()

            # Continue if date is out of range
            if date < start[0:10] or date > stop[0:10]:
                continue

            # Collect information
            tools     = element.element('tools',0)
            ntools    = tools.elements()
            sum_calls = 0
            sum_wall  = 0.0
            sum_cpu   = 0.0
            sum_gCO2e = 0.0
            for k in range(ntools):

                # Extract information
                tool  = tools.element(k)
                name  = tool.name()
                calls = int(tool.attribute('calls'))
                wall  = float(tool.attribute('wall'))
                cpu   = float(tool.attribute('cpu'))
                gCO2e = float(tool.attribute('gCO2e'))
                
                # Update sums
                sum_calls += calls
                sum_wall  += wall
                sum_cpu   += cpu
                sum_gCO2e += gCO2e

                # Update tools data
                if name in statistics['tools']:
                    statistics['tools'][name]['calls'] += calls
                    statistics['tools'][name]['wall']  += wall
                    statistics['tools'][name]['cpu']   += cpu
                    statistics['tools'][name]['gCO2e'] += gCO2e
                else:
                    entry = {'name': name, 'calls': calls, 'wall': wall,
                             'cpu': cpu, 'gCO2e': gCO2e}
                    statistics['tools'].append(entry)

            # Update totals
            total_calls += sum_calls
            total_wall  += sum_wall
            total_cpu   += sum_cpu
            total_gCO2e += sum_gCO2e

            # Set entry
            entry = {'date': date, 'calls': sum_calls, 'wall': sum_wall,
                     'cpu': sum_cpu, 'gCO2e': sum_gCO2e}
            statistics['daily'].append(entry)

        # Update total statistics
        statistics['calls'] = total_calls
        statistics['wall']  = total_wall
        statistics['cpu']   = total_cpu
        statistics['gCO2e'] = total_gCO2e

        # Return statistics
        return statistics

    # Format time
    def _format_time(self, seconds):
        """
        Format time

        Parameters
        ----------
        seconds : float
            Time in seconds
        """
        # Format according to precision
        if seconds < 60.0:
            format = '%.3f seconds' % (seconds)
        elif seconds < 3600.0:
            format = '%.3f minutes' % (seconds/60.0)
        else:
            format = '%.3f hours' % (seconds/3600.0)

        # Return format
        return format

    # Format carbon footprint
    def _format_footprint(self, gCO2e):
        """
        Format carbon footprint

        Parameters
        ----------
        gCO2e : float
            Carbon footprint in gCO2e
        """
        # Format according to precision
        if gCO2e < 1000.0:
            format = '%.3f g CO2e' % (gCO2e)
        elif seconds < 1.0e6:
            format = '%.3f kg CO2e' % (gCO2e/1000.0)
        else:
            format = '%.3f t CO2e' % (gCO2e/1.0e6)

        # Return format
        return format

    # Log global statistics
    def _global_statistics(self, statistics):
        """
        Log global statistics

        Parameters
        ----------
        statistics : dict
            Statistics dictionary
        """
        # Log header
        self._log_header1(gammalib.TERSE, 'Global statistics')

        # Derive information
        tstart   = gammalib.GTime(statistics['use_start'])
        tstop    = gammalib.GTime(statistics['use_stop'])
        duration = tstop - tstart
        if statistics['cpu'] != 0.0:
            ci_cpu = statistics['gCO2e']/(statistics['cpu']/3600.0)
        else:
            ci_cpu = 0.0
        if duration != 0.0:
            fp_dur = statistics['gCO2e']/(duration * gammalib.sec2day)
        else:
            fp_dur = 0.0
        fp_yr      = fp_dur * 365.25
        dur_str    = self._format_time(duration)
        if statistics['wall'] != 0.0:
            load_str = '%.1f %%' % (statistics['cpu']/statistics['wall']*100.0)
        else:
            load_str = 'undefined'
        ci_cpu_str = self._format_footprint(ci_cpu) + ' / CPU hour'
        fp_dur_str = self._format_footprint(fp_dur) + ' / day'
        fp_yr_str  = self._format_footprint(fp_yr) + ' / year'
        wall_str   = self._format_time(statistics['wall'])
        cpu_str    = self._format_time(statistics['cpu'])
        gCO2e_str  = self._format_footprint(statistics['gCO2e'])
        date_str   = '%s - %s' % (statistics['start'], statistics['stop'])
        use_str    = '%s - %s' % (statistics['use_start'], statistics['use_stop'])

        # Log report information
        self._log_value(gammalib.NORMAL, 'Creation date', statistics['creation'])
        self._log_value(gammalib.NORMAL, 'Last statistics update', statistics['modified'])
        self._log_value(gammalib.NORMAL, 'Statistics date interval', date_str)
        self._log_value(gammalib.NORMAL, 'Used date interval', use_str)
        self._log_value(gammalib.NORMAL, 'Duration of used interval', dur_str)
        self._log_value(gammalib.NORMAL, 'Total number of ctool runs', statistics['calls'])
        self._log_value(gammalib.NORMAL, 'Total wall clock time', wall_str)
        self._log_value(gammalib.NORMAL, 'Total CPU time', cpu_str)
        self._log_value(gammalib.NORMAL, 'Average CPU load', load_str)
        self._log_value(gammalib.NORMAL, 'Total carbon footprint', gCO2e_str)
        self._log_value(gammalib.NORMAL, 'Average carbon intensity', ci_cpu_str)
        self._log_value(gammalib.NORMAL, 'Average daily footprint', fp_dur_str)
        self._log_value(gammalib.NORMAL, 'Expected annual footprint', fp_yr_str)

        # Return
        return

    # Log daily statistics
    def _daily_statistics(self, statistics):
        """
        Log daily statistics

        Parameters
        ----------
        statistics : dict
            Statistics dictionary
        """
        # Log header
        self._log_header1(gammalib.TERSE, 'Daily statistics')

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Carbon footprint')

        # Loop over daily entries
        for entry in statistics['daily']:
            self._log_value(gammalib.NORMAL, entry['date'], self._format_footprint(entry['gCO2e']))

        # Log daily run statistics
        self._log_header3(gammalib.NORMAL, 'ctools or cscript calls')

        # Loop over daily entries
        for entry in statistics['daily']:
            self._log_value(gammalib.NORMAL, entry['date'], entry['calls'])

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Used wall clock time')

        # Loop over daily entries
        for entry in statistics['daily']:
            self._log_value(gammalib.NORMAL, entry['date'], self._format_time(entry['wall']))

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Used CPU time')

        # Loop over daily entries
        for entry in statistics['daily']:
            self._log_value(gammalib.NORMAL, entry['date'], self._format_time(entry['cpu']))

        # Return
        return

    # Log tools statistics
    def _tools_statistics(self, statistics):
        """
        Log tools statistics

        Parameters
        ----------
        statistics : dict
            Statistics dictionary
        """
        # Get Python version information
        req_version = (2,4)
        cur_version = sys.version_info

        # Log header
        self._log_header1(gammalib.TERSE, 'ctools and cscripts statistics')

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Carbon footprint')

        # Optionally sort list
        if cur_version > req_version:
            sorted_entries = sorted(statistics['tools'], key=lambda d: d['gCO2e'], reverse=True)
        else:
            sorted_entries = statistics['tools']

        # Loop over entries
        for i, entry in enumerate(sorted_entries):
            if i < 10:
                level = gammalib.NORMAL
            else:
                level = gammalib.EXPLICIT
            self._log_value(level, entry['name'], self._format_footprint(entry['gCO2e']))
        if len(sorted_entries) > 9 and self['chatter'].integer() < 3:
            self._log_string(gammalib.NORMAL, ' ... (list truncated after 10 entries) ...')

        # Log daily run statistics
        self._log_header3(gammalib.NORMAL, 'ctools or cscript calls')

        # Optionally sort list
        if cur_version > req_version:
            sorted_entries = sorted(statistics['tools'], key=lambda d: d['calls'], reverse=True)
        else:
            sorted_entries = statistics['tools']

        # Loop over entries
        for i, entry in enumerate(sorted_entries):
            if i < 10:
                level = gammalib.NORMAL
            else:
                level = gammalib.EXPLICIT
            self._log_value(level, entry['name'], entry['calls'])
        if len(sorted_entries) > 9 and self['chatter'].integer() < 3:
            self._log_string(gammalib.NORMAL, ' ... (list truncated after 10 entries) ...')

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Used wall clock time')

        # Optionally sort list
        if cur_version > req_version:
            sorted_entries = sorted(statistics['tools'], key=lambda d: d['wall'], reverse=True)
        else:
            sorted_entries = statistics['tools']

        # Loop over entries
        for i, entry in enumerate(sorted_entries):
            if i < 10:
                level = gammalib.NORMAL
            else:
                level = gammalib.EXPLICIT
            self._log_value(level, entry['name'], self._format_time(entry['wall']))
        if len(sorted_entries) > 9 and self['chatter'].integer() < 3:
            self._log_string(gammalib.NORMAL, ' ... (list truncated after 10 entries) ...')

        # Log daily carbon footprint statistics
        self._log_header3(gammalib.NORMAL, 'Used CPU time')

        # Optionally sort list
        if cur_version > req_version:
            sorted_entries = sorted(statistics['tools'], key=lambda d: d['cpu'], reverse=True)
        else:
            sorted_entries = statistics['tools']

        # Loop over entries
        for i, entry in enumerate(sorted_entries):
            if i < 10:
                level = gammalib.NORMAL
            else:
                level = gammalib.EXPLICIT
            self._log_value(level, entry['name'], self._format_time(entry['cpu']))
        if len(sorted_entries) > 9 and self['chatter'].integer() < 3:
            self._log_string(gammalib.NORMAL, ' ... (list truncated after 10 entries) ...')
            

        # Return
        return

    # Create figure
    def _create_figure(self, outfile, statistics):
        """
        Create figure

        Parameters
        ----------
        outfile : str
            Figure file name
        statistics : dict
            Statistics dictionary
        """
        # Optionally use matplotlib to create a figure
        try:
        
            # Import matplotlib
            import matplotlib.pyplot   as plt
            import matplotlib.gridspec as gridspec

            # Create figure
            fig = plt.figure(figsize=(15.0,7.8))

            # Create title and subtitle
            title    = 'ctools carbon footprint report'
            subtitle = r'Dates: %s - %s' % \
                       (statistics['use_start'], statistics['use_stop'])
            fig.suptitle(title, fontsize=16)
            fig.text(0.5, 0.93, subtitle, fontsize=12, ha='center')

            # Set plot margins
            fig.subplots_adjust(left=0.07, bottom=0.07, right=0.97, top=0.88, wspace=0.3, hspace=0.35)

            # Divide figure
            gs1 = gridspec.GridSpec(3,3)
            ax1 = fig.add_subplot(gs1[0,0:2])
            ax2 = fig.add_subplot(gs1[0,2])

            # Plot daily footprint
            self._plot_daily(ax1, statistics)

            #plt.show()

            # Save figure
            fig.savefig(outfile, dpi=300)

            # Log file creation
            self._log_value(gammalib.NORMAL, 'Graphics file', outfile)

        # Catch exceptions
        except (ImportError, RuntimeError):

            # Log file creation
            self._log_value(gammalib.NORMAL, 'Graphics file', 'matplotlib not available')

        # Return
        return

    # Plot daily footprint
    def _plot_daily(self, ax, statistics):
        """
        Plot daily footprint

        Parameters
        ----------
        ax : pyplot
            Plotting frame
        statistics : dict
            Statistics dictionary
        """
        # Create bar data
        days = [i              for i, _  in enumerate(statistics['daily'])]
        data = [entry['gCO2e'] for entry in statistics['daily']]

        # Plot bar data
        ax.bar(days, data, 1.0, bottom=0.0, color='red')

        # Set labels
        ax.set_title('Footprint')
        ax.set_xlabel('Days since %s' % statistics['use_start'][0:10])
        ax.set_ylabel(r'g CO$_2$e')
        
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
        self._log_header1(gammalib.TERSE, 'Load statistics data')

        # Load statistics
        self._statistics = self._load_statistics()

        # Log statistics
        self._global_statistics(self._statistics)
        self._daily_statistics(self._statistics)
        self._tools_statistics(self._statistics)

        # Return
        return

    def save(self):
        """
        Save something
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save graphics')

        # Continue only if filename is valid
        if self['outfile'].is_valid():

            # Get outfile
            outfile = self['outfile'].filename()

            # Create figure
            self._create_figure(outfile.url(), self._statistics)

        # ... signal that no graphics file was specified
        else:
            self._log_value(gammalib.NORMAL, 'Graphics file', 'not specified')

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csfootprint(sys.argv)

    # Execute application
    app.execute()
