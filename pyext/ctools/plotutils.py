# ==========================================================================
# This script provides a number of functions that are useful for handling
# CTA observations.
#
# Copyright (C) 2014 Rolf Buehler
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
import ctools
import gammalib
try:
        import matplotlib.pyplot as plt
        has_matplotlib = True
except:
        has_matplotlib = False
try:
        import aplpy
        has_aplpy = True
except:
        has_aplpy = False


# ============= #
# Plot skymaps  #
# ============= #
def plotskymaps(modelmap, countmap):
    """
    Plot the model map.

    Parameters:    
     modname - model input name
     evtname - events input name
    Keywords:
     None
    """
    # Load model map
    fig = plt.figure()
    f1 = aplpy.FITSFigure(countmap, figure=fig, subplot=[0.1,0.1,0.35,0.8])
    f1.set_tick_labels_font(size='small')
    f1.set_axis_labels_font(size='small')
    f1.show_colorscale()

    f1.tick_labels.set_yformat('dd:mm:ss')
    f1.tick_labels.set_xformat('hh:mm:ss')
    f1.axis_labels.set_xtext('Right Ascension (J2000)')
    f1.axis_labels.set_ytext('Declination (J2000)')
    f1.ticks.set_length(10.5, minor_factor=0.5)

    f2 = aplpy.FITSFigure(modelmap, figure=fig, subplot=[0.5,0.1,0.35,0.8])
    f2.set_tick_labels_font(size='small')
    f2.set_axis_labels_font(size='small')
    f2.show_colorscale()

    f2.hide_yaxis_label()
    f2.hide_ytick_labels()

    fig.canvas.draw()
    fig.show()
    
