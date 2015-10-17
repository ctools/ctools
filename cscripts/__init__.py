# ==========================================================================
# cscripts Python module
#
# Copyright (C) 2015 Juergen Knoedlseder
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

# List all scripts in this module
__all__ = [
    "cscaldb",
    "cshessobs",
    "csinfo",
    "cslightcrv",
    "csobsdef",
    "cspull",
    "csresmap",
    "cssens",
    "csspec",
    "cstsdist",
    "obsutils",
    "test",
]

# Import
from .cscaldb    import cscaldb
from .cshessobs  import cshessobs
from .csinfo     import csinfo
from .cslightcrv import cslightcrv
from .csobsdef   import csobsdef
from .cspull     import cspull
from .csresmap   import csresmap
from .cssens     import cssens
from .csspec     import csspec
from .cstsdist   import cstsdist
from . import obsutils

# Add test function
def test():
    """
    Run cscripts tests.
    """
    print("No cscripts unit tests have been implemented so far.")
