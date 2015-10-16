# ==========================================================================
# This script tests the cscripts.
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
import os
from tempfile import mkdtemp
import unittest


# ================== #
# Analysis pipelines #
# ================== #
class scripts(unittest.TestCase):
    """
    Test the cscripts.
    """
    # Setup of test class
    def setUp(self):
        gammalib_dir    = os.environ['GAMMALIB']
        self.model_name = os.path.join(gammalib_dir, "share/models/crab.xml")
        self.workdir    = mkdtemp()

    # Test cscaldb
    def test_cscaldb(self):
        """
        Test cscaldb.
        """
        # Test cscaldb
        #caldb = cscaldb()
        #caldb["debug"] = True
        #caldb.execute()

        # TODO: add asserts
        assert True
