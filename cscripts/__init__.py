"""cscripts --- Python scripts from the ctools package.

http://cta.irap.omp.eu/ctools/
"""


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


def test():
    """Run ctools tests.
    """
    import unittest
    from .tests.test_unbinned import TestUnbinnedAnalysis
    from .tests.test_binned import TestBinnedAnalysis
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestUnbinnedAnalysis)
    # unittest.TextTestRunner(verbosity=2).run(suite)
    test_suite = unittest.TestSuite()

    # Test cases must be listed manually here ...
    # we currently don't use unittest for test collection
    test_suite.addTest(unittest.makeSuite(TestUnbinnedAnalysis))
    test_suite.addTest(unittest.makeSuite(TestBinnedAnalysis))

    unittest.TextTestRunner(verbosity=2).run(test_suite)
