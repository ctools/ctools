# import gammalib first
import gammalib

# Import sub modules
from .tools import *

# Remove swig register stuff
bad_entries = [entry for entry in list(locals())
               if entry.endswith('_swigregister')]
for entry in bad_entries:
    del locals()[entry]

# Clean-up
del tools
del gammalib
del bad_entries
del entry

# Add test function
def test(verbosity=2):
    """
    Run ctools tests.
    """
    # Imports
    import unittest
    from .tests.unbinned import pipelines

    # Create a test suite
    test_suite = unittest.TestSuite()

    # Test cases must be listed manually here ...
    # We currently don't use unittest for test collection
    test_suite.addTest(unittest.makeSuite(pipelines))

    # Create a test runner and run the test
    test_runner = unittest.TextTestRunner(verbosity=verbosity)
    test_runner.run(test_suite)
