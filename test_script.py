import unittest
import numpy as np
import pyUniDOE as pydoe

class TestpyUniDOE(unittest.TestCase):
  """
  pyUniDOE test class.
  """

  def setUp(self):
    """ Initialise test suite - no-op. """
    pass

  def tearDown(self):
    """ Clean up test suite - no-op. """
    pass

  def test_DesignEval(self):
    """ Test DesignEval. """
    x = [[1, 2],
         [3, 3],
         [2, 1]]
    self.assertTrue(np.abs(pydoe.DesignEval(x,crit="CD2") - 0.029578) <0.001)

if __name__ == '__main__':
    unittest.main()
