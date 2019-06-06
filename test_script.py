import unittest
import numpy as np
import pyunidoe as pydoe

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

  def test_design_eval(self):
    """ Test DesignEval. """
    x = np.array([[1, 2],
                  [3, 3],
                  [2, 1]])
    self.assertTrue(np.abs(pydoe.design_eval(
      x, crit="CD2") - 0.029578) < 0.001)

if __name__ == '__main__':
    unittest.main()
