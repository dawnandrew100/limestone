
from __future__ import annotations
import unittest
from limestone import editdistance

class TestNeedlemanWunsch(unittest.TestCase):
        
    def test_align(self):
        nwalignment = editdistance.needlemanWunsch.align("BA", "ABA")
        hsalignment = editdistance.hirschbergDist.align("BA","ABA")
        self.assertEqual(nwalignment, hsalignment)

if __name__ == '__main__':
    unittest.main()
