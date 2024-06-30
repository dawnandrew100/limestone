
from __future__ import annotations
import unittest
from limestone import editdistance

class TestHirschberg(unittest.TestCase):
        
    def test_align1(self):
        nwalignment = editdistance.needlemanWunsch.align("BA", "ABA")
        hsalignment = editdistance.hirschbergDist.align("BA","ABA")
        self.assertEqual(nwalignment, hsalignment)
    
    def test_align2(self):
        nwalignment = editdistance.needlemanWunsch.align("ATCG","ATCG")
        hsalignment = editdistance.hirschbergDist.align("ATCG","ATCG")
        self.assertEqual(nwalignment, hsalignment)

    def test_align3(self):
        nwalignment = editdistance.needlemanWunsch.align("ACTG","FHYU")
        hsalignment = editdistance.hirschbergDist.align("ACTG","FHYU")
        self.assertEqual(nwalignment, hsalignment)

    def test_align4(self):
        nwalignment = editdistance.needlemanWunsch.align("AGTACGCA","TATGC")
        hsalignment = editdistance.hirschbergDist.align("AGTACGCA","TATGC")
        self.assertEqual(nwalignment, hsalignment)

if __name__ == '__main__':
    unittest.main()
