from __future__ import annotations
import unittest
from limestone.editdistance import needlemanWunsch
from limestone.editdistance import hirschbergDist

class TestHirschberg(unittest.TestCase):
        
    def test_align1(self):
        nwalignment = needlemanWunsch.align("BA", "ABA")
        hsalignment = hirschbergDist.align("BA","ABA")
        self.assertEqual(nwalignment, hsalignment)
    
    def test_align2(self):
        nwalignment = needlemanWunsch.align("ATCG","ATCG")
        hsalignment = hirschbergDist.align("ATCG","ATCG")
        self.assertEqual(nwalignment, hsalignment)

    def test_align3(self):
        nwalignment = needlemanWunsch.align("ACTG","FHYU")
        hsalignment = hirschbergDist.align("ACTG","FHYU")
        self.assertEqual(nwalignment, hsalignment)

    def test_align4(self):
        nwalignment = needlemanWunsch.align("AGTACGCA","TATGC")
        hsalignment = hirschbergDist.align("AGTACGCA","TATGC")
        self.assertEqual(nwalignment, hsalignment)

    def test_align5(self):
        nwalignment = needlemanWunsch.align("TACOBELL","TABELL")
        hsalignment = hirschbergDist.align("TACOBELL","TABELL")
        self.assertEqual(nwalignment, hsalignment)

    def test_align6(self):
        nwalignment = needlemanWunsch.align("FASTFOAODKING","FASTAKING")
        hsalignment = hirschbergDist.align("FASTFOAODKING","FASTAKING")
        self.assertEqual(nwalignment, hsalignment)


if __name__ == '__main__':
    unittest.main()
