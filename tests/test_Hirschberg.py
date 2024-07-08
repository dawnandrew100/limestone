from __future__ import annotations
import unittest
from limestone.editdistance import needleman_wunsch
from limestone.editdistance import hirschberg

class TestHirschberg(unittest.TestCase):
        
    def test_align1(self):
        nwalignment = needleman_wunsch.align("BA", "ABA")
        hsalignment = hirschberg.align("BA","ABA")
        self.assertEqual(nwalignment, hsalignment)
    
    def test_align2(self):
        nwalignment = needleman_wunsch.align("ATCG","ATCG")
        hsalignment = hirschberg.align("ATCG","ATCG")
        self.assertEqual(nwalignment, hsalignment)

    def test_align3(self):
        nwalignment = needleman_wunsch.align("ACTG","FHYU")
        hsalignment = hirschberg.align("ACTG","FHYU")
        self.assertEqual(nwalignment, hsalignment)

    def test_align4(self):
        nwalignment = needleman_wunsch.align("AGTACGCA","TATGC")
        hsalignment = hirschberg.align("AGTACGCA","TATGC")
        self.assertEqual(nwalignment, hsalignment)

    def test_align5(self):
        nwalignment = needleman_wunsch.align("TACOBELL","TABELL")
        hsalignment = hirschberg.align("TACOBELL","TABELL")
        self.assertEqual(nwalignment, hsalignment)

    def test_align6(self):
        nwalignment = needleman_wunsch.align("FASTFOAODKING","FASTAKING")
        hsalignment = hirschberg.align("FASTFOAODKING","FASTAKING")
        self.assertEqual(nwalignment, hsalignment)


if __name__ == '__main__':
    unittest.main()
