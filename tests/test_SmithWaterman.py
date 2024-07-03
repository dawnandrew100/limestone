from __future__ import annotations
import unittest
from limestone.editdistance import smithWaterman

class TestWatermanSmithBayer(unittest.TestCase):
    
    def test_distance_diff(self):
        dist = smithWaterman.distance("ACTG","FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = smithWaterman.similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = smithWaterman.normalized_distance("ACTG","FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = smithWaterman.normalized_similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)
    
    def test_distance_sim(self):
        dist = smithWaterman.distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = smithWaterman.similarity("ACTG","ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = smithWaterman.normalized_distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)
 
    def test_norm_distance1(self):
        dist = smithWaterman.normalized_distance("ACTG","ACTF")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = smithWaterman.normalized_distance("ACTG","ACFF")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = smithWaterman.normalized_distance("ACTG","AFFF")
        self.assertEqual(dist, 0.75)
    
    def test_norm_similarity_sim(self):
        sim = smithWaterman.normalized_similarity("ACTG","ACTG")
        self.assertEqual(sim, 1.0)
 
    def test_norm_similarity1(self):
        dist = smithWaterman.normalized_similarity("ACTG","ACTF")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = smithWaterman.normalized_similarity("ACTG","ACFF")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = smithWaterman.normalized_similarity("ACTG","AFFF")
        self.assertEqual(dist, 0.25) 

    def test_align(self):
        alignment = smithWaterman.align("BA", "ABA")
        self.assertEqual(alignment, "BA\nBA")
    
    def test_align2(self):
        alignment = smithWaterman.align("AGTACGCA","TATGC")
        self.assertEqual(alignment, "TACGC\nTATGC")

if __name__ == '__main__':
    unittest.main()
