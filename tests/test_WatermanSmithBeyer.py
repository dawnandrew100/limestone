from __future__ import annotations
import unittest
from limestone.editdistance import watermanSmithBeyer

class TestWatermanSmithBayer(unittest.TestCase):
    
    def test_distance_diff(self):
        dist = watermanSmithBeyer.distance("ACTG","FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = watermanSmithBeyer.similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = watermanSmithBeyer.normalized_distance("ACTG","FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = watermanSmithBeyer.normalized_similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)
    
    def test_distance_sim(self):
        dist = watermanSmithBeyer.distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = watermanSmithBeyer.similarity("ACTG","ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = watermanSmithBeyer.normalized_distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = watermanSmithBeyer.normalized_similarity("ACTG","ACTG")
        self.assertEqual(sim, 1.0)
 
    def test_norm_distance1(self):
        dist = watermanSmithBeyer.normalized_distance("ACTG","AATG")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = watermanSmithBeyer.normalized_distance("ACTG","AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = watermanSmithBeyer.normalized_distance("ACTG","AAAA")
        self.assertEqual(dist, 0.75)
 
    def test_norm_similarity1(self):
        dist = watermanSmithBeyer.normalized_similarity("ACTG","AATG")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = watermanSmithBeyer.normalized_similarity("ACTG","AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = watermanSmithBeyer.normalized_similarity("ACTG","AAAA")
        self.assertEqual(dist, 0.25)

    def test_align(self):
        alignment = watermanSmithBeyer.align("BA", "ABA")
        self.assertEqual(alignment, "-BA\nABA")

if __name__ == '__main__':
    unittest.main()
