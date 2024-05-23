from __future__ import annotations
import unittest
from limestone import textdistance

class TestWatermanSmithBayer(unittest.TestCase):
    
    def test_distance_diff(self):
        dist = textdistance.smithWaterman.distance("ACTG","FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = textdistance.smithWaterman.similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = textdistance.smithWaterman.normalized_distance("ACTG","FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = textdistance.smithWaterman.normalized_similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)
    
    def test_distance_sim(self):
        dist = textdistance.smithWaterman.distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = textdistance.smithWaterman.similarity("ACTG","ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = textdistance.smithWaterman.normalized_distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = textdistance.smithWaterman.normalized_similarity("ACTG","ACTG")
        self.assertEqual(sim, 1.0)
    
    def test_align(self):
        alignment = textdistance.smithWaterman.align("BA", "ABA")
        self.assertEqual(alignment, "BA\nBA")

if __name__ == '__main__':
    unittest.main()
