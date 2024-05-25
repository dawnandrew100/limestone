from __future__ import annotations
import unittest
from limestone import editdistance

class TestLevenshtein(unittest.TestCase):
    
    def test_distance_diff(self):
        dist = editdistance.levenshteinDist.distance("ACTG","FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = editdistance.levenshteinDist.similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = editdistance.levenshteinDist.normalized_distance("ACTG","FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = editdistance.levenshteinDist.normalized_similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)
    
    def test_distance_sim(self):
        dist = editdistance.levenshteinDist.distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = editdistance.levenshteinDist.similarity("ACTG","ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = editdistance.levenshteinDist.normalized_distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = editdistance.levenshteinDist.normalized_similarity("ACTG","ACTG")
        self.assertEqual(sim, 1.0)
    
    def test_align(self):
        alignment = editdistance.levenshteinDist.align("BA", "ABA")
        self.assertEqual(alignment, "-BA\nABA")

if __name__ == '__main__':
    unittest.main()
