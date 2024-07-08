from __future__ import annotations
import unittest
from limestone.editdistance import levenshtein 

class TestLevenshtein(unittest.TestCase):
    def test_distance_diff(self):
        dist = levenshtein.distance("ACTG", "FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = levenshtein.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = levenshtein.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = levenshtein.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = levenshtein.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = levenshtein.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = levenshtein.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = levenshtein.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance1(self):
        dist = levenshtein.normalized_distance("ACTG", "AATG")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = levenshtein.normalized_distance("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = levenshtein.normalized_distance("ACTG", "AAAA")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity1(self):
        dist = levenshtein.normalized_similarity("ACTG", "AATG")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = levenshtein.normalized_similarity("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = levenshtein.normalized_similarity("ACTG", "AAAA")
        self.assertEqual(dist, 0.25)

    def test_align(self):
        alignment = levenshtein.align("BA", "ABA")
        self.assertEqual(alignment, "-BA\nABA")

if __name__ == '__main__':
    unittest.main()
