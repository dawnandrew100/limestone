
from __future__ import annotations
import unittest
from goombay import lowrance_wagner 

class TestLevenshtein(unittest.TestCase):
    def test_distance_diff(self):
        dist = lowrance_wagner.distance("ACTG", "FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = lowrance_wagner.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = lowrance_wagner.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = lowrance_wagner.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = lowrance_wagner.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = lowrance_wagner.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = lowrance_wagner.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = lowrance_wagner.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance1(self):
        dist = lowrance_wagner.normalized_distance("ACTG", "AATG")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = lowrance_wagner.normalized_distance("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = lowrance_wagner.normalized_distance("ACTG", "AAAA")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity1(self):
        dist = lowrance_wagner.normalized_similarity("ACTG", "AATG")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = lowrance_wagner.normalized_similarity("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = lowrance_wagner.normalized_similarity("ACTG", "AAAA")
        self.assertEqual(dist, 0.25)

    def test_align1(self):
        alignment = lowrance_wagner.align("BA", "ABA")
        self.assertEqual(alignment, "-BA\nABA")

    def test_align2(self):
        alignment = lowrance_wagner.align("ACTG", "ATCG")
        self.assertEqual(alignment, "ACTG\nACTG")

    def test_align3(self):
        alignment = lowrance_wagner.align("ACTGGTAC", "ATCGGATC")
        self.assertEqual(alignment, "ACTGGTAC\nACTGGTAC")

if __name__ == '__main__':
    unittest.main()
