from __future__ import annotations
import unittest
from goombay import needleman_wunsch

class TestNeedlemanWunsch(unittest.TestCase):
    def test_distance_diff(self):
        dist = needleman_wunsch.distance("ACTG", "FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = needleman_wunsch.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = needleman_wunsch.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = needleman_wunsch.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = needleman_wunsch.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = needleman_wunsch.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = needleman_wunsch.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = needleman_wunsch.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance1(self):
        dist = needleman_wunsch.normalized_distance("ACTG", "AATG")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = needleman_wunsch.normalized_distance("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = needleman_wunsch.normalized_distance("ACTG", "AAAA")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity1(self):
        dist = needleman_wunsch.normalized_similarity("ACTG", "AATG")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = needleman_wunsch.normalized_similarity("ACTG", "AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = needleman_wunsch.normalized_similarity("ACTG", "AAAA")
        self.assertEqual(dist, 0.25)

    def test_align(self):
        alignment = needleman_wunsch.align("BA", "ABA")
        self.assertEqual(alignment, "-BA\nABA")

if __name__ == '__main__':
    unittest.main()
