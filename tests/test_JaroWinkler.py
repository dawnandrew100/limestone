from __future__ import annotations
import unittest
from limestone import jaro_winkler

class TestJaroWinkler(unittest.TestCase):
    def test_distance_diff(self):
        dist = jaro_winkler.distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_similarity_diff(self):
        sim = jaro_winkler.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = jaro_winkler.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = jaro_winkler.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = jaro_winkler.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = jaro_winkler.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance_sim(self):
        dist = jaro_winkler.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = jaro_winkler.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance1(self):
        dist = jaro_winkler.normalized_distance("ACTG", "AATG")
        self.assertEqual(dist, 0.15)

    def test_norm_distance2(self):
        dist = jaro_winkler.normalized_distance("ACTG", "AAAG")
        self.assertEqual(dist, 0.3)

    def test_norm_distance3(self):
        dist = jaro_winkler.normalized_distance("ACTG", "AAAA")
        self.assertEqual(dist, 0.45)

    def test_norm_similarity1(self):
        dist = jaro_winkler.normalized_similarity("ACTG", "AATG")
        self.assertEqual(dist, 0.85)

    def test_norm_similarity2(self):
        dist = jaro_winkler.normalized_similarity("ACTG", "AAAG")
        self.assertEqual(dist, 0.7)

    def test_norm_similarity3(self):
        dist = jaro_winkler.normalized_similarity("ACTG", "AAAA")
        self.assertEqual(dist, 0.55)

if __name__ == '__main__':
    unittest.main()
