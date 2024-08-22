from __future__ import annotations
import unittest
from limestone import jaro

class TestJaro(unittest.TestCase):
    def test_distance_diff(self):
        dist = jaro.distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_similarity_diff(self):
        sim = jaro.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = jaro.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = jaro.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = jaro.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = jaro.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance_sim(self):
        dist = jaro.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = jaro.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_distance1(self):
        dist = jaro.normalized_distance("ACTG", "AATG")
        self.assertEqual(dist, 0.17)

    def test_norm_distance2(self):
        dist = jaro.normalized_distance("ACTG", "AAAG")
        self.assertEqual(dist, 0.33)

    def test_norm_distance3(self):
        dist = jaro.normalized_distance("ACTG", "AAAA")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity1(self):
        dist = jaro.normalized_similarity("ACTG", "AATG")
        self.assertEqual(dist, 0.83)

    def test_norm_similarity2(self):
        dist = jaro.normalized_similarity("ACTG", "AAAG")
        self.assertEqual(dist, 0.67)

    def test_norm_similarity3(self):
        dist = jaro.normalized_similarity("ACTG", "AAAA")
        self.assertEqual(dist, 0.5)

if __name__ == '__main__':
    unittest.main()
