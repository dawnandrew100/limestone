from __future__ import annotations
import unittest
from limestone import longest_common_subsequence

class TestLCS(unittest.TestCase):
    def test_distance_diff1(self):
        dist = longest_common_subsequence.distance("ACTG", "FHYU")
        self.assertEqual(dist, 4.0)

    def test_distance_diff2(self):
        dist = longest_common_subsequence.distance("ACTG", "AHYG")
        self.assertEqual(dist, 2.0)

    def test_distance_diff3(self):
        dist = longest_common_subsequence.distance("ACTG", "ACYU")
        self.assertEqual(dist, 2.0)

    def test_distance_diff4(self):
        dist = longest_common_subsequence.distance("ACTG", "AHYU")
        self.assertEqual(dist, 3.0)

    def test_similarity_diff(self):
        sim = longest_common_subsequence.similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = longest_common_subsequence.normalized_distance("ACTG", "FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = longest_common_subsequence.normalized_similarity("ACTG", "FHYU")
        self.assertEqual(sim, 0.0)

    def test_distance_sim(self):
        dist = longest_common_subsequence.distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = longest_common_subsequence.similarity("ACTG", "ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = longest_common_subsequence.normalized_distance("ACTG", "ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_distance1(self):
        dist = longest_common_subsequence.normalized_distance("ACTG", "ACTF")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = longest_common_subsequence.normalized_distance("ACTG", "ACFF")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = longest_common_subsequence.normalized_distance("ACTG", "AFFF")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity_sim(self):
        sim = longest_common_subsequence.normalized_similarity("ACTG", "ACTG")
        self.assertEqual(sim, 1.0)

    def test_norm_similarity1(self):
        dist = longest_common_subsequence.normalized_similarity("ACTG", "ACTF")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = longest_common_subsequence.normalized_similarity("ACTG", "ACFF")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = longest_common_subsequence.normalized_similarity("ACTG", "AFFF")
        self.assertEqual(dist, 0.25) 

    def test_align(self):
        alignment = longest_common_subsequence.align("BA", "ABA")
        self.assertEqual(alignment, "BA")

    def test_align2(self):
        alignment = longest_common_subsequence.align("AGTACGCA", "TATGC")
        self.assertEqual(alignment, "ATGC")
    
    def test_align3(self):
        alignment = longest_common_subsequence.align("ATCGATCGATGCTAGCTA", "ATGATCGAGCTA")
        self.assertEqual(alignment, "ATGATCGAGCTA")
    
    def test_align4(self):
        alignment = longest_common_subsequence.align("ATCGATCGATGCTAGCTA", "ATGATCGAGCTAVJONVJNV")
        self.assertEqual(alignment, "ATGATCGAGCTA")

if __name__ == '__main__':
    unittest.main()
