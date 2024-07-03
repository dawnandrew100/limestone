from __future__ import annotations
import unittest
from limestone.editdistance import hammingDist

class TestHamming(unittest.TestCase):    
    def test_distance_diff(self):
        dist = hammingDist.distance("ACTG","FHYU")
        self.assertEqual(dist, 4.0)

    def test_similarity_diff(self):
        sim = hammingDist.similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)

    def test_norm_distance_diff(self):
        dist = hammingDist.normalized_distance("ACTG","FHYU")
        self.assertEqual(dist, 1.0)

    def test_norm_similarity_diff(self):
        sim = hammingDist.normalized_similarity("ACTG","FHYU")
        self.assertEqual(sim, 0.0)
    
    def test_distance_sim(self):
        dist = hammingDist.distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_similarity_sim(self):
        sim = hammingDist.similarity("ACTG","ACTG")
        self.assertEqual(sim, 4.0)

    def test_norm_distance_sim(self):
        dist = hammingDist.normalized_distance("ACTG","ACTG")
        self.assertEqual(dist, 0.0)

    def test_norm_similarity_sim(self):
        sim = hammingDist.normalized_similarity("ACTG","ACTG")
        self.assertEqual(sim, 1.0)
    
    def test_norm_distance1(self):
        dist = hammingDist.normalized_distance("ACTG","AATG")
        self.assertEqual(dist, 0.25)

    def test_norm_distance2(self):
        dist = hammingDist.normalized_distance("ACTG","AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_distance3(self):
        dist = hammingDist.normalized_distance("ACTG","AAAA")
        self.assertEqual(dist, 0.75)
 
    def test_norm_similarity1(self):
        dist = hammingDist.normalized_similarity("ACTG","AATG")
        self.assertEqual(dist, 0.75)

    def test_norm_similarity2(self):
        dist = hammingDist.normalized_similarity("ACTG","AAAG")
        self.assertEqual(dist, 0.5)

    def test_norm_similarity3(self):
        dist = hammingDist.normalized_similarity("ACTG","AAAA")
        self.assertEqual(dist, 0.25)

    def test_diff_len(self):
        dist = hammingDist.distance("ACTG","AATGA")
        self.assertEqual(dist,2.0)

    def test_diff_len2(self):
        dist = hammingDist.distance("AATGA","ACTG")
        self.assertEqual(dist,2.0)

    def test_binary_diff(self):
        dist = hammingDist.binary_distance_array("ACTG","AATG")
        ans = [1,0,1,1]
        self.assertEqual(dist, ans)

    def test_binary_sim(self):
        dist = hammingDist.binary_similarity_array("ACTG","AATG")
        ans = [0,1,0,0]
        self.assertEqual(dist, ans)

    def test_align1(self):
        dist = hammingDist.align("ACTG","ATGA")
        ans = f"ACTG\nATGA"
        self.assertEqual(dist,ans)

    def test_align2(self):
        dist = hammingDist.align("ACTGAA","ATGA")
        ans = f"ACTGAA\nATGA"
        self.assertEqual(dist,ans)

if __name__ == '__main__':
    unittest.main()
