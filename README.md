# Limestone
This project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Smith-Waterman, Wagner-Fischer, and Waterman-Smith-Beyer algorithms. 

***Please ensure that numpy is installed so that this project can work correctly***

# Implementation

**Below is a table of the features of each algorithm.**

| Algorithm                    | Alignment | Matrices | Distance/Similarity/Normalized |
| ------------------           | --------- | -------- | ------------------------------ |
|Needleman-Wunsch              |    [x]    |    [x]   |               [x]              |
|Smith-Waterman                |    [x]    |    [x]   |               [x]              |
|Waterman-Smith-Beyer          |    [x]    |    [x]   |               [x]              |
|Wagner-Fischer                |    [x]    |    [x]   |               [x]              |
|Lowrance-Wagner               |    [x]    |    [x]   |               [x]              |
|Hamming                       |    [x]    |    [ ]   |               [x]              |
|Hirschberg                    |    [x]    |    [ ]   |               [ ]              |
|Longest Common Subsequence    |    [x]    |    [x]   |               [x]              |
|Shortest Common Supersequence |    [x]    |    [x]   |               [x]              |

## Algorithms Explained
[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)

[Smith-Waterman ](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)

[Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer)

[Wagner-Fischer](https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm) <- Levenshtein distance

[Lowrance-Wagner](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2819-0) <- Damerau–Levenshtein distance (Levenshtein distance plus adjacent swapping)

[Hamming](https://en.wikipedia.org/wiki/Hamming_distance)

[Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)

[Longest Common Subsequence](https://en.wikipedia.org/wiki/Longest_common_subsequence)

[Shortest Common Supersequence](https://en.wikipedia.org/wiki/Shortest_common_supersequence)

# Code Examples

**Hamming Distance**
```python
from limestone.editdistance import hammingDist

qs = "AFTG"
ss = "ACTG"

print(hammingDist.distance(qs, ss))
# 1
print(hammingDist.similarity(qs, ss))
# 3 
print(hammingDist.binary_distance_array(qs, ss))
# [1,0,1,1]
print(hammingDist.binary_similarity_array(qs, ss))
# [0,1,0,0]
print(hammingDist.normalized_distance(qs, ss))
# 0.25
print(hammingDist.normalized_similarity(qs, ss))
# 0.75
```

**Needleman-Wunsch**
```python
from limestone.editdistance import needlemanWunsch

print(needlemanWunsch.distance("ACTG","FHYU"))
# 4
print(needlemanWunsch.distance("ACTG","ACTG"))
# 0
print(needlemanWunsch.similarity("ACTG","FHYU"))
# 0
print(needlemanWunsch.similarity("ACTG","ACTG"))
# 4
print(needlemanWunsch.normalized_distance("ACTG","AATG"))
#0.25
print(needlemanWunsch.normalized_similarity("ACTG","AATG"))
#0.75
print(needlemanWunsch.align("BA","ABA"))
#-BA
#ABA
print(needlemanWunsch.matrix("AFTG","ACTG"))
[[0. 2. 4. 6. 8.]
 [2. 0. 2. 4. 6.]
 [4. 2. 1. 3. 5.]
 [6. 4. 3. 1. 3.]
 [8. 6. 5. 3. 1.]]
 ```

# Work In Progress

-- To be continued

# Caveats

Due to the recursive nature of the Hirschberg algorithm, if a distance score or matrix is needed it is best to use the Needleman-Wunsch algorithm instead.

Note that due to the fact that the Hamming distance does not allow for substitutions, insertions, or deletions, the "aligned sequence" that is returned is just the original sequences in a formatted string.

My Waterman-Smith-Beyer implementation does not always align with that of [Freiburg University](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer), the site I've been using for alignment validation.
It is possible that their implementation has an issue and not mine but I wanted to mention this here and provide the link to my [StackOverflow](https://bioinformatics.stackexchange.com/questions/22683/waterman-smith-beyer-implementation-in-python) question for the sake of posterity.

During the beginning of this project I thought that the Levenshtein distance was an algorithm, but it is the end result that is being calculated with an approach such as Wagner-Fischer which uses Needleman-Wunsch-esque matrices to calculate the Levenshtein distance.
Thusly, the Levenshtein distance implementation has been switched with the Wagner-Fischer algorithm.
Damerau-Levenshtein distance is found using the Lowrance-Wagner algorithm.
