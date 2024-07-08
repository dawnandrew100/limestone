# Limestone
This project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Smith-Waterman, and Waterman-Smith-Beyer algorithms. 

***Please ensure that numpy is installed so that this project can work correctly***

# Implementation

**Below is a table of the features of each algorithm.**

| Algorithm          | Alignment | Matrices | Distance/Similarity/Normalized |
| ------------------ | --------- | -------- | ------------------------------ |
|Needleman-Wunsch    |    [x]    |    [x]   |               [x]              |
|Smith-Waterman      |    [x]    |    [x]   |               [x]              |
|Waterman-Smith-Beyer|    [x]    |    [x]   |               [x]              |
|Levenshtein         |    [x]    |    [x]   |               [x]              |
|Hamming             |    [x]    |    [ ]   |               [x]              |
|Hirschberg          |    [x]    |    [ ]   |               [ ]              |

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

My Waterman-Smith-Beyer implementation does not align with that of [Freiburg University](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer), the site I've been using for alignment validation.
It is possible that their implementation has an issue and not mine but I wanted to mention this here and provide the link to my [StackOverflow](https://bioinformatics.stackexchange.com/questions/22683/waterman-smith-beyer-implementation-in-python) question for the sake of posterity.
