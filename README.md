# Limestone
This project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Smith-Waterman, Wagner-Fischer, Waterman-Smith-Beyer, Wagner-Fischer, Lowrance-Wagner, Longest Common Subsequence, and Shortest Common Supersequence algorithms. 

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
|Jaro                          |    [ ]    |    [x]   |               [x]              |
|Jaro Winkler                  |    [ ]    |    [x]   |               [x]              |
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

[Jaro & Jaro-Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance)

[Longest Common Subsequence](https://en.wikipedia.org/wiki/Longest_common_subsequence)

[Shortest Common Supersequence](https://en.wikipedia.org/wiki/Shortest_common_supersequence)

# Code Examples

**Hamming Distance**
```python
from limestone import hamming

qs = "AFTG"
ss = "ACTG"

print(hamming.distance(qs, ss))
# 1
print(hamming.similarity(qs, ss))
# 3
print(hamming.binary_distance_array(qs, ss))
# [0,1,0,0]
print(hamming.binary_similarity_array(qs, ss))
# [1,0,1,1]
print(hamming.normalized_distance(qs, ss))
# 0.25
print(hamming.normalized_similarity(qs, ss))
# 0.75
```

**Needleman-Wunsch**
```python
from limestone import needleman_wunsch

print(needleman_wunsch.distance("ACTG","FHYU"))
# 4
print(needleman_wunsch.distance("ACTG","ACTG"))
# 0
print(needleman_wunsch.similarity("ACTG","FHYU"))
# 0
print(needleman_wunsch.similarity("ACTG","ACTG"))
# 4
print(needleman_wunsch.normalized_distance("ACTG","AATG"))
#0.25
print(needleman_wunsch.normalized_similarity("ACTG","AATG"))
#0.75
print(needleman_wunsch.align("BA","ABA"))
#-BA
#ABA
print(needleman_wunsch.matrix("AFTG","ACTG"))
[[0. 2. 4. 6. 8.]
 [2. 0. 2. 4. 6.]
 [4. 2. 1. 3. 5.]
 [6. 4. 3. 1. 3.]
 [8. 6. 5. 3. 1.]]
 ```

# Work In Progress

Jaro and Jaro-Winkler algorithms.
Importing and parsing FASTA, FASTQ, and PDB files.

# Caveats

Due to the recursive nature of the Hirschberg algorithm, if a distance score or matrix is needed it is best to use the Needleman-Wunsch algorithm instead.

Note that due to the fact that the Hamming distance does not allow for insertions, or deletions, the "aligned sequence" that is returned is just the original sequences in a formatted string. 
This is due to the fact that actually aligning the two sequences using this algorithm would just lead to two lines of the query sequence. 
It should also be noted that the Hamming distance is intended to only be used with sequences of the same length. 
To compensate for strings of differing lengths, my algorithm adds 1 extra point to the distance for every additional letter in the longer sequence since this can be seen as "swapping" the empty space for a letter or vice versa. However, any distance obtained this way **will not reflect an accurate Hamming distance**.

My Waterman-Smith-Beyer implementation does not always align with that of [Freiburg University](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer), the site I've been using for alignment validation.
It is possible that their implementation has an issue and not mine but I wanted to mention this here and provide the link to my [StackOverflow](https://bioinformatics.stackexchange.com/questions/22683/waterman-smith-beyer-implementation-in-python) question for the sake of posterity.

During the beginning of this project I thought that the Levenshtein distance was an algorithm, but it is the end result that is being calculated with an approach such as Wagner-Fischer which uses Needleman-Wunsch-esque matrices to calculate the Levenshtein distance.
Thusly, the Levenshtein distance implementation has been switched with the Wagner-Fischer algorithm.
Damerau-Levenshtein distance is found using the Lowrance-Wagner algorithm.
