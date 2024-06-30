# Limestone
Sequence aligner that can produce scoring matrices for Needleman-Wunsch, Smith-Waterman, and Waterman-Smith-Beyer algorithms. 
This project can also produce distance, similarity, normalised distance, and normalised similarity scores for all alignment algorithms with full scoring matrices.

Due to the recursive nature of the Hirschberg algorithm, if a distance score or matrix is needed it is best to use the Needleman-Wunsch algorithm instead.

This project uses the numpy library in order to perform matrix manipulation that all of the underlining scores and matrices are derived from. 
