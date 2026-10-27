# BINF6250

1. Read the input sequences from the FASTA file.
       a. Each header line starts with '>', and we keep only the last word as the ID (DM1, P10, ...).
       b. Sequence lines get joined together until the next header shows up.
       c. The last sequence is saved after the loop, since no header comes after it.

2. Store the sequences in a dictionary mapping each ID to its sequence string.

3. Compute the distance matrix:
       a. Start with an n x n matrix of zeros (the diagonal stays 0, since a sequence is 0 distance from itself).
       b. For each pair of sequences (i, j) in the upper triangle only (j starts at i + 1):
           Align the two sequences with our week 5 Smith-Waterman (cal_score from utils.py),
               passing in the notebook's scoring values (match 2, mismatch -1, gap -1).
           Find the highest scoring cell anywhere in the matrix (not the corner, since it's a local alignment).
           Trace back from that cell to get the two aligned sequences (traceback from utils.py).
           Compute the normalized Hamming distance between the aligned sequences.
           Store the same distance at (i, j) and (j, i), since distance is symmetric.

4. Make a list of labels, one per sequence. At the start each label is just the sequence ID,
   but as nodes get joined, each label becomes the Newick text for its whole subtree.

5. While n > 2:
       Compute the row totals (r_i) and the Q-matrix from the current distance matrix.
       Find the pair of nodes (i, j) with the smallest Q value (on a tie, the first pair found wins).
       Calculate the branch lengths from i and j to the new internal node k (d_i and d_j),
           before touching the matrix, since this step needs the old values.
       If a branch length comes out negative, set it to 0.
       Make the label for node k: (label_i:d_i,label_j:d_j), with branch lengths rounded to 3 decimals in the text only.
       Make a new matrix of size n - 1.
       Copy the distances between the remaining nodes into the new matrix.
       Calculate the distances between k and every remaining node, and put them in the last row and column.
       Replace the old distance matrix with the new one.
       Update the label list in the same order: remaining labels first, then k's label at the end.
       Update n -> n - 1.

6. When 2 nodes are left, join them with one branch of length d(0,1), putting the whole length
   on the first node and 0 on the second, then close the string with ';'.
   (This closing point is only where the Newick string ends, it's not a real biological root.)

7. Plot the phylogenetic tree from the Newick string.