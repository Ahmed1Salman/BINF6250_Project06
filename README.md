# BINF6250
Pseudo Code
```
def read_fasta (file name)
  read fasta file line by line
  initialize sequence list and dictionary
  if it is the first line, store as header
  loop through lines
  if line starts with '>', join the sequence list and store sequences and header as key value pair, then store current line as new header and reset the sequence list.
  else store line in sequence list
  after loop, join the sequence list one last time and store header and sequence as key value pair in dictionary
  return dictionary

def smith_waterman (both sequences, scoring algorithm)
  initialize score matrix with rows and columns equal to the length of the sequences
  loop through sequences and calculate score for each position
  return maximum score normalized by dividing over the maximum possible score

def build_distance_matrix (sequence dictionary)
  initialize sequence matrix with rows and columns equal to number of sequences in dictionary
  store dictionary keys as ids
  loop through the matrix and add smith waterman scores for every position based on sequences compared
  return sequence matrix, ids

def neighbor_joining (distance matrix, id list)
  Copy matrix to avoid destruction of original, increases the memory usage heavily, so remove if needed.
  Extract node names from labels list by removing extra information
  Loop through the nodes while the matrix shape is more than 2x2
  Sum rows from original matrix
  Calculate first D* iteration
  Replace the zeros in the diagonal with a very large number to prevent the min function from picking the diagonal.
  Min function to pick lowest value in matrix
  Calculate delta between nodes, and then limb lengths
  Calculate the distances between the new node created above with the rest of the matrix then add the new distances to the matrix
  Delete the nodes that have already been merged from the matrix
  Append the new name to the list, then remove the old names from the list based on index
  Subtract 1 from the number of nodes and repeat the process for all iterations
  Produce the final string that will be returned

```
