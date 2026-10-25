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



```
