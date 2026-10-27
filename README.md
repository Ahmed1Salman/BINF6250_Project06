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
Citations:
https://primerdigital.com/tools/phylo.html
This website was used to test our code and see how it compares to a publicly available neighbor joining algorithm. I fed it the data file and compared our produced tree to this. There were some differences that I couldn't explain myself but the result seemed extremely close that we decided it was good enough

AI usage:
Ahmed Salman - Claude AI (Sonnet 5.5) was used to understand how to use 3 specific functions: `np.fill_diagonal`, `np.stack` and the v and h versions, and lastly `array.reshape`. Prompts used are:
"How to prevent `array.min` from selecting anything in the diagonal of a symmetrical numpy array in python?". This lead to `np.fill_diagonal`. 
"How to use `np.stack` to add a row or column to an existing array?". This lead to `np.vstack` and `np.hstack`. `np.stack` only takes the arrays you feed it without altering it and so it requires the inputs to already be 2d arrays if you are planning on adding them as columns. `np.vstack` automatically converts 1d arrays to 2d arrays and works on them, allowing me to shorten my code and "skip" a step.
For `array.reshape` I simply pasted the error I kept receiving after trying to run the code. Claude suggested using this which worked similarily to `np.vstack` and `np.hstack` but without adding anything, simply converting the current 1d array to a 2d array.
