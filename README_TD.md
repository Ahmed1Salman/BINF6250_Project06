# Introduction
In this project, we implemented **neighbor joining (NJ)**, a **distance-based phylogeny** method, in Python, and used it to build an **unrooted tree** of 20 HIV-1 reverse transcriptase (RT) protein sequences.

A phylogeny is a hypothesis about how a set of sequences descended from common ancestors. Distance-based methods first squeeze every pair of sequences into one number (how different they are), then build a tree that fits those numbers. The simplest distance method, UPGMA, always merges the two closest clusters, which assumes every lineage evolves at the same rate (a molecular clock). NJ does not assume a clock. Instead of joining the closest pair, it corrects each distance by how far each node is from everyone else, and joins the pair with the smallest corrected value (Saitou & Nei, 1987; Studier & Keppler, 1988).

Our program does this in three stages:
1. **Distance table:** read the sequences, score every pair with Smith-Waterman local alignment (adapted from our Project 05 code), and convert those similarity scores into distances.
2. **Build the tree:** run the NJ loop, joining one pair per round and shrinking the distance table by one each time, until two nodes remain.
3. **Newick and drawing:** write the finished tree as a Newick string and draw it with Biopython and Matplotlib.


# Documentation and Results
### Requirements
- Python 
- Jupyter Notebook
- numpy, for the score and distance matrices
- biopython, for reading the Newick string and drawing the tree (`Bio.Phylo`)
- matplotlib, for displaying the tree

The starter notebook suggests this environment:
```bash
micromamba create -n module06 -c bioconda -c conda-forge python numpy matplotlib biopython
```

### Usage
Open `project06.ipynb` in Jupyter and run all cells in order, top to bottom. `lafayette_SARS_RT.fasta` must be in the same folder as the notebook.\
The notebook is split into:\
Imports: `typing`, `StringIO`, `matplotlib.pyplot`, `Bio.Phylo`, `numpy`\
Stage 1 functions: `read_fasta()`, `smith_waterman()`, `build_distance_matrix()`\
Stage 2 function: `neighbor_joining()`\
Stage 3 functions: our Newick helper, called at the end of `neighbor_joining()`, and `plot_tree()` (provided in the starter notebook)\
Main: reads the sequences, builds the distance matrix, runs NJ, and plots the tree

### Project structure
```
project06/
├── project06.ipynb
├── lafayette_SARS_RT.fasta
├── README.md
```

### Data
`lafayette_SARS_RT.fasta` holds **20 HIV-1 reverse transcriptase protein sequences** (despite "SARS" in the file name). Each header ends with a short label, which we use as the sequence name: DM1–DM3, JT1–JT2, and P1–P15. Example header:
```
>gi|24209940|gb|AAN41462.1| reverse transcriptase, partial [Human immunodeficiency virus 1] DM1
```

### Settings
These are the parameters of `smith_waterman()`, as named in the Project 06 starter notebook:\
Setting- Default (Meaning)\
`sequence1`- none (first sequence)\
`sequence2`- none (second sequence)\
`match_score`- `2.0` (added on the diagonal when the two letters are the same)\
`mismatch_penalty`- `-1.0` (added on the diagonal when the two letters are different)\
`gap_penalty`- `-1.0` (added for an up or left move)

### Distance transformation
Smith-Waterman gives a **similarity** (bigger means more alike), but NJ needs a **distance** (bigger means more different, 0 for identical). Raw local alignment scores also grow with sequence length.

Checks against the properties of a distance:
- **Identity:** 
- **Symmetry:**
- **Non-negativity:** 
- **Triangle inequality:** 

### Output
`build_distance_matrix()` returns the 20 × 20 distance matrix and the list of sequence names.\
`neighbor_joining()` returns the tree as a Newick string.\
`plot_tree()` draws the tree.

### Note about the materials: 
Sequence count: the build_distance_matrix() docstring expects len(ids) == 19 and the neighbor_joining() docstring expects 18 commas, but the FASTA file has 20 sequences (a binary tree with 20 leaves has 19 commas)\
Drawing library: the notebook's background says the tree is visualized with ete3, but plot_tree() uses Biopython's Phylo.

# Pseudocode
## Stage 1: Distance table
1. Read the sequences
```
def read_fasta(filename):
    #1. Set up an empty dictionary {name: sequence} and "no current name yet"
    sequences = empty dictionary
    current_name = nothing

    #2. Read the file one line at a time
    for each line in the file:
        # remove the newline and spaces at both ends
        line = line stripped
        # blank line: nothing to read, skip it
        if line is empty: skip it
        # header line: the name is the LAST word (e.g. "DM1", "P1")
        else if line starts with ">":
            current_name = last word of line
            # start an empty sequence for this name
            sequences[current_name] = ""
        # sequence line: glue it onto the current sequence (each sequence spans 4 lines)
        else:
            sequences[current_name] = sequences[current_name] + line

    #3. Give back the dictionary
    return sequences
```
2. Score two sequences (adapted from Project 05: same fill loop and cell logic as `cal_score()`; traceback, printing, and random tie-break removed because only the score is needed)
a. Raw similarity (number)
```
def raw_smith_waterman(seq1, seq2, match_score=2.0, mismatch_penalty=-1.0, gap_penalty=-1.0):
    #1. Set up the score grid; +1 adds the zero border row and column
    score_matrix = zeros of size (len(seq1)+1) x (len(seq2)+1)

    #2. Walk through each cell, row by row, top to bottom, left to right
    for i from 1 to len(seq1):
        for j from 1 to len(seq2):
            # look up the cell's two letters (i-1, j-1 because of the border)
            if seq1[i-1] == seq2[j-1]: letter_score = match_score
            else:                      letter_score = mismatch_penalty
            # the four options, keep the biggest
            # option 1: up-left box + match/mismatch score
            diag = score_matrix[i-1][j-1] + letter_score
            # option 2: box above + gap penalty
            up   = score_matrix[i-1][j]   + gap_penalty
            # option 3: box to the left + gap penalty
            left = score_matrix[i][j-1]   + gap_penalty
            # option 4: 0 = start over (this is what makes it LOCAL)
            score_matrix[i][j] = max(0, diag, up, left)

    #3. The score is the BIGGEST cell anywhere in the grid, tie-break do not matter here
    return the biggest number in score_matrix
```
b. Normalized similarity (%) 
*Note: This is a simple normalization design choice I had picked after compared it with AI provided Feng & Doolittle (1987) approach's code. See notebook for detail but as for reasoning, I added check_distance_properties() using the criterias in the slide (refered in ### Distance transformation section of this README file). The result is saying Feng & Doolittle has Triangle inequality violations where the simple version doesn't so we went with that one. 
```
def smith_waterman(sequence1, sequence2, match_score=2.0, mismatch_penalty=-1.0, gap_penalty=-1.0):
    #1. Raw score for the pair
    s_ab = raw_smith_waterman(sequence1, sequence2, ...)      

    #2. Each sequence's perfect score (aligned with itself)
    s_aa = raw_smith_waterman(sequence1, sequence1, ...)      
    s_bb = raw_smith_waterman(sequence2, sequence2, ...)      

    #3.Compare to perfect -> 1.0 for identical sequences
    return s_ab / sqrt(s_aa * s_bb)                           
```
3. Build the distance matrix
```
def build_distance_matrix(sequences, match_score=2.0, mismatch_penalty=-1.0, gap_penalty=-1.0):
    #1. List the names, count them, and make an N x N table of 0s
    ids = list of names
    N = number of ids
    distance = N x N zeros            # the diagonal stays 0 (identity)

    #2. Every pair once (j starts after i, because the table is symmetric)
    for i from 0 to N-1:
        for j from i+1 to N-1:
            # how alike: smith_waterman does STEP 1 + STEP 2
            similarity = smith_waterman(seq[ids[i]], seq[ids[j]], ...)    

            # flip "how alike" into "how different"
            d = 1 - similarity                                           

            # guard: a distance is never negative
            if d < 0: d = 0

            # write it in BOTH halves (symmetry)
            distance[i][j] = d
            distance[j][i] = d

    #3. Give back the table and the names, in the same order as the rows
    return distance, ids
```
## Stage 2: Build the tree
```
def neighbor_joining(distance_matrix, labels, seed=None):
    
    #1. Input: 
    # Make a copy of the distance table so the original stays the same
    D = copy of distance_matrix
    # Set up a random picker with a seed, so ties are broken the same way every run
    rng = random number maker using seed

    #2. Initialization:
    # Start with every sequence as its own leaf
    # Keep the names in the same order as the rows of D
    active = copy of labels
    # Start an empty tree
    # Each new node we create will remember the two nodes it joined and their branch lengths
    T = empty dictionary
    # Count how many nodes are left
    n = number of active
    
    #3. Iteration: 
    # Keep joining pairs until only 2 nodes are left
    while n > 2:

        # For each node, add up its distances to everyone else
        r[i] = sum of row i of D, for every i

        # Calculate Q for every pair (this finds who is close to each other, not just close to everyone)
        Q[i][j] = (n - 2) * D[i][j] - r[i] - r[j], for every pair i < j

        # Find the smallest Q
        best_Q = smallest Q value

        # Find every pair that has that smallest Q (use "close to", because decimals are not stored exactly)
        tied_pairs = every (i, j) where Q[i][j] is close to best_Q

        # If there is a tie, pick one pair at random
        (i, j) = random pick from tied_pairs using rng

        # Work out how long the branch is from the new node to i, and to j
        d_i = (D[i][j] + (r[i] - r[j]) / (n - 2)) / 2
        d_j = D[i][j] - d_i

        # A negative branch = 0
        if d_i < 0: d_i = 0
        if d_j < 0: d_j = 0

        # Make a new node k that joins i and j, and save its two branches in the tree
        k = new unique name, e.g. "node1"
        T[k] = [(active[i], d_i), (active[j], d_j)]

        # Work out the distance from the new node k to every other node x
        for each x that is not i or j:
            d(k,x) = (D[i][x] + D[j][x] - D[i][j]) / 2

        # Remove i and j from the table, and add k as the last row and column
        D = D without rows/cols i and j, plus d(k,x) as the new last row/col

        # Do the same to the list of names, so names and rows still match
        active = active without i and j, plus k at the end

        # One fewer node is left now
        n = n - 1

    #3. Termination: Only 2 nodes are left: connect them with one last branch and write the tree as text
    return to_newick(T, active[0], active[1], D[0][1])
```
## Stage 3: Newick and drawing
*Note: Trang's version of pseudo code and code is recursive. Do not know which one we ended up using but I put some reasoning in my reflection in term of why I think Ahmed's verison is better. 
```
def node_string(T, node):

    # If the node is a leaf (not in the tree dictionary), just write its name
    if node is not in T:
        return node

    # Otherwise, write each of its two branches as "name:length"
    # If the other end is also a new node, write that one first (the function calls itself)
    parts = empty list
    for each (other_node, length) in T[node]:
        add  node_string(T, other_node) + ":" + length  to parts

    # Put the two parts in brackets, separated by a comma
    return "(" + parts joined with "," + ")"


def to_newick(T, i, j, d_ij):

    # Make sure i is the new node (swap them if i is a leaf and j is not)
    if i is not in T and j is in T:
        swap i and j

    # If both are leaves, there were only 2 sequences, so write them as a simple pair
    if i is not in T:
        return "(" + i + ":" + d_ij + "," + j + ":0);"

    # Write the two branches of node i
    parts = each branch of T[i] written as  node_string(T, other_node) + ":" + length

    # Add the last branch to j as a third item
    add  node_string(T, j) + ":" + d_ij  to parts

    # Put everything in brackets and end with ";"
    return "(" + parts joined with "," + ");"

def plot_tree(newick_tree: str): #Provided - no changes 
```


# Successes

# Struggles

**Note on overall limitations of our program:**
- **Scoring is not built for protein.** As we noted in Project 05, our Smith-Waterman gives every match the same +2 and every mismatch the same −1. Real amino acids are not equally different: some substitutions (like two similar amino acids) are common and mild, others are rare and drastic. A substitution matrix such as BLOSUM62 would score these differently. We kept the uniform scoring because the starter notebook's defaults use it, and because our self-score normalization would still work if we switched.

**Things we still don't understand:**

# Personal Reflections
## Group Leader

## Other members
Trang Do: did not see Trivedi version yet, will update and talk about what I learnt from his later. 
Meeting as a group really helped. We aligned on the math and went through the additional materials together, and we kept communication open all week so we could help each other out.

I really liked that our individual versions turned out different this time. I went in assuming recursion was the way to build the Newick string, because of the tree’s shape, even though we’d discussed in Project 04’s peer review that Python might discourage recursion for a reason. Ahmed went iterative and got the exact same tree. Recursion is fine for our 20 sequences, but his approach avoids Python’s recursion limit on much bigger or lopsided trees, so I think his is the safer choice.

I also learned to test my choices instead of just picking them. I was curious about the Feng & Doolittle (1987) approach the assignment points to, so I had Claude write that version and ran both through the distance rules from the Module 06 slides (check_distance_properties()). Feng & Doolittle had 114 triangle-inequality violations and the simple version had none, so we went with the simple one. Both still gave the same tree, so the choice didn’t change our result.
# Generative AI Appendix
**Generative AI:**\
Trang Do: I used Claude Opus 5.5 (Anthropic) to help me understand the algorithms after reading the slides and the starter notebook. To test my understanding, I had it generate a small example (a 4-sequence toy dataset) that I worked through by hand: the Smith-Waterman grid, the distance table, and each round of neighbor joining. I also had it explain the Newick format and how to build a Newick string from the finished tree.
For similarity normalization: I tried two approaches: a simple normalization (score ÷ √(self-score × self-score), then distance = 1 − similarity), and the Feng & Doolittle (1987) approach the assignment points to (correcting for self-scores and shuffled sequences, then −ln). I had Claude wrote code Feng & Doolittle version and then provided the properties of a distance check from our Module 6 slides and had it wrote check_distance_properties(), comapre and then chose (detailed above). 
**Citations**\


