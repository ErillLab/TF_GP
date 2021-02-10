#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 09:41:42 2020

@author: elia
"""

# IMPORTS AND CONFIG ----------------------------------------------------------

from objects.organism_factory import OrganismFactory
import json
import numpy as np


def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content

    Args:
        filename: Name of the json file to read

    Returns:
        Dictionary with the json file info
    """

    with open(filename) as json_content:
        return json.load(json_content)


JSON_CONFIG_FILENAME = "config.json"
config = read_json_file(JSON_CONFIG_FILENAME)

configOrganism = config["organism"]
configOrganismFactory = config["organismFactory"]
configConnector = config["connector"]
configPssm = config["pssm"]


organism_factory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
)


# TEST GET PLACEMENT ----------------------------------------------------------

# Example DNA sequence
example_dna_seq = "AGAAGGAATGACCCAGGCCGGGCACGATACTAAACACAAAAAGACAAAGCAATATGCACTGGTAACAGCAAGTGTACGTAACACGTCAACCCAATTAAGCGACCGAAATCGTAAATTCTGATAGAGTCTATGGTGGAGGCATATCCAAGATCTAACCAATATACCACAAGCAGAGTCAGATTAATGTTCCGAACGGAAGA"
example_dna_seq = example_dna_seq.lower()

example_dna_seq2 = "AGAAGGAATGACCCTTGTCTTGGCACGATACTAAACACAAAAAGACAAAGCAATATGCACTGGTAACAGCAAGTGTACGTAACACGTCAACCCAATTAAGCGACCGAAATCGTAAATTCTGATAGAGTCTATGGTGGAGGCATATCCAAGATCTAACCAATATACCACAAGCAGAGTCAGATTAATGTTCCGAACGGAAGA"
example_dna_seq2 = example_dna_seq2.lower()

example_dna_seq3 = "AGAAGGAATGACCCTTGCTTGGCACGATACTAAACACAAAAAGACAAAGCAATATGCACTGGTAACAGCAAGTGTACGTAACACGTCAACCCAATTAAGCGACCGAAATCGTAAATTCTGATAGAGTCTATGGTGGAGGCATATCCAAGATCTAACCAATATACCACAAGCAGAGTCAGATTAATGTTCCGAACGGAAGA"
example_dna_seq3 = example_dna_seq3.lower()

# Import example organisms
org_list = organism_factory.import_organisms("example_organisms.json")
org_3_nodes, org_3_nodes2, org_5_nodes, org_9_nodes = org_list
print("\n")
org_3_nodes.print()
print("\n")
org_5_nodes.print()
print("\n")
org_9_nodes.print()

# Or generate a new example organism
new_example_organism = organism_factory.get_organism()
new_example_organism.print()



# Chose the organism to use as example in the demo

example_organism = org_3_nodes
example_organism2 = org_3_nodes2



# Get best placement
#results = example_organism.get_placement(example_dna_seq)
results = example_organism.get_placement(example_dna_seq, print_out = True)

results = example_organism2.get_placement(example_dna_seq2, print_out = True)

results = example_organism2.get_placement(example_dna_seq3, print_out = True)

example_organism2.get_seq_fitness(example_dna_seq3, print_out = True)

print(example_organism2.print_result(example_dna_seq3))

results["energy"]
results["recognizers_scores"]
results["connectors_scores"]






"""



# WALK THROUGH THE CODE OF get_placement --------------------------------------


# Initialize the two matrices

# Number of rows
m = example_organism.sum_pssm_lengths()
# Number of columns
n = len(example_dna_seq)

# Initialize matrix of scores
scores_matrix = np.full((m+1, n+1), -1 * np.inf)
scores_matrix[0,:] = 0

# Initialize matrix of pointers
pointers_matrix = np.full((2, m+1, n+1), None)

# Fill the matrices
for i in range(1, m + 1):
    
    # Diagonal scores over row i
    for j in range(1, n + 1):
        
        diag_score = example_organism.get_diag_score(pointers_matrix, i, j, example_dna_seq)
        scores_matrix[i,j] = scores_matrix[i-1, j-1] + diag_score
        # Annotate "where we came from" in the pointers_matrix
        pointers_matrix[0][i,j] = i - 1  # row idx of the origin
        pointers_matrix[1][i,j] = j - 1  # column idx of the origin
        
    
    # Horizontal scores over row i
    # (only in rows at the interface with the next PSSM)
    if example_organism.is_last(i) and i != m:
        
        # Scores and pointers from horizontal moves are temporarily stored in
        # the following arrays. They will be written altogether at the end of
        # the loop over j, to avoid reading them as starting scores for other
        # horizontal moves at a later cycles in the for loop over j
        tmp_gap_scores = scores_matrix[i,:].copy()  # vector of length n+1
        tmp_gap_pointers = pointers_matrix[:,i,:].copy()  # 2 x (n+1) matrix
        
        for j in range(1, n + 1):
            
            # Compute all the possible values from all the possible
            # horizontal moves (gaps) that land on [i,j]
            for start in range(j):
                gap_size = j - start
                pssm_idx = example_organism.row_to_pssm[i][0]
                gap_score = example_organism.get_gap_score(pssm_idx, gap_size, n)
                candidate_score = scores_matrix[i, start] + gap_score
                
                if candidate_score >= tmp_gap_scores[j]:
                    # Store the score
                    tmp_gap_scores[j] = candidate_score
                    # Annotate "where we came from" in the tmp_gap_pointers
                    tmp_gap_pointers[0,j] = i  # row idx of the origin
                    tmp_gap_pointers[1,j] = start  # column idx of the origin
        
        
        # Update the original matrices
        scores_matrix[i,:] = tmp_gap_scores
        pointers_matrix[:,i,:] = tmp_gap_pointers
    
# Get best binding energy
last_row = scores_matrix[-1,:]
best = max(last_row)

# BACKTRACKING

# Position of best (where backtracking starts from)
best_i = m  # it always comes from the last row by definition
best_j = int(np.where(last_row == best)[0])  # column of best

# Traverse back the matrix from the best element in the last row
# and store the alignment path
alignment_path = []
alignment_path = example_organism.traverse_matrix(pointers_matrix, best_i, best_j, alignment_path)
alignment_path.reverse()  # Top-down instead of bottom-up

# Get scores and positions of all the nodes of the organism
node_scores, node_positions, cols_of_0_gaps = example_organism.get_node_positions_and_energies(
    alignment_path, scores_matrix, pointers_matrix, example_dna_seq
)

# Print placement
example_organism.print_placement(node_positions, node_scores, cols_of_0_gaps, example_dna_seq)

# Split node-scores in recognizers-scores and connectors-scores
node_scores = node_scores[1:]  # Remove token node
recognizers_scores = []
connectors_scores = []
for i in range(len(node_scores)):
    if i % 2 == 0:
        recognizers_scores.append(node_scores[i])
    else:
        connectors_scores.append(node_scores[i])

# Return output dictionary
output_dict = {"energy": best,
               "recognizers_scores": recognizers_scores,
               "connectors_scores": connectors_scores}
    


"""
















