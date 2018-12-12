# test script to compute Katz similarity
import os

file_name_Tx1 = './sampleTaxonomy/TxX_edge_list.txt'
file_name_Tx2 = './sampleTaxonomy/TxC_edge_list.txt'
output_file_name =  './Katz_similarity_output.txt'
output_file_name_importance = './Katz_similarity_output_VI.txt'

alpha = 0.8;
num_vertices = 9;
# this is for the example taxonomy pair, you'll need to know num_vertices for the two taxonomies you're comparing, this will be the size of the union of the vertex set for both graphs
gamma = 1.0/(num_vertices*num_vertices);

# for Katz Similarity for Graphs
method_id = 1


if(method_id==1):
    cmd = './KatzDAGSimilarity ' + file_name_Tx1 + ' ' + file_name_Tx2 + ' ' + str(alpha) + ' ' + str(gamma) + ' ' + str(num_vertices) + ' ' + str(method_id) + ' ' + output_file_name + ' ' + output_file_name_importance

# for Grouped Katz Similarity for Graphs
if(method_id == 2):
    num_groups = 3
    cmd = './KatzDAGSimilarity ' + file_name_Tx1 + ' ' + file_name_Tx2 + ' ' + str(alpha) + ' ' + str(gamma) + ' ' + str(num_vertices) + ' ' + str(method_id) + ' ' + output_file_name + ' ' + str(num_groups) + ' ' + output_file_name_importance

# for Katz Index Similarity for Graphs
if(method_id == 3):
    cmd = './KatzDAGSimilarity ' + file_name_Tx1 + ' ' + file_name_Tx2 + ' ' + str(alpha) + ' ' + str(gamma) + ' ' + str(num_vertices) + ' ' + str(method_id) + ' ' + output_file_name + ' ' + output_file_name_importance

# call the Katz Similarity executable
os.system(cmd)

# ------------------- Read the output -------------------------
import numpy as np
c = np.loadtxt(output_file_name)
graph_similarity = c[0]
vector_difference_norm = c[1]
print 'The similarity between the graphs by method ' + str(method_id) + ' is ' + str(graph_similarity) + ' and '
print 'the corresponding vector difference norm is ' + str(vector_difference_norm)

VI = np.loadtxt(output_file_name_importance)

# Compute the top 3 vertices that changed the most and their corresponding importances
ix = np.argsort(VI)

print 'The top 3 vertices in terms of their importance to the graph similairity are: '
print 'Vertex id ' + str(ix[-1] + 1) + ', importance: ' + str(VI[ix[-1]])
print 'Vertex id ' + str(ix[-2] + 1) + ', importance: ' + str(VI[ix[-2]])
print 'Vertex id ' + str(ix[-3] + 1) + ', importance: ' + str(VI[ix[-3]])