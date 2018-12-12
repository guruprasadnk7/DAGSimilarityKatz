-------- Format of the input files -----------------
The Katz Similarity code expects the two Directed Acyclic Graphs (DAGs) being compared to be input as edge list files. The vertices are expected to be indexed starting from 1 to num_vertices, where num_vertices is the number of vertices in the union of the vertex sets of the 2 DAGs. Ensure that the vertex numbering is consistent between the two graphs. i.e vertex id v corresponds to the same vertex in DAG 1 as in DAG 2.
For every edge (u->v) in a DAG, the edge list file will have a line 'u v'. Every line of the edge list fie corresponds to one edge in the graph.

The sampleTaxonomy folder has the edge lists for 3 taxonomies from the caricature in figure 1 in the paper. The vertices are indexed as follows -

(Vertex ID, Vertex name)
(1, Animals)
(2, Mammals)
(3, Felines)
(4, Tigers)
(5, Domestic Cats)
(6, Bovines)
(7, Reptiles)
(8, Lizards)
(9, Snakes)

------------------- Running the code ----------------
Compile the C code like so,
gcc ./KatzDAGSimilarity.c -lm -o ./KatzDAGSimilarity

Once the executable is made, you'll need to pass it the appropriate arguments. Look at the test_script.py python script that lists the arguments for each of the 3 cases (i.e Katz Similarity for Graphs, Grouped Katz Similarity for Graphs,Katz Index Similarity for Graphs).


------------------- Citation ---------------------
If you make use of this code, please cite the following paper -

Nayak, G., Dutta, S., Ajwani, D., Nicholson, P., & Sala, A. "Automated assessment of knowledge hierarchy evolution: comparing directed acyclic graphs." Information Retrieval Journal