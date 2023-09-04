# ASI-GRAPH
Algorithm (adapted in Python) to find approximal sub-graph isomorphism between two graphs

The degree of freedom can be split into three categories:
Mismatch of labels
An unequal amount of edges (more or less than the reference graph)
Inequal amount of node (more or less than the reference graph) - NOT IMPLEMENTED

This algorithm uses brute force to calculate all match possibilities respecting the penalty imposed by the user.

The implementation works for multi-directed graphs, directed graphs, and undirected graphs.

---Utilisation---
compare2(g,Motif,nb_miss,arr_plus,arr_moins,range_diff)

g=Main graph
Motif=graph to mach with subgraphs of g
nb_miss= maximal number of mismatches allowed in results of approximal isomorphism
arr_plus= maximal number of additional edges allowed in results of approximal isomorphism
arr_moins= maximal number of edges allowed to be missed in results of approximal isomorphism
range_diff=maximal number of mismatches of range (ie. long-range or short-range)

return:
A list of edge match (e1,e2) with e1 belonging to g and e2 to Motif
