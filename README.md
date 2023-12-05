# r-admissability
A Boost-based poly-time algorithm to approximate the r-admissability of a (sparse) graph

# Theoretical Background 
R-admissability is a notion in the theory of (sparse) graphs that can be used for approximating a lower bound for the weak-coloring number, a property of sparse graphs that is very important for many real-world algorithms. 
The r-admissability of a graph adm<sub>r</sub>(G) is an ordering of the graph’s vertices such that the greatest of the r-admissable sets of the vertices is minimal under all orderings. The r-admissable set of a vertex v are all those vertices that are smaller than v w.r.t. the ordering and reachable from v via disjoint paths of length at most r and with inner vertices higher in the order than v.

# The Algorithm
The algorithm builds the order from right to left always putting a vertex with an r-admissable set of minimum size. The number of paths starting at a vertex v is determined by Depth-First-Search (DFS) in a greedy manner, i.e. always choosing a path where the first inner vertex (a vertex greater than v in the ordering) comes first in the neighbor-list of v and is not part of a path yet.

# COMPLEXITY
Let n be the number of vertices in the graph and m be the number of its edges.
The complexity of the algorithm is O(n<sup>3</sup>m), since it takes O(n) steps to build the order, where each of them consists of comparing the r-admissable sets of O(n) vertices. The calculation of each set takes O(nm) steps since DFS is performed at most m times.

# INPUT AND PREREQUISITES
The algorithm makes use of the Boost library as well as c++14 features.
The input is expected to be either of file format „csv“ or „txtg“ and an edge list where one edge spans a line and the endpoints are separated by a tab. 
