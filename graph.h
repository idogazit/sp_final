/*
* GRAPH Summary:
*
* A module that implements a graph class
*
*
* build_graph_A			- Creates a new graph matrix A
* kill_graph			- Terminates all memory allocations related to the graph
* compute_M				- Computes the sum of the degrees of nodes in the graph
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef _GRAPH_H
#define _GRAPH_H

/*
* Type used for containing all relevant information regarding to a graph
*/
typedef struct {
	int** mat_A;
	int* vec_k;
	int m;
	int num_nodes;
	double norm;
	double* tmp_vec;

}Graph_A;


void build_graph_A(Graph_A* graph, FILE* file);
void kill_graph(Graph_A* graph);
int compute_M(int* vec_K, int size_K);

#endif