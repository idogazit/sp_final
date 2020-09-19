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
	int max_num_iters;

}Graph_A;


/*
* This function gets an address for a graph, reads the relevant info from the file
* and keeps it in the struct vars
*
* @param: graph - the target graph address, file - the input file address
*/
void build_graph_A(Graph_A* graph, FILE* file);


/*
* This function gets an address for a graph and frees all allocated memory of it
*
* @param : graph - the target graph address
*/
void kill_graph(Graph_A* graph);


/*
* This function gets a vector of nodes degrees and it's dimension
* and it returns the sum of the vector (value M)
*
* @param : vec_K - the target vector address, dim - the vector's size
* @return: sum of vec_K
*/
int compute_M(int* vec_K, int dim);

#endif