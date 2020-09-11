#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef _GRAPH_H
#define _GRAPH_H


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