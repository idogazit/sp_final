#include <stdio.h>
#include <stdlib.h>
#include "graph.h"


void kill_graph(Graph_A* graph) {
	int** p_node;
	free(graph->tmp_vec);
	for (p_node = graph->mat_A; p_node < graph->mat_A + graph->num_nodes; p_node++) {
		free(*p_node);
	}
	free(graph->mat_A);

	free(graph->vec_k);
}

int compute_M(int* vec_K, int dim) {
	int m = 0, * p;
	for (p = vec_K; p < &vec_K[dim]; p++) {
		m += *p;
	}
	return m;
}

void build_graph_A(Graph_A* graph, FILE* file) {
	int dim[1], ** pi, * deg;
	int n;

	n = fread(dim, sizeof(int), 1, file);
	if (n != 1)
	{
		printf("Error in reading from input file\n");
		exit(-1);
	}
	graph->num_nodes = *dim;
	graph->mat_A = (int**)malloc(*dim * sizeof(int*));
	graph->vec_k = (int*)malloc(*dim * sizeof(int));
	graph->tmp_vec = (double*)malloc(*dim * sizeof(double));
	for (pi = graph->mat_A, deg = graph->vec_k; pi < &(graph->mat_A[*dim]); pi++, deg++) {
		n = fread(deg, sizeof(int), 1, file);
		if (n != 1)
		{
			printf("Error in reading from input file\n");
			exit(-1);
		}
		*pi = (int*)malloc(*deg * sizeof(int));
		n = fread(*pi, sizeof(int), *deg, file);
		if (n != *deg)
		{
			printf("Error in reading from input file\n");
			exit(-1);
		}
	}
	graph->m = compute_M(graph->vec_k, graph->num_nodes);
	graph->max_num_iters = 10000 * graph->num_nodes;
}