#include <stdio.h>
#include <stdlib.h>

typedef struct {
	int** mat_A;
	int* vec_k;
	int m;
	int num_nodes;
	double norm;
	double* tmp_vec;

}Graph_A;

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
	/*	int f*/

	/*f =*/ fread(dim, sizeof(int), 1, file);
	/*check f==1 */
	graph->num_nodes = *dim;
	graph->mat_A = (int**)calloc(*dim, sizeof(int*));
	graph->vec_k = (int*)calloc(*dim, sizeof(int));
	graph->tmp_vec = (double*)calloc(*dim, sizeof(double));
	for (pi = graph->mat_A, deg = graph->vec_k; pi < &(graph->mat_A[*dim]); pi++, deg++) {
		/*	f =*/ fread(deg, sizeof(int), 1, file);
		/*check f==1 */
		*pi = (int*)calloc(*deg, sizeof(int));
		/*	f =*/ fread(*pi, sizeof(int), *deg, file);
		/*check f==*dim */
	}
	graph->m = compute_M(graph->vec_k, graph->num_nodes);

}