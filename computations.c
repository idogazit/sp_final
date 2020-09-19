#include "partition.h"
#include "computations.h"
#define Epsilon 0.00001

void make_vec_of_1s(double* vec, int dim) {
	double* p;
	for (p = vec; p < &vec[dim]; p++) {
		*p = 1.0;
	}
}

void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node, int i) {
	int* p_grp, m = graph->m, dim = group->size_g;
	double* p_Bg, f_g_node = 0.0, tmp = (double)graph->vec_k[node] / (double)m;


	build_row_g(row_Bg_hat, graph, group, node);

	for (p_grp = group->arr_g, p_Bg = row_Bg_hat; p_grp < group->arr_g + dim; p_Bg++, p_grp++) {
		*p_Bg -= tmp * (double)graph->vec_k[*p_grp];
		f_g_node += *p_Bg;
	}
	row_Bg_hat[i] -= f_g_node;

}

void build_row_g(double* row_g, Graph_A* graph, Group* group, int node) {
	int size = group->size_g;
	int deg = graph->vec_k[node];
	int* p_node_ngbrs = graph->mat_A[node], * p_grp = group->arr_g;
	double* p_rowg = row_g;


	for (; p_rowg < row_g + size; p_rowg++) {
		*p_rowg = 0.0;
	}
	p_rowg = row_g;
	while ((p_grp < group->arr_g + size) && (p_node_ngbrs < graph->mat_A[node] + deg)) {
		if (*p_grp < *p_node_ngbrs) {
			p_rowg++;
			p_grp++;
		}
		else if (*p_grp > * p_node_ngbrs) {
			p_node_ngbrs++;
		}
		else {
			*p_rowg = 1.0;
			p_rowg++;
			p_grp++;
			p_node_ngbrs++;
		}
	}
}



double compute_graph_norm(Graph_A* graph, Group* trivial_group) {
	int size = graph->num_nodes;
	int* p_grp, i = 0;
	double max = 0.0, * p_row_Bgh, * p_col;
	double* col_sum = (double*)calloc(size, sizeof(double));


	for (p_grp = trivial_group->arr_g; p_grp < trivial_group->arr_g + size; p_grp++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, trivial_group, *p_grp, i++);
		for (p_row_Bgh = graph->tmp_vec, p_col = col_sum; p_col < col_sum + size; p_row_Bgh++, p_col++) {
			*p_col += fabs(*p_row_Bgh);
		}
	}
	for (p_col = col_sum; p_col < &col_sum[size]; p_col++) {
		if (max < *p_col) { max = *p_col; }
	}
	free(col_sum);
	return max;
}

double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group) {
	int size = group->size_g;
	int* p_grp, i = 0;
	double val = 0.0, * p_vec1, * p_vec2, * p_Bg;
	for (p_vec1 = vec, p_grp = group->arr_g; p_grp < group->arr_g + size; p_vec1++, p_grp++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, group, *p_grp, i++);
		for (p_vec2 = vec, p_Bg = graph->tmp_vec; p_vec2 < vec + size; p_Bg++, p_vec2++) {
			val += (double)*p_vec1 * *p_vec2 * *p_Bg;;
		}
	}
	return val;
}

double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A* graph, Group* group) {
	double val, col_m_row;
	col_m_row = row_multiply_col(leading_eigenvec, leading_eigenvec, group->size_g);
	if (col_m_row == 0)
	{
		printf("Error in computing leading eigen value\n");
		exit(-1);
	}
	val = compute_vec_BgH_vec(leading_eigenvec, graph, group);
	val /= col_m_row;

	return val;
}


void compute_vec_s_on_eigen_vec(double* eigenvec, int dim) {
	double* q;
	for (q = eigenvec; q < &eigenvec[dim]; q++) {
		if (*q > Epsilon) { *q = 1.0; }
		else { *q = -1.0; }
	}
}
void devide_according_to_s(Devision* devision, Group* group, double* vec_s) {
	int* p_dev1, * p_dev2, * p_g, dim = group->size_g;
	double* p_s;
	devision->group1.size_g = 0;
	devision->group2.size_g = 0;


	for (p_s = vec_s; p_s < &vec_s[dim]; p_s++) {
		if (*p_s > 0.0) {
			devision->group1.size_g++;
		}
		else {
			devision->group2.size_g++;
		}
	}

	devision->group1.arr_g = (int*)malloc(devision->group1.size_g * sizeof(int));
	devision->group2.arr_g = (int*)malloc(devision->group2.size_g * sizeof(int));

	p_dev1 = devision->group1.arr_g;
	p_dev2 = devision->group2.arr_g;

	for (p_g = group->arr_g, p_s = vec_s; p_g < &group->arr_g[dim]; p_s++, p_g++) {
		if (*p_s > Epsilon) {
			*p_dev1 = *p_g;
			p_dev1++;
		}
		else {
			*p_dev2 = *p_g;
			p_dev2++;
		}
	}
}


void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group) {
	double* vec0, * curr_vec;
	int size = group->size_g;
	int i = 0;
	vec0 = (double*)malloc(size * sizeof(double));
	curr_vec = (double*)malloc(size * sizeof(double));


	generate_rand_vec0(vec0, size);
	copy_vector(vec0, curr_vec, size);
	generate_next_vec(eigenvec, curr_vec, graph, group);
	while ((vectors_difference_is_small(curr_vec, eigenvec, size) == 0)) {
		if (i++ == graph->max_num_iters)
		{
			printf("Infinite loop accured\n");
			exit(-1);
		}
		copy_vector(eigenvec, curr_vec, size);
		generate_next_vec(eigenvec, curr_vec, graph, group);
	}
	free(vec0);
	free(curr_vec);
}
void generate_rand_vec0(double* vec0, int dim) {

	double* p;

	for (p = vec0; p < &vec0[dim]; p++) {
		*p = (double)(rand() % 5000) + 1;
	}
}
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group) {
	int size = group->size_g;
	int i;
	double* p_next_vec;
	double norm = 0.0;
	for (i = 0; i < size; i++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, group, group->arr_g[i], i);

		graph->tmp_vec[i] += graph->norm;
		next_vec[i] = row_multiply_col(graph->tmp_vec, curr_vec, size);
		norm += next_vec[i] * next_vec[i];
	}
	norm = sqrt(norm);
	if (norm == 0)
	{
		printf("Error in generating next vector\n");
		exit(-1);
	}
	for (p_next_vec = next_vec; p_next_vec < &next_vec[size]; p_next_vec++) {
		*p_next_vec /= norm;
	}

}
double row_multiply_col(double* row, double* col, int dim) {

	double* p_row, * p_col;
	double sum = 0;
	p_col = col;

	for (p_row = row; p_row < &row[dim]; p_row++) {
		sum += *p_row * *p_col;
		p_col++;
	}

	return sum;
}
void copy_vector(double* org_vec, double* new_vec, int dim) {
	double* org_p, * new_p;
	new_p = new_vec;
	for (org_p = org_vec; org_p < &org_vec[dim]; org_p++) {
		*new_p = *org_p;
		new_p++;
	}
}

int vectors_difference_is_small(double* vector1, double* vector2, int dim) {

	double* p1, * p2;
	p2 = vector2;

	for (p1 = vector1; p1 < &vector1[dim]; p1++) {
		if (fabs(*p1 - *p2) >= Epsilon) return 0;
		p2++;
	}

	return 1;
}
void reset_unmoved(int* unmoved, int len)
{
	int* p;
	for (p = unmoved; p < &unmoved[len]; p++)
	{
		*p = 0;
	}
}

void compute_vec_Bg_ii(double* vec_Bg_ii, Graph_A* graph, Group* group) {
	double* p_vec;
	int* p_grp;

	for (p_vec = vec_Bg_ii, p_grp = group->arr_g; p_vec < vec_Bg_ii + group->size_g; p_vec++, p_grp++) {
		*p_vec = (double)(graph->vec_k[*p_grp] * graph->vec_k[*p_grp]) / (double)graph->m;
	}
}

void compute_row_Bg(double* row_Bg, Graph_A* graph, Group* group, int node) {
	int* p_grp, m = graph->m, dim = group->size_g;
	double* p_Bg, tmp = (double)graph->vec_k[node] / (double)m;

	build_row_g(row_Bg, graph, group, node);

	for (p_grp = group->arr_g, p_Bg = row_Bg; p_grp < group->arr_g + dim; p_Bg++, p_grp++) {
		*p_Bg -= tmp * (double)graph->vec_k[*p_grp];
	}
}
