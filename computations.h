#include "partition.h"

#ifndef _COMPUTATIONS_H
#define _COMPUTATIONS_H



int vectors_difference_is_small(double* vector1, double* vector2, int dim);
void devide_according_to_s(Devision* devision, Group* group, double* vec_s);
double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group);
void build_row_g(double* row_g, Graph_A* graph, Group* group, int node);
void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node, int i);
void compute_row_Bg(double* row_Bg, Graph_A* graph, Group* group, int node);
double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A*, Group* group);
void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group);
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group);
double compute_graph_norm(Graph_A* graph, Group* trivial_group);
void compute_vec_rowBgH_s(double* vec_rowBgH_s, Graph_A* graph, Group* group, double* vec_s);
void compute_vec_Bg_ii(double* vec_Bg_ii, Graph_A* graph, Group* group);

/*General functions: */
void compute_vec_s_on_eigen_vec(double* eigenvec, int dim);
void make_vec_of_1s(double* vec, int dim);
void generate_rand_vec0(double* vec0, int dim);
double row_multiply_col(double* row, double* col, int dim);
void copy_vector(double* org_vec, double* new_vec, int dim);
void reset_unmoved(int* unmoved, int len);


#endif