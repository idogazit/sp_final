/*
* COMPUTATIONS Summary:
*
* A module that holds all the low-level and simple computations needed in the main algorithms
*
*
* vectors_difference_is_small		- This function gets two vectors and checks if the metrics between them is smaller than epsilon
* devide_according_to_s				- This function gets a group and a vector of +1 OR -1 and divides the group into two groups according to the vector
* compute_vec_BgH_vec				- This function computes the product of a vector multiples the modularity function from both sides
* build_row_g						- This function calculates the a current row of the neighbors matrix according to the graph and the current group
* compute_row_Bg_hat				- This function calculates the a current row of the modularity matrix according to the graph and the current group
* compute_leading_eigenvalue		- This function computes the leading eigenvalue according to the leading eigenvector of the modularity matrix
* compute_leading_eigenvec			- This function computes the leading eigenvector of the modularity matrix
* generate_next_vec					- This function generates the next vector in the power iteration
* compute_graph_norm				- This function computes the norm of the modularity matrix of the entire graph
* compute_vec_Bg_ii					- This function computes the vector of the i,i index of the modularity matrix accoring to the current group
* compute_vec_s_on_eigen_vec		- This function computes vector of +1 OR -1 according to the index sign in the eigenvector
* make_vec_of_1s					- This function makes a vector of 1s
* generate_rand_vec0				- This function generates a random vector
* row_multiply_col					- This function calculates the product of a row vector multiply a column vector
* copy_vector						- This function copies the info of a double vector to a new vector
* reset_unmoved						- Resets the set of unmoved nodes to be zeros
*/

#include "partition.h"

#ifndef _COMPUTATIONS_H
#define _COMPUTATIONS_H


/*
 * This function gets two vectors and checks if the metrics between them is smaller than epsilon
 *
 * @param : vector1 - first vector, vector2 - second vector, dim - size of the vectors
 * @return : returns 1 if the metrics between them is smaller than Epsilon, else returns 0
 */
int vectors_difference_is_small(double* vector1, double* vector2, int dim);


/*
 * This function gets a group and a vector of +1 OR -1 and divides the group into two groups according to the vector
 *
 * @param : devision - the division of the group address,
 * 			group - the current relevant group of nodes address
 * 			vec_s - the address of the vector that the group gets divided according to it
 */
void devide_according_to_s(Devision* devision, Group* group, double* vec_s);


/*
 * This function computes the product of a vector multiples the modularity function from both sides
 *
 * @param : vec - the address of the vector that gets multiply, graph - the entire graph info address
 * 			group - the current relevant group of nodes address
 * @return : returns the product of the multiplication
 */
double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group);


/*
 * This function calculates the a current row of the neighbors matrix according to the graph and the current group
 *
 * @param : row_g - the address in which the row will be written to,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 * 			node - the node in the group that according to it the row is calculated
 */
void build_row_g(double* row_g, Graph_A* graph, Group* group, int node);


/*
This function calculates the current row of the modularity (after normalization) matrix according to the graph and the current group
 *
 * @param : row_Bg_hat - the address in which the row will be written to,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 * 			node - the node in the group that according to it the row is calculated
 * 			i - the number of the row in the modularity matrix
 */
void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node, int i);


/*
 * This function computes the leading eigenvalue according to the leading eigenvector of the modularity matrix
 *
 * @param : leading_eigenvec - the leading eigenvector address,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 * 	@return : the leading eigenvalue
 */
double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A*, Group* group);


/*
 * This function computes the leading eigenvector of the modularity matrix
 *
 * @param : eigenvec - the address in which the vector will be written to,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 */
void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group);


/*
 * This function generates the next vector in the power iteration
 *
 * @param : next_vec - the address in which the vector will be written to,
 * 			curr_vec - the address of the current vector in the power iteration calculation,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 */
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group);


/*
 * This function computes the norm of the modularity matrix of the entire graph
 *
 * @param : graph - the entire graph info address,
 * 			trivial_group - the address of a group that contains all the nodes in the graph
 * @return : the norm of the modularity matrix of the entire graph
 */
double compute_graph_norm(Graph_A* graph, Group* trivial_group);


/*
 * This function computes the vector of the i,i index of the modularity matrix accoring to the current group
 *
 * @param : vec_Bg_ii - the address in which the vector will be written to,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 */
void compute_vec_Bg_ii(double* vec_Bg_ii, Graph_A* graph, Group* group);


/*
 * This function computes vector of +1 OR -1 according to the index sign in the eigenvector
 *
 * @param : eigenvec - the eigenvector address whih according to is vector S is calculated AND
 *  				   the address that the vector will be written to,
 *  		dim - the size of the vector
 */
void compute_vec_s_on_eigen_vec(double* eigenvec, int dim);


/*
 * This function makes a vector of 1s
 *
 * @param : vec - the address in which the vector will be written to,
 *   		dim - the size of the vector
 */
void make_vec_of_1s(double* vec, int dim);


/*
 * This function generates a random vector
 *
 * @param : vec0 - the address in which the vector will be written to,
 *   		dim - the size of the vector
 */
void generate_rand_vec0(double* vec0, int dim);


/*
 * This function calculates the product of a row vector multiply a column vector
 *
 * @param : row - the row vector address,
 * 			col - the column vector address,
 * 			dim - the size of the vectors
 * @return : the product of the two vectors
 */
double row_multiply_col(double* row, double* col, int dim);


/*
 * This function copies the info of a double vector to a new vector
 *
 * @param : org_vec - the original vector address,
 * 			new_vec - the new vector address,
 * 			dim - the size of the vectors
 */
void copy_vector(double* org_vec, double* new_vec, int dim);


/*
 * This function resets an array of integers to array of 0s
 *
 * @param : unmoved - the address of the array,
 * 			len - the size of the array
 */
void reset_unmoved(int* unmoved, int len);


/*
This function calculates the current row of the modularity matrix according to the graph and the current group
 *
 * @param : row_Bg - the address in which the row will be written to,
 * 			graph - the entire graph info address,
 * 			group - the current relevant group of nodes address
 * 			node - the node in the group that according to it the row is calculated
 */
void compute_row_Bg(double* row_Bg, Graph_A* graph, Group* group, int node);


#endif