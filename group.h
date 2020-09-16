/*
* GROUP Summary:
*
* A module that implements a group class
*
*
* copy_group			- Copies of a group
* trivial_group			- Creates a the trivial nodes group: {0,1,...,n-1}
* output_groups			- Writes a set of groups into output file
* copy_vector_int		- Copies a vector
*
*/
#include "graph.h"

#ifndef _GROUP_H
#define _GROUP_H

/*
* Type used for containing all relevant information regarding to a group
*/
typedef struct {
	int* arr_g;
	int size_g;
}Group;

void copy_group(Group* org_group, Group* new_group);
void trivial_group(Group* triv_group, Graph_A* graph);
int output_groups(Group* O, int num_groups, char* out_file);
void copy_vector_int(int* org_vec, int* new_vec, int dim);

#endif