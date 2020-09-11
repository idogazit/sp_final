#include "graph.h"

#ifndef _GROUP_H
#define _GROUP_H

typedef struct {
	int* arr_g;
	int size_g;
}Group;

void copy_group(Group* org_group, Group* new_group);
void trivial_group(Group* triv_group, Graph_A* graph);
void output_groups(Group* O, int num_groups, char* out_file);
void copy_vector_int(int* org_vec, int* new_vec, int dim);

#endif