#include "graph.h"
#include "group.h"


void copy_vector_int(int* org_vec, int* new_vec, int dim) {
	int* org_p, * new_p;
	new_p = new_vec;
	for (org_p = org_vec; org_p < &org_vec[dim]; org_p++) {
		*new_p = *org_p;
		new_p++;
	}
}

void copy_group(Group* org_group, Group* new_group) {
	new_group->size_g = org_group->size_g;
	new_group->arr_g = (int*)malloc(new_group->size_g * sizeof(int));
	copy_vector_int(org_group->arr_g, new_group->arr_g, new_group->size_g);
}

void trivial_group(Group* triv_group, Graph_A* graph)
{
	int i;

	triv_group->size_g = graph->num_nodes;
	triv_group->arr_g = (int*)malloc(graph->num_nodes * sizeof(int));

	for (i = 0; i < triv_group->size_g; i++)
	{
		triv_group->arr_g[i] = i;
	}
}

int output_groups(Group* O, int num_groups, char* out_file)
{
	FILE* out_div;
	int n;
	Group* pointer;
	Group g;
	out_div = fopen(out_file, "w");

	n = fwrite(&num_groups, sizeof(int), 1, out_div);
	if (n != 1)
	{
		printf("Error writing to output file\n");
		return -1;
	}

	for (pointer = &O[0]; pointer < &O[num_groups]; pointer++)
	{
		g = *pointer;
		n = fwrite(&(g.size_g), sizeof(int), 1, out_div);
		if (n != 1)
		{
			printf("Error writing to output file\n");
			return -1;
		}
		n = fwrite(g.arr_g, sizeof(int), g.size_g, out_div);
		if (n != g.size_g)
		{
			printf("Error writing to output file\n");
			return -1;
		}
	}

	fclose(out_div);

	return 0;
}

