#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

typedef struct {
	int** mat_A;
	int* vec_k;
	int m;
	int num_nodes;
}Graph_A;

typedef struct {
	int* arr_g;
	int size_g;
	int* tmp_vec;
}Group;

typedef struct {
	Group group1, group2;
}Devision;

double Epsilon = 0.00001;

Group* alg3(Graph_A* graph);
Group trivial_group(Graph_A* graph);
void outout_groups(Group* O, int num_groups, char* out_file);

Group* alg3(Graph_A* graph)
{
	Group triv_g = trivial_group(graph);
	Group* P = calloc(graph->num_nodes, sizeof(int));
	Group* O;
	Group g;
	Devision temp;
	int size_P = 1, last_ind_P = 1;
	int size_O = 0, last_ind_O = 0;
	P[0] = triv_g;
	O = calloc(graph->num_nodes, sizeof(Group));
	while (size_P > 0)
	{
		g = P[last_ind_P];
		--last_ind_P;
		--size_P;
		
		/*here comes algorithm 2 into part*/
		//temp = alg_2(g);

		if (temp.group1.size_g == 0)
		{
			/*add group2 into O*/
			O[last_ind_O] = temp.group2;
			last_ind_O++;
			size_O++;
		}
		if (temp.group2.size_g == 0)
		{
			/* add group1 into O*/
			O[last_ind_O] = temp.group1;
			last_ind_O++;
			size_O++;
		}
		else
		{
			if (temp.group1.size_g == 1)
			{
				/*add group 1 into O*/
				O[last_ind_O] = temp.group1;
				last_ind_O++;
				size_O++;
			}
			else
			{
				/*add group 1 into P*/
				P[last_ind_P] = temp.group1;
				last_ind_P++;
				size_P++;
			}
			if (temp.group2.size_g == 1)
			{
				/*add group 2 into O*/
				O[last_ind_O] = temp.group2;
				last_ind_O++;
				size_O++;
			}
			else
			{
				/*add group 2 into P*/
				P[last_ind_P] = temp.group2;
				last_ind_P++;
				size_P++;
			}
		}
	}
	return O;
}
Group trivial_group(Graph_A* graph)
{
	Group full;
	int* arr;
	int i;
	full.size_g = graph->num_nodes;
	arr = calloc(full.size_g, sizeof(int));
	full.tmp_vec = calloc(full.size_g, sizeof(int));
	for (i = 0; i < full.size_g; i++)
		arr[i] = i;
	full.arr_g = arr;
	return full;
}
void outout_groups(Group* O, int num_groups, char* out_file)
{
	FILE* out_div;
	int n;
	Group* pointer;
	Group g;
	out_div = fopen(out_file, "w");
	n = fwrite(&num_groups, sizeof(int), 1, out_div);
	if (n != 1)
	{
		printf("Error writing to output file");
	}
	for (pointer = &O[0]; pointer < &O[num_groups]; pointer++)
	{
		g = *pointer;
		n = fwrite(&g.size_g, sizeof(int), 1, out_div);
		if (n != 1)
		{
			printf("Error writing to output file");
		}
		n = fwrite(&g.arr_g, sizeof(int), g.size_g, out_div);
		if (n != g.size_g)
		{
			printf("Error writing to output file");
		}
	}
	fclose(out_div);

}
int main(int argc, char* argv[]) {
	/*Graph_A graph;
	int** P;
	int* p[6];
	graph.num_nodes = 6;
	graph.vec_k = NULL;
	graph.m = 0;
	graph.mat_A = NULL;
	for (int i = 0; i < graph.num_nodes; i++)
	{
		p[i] = i + 1;
	}
	for (int i = 0; i < graph.num_nodes; i++)
	{
		printf("%d",p[i]);
	}*/
}