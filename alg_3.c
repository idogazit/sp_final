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

int alg3(Graph_A* graph, Group* O);
Group trivial_group(Graph_A* graph);
void output_groups(Group* O, int num_groups, char* out_file);

int alg3(Graph_A* graph, Group* O)
{
	Group triv_g = trivial_group(graph);
	Group* P = calloc(graph->num_nodes, sizeof(Group));
	Group g;
	Devision temp;
	int size_P = 1;
	int size_O = 0;
	P[0] = triv_g;

	while (size_P > 0)
	{
		g = P[size_P];
		--size_P;
		
		/*here comes algorithm 2 into part*/
		//temp = alg_2(g);

		if (temp.group1.size_g == 0 || temp.group2.size_g == 0)
		{
			/*add group g into O*/
			O[size_O] = g;
			size_O++;
		}
		else
		{
			if (temp.group1.size_g == 1)
			{
				/*add group 1 into O*/
				O[size_O] = temp.group1;
				size_O++;
			}
			else
			{
				/*add group 1 into P*/
				P[size_P] = temp.group1;
				size_P++;
			}
			if (temp.group2.size_g == 1)
			{
				/*add group 2 into O*/
				O[size_O] = temp.group2;
				size_O++;
			}
			else
			{
				/*add group 2 into P*/
				P[size_P] = temp.group2;
				size_P++;
			}
		}
	}
	
	free(triv_g.arr_g);
	free(triv_g.tmp_vec);
	free(g.arr_g);
	free(g.tmp_vec);
	free(P);
	
	return size_O;
}
Group trivial_group(Graph_A* graph)
{
	Group full;
	int i;
	
	full.size_g = graph->num_nodes;
	full.arr_g = calloc(full.size_g, sizeof(int));
	full.tmp_vec = calloc(full.size_g, sizeof(int));

	for (i = 0; i < full.size_g; i++)
	{
		full.arr_g[i] = i;
	}

	return full;
}
void output_groups(Group* O, int num_groups, char* out_file)
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
	}
	
	for (pointer = &O[0]; pointer < &O[num_groups]; pointer++)
	{
		g = *pointer;
		n = fwrite(&(g.size_g), sizeof(int), 1, out_div);
		if (n != 1)
		{
			printf("Error writing to output file\n");
		}
		n = fwrite(&(g.arr_g), sizeof(int), g.size_g, out_div);
		if (n != g.size_g)
		{
			printf("Error writing to output file\n");
		}
	}
	
	fclose(out_div);
}
int main(int argc, char* argv[]) {
	Graph_A* graph;
	Group* final_div;
	int num_groups;
	FILE* input_file;

	input_file = fopen(argv[1], "r"); /*opens input file and creates graph A*/
	graph = build_graph_A(input_file);
	fclose(input_file);

	num_groups = alg3(graph, &final_div); /*calls for algorithm 3*/

	output_groups(final_div, num_groups, argv[2]); /*writes final groups division into output file*/
	
	return 0;
}