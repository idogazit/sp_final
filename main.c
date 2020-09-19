
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "algorithms.h"

int main(int argc, char* argv[]) {
	Graph_A graph;
	Partition final_partition;
	FILE* input_file;
	int n;

	
	if (argc != 3)
	{
		printf("Error in inputs\n");
		return -1;
	}
	
	input_file = fopen(argv[1], "r");
	if (input_file == NULL)
	{
		printf("Error in opening file\n");
		return -1;
	}
	
	build_graph_A(&graph,input_file);
	fclose(input_file);
	
	final_partition.groups = (Group*)malloc(graph.num_nodes * sizeof(Group));
	final_partition.num_of_groups = 0;
	
	if(graph.m == 0)
	{
		trivial_partition(&final_partition,graph.num_nodes);
	}
	else
	{
		alg3(&graph, &final_partition);
	}
	
	n = output_groups(final_partition.groups, final_partition.num_of_groups, argv[2]);
	if (n != 0)
	{
		printf("Error in writing to destination file\n");
		return -1;
	}

	kill_partition(&final_partition);
	kill_graph(&graph);
	
	return 0;
}