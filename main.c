/*
 * alg_2.c
 *
 *  Created on: 24 ���� 2020
 *      Author: USER
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "algorithms.h"

int main(int argc, char* argv[]) {
	Graph_A graph;
	Partition final_partition;
	FILE* input_file;
	clock_t start, end;
	start = clock();
	if (argc != 3)
	{
		printf("Error in inputs\n");
		return -1;
	}
	
	input_file = fopen(argv[1], "r");
	build_graph_A(&graph,input_file);
	fclose(input_file);
	final_partition.groups =(Group*)calloc(graph.num_nodes, sizeof(Group));
	final_partition.num_of_groups = 0;
	if(graph.m == 0){
		trivial_partition(&final_partition,graph.num_nodes);
	}
	else{
		alg3(&graph, &final_partition);
	}
	output_groups(final_partition.groups, final_partition.num_of_groups, argv[2]);

	kill_partition(&final_partition);
	kill_graph(&graph);
	end = clock();
	printf("Execution time is %f seconds\n",
				((double)(end-start) / CLOCKS_PER_SEC));
	/*print_output(argv[2]);*/
	return 0;
}