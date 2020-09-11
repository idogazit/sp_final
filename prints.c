void print_graph(Graph_A* graph) /*debugging function*/
{
	int i = 0, j = 0;
	printf("=========Graph========\n");
	printf("M = %d\n", graph->m);
	printf("Number of nodes = %d\n", graph->num_nodes);
	printf("Norm = %f\n", graph->norm);
	printf("Vec K = [");
	for (i = 0; i < graph->num_nodes; i++)
	{
		printf("%d,", graph->vec_k[i]);
	}
	printf("]\n");

	printf("Mat: \n");
	for (i = 0; i < graph->num_nodes; i++)
	{
		printf("Node #%d neighbors: ", i);
		for (j = 0; j < graph->vec_k[i]; j++)
		{
			printf("%d,", graph->mat_A[i][j]);
		}
		printf("\n");
	}
}
void print_devision(Devision d) /*debugging function*/
{
	printf("$$$$$$$$$$$ Devision: $$$$$$$$$$$\n");
	printf("First Group:\n");
	print_group(d.group1);
	printf("Second Group:\n");
	print_group(d.group2);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

}
void print_group(Group g) /*debugging function*/
{
	int i = 0;
	printf("@@@@@@@@@@ Group: @@@@@@@@@@\n");

	printf("Group size = %d\n", g.size_g);

	printf("Group nodes: {");
	for (i = 0; i < g.size_g; i++)
	{
		printf("%d, ", g.arr_g[i]);
	}
	printf("}\n");
	/*
	printf("Group junk array: {");
	for (i = 0; i < g.size_g; i++)
	{
		printf("%f, ", g.tmp_vec[i]);
	}
	printf("}\n");
	*/
	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
}
void print_groups(Group* gs, int num_gs)/*debugging function*/
{
	int i;
	for (i = 0; i < num_gs; i++)
		print_group(gs[i]);
}
void print_mat_g(Graph_A* graph, Group* group) {
	int size = group->size_g;
	int i, j;
	printf("\n!!!!!!!!!!!!!mat_g!!!!!!!!!!!!!!!\n");
	for (i = 0; i < size; i++) {
		build_row_g(graph->tmp_vec, graph, group, group->arr_g[i]);
		for (j = 0; j < size; j++) {
			printf("%.0f, ", graph->tmp_vec[j]);
		}
		printf("\n");
	}
	printf("\n!!!!!!!!!!!!!!!!!!!\n\n");
}
void print_mat_Bgh(Graph_A* graph, Group* group) {
	int size = group->size_g;
	int i, j;
	printf("\n!!!!!!!!!!!!!mat_Bgh!!!!!!!!!!!!!!!\n");
	for (i = 0; i < size; i++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, group, group->arr_g[i], i);
		for (j = 0; j < size; j++) {
			printf("%.2f,	", graph->tmp_vec[j]);
		}
		printf("\n");
	}
	printf("\n!!!!!!!!!!!!!!!!!!!\n\n");
}
void print_output(char* out_file)
{
	FILE* output;
	int n;
	int num_groups, i, j;
	int num_nodes;
	int* nodes;

	output = fopen(out_file, "r");

	n = fread(&num_groups, sizeof(int), 1, output);
	assert(n == 1);
	printf("Number of groups in division: %d\n", num_groups);

	for (i = 0; i < num_groups; i++)
	{
		n = fread(&num_nodes, sizeof(int), 1, output);
		printf("Number of nodes in Group: %d\n", num_nodes);
		assert(n == 1);

		nodes = (int*)calloc(num_nodes, sizeof(int));
		n = fread(nodes, sizeof(int), num_nodes, output);
		assert(n == num_nodes);

		printf("Group nodes: ");
		for (j = 0; j < num_nodes; j++)
		{
			printf("%d, ", nodes[j]);
		}
		printf("\n");

		free(nodes);
	}

	fclose(output);
}