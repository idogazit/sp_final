/*
 * alg_2.c
 *
 *  Created on: 24 ���� 2020
 *      Author: USER
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <prints.h>

typedef struct {
	int** mat_A;
	int* vec_k;
	int m;
	int num_nodes;
	double norm;
	double* tmp_vec;

}Graph_A;

typedef struct {
	int* arr_g;
	int size_g;
}Group;

typedef struct {
	Group group1, group2;
}Devision;

typedef struct {
	Group* groups;
	int num_of_groups;
}Partition;

double Epsilon = 0.00001;

int compute_M(int* vec_K, int size_K);
void devide_group_into_two(Devision* devision, Group* group, Graph_A* graph);
void generate_rand_vec0(double* vec0, int dim);
double row_multiply_col(double* row, double* col, int dim);
void copy_vector(double* org_vec, double* new_vec, int dim);
int vectors_difference_is_small(double* vector1, double* vector2, int dim);
void make_vec_of_1s(double* vec, int dim);
void devide_according_to_s(Devision* devision, Group* group, double* vec_s);
void build_graph_A(Graph_A* graph, FILE* file);
double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group);
void trivial_partition(Partition* partition,int num_of_nodes);
void build_row_g(double* row_g, Graph_A* graph, Group* group, int node);
void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node,int i);
void compute_row_Bg(double* row_Bg, Graph_A* graph, Group* group, int node);
double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A*, Group* group);
void compute_vec_s_on_eigen_vec(double* eigenvec, int dim);
void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group);
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group);
double compute_graph_norm(Graph_A* graph, Group* trivial_group);
void kill_graph(Graph_A* graph);

void alg3(Graph_A* graph, Partition* O);
void trivial_group(Group* triv_group,Graph_A* graph);
void output_groups(Group* O, int num_groups, char* out_file);

void alg_4(double* vec_s, Graph_A* graph, Group* group);
int max_ind(double* arr, int len);
void reset_unmoved(int* unmoved, int len);
void compute_vec_rowBgH_s(double* vec_rowBgH_s, Graph_A* graph, Group* group, double* vec_s);
void compute_vec_Bg_ii(double* vec_Bg_ii, Graph_A* graph, Group* group);
/*debugging functions: */
void print_graph(Graph_A* graph);
void print_devision(Devision d);
void print_group(Group g);
void print_groups(Group* gs, int num_gs);
void print_mat_g(Graph_A* graph,Group* group);
void print_mat_Bgh(Graph_A* graph,Group* group);
void print_output(char* out_file);


void push_partition(Partition* partition, Group* group);
void pop_partition(Group* pop_group,Partition* partition);
void kill_partition(Partition* partition);
void copy_group(Group* org_group, Group* new_group);

void make_vec_of_1s(double* vec, int dim) {
	double* p;
	for (p = vec; p < &vec[dim]; p++) {
		*p = 1.0;
	}
}
int compute_M(int* vec_K, int dim) {
	int m = 0, * p;
	for (p = vec_K; p < &vec_K[dim]; p++) {
		m += *p;
	}
	return m;
}
void trivial_partition(Partition* partition,int num_of_nodes){
	int *i;
	for(i=0;i<num_of_nodes;i++){
		partition->groups[i].arr_g =(int*)malloc(sizeof(int));
		partition->groups[i].arr_g[0] = i;
		partition->groups[i].size_g = 1;
	}
	partition->num_of_groups = num_of_nodes;
}


void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node,int i) {
	int* p_grp, m = graph->m, dim = group->size_g;
	double* p_Bg, f_g_node = 0.0, tmp =(double)graph->vec_k[node] / (double)m;


	build_row_g(row_Bg_hat, graph, group, node);

	for (p_grp = group->arr_g, p_Bg = row_Bg_hat; p_grp < group->arr_g + dim; p_Bg++, p_grp++) {
		*p_Bg -= tmp *(double)graph->vec_k[*p_grp];
		f_g_node += *p_Bg;
	}
	row_Bg_hat[i] -= f_g_node;

}

void build_row_g(double* row_g, Graph_A* graph, Group* group, int node) {
	int size = group->size_g;
	int deg = graph->vec_k[node];
	int* p_node_ngbrs = graph->mat_A[node], * p_grp = group->arr_g;
	double* p_rowg = row_g;


	for (; p_rowg < row_g + size; p_rowg++) {
		*p_rowg = 0.0;
	}
	p_rowg = row_g;
	while ((p_grp < group->arr_g + size) && (p_node_ngbrs < graph->mat_A[node] + deg)) {
		if (*p_grp < *p_node_ngbrs) {
			p_rowg++;
			p_grp++;
		}
		else if (*p_grp > *p_node_ngbrs) {
			p_node_ngbrs++;
		}
		else {
			*p_rowg = 1.0;
			p_rowg++;
			p_grp++;
			p_node_ngbrs++;
		}
	}
}

 void build_graph_A(Graph_A* graph,FILE* file) {
	int dim[1], ** pi, * deg;
	/*	int f*/

	/*f =*/ fread(dim, sizeof(int), 1, file);
	/*check f==1 */
	graph->num_nodes = *dim;
	graph->mat_A = (int**)calloc(*dim, sizeof(int*));
	graph->vec_k = (int*)calloc(*dim, sizeof(int));
	graph->tmp_vec = (double*)calloc(*dim,sizeof(double));
	for (pi = graph->mat_A, deg = graph->vec_k; pi < &(graph->mat_A[*dim]); pi++, deg++) {
		/*	f =*/ fread(deg, sizeof(int), 1, file);
		/*check f==1 */
		*pi = (int*)calloc(*deg, sizeof(int));
		/*	f =*/ fread(*pi, sizeof(int), *deg, file);
		/*check f==*dim */
	}
	graph->m = compute_M(graph->vec_k, graph->num_nodes);

}

double compute_graph_norm(Graph_A* graph, Group* trivial_group) {
	int size = graph->num_nodes;
	int* p_grp, i=0;
	double max = 0.0, * p_row_Bgh, * p_col;
	double* col_sum = (double*)calloc(size, sizeof(double));


	for (p_grp = trivial_group->arr_g; p_grp < trivial_group->arr_g + size; p_grp++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, trivial_group, *p_grp,i++);
		for (p_row_Bgh = graph->tmp_vec, p_col = col_sum; p_col < col_sum + size; p_row_Bgh++, p_col++) {
			*p_col += fabs(*p_row_Bgh);
		}
	}
	for (p_col = col_sum; p_col < &col_sum[size]; p_col++) {
		if (max < *p_col) { max = *p_col; }
	}
	free(col_sum);
	return max;
}

double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group) {
	int size = group->size_g;
	int* p_grp, i=0;
	double val = 0.0, * p_vec1, * p_vec2, * p_Bg;
	for (p_vec1 = vec, p_grp = group->arr_g; p_grp < group->arr_g + size; p_vec1++, p_grp++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, group, *p_grp,i++);
		for (p_vec2 = vec, p_Bg = graph->tmp_vec; p_vec2 < vec + size; p_Bg++, p_vec2++) {
			val += (double)*p_vec1 * *p_vec2 * *p_Bg;;
		}
	}
	return val;
}

double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A* graph, Group* group) {
	double val;
	val = compute_vec_BgH_vec(leading_eigenvec, graph, group);
	val /= row_multiply_col(leading_eigenvec, leading_eigenvec, group->size_g);

	return val;
}


void compute_vec_s_on_eigen_vec(double* eigenvec, int dim) {
	double* q;
	for (q = eigenvec; q < &eigenvec[dim]; q++) {
		if (*q > Epsilon) { *q = 1.0; }
		else { *q = -1.0; }
	}
}
void devide_according_to_s(Devision* devision, Group* group, double* vec_s) {
	int* p_dev1, * p_dev2, * p_g, dim = group->size_g;
	double* p_s;
	devision->group1.size_g = 0;
	devision->group2.size_g = 0;


	for (p_s = vec_s; p_s < &vec_s[dim]; p_s++) {
		if (*p_s > 0.0) {
			devision->group1.size_g++;
		}
		else {
			devision->group2.size_g++;
		}
	}

	devision->group1.arr_g = (int*)calloc(devision->group1.size_g, sizeof(int));
	devision->group2.arr_g = (int*)calloc(devision->group2.size_g, sizeof(int));

	p_dev1 = devision->group1.arr_g;
	p_dev2 = devision->group2.arr_g;

	for (p_g = group->arr_g, p_s = vec_s; p_g < &group->arr_g[dim]; p_s++, p_g++) {
		if (*p_s > Epsilon) {
			*p_dev1 = *p_g;
			p_dev1++;
		}
		else {
			*p_dev2 = *p_g;
			p_dev2++;
		}
	}
}
void devide_group_into_two(Devision* devision, Group* group, Graph_A* graph) {
	double* vector, val_eigen, s_BgH_s;
	vector = (double*)calloc(group->size_g, sizeof(double));

	compute_leading_eigenvec(vector, graph, group);

		
	val_eigen = compute_leading_eigenvalue(vector, graph, group);
		
	if (val_eigen <= 0.0) {
		make_vec_of_1s(vector, group->size_g);
	}
	else {
		compute_vec_s_on_eigen_vec(vector, group->size_g);
		s_BgH_s = compute_vec_BgH_vec(vector, graph, group);
		if (s_BgH_s <= 0.0) {
			make_vec_of_1s(vector, group->size_g);
		}
		else
		{
			alg_4(vector, graph, group);
		}
	}
	devide_according_to_s(devision, group, vector);

	free(vector);
}

void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group) {
	double* vec0, * curr_vec;
	int size = group->size_g;
	int i;
	vec0 = (double*)calloc(size, sizeof(double));
	curr_vec = (double*)calloc(size, sizeof(double));


	generate_rand_vec0(vec0, size);
	copy_vector(vec0, curr_vec, size);
	generate_next_vec(eigenvec, curr_vec, graph, group);
	while ((vectors_difference_is_small(curr_vec, eigenvec, size) == 0)) {
		++i;
		copy_vector(eigenvec, curr_vec, size);
		generate_next_vec(eigenvec, curr_vec, graph, group);
	}
	/*printf("\n\n number of power iterations: %d\n",i);
	printf("group size: %d\n\n",group->size_g);*/
	free(vec0);
	free(curr_vec);
}
void generate_rand_vec0(double* vec0, int dim) {

	double* p;

	for (p = vec0; p < &vec0[dim]; p++) {
		*p = (double)(rand() % 5000) + 1;
	}
}
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group) {
	int size = group->size_g;
	int i;
	double* p_next_vec;
	double norm = 0.0;
	for (i = 0; i < size; i++) {
		compute_row_Bg_hat(graph->tmp_vec, graph, group, group->arr_g[i],i);

		graph->tmp_vec[i] += graph->norm;
		next_vec[i] = row_multiply_col(graph->tmp_vec, curr_vec, size);
		norm += next_vec[i] * next_vec[i];
	}
	norm = sqrt(norm);
	assert(norm != 0);
	for (p_next_vec = next_vec; p_next_vec < &next_vec[size]; p_next_vec++) {
		*p_next_vec /= norm;
	}

}
double row_multiply_col(double* row, double* col, int dim) {

	double* p_row, * p_col;
	double sum = 0;
	p_col = col;

	for (p_row = row; p_row < &row[dim]; p_row++) {
		sum += *p_row * *p_col;
		p_col++;
	}

	return sum;
}
void copy_vector(double* org_vec, double* new_vec, int dim) {
	double* org_p, * new_p;
	new_p = new_vec;
	for (org_p = org_vec; org_p < &org_vec[dim]; org_p++) {
		*new_p = *org_p;
		new_p++;
	}
}
void copy_vector_int(int* org_vec, int* new_vec, int dim) {
	int* org_p, * new_p;
	new_p = new_vec;
	for (org_p = org_vec; org_p < &org_vec[dim]; org_p++) {
		*new_p = *org_p;
		new_p++;
	}
}
int vectors_difference_is_small(double* vector1, double* vector2, int dim) {

	double* p1, * p2;
	p2 = vector2;

	for (p1 = vector1; p1 < &vector1[dim]; p1++) {
		if (fabs(*p1 - *p2) >= Epsilon) return 0;
		p2++;
	}

	return 1;
}

void alg3(Graph_A* graph, Partition* O)
{
	Group triv_g ;
	Partition P;
	Group grp;
	Devision temp_dev;
	P.groups = (Group*)calloc(graph->num_nodes, sizeof(Group));
	P.num_of_groups = 0;

	trivial_group(&triv_g, graph);
	graph->norm = compute_graph_norm(graph,&triv_g);
	push_partition(&P,&triv_g);
	while (P.num_of_groups > 0)
	{

		pop_partition(&grp, &P);
		
		/*here comes algorithm 2 into part*/
		devide_group_into_two(&temp_dev, &grp, graph);/*********/
		
		
		if (temp_dev.group1.size_g == 0 || temp_dev.group2.size_g == 0)
		{
			/*add group grp into O*/
			
			push_partition(O,&grp);
		}
		else
		{
			if (temp_dev.group1.size_g == 1)
			{
				/*add group 1 into O*/
				printf("\n\nIn if for group 1\n\n");

				push_partition(O,&(temp_dev.group1));
			}
			else
			{
				/*add group 1 into P*/

				push_partition(&P,&(temp_dev.group1));
			}
			if (temp_dev.group2.size_g == 1)
			{
				/*add group 2 into O*/
				printf("\n\nIn if for group 2\n\n");
				push_partition(O,&(temp_dev.group2));
			}
			else
			{
				/*add group 2 into P*/

				push_partition(&P,&(temp_dev.group2));
			}
		}
		free(temp_dev.group1.arr_g);
		free(temp_dev.group2.arr_g);
		free(grp.arr_g);

	}
	
	free(triv_g.arr_g);
	kill_partition(&P);
}
void trivial_group(Group* triv_group,Graph_A* graph)
{
	int i;

	triv_group->size_g = graph->num_nodes;
	triv_group->arr_g = calloc(graph->num_nodes, sizeof(int));

	for (i = 0; i < triv_group->size_g; i++)
	{
		triv_group->arr_g[i] = i;
	}
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
		n = fwrite(g.arr_g, sizeof(int), g.size_g, out_div);
		if (n != g.size_g)
		{
			printf("Error writing to output file\n");
		}
	}

	fclose(out_div);
}
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
void push_partition(Partition* partition, Group* group){
	Group new_group;
	copy_group(group,&new_group);
	partition->groups[(partition->num_of_groups)++] = new_group;
}
void pop_partition(Group* pop_group,Partition* partition){
	copy_group(&(partition->groups[partition->num_of_groups - 1]),pop_group);
	free(partition->groups[--(partition->num_of_groups)].arr_g);
}
void kill_partition(Partition* partition){
	Group *p_groups;
	for(p_groups = partition->groups; p_groups < partition->groups + partition->num_of_groups; p_groups++){
		free(p_groups->arr_g);
	}
	free(partition->groups);
}
void kill_graph(Graph_A* graph){
	int **p_node;
	free(graph->tmp_vec);
	for(p_node=graph->mat_A;p_node < graph->mat_A + graph->num_nodes; p_node++){
		free(*p_node);
	}
	free(graph->mat_A);

	free(graph->vec_k);
}
void copy_group(Group* org_group, Group* new_group){
		new_group->size_g = org_group->size_g;
		new_group->arr_g = (int*)calloc(new_group->size_g,sizeof(int));
		copy_vector_int(org_group->arr_g,new_group->arr_g,new_group->size_g);
}
void print_mat_g(Graph_A* graph,Group* group){
	int size = group->size_g;
	int i,j;
	printf("\n!!!!!!!!!!!!!mat_g!!!!!!!!!!!!!!!\n");
	for(i=0;i<size;i++){
		build_row_g(graph->tmp_vec,graph,group,group->arr_g[i]);
		for(j=0;j<size;j++){
			printf("%.0f, ",graph->tmp_vec[j]);
		}
		printf("\n");
	}
	printf("\n!!!!!!!!!!!!!!!!!!!\n\n");
}
void print_mat_Bgh(Graph_A* graph,Group* group){
	int size = group->size_g;
	int i,j;
	printf("\n!!!!!!!!!!!!!mat_Bgh!!!!!!!!!!!!!!!\n");
	for(i=0;i<size;i++){
		compute_row_Bg_hat(graph->tmp_vec,graph,group,group->arr_g[i],i);
		for(j=0;j<size;j++){
			printf("%.2f,	",graph->tmp_vec[j]);
		}
		printf("\n");
	}
	printf("\n!!!!!!!!!!!!!!!!!!!\n\n");
}

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
int max_ind(double* arr, int len)
{
	int max = 0, i = 0;
	
	for (i = 0; i < len; i++)
	{
		if (arr[i] > arr[max])
		{
			max = i;
		}
	}
	return max;
}
void reset_unmoved(int* unmoved, int len)
{
	int i;
	for (i = 0; i < len; i++)
	{
		unmoved[i] = 0;
	}
}
void compute_vec_rowBgH_s(double* vec_rowBgH_s, Graph_A* graph, Group* group, double* vec_s){
	int i;

	for(i=0;i<group->size_g;i++){
		compute_row_Bg_hat(graph->tmp_vec, graph, group,group->arr_g[i],i);
		vec_rowBgH_s[i] = row_multiply_col(graph->tmp_vec, vec_s,group->size_g);
	}
}
void compute_vec_Bg_ii(double* vec_Bg_ii, Graph_A* graph, Group* group){
	double *p_vec;
	int *p_grp;

	for(p_vec=vec_Bg_ii,p_grp=group->arr_g;p_vec<vec_Bg_ii + group->size_g;p_vec++,p_grp++){
		*p_vec = (double)(graph->vec_k[*p_grp] * graph->vec_k[*p_grp]) / (double)graph->m;
	}
}
void compute_row_Bg(double* row_Bg, Graph_A* graph, Group* group, int node) {
	int* p_grp, m = graph->m, dim = group->size_g;
	double* p_Bg, tmp =(double)graph->vec_k[node] / (double)m;

	build_row_g(row_Bg, graph, group, node);

	for (p_grp = group->arr_g, p_Bg = row_Bg; p_grp < group->arr_g + dim; p_Bg++, p_grp++) {
		*p_Bg -= tmp *(double)graph->vec_k[*p_grp];
	}
}
void alg_4(double* vec_s, Graph_A* graph, Group* group)
{
	int k, i, j = 0;
	int score_max_ind, improve_max_ind;
	int* unmoved;
	double* score;
	double* improve;
	int* indices;
	double /*Q0,*/ delta_Q = 0;
	int len = group->size_g;

	double rowBg_s, *vec_Bg_ii ;

	unmoved = (int*)calloc(len, sizeof(int));
	indices = (int*)calloc(len, sizeof(int));
	score = (double*)calloc(len, sizeof(double));
	improve = (double*)calloc(len, sizeof(double));

	vec_Bg_ii = (double*)calloc(len, sizeof(double));

	compute_vec_Bg_ii(vec_Bg_ii,graph,group);


	do 
	{
		/*step 1*/
		reset_unmoved(unmoved, len);
		improve_max_ind = 0;
		/*step 2*/
		for (i = 0; i < len; i++)
		{
			/*a*/
			/*Q0 = compute_vec_BgH_vec(vec_s, graph, group);
			*/score_max_ind = -1;
			/*b*/
			for (k = 0; k < len; k++)
			{
				if (unmoved[k] == 0)
				{
					compute_row_Bg(graph->tmp_vec,graph,group,group->arr_g[k]);
					rowBg_s = row_multiply_col(graph->tmp_vec,vec_s,len);
					score[k] = 4 * (vec_Bg_ii[k] - (vec_s[k] * rowBg_s));

					/*vec_s[k] *= -1;
					score[k] = compute_vec_BgH_vec(vec_s,graph,group) - Q0;
					vec_s[k] *= -1;
					*/if(score_max_ind == -1){
						score_max_ind = k;
					}
					if(score[score_max_ind] < score[k]){
						score_max_ind = k;
					}
				}
			}

			/*c*/
			/*score_max_ind = max_ind(score, len);
			*/
			/*d*/
			vec_s[score_max_ind] *= -1;

			/*e*/
			indices[i] = score_max_ind;
			
			/*f*/
			if (i == 0)
			{
				improve[i] = score[score_max_ind];
			}
			else
			{
				improve[i] = improve[i - 1] + score[score_max_ind];
			}
			if(improve[improve_max_ind] < improve[i]){
				improve_max_ind = i;
			}
			/*g*/
			unmoved[score_max_ind] = -1;
		}

		/*step 3*/
	/*	improve_max_ind = max_ind(improve, len);
	*/
		if (improve[len - 1] == 0 && improve[improve_max_ind] == 0)
		{
			improve_max_ind = len - 1;
		}

		/*step 4*/
		for (i = len - 1; i >= improve_max_ind + 1; i--)
		{
			j = indices[i];
			vec_s[j] *= -1;
		}

		/*step 5*/
		if (improve_max_ind == len - 1)
		{
			delta_Q = 0;
		}
		else
		{
			delta_Q = improve[improve_max_ind];
		}
	} while (delta_Q > 0);
	
	free(unmoved);
	free(indices);
	free(score);
	free(improve);

	free(vec_Bg_ii);
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
