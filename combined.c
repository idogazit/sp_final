/*
 * alg_2.c
 *
 *  Created on: 24 баев 2020
 *      Author: USER
 */
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
	double norm;
}Graph_A;

typedef struct {
	int* arr_g;
	int size_g;
	double* tmp_vec;
}Group;

typedef struct {
	Group group1, group2;
}Devision;

double Epsilon = 0.00001;

int compute_M(int* vec_K, int size_K);
Devision devide_group_into_two(Group* group, Graph_A* graph);
void generate_rand_vec0(double* vec0, int dim);
double row_multiply_col(double* row, double* col, int dim);
void copy_vector(double* prev_bk, double* curr_vec, int dim);
int vectors_difference_is_small(double* vector1, double* vector2, int dim);
void make_vec_of_1s(double* vec, int dim);
void devide_according_to_s(Devision* devision, Group* group, double* vec_s);
Graph_A build_graph_A(FILE* file);

void build_row_g(double* row_g, Graph_A* graph, Group* group, int node);
void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node);
double compute_leading_eigenvalue(double* leading_eigenvec, Graph_A*, Group* group);
void compute_vec_s_on_eigen_vec(double* eigenvec, int dim);
void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group);
void generate_next_vec(double* next_vec, double* curr_vec, Graph_A* graph, Group* group);
double compute_graph_norm(Graph_A* graph);

int alg3(Graph_A* graph, Group* O);
Group trivial_group(Graph_A* graph);
void output_groups(Group* O, int num_groups, char* out_file);

/*debugging functions: */
void print_graph(Graph_A graph);
void print_devision(Devision d);
void print_group(Group g);
void print_groups(Group* gs, int num_gs);

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

void compute_row_Bg_hat(double* row_Bg_hat, Graph_A* graph, Group* group, int node) {
	int* p_grp, m = graph->m, dim = group->size_g;
	double* p_Bg, f_g_node = 0.0, tmp = ((double)graph->vec_k[node]) / ((double)m);


	build_row_g(row_Bg_hat, graph, group, node);

	for (p_grp = group->arr_g, p_Bg = row_Bg_hat; p_grp < group->arr_g + dim; p_Bg++, p_grp++) {
		*p_Bg -= tmp * graph->vec_k[*p_grp];
		f_g_node += *p_Bg;
	}
	row_Bg_hat[node] -= f_g_node;
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
		else if (*p_grp > * p_node_ngbrs) {
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
Graph_A build_graph_A(FILE* file) {
	int dim[1], ** pi, * deg;
	/*	int f*/
	Graph_A graph;
	/*f =*/ fread(dim, sizeof(int), 1, file);
	/*check f==1 */
	graph.num_nodes = *dim;
	graph.mat_A = (int**)calloc(*dim, sizeof(int*));
	graph.vec_k = (int*)calloc(*dim, sizeof(int));
	for (pi = graph.mat_A, deg = graph.vec_k; pi < &graph.mat_A[*dim]; pi++, deg++) {
		/*	f =*/ fread(deg, sizeof(int), 1, file);
		/*check f==1 */
		*pi = (int*)calloc(*deg, sizeof(int));
		/*	f =*/ fread(*pi, sizeof(int), *deg, file);
		/*check f==*dim */
	}
	graph.m = compute_M(graph.vec_k, graph.num_nodes);
	graph.norm = compute_graph_norm(&graph);
	return graph;
}

double compute_graph_norm(Graph_A* graph) {
	int size = graph->num_nodes;
	Group group;
	int* p_grp;
	double max = 0, * p_row_Bgh, * p_col;
	double* col_sum = (double*)calloc(size, sizeof(double));
	group.size_g = size;
	group.tmp_vec = (double*)malloc(size * sizeof(double));
	group.arr_g = (int*)calloc(size, sizeof(int));

	for (p_grp = group.arr_g; p_grp < group.arr_g + size; p_grp++) {
		for (p_row_Bgh = group.tmp_vec, p_col = col_sum; p_col < col_sum + size; p_row_Bgh++, p_col++) {
			compute_row_Bg_hat(group.tmp_vec, graph, &group, *p_grp);
			*p_col += *p_row_Bgh;
		}
	}
	for (p_col = col_sum; p_col < &col_sum[size]; p_col++) {
		if (max < *p_col) { max = *p_col; }
	}
	free(col_sum);
	free(group.arr_g);
	free(group.tmp_vec);
	return max;
}

double compute_vec_BgH_vec(double* vec, Graph_A* graph, Group* group) {
	int size = group->size_g;
	int* p_grp;
	double val = 0.0, * p_vec1, * p_vec2, * p_Bg;
	for (p_vec1 = vec, p_grp = group->arr_g; p_grp < group->arr_g + size; p_vec1++, p_grp++) {
		compute_row_Bg_hat(group->tmp_vec, graph, group, *p_grp);
		for (p_vec2 = vec, p_Bg = group->tmp_vec; p_vec2 < vec + size; p_Bg++, p_vec2++) {
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
	devision->group1.tmp_vec = (double*)calloc(devision->group1.size_g, sizeof(double));
	devision->group2.tmp_vec = (double*)calloc(devision->group2.size_g, sizeof(double));

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
Devision devide_group_into_two(Group* group, Graph_A* graph) {
	double* vector, val_eigen, s_Bg_s;
	Devision devision;


	vector = (double*)calloc(group->size_g, sizeof(double));


	/*	printf("before compute_leading_eigenvec\n");*/
	compute_leading_eigenvec(vector, graph, group);
	/*	printf("after compute_leading_eigenvec\n");
		printf("eigen vector =( ");
		for(i=0;i<group->size_g;i++){
			printf("%f, ",vector[i]);}
		printf(")\n");
		printf("before compute_leading_eigenvalue\n");*/
	val_eigen = compute_leading_eigenvalue(vector, graph, group);
	/*	printf("after compute_leading_eigenvalue\n");
		printf("eigen value = %f\n",val_eigen);
	*/	if (val_eigen <= 0.0) {
	/*		printf("val_eigen<=0.0");
			printf("before make_vec_of_1\n");
	*/		make_vec_of_1s(vector, group->size_g);
	/*		printf("after make_vec_of_1s\n");*/
	}
	else {
		/*		printf("before compute_vec_s\n");
		*/		compute_vec_s_on_eigen_vec(vector, group->size_g);
		/*		printf("after compute_vec_s\n");
				printf("vec_s:\n(");
				for(i=0;i<group->size_g;i++){
					printf("%f, ",vector[i]);
				}
				printf("\n");
				printf("before compute_vec_BgH_vec\n");
		*/		s_Bg_s = compute_vec_BgH_vec(vector, graph, group);
		/*		printf("after compute_vec_BgH_vec\n");
				printf("s_Bg_s = %f\n",s_Bg_s);


		*/		if (s_Bg_s <= 0.0) {
		/*			printf("s_Bg_s<=0.0");
		*/			make_vec_of_1s(vector, group->size_g);
		}
	}
	/*	printf("before devide_according_to_s\n");
	*/	devide_according_to_s(&devision, group, vector);
	/*	printf("after devide_according_to_s\n");
	*/
	free(vector);
	return devision;
}

void compute_leading_eigenvec(double* eigenvec, Graph_A* graph, Group* group) {
	double* vec0, * curr_vec;
	int size = group->size_g;

	vec0 = (double*)calloc(size, sizeof(double));
	curr_vec = (double*)calloc(size, sizeof(double));


	generate_rand_vec0(vec0, size);
	copy_vector(vec0, curr_vec, size);
	/*	printf("before generate_next_vec\n");
	*/	generate_next_vec(eigenvec, curr_vec, graph, group);
	/*	printf("after generate_next_vec\n");
	*/	while (vectors_difference_is_small(curr_vec, eigenvec, size) == 0) {
		copy_vector(eigenvec, curr_vec, size);
		generate_next_vec(eigenvec, curr_vec, graph, group);
	}
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
	double norm = 0;

	for (i = 0; i < size; i++) {
		compute_row_Bg_hat(group->tmp_vec, graph, group, group->arr_g[i]);
		group->tmp_vec[i] += graph->norm;
		next_vec[i] = row_multiply_col(group->tmp_vec, curr_vec, size);
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
void copy_vector(double* prev_bk, double* curr_vec, int dim) {
	double* prev_p, * curr_p;
	curr_p = curr_vec;
	for (prev_p = prev_bk; prev_p < &prev_bk[dim]; prev_p++) {
		*curr_p = *prev_p;
		curr_p++;
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
/*
int main (int argc, char* argv[]){
	int** mat_A = (int**)calloc(6,sizeof(int*));
	int i;
	Devision devision;
	Group group;
	Graph_A graph;
	int vec_k[]={2,2,2,2,3,3};
	group.size_g = 3;
	group.arr_g =(int*) calloc(3,sizeof(int));
	group.arr_g[0] = 1;
	group.arr_g[1] = 2;
	group.arr_g[2] = 4;

	group.tmp_vec= (double*) calloc(3,sizeof(double));


	mat_A[0] = (int*)calloc(2,sizeof(int));
	mat_A[1] = (int*)calloc(2,sizeof(int));
	mat_A[2] = (int*)calloc(2,sizeof(int));
	mat_A[3] = (int*)calloc(2,sizeof(int));
	mat_A[4] = (int*)calloc(3,sizeof(int));
	mat_A[5] = (int*)calloc(3,sizeof(int));

	mat_A[0][0]=3;
	mat_A[0][1]=5;

	mat_A[1][0]=2;
	mat_A[1][1]=4;

	mat_A[2][0]=1;
	mat_A[2][1]=4;

	mat_A[3][0]=0;
	mat_A[3][1]=5;


	mat_A[4][0]=1;
	mat_A[4][1]=2;
	mat_A[4][2]=5;

	mat_A[5][0]=0;
	mat_A[5][1]=3;
	mat_A[5][2]=4;




	graph.mat_A = mat_A;
	graph.vec_k=vec_k;
	graph.num_nodes=6;
	graph.m=compute_M(vec_k,6);
	printf("before compute_graph_norm\n");
	graph.norm = compute_graph_norm(&graph);
	printf("after compute_graph_norm\n");
	printf("graph norm = %f\n",graph.norm);

	devision = devide_group_into_two(&group,&graph);

	int i;
	FILE* input_file = fopen(argv[1], "r");
	Devision devision;
	Graph_A graph = build_graph_A(input_file);
	Group group;
	group.size_g = graph.num_nodes;
	group.arr_g = (int*)calloc(graph.num_nodes,sizeof(int));
	for(i=0;i<group.size_g;i++){
		group.arr_g[i] = i;
	}
	devision = devide_group_into_two(&group,&graph);
	printf("size of group 1: %d\n",devision.group1.size_g);
	printf("group 1: ");
	for(i=0;i<devision.group1.size_g;i++){
		printf("%d ",devision.group1.arr_g[i]);
	}
	printf("\n");
	printf("size of group 2: %d\n",devision.group2.size_g);
		printf("group 2: ");
		for(i=0;i<devision.group2.size_g;i++){
			printf("%d ",devision.group2.arr_g[i]);
		}

	double** mat = (double**)calloc(6,sizeof(double*));
	double eigen_vec[6];
	int i;
	for(i =0;i<6;i++){
		mat[i]=(double*)calloc(6,sizeof(double));
		mat[i][i]=i+0.1;

	}
	compute_leading_eigenvec(eigen_vec,mat,6);

	for(i=0;i<6;i++){
		printf("%f",eigen_vec[i]);
	}

	fclose(input_file);

	printf("\n %d, %s\n",argc,argv[0]);


	return 0;
}
*/
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
		printf("\n\nP:\n\n");/*debugging prints*/
		print_groups(P, size_P);
		printf("\n\nO:\n\n");
		print_groups(O, size_O);

		g = P[size_P - 1];
		--size_P;
		
		printf("\n\nCurrent g:\n\n");/*debugging print*/
		print_group(g);

		/*here comes algorithm 2 into part*/
		temp = devide_group_into_two(&g, graph);
		
		print_devision(temp); /*debugging print*/
		
		if (temp.group1.size_g == 0 || temp.group2.size_g == 0)
		{
			/*add group g into O*/
			
			/*O[size_O] = g;
			size_O++;*/
			*O = g;
			O++;
		}
		else
		{
			if (temp.group1.size_g == 1)
			{
				/*add group 1 into O*/
				
				/*O[size_O] = temp.group1;
				size_O++;*/

				*O = temp.group1;
				O++;
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
				
				/*O[size_O] = temp.group2;
				size_O++;*/

				*O = temp.group1;
				O++;
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
	full.tmp_vec = calloc(full.size_g, sizeof(double));

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
void print_graph(Graph_A graph) /*debugging function*/
{
	int i = 0, j = 0;
	printf("=========Graph========\n");
	printf("M = %d\n", graph.m);
	printf("Number of nodes = %d\n", graph.num_nodes);
	printf("Norm = %f\n", graph.norm);
	printf("Vec K = [");
	for (i = 0; i < graph.num_nodes; i++)
	{
		printf("%d,", graph.vec_k[i]);
	}
	printf("]\n");

	printf("Mat: \n");
	for (i = 0; i < graph.num_nodes; i++)
	{
		printf("Node #%d neighbors: ", i);
		for (j = 0; j < graph.vec_k[i]; j++)
		{
			printf("%d,", graph.mat_A[i][j]);
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


int main(int argc, char* argv[]) {
	Graph_A graph;
	Group* final_div = malloc(sizeof(Group*));
	int num_groups;
	FILE* input_file;

	if (argc != 3)
	{
		printf("Error in inputs\n");
		return -1;
	}
	
	input_file = fopen(argv[1], "r"); /*opens input file and creates graph A*/
	graph = build_graph_A(input_file);
	fclose(input_file);
	
	final_div->arr_g = calloc(graph.num_nodes, sizeof(Group));
	final_div->tmp_vec = calloc(graph.num_nodes, sizeof(Group));
	final_div->size_g = 0;

	num_groups = alg3(&graph, final_div); /*calls for algorithm 3*/

	output_groups(final_div, num_groups, argv[2]); /*writes final groups division into output file*/
	
	free(final_div->arr_g);
	free(final_div->tmp_vec);
	
	return 0;
}

