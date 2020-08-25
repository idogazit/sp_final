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

typedef struct{
	int** mat_A;
	int* vec_k;
	int m;
	int num_nodes;
}Graph_A;

typedef struct{
	int* arr_g;
	int size_g;
}Group;

typedef struct{
	Group group1,group2;
}Devision;

double Epsilon = 0.00001;

int compute_M(int* vec_K, int size_K);
void compute_mat_Bg(double** mat_Bg, Graph_A* graph, Group* group);
void compute_leading_eigenvec(double* eigenvec,double** mat, int dim );
double compute_leading_eigenvalue(double* leading_eigenvec, double** mat, int dim);
void compute_vec_s(double* vec_s, double* eigenvec, int dim);
void compute_vec_f_g(double* vec_f_g, double** mat_Bg,int size_Bg);
void compute_mat_Bg_hat(double** mat_Bg_hat,Graph_A* graph, Group* group);
void build_mat_g(double** mat_g, Graph_A* graph, Group* group);
Devision devide_group_into_two(Group* group, Graph_A* graph);
double compute_mat_norm(double** mat, int dim);
void generate_rand_vec0(double* vec0, int dim);
void generate_next_vec(double* next_vec, double* curr_vec, double** mat, int dim);
double row_multiply_col(double* row, double* col, int dim);
void copy_vector(double* prev_bk, double* curr_vec,int dim);
int vectors_difference_is_small(double* vector1, double* vector2, int dim);
void kill_mat_double(double** mat, int dim);
void kill_mat_int(int** mat, int dim);
void compute_mat_multiply_vec(double* vec_comp, double** mat, double* vec, int dim);
int compute_vec1_mult_mat_mult_vec2(double* vec1, double** mat, double* vec2, int dim);
void make_vec_of_1s(double* vec,int dim);
void devide_according_to_s(Devision* devision, Group* group, double* vec_s);
double** build_mat_double(int dim);
Graph_A build_graph_A(FILE* file);
void change_mat_to_mat_tag(double** mat, double mat_norm, int dim);
void change_mat_tag_to_mat(double** mat, double mat_norm, int dim);


void make_vec_of_1s(double* vec,int dim){
	double *p;
	for(p=vec; p<&vec[dim]; p++){
		*p = 1.0;
	}
}
int compute_M(int* vec_K, int dim){
	int m=0,*p;
	for(p=vec_K; p<&vec_K[dim]; p++){
		m+=*p;
	}
	return m;
}
void compute_mat_Bg_hat(double** mat_Bg_hat,Graph_A* graph, Group* group){
	int i,dim=group->size_g;
	double *vec_f_g;
	vec_f_g = (double*)calloc(dim,sizeof(double));
	printf("before compute_mat_Bg\n");
	compute_mat_Bg(mat_Bg_hat,graph, group);
	printf("after compute_mat_Bg\n");
	printf("before compute_vec_f_g\n");
	compute_vec_f_g(vec_f_g,mat_Bg_hat,dim);
	printf("after compute_vec_f_g\n");
	for(i=0;i<dim;i++){
		mat_Bg_hat[i][i] -= vec_f_g[i];
	}
	free(vec_f_g);
}
void compute_mat_Bg(double** mat_Bg, Graph_A* graph, Group* group){
	int *pi_grp, *pj_grp ,m,dim=group->size_g;
	double **pi_Bg, *pj_Bg, tmp;
	int i,j;
	m = graph->m;
	printf("M = %d\n",m);
	printf("before build_mat_g\n");
	build_mat_g(mat_Bg,graph,group);
	printf("after build_mat_g\n");
	/*printf("vector k_g: (");*/
	for(pi_grp=group->arr_g, pi_Bg=mat_Bg; pi_grp< group->arr_g + dim; pi_Bg++, pi_grp++){
		/*printf("%d, ",graph->vec_k[*pi_grp]);*/
		for(pj_grp=group->arr_g, pj_Bg=*pi_Bg; pj_grp< group->arr_g + dim; pj_Bg++, pj_grp++){
			tmp =(double)((graph->vec_k[*pi_grp]) * (graph->vec_k[*pj_grp]))/m;
			printf("tmp = %f/n",tmp);
			*pj_Bg -= tmp;

		}
	}
	printf("mat_Bg:\n");
		for(i=0;i<group->size_g;i++){
			for(j=0;j<group->size_g;j++){
				printf("%f ",mat_Bg[i][j]);
			}
			printf("\n");
		}
}
void build_mat_g(double** mat_g, Graph_A* graph, Group* group){
	int i,j;
	int size_g=group->size_g;
	for(i=0;i<size_g;i++){
		int deg = graph->vec_k[group->arr_g[i]];
		int j=0, d=0;
		while((j <  size_g ) && (d < deg)){
				if		(group->arr_g[j] < graph->mat_A[group->arr_g[i]][d])	{j++;}
				else if	(group->arr_g[j] > graph->mat_A[group->arr_g[i]][d])	{d++;}
				else	{
						mat_g[i][j] = 1;
						j++;
						d++;
				}
		}
	}
	printf("mat_g:\n");
		for(i=0;i<group->size_g;i++){
			for(j=0;j<group->size_g;j++){
				printf("%f ",mat_g[i][j]);
			}
			printf("\n");
		}
}
Graph_A build_graph_A(FILE* file){
	int dim[1] ,**pi, *deg;
/*	int f*/
	Graph_A graph;
	/*f =*/ fread(dim,sizeof(int),1,file);
	/*check f==1 */
	graph.num_nodes = *dim;
	graph.mat_A = (int**)calloc(*dim,sizeof(int*));
	graph.vec_k = (int*)calloc(*dim,sizeof(int));
	for(pi=graph.mat_A, deg=graph.vec_k; pi<&graph.mat_A[*dim]; pi++, deg++){
	/*	f =*/ fread(deg,sizeof(int),1,file);
		/*check f==1 */
		*pi = (int*)calloc(*deg,sizeof(int));
	/*	f =*/ fread(*pi,sizeof(int),*deg,file);
		/*check f==*dim */
	}
	return graph;
}
void compute_vec_f_g(double* vec_f_g, double** mat_Bg,int size_Bg){
	double **pi_mat, *pj_mat, *p_f;
	for(pi_mat=mat_Bg, p_f=vec_f_g; p_f < vec_f_g + size_Bg; pi_mat++, p_f++){
		for(pj_mat=*pi_mat; pj_mat < *pi_mat + size_Bg; pj_mat++){
			*p_f += *pj_mat;
		}
	}
}
double compute_mat_norm(double** mat, int dim){
	double max=0, *q ,**pi, *pj, *col_sum;
	col_sum = (double*)calloc(dim,sizeof(double));

	for(pi=mat; pi<&mat[dim]; pi++){
		for(pj=*pi, q=col_sum; q<&col_sum[dim]; pj++,q++){
			*q += fabs(*pj);
		}
	}
	for(q=col_sum; q<&col_sum[dim]; q++){
		if(max<*q){max=*q;}
	}
	free(col_sum);
	return max;
}
void compute_mat_multiply_vec(double* vec_comp, double** mat, double* vec, int dim){
	double *p_vec, **p_mat;
	for(p_vec=vec_comp, p_mat=mat;
				p_vec<&vec_comp[dim]; p_vec++, p_mat++){
			*p_vec = row_multiply_col(*p_mat,vec,dim);
		}
}

int compute_vec1_mult_mat_mult_vec2(double* vec1, double** mat, double* vec2, int dim){
	double val, *tmp_vec;
	tmp_vec=(double*)calloc(dim,sizeof(double));
	compute_mat_multiply_vec(tmp_vec,mat,vec2,dim);
	val = row_multiply_col(vec1,tmp_vec,dim);
	free(tmp_vec);
	return val;
}
double compute_leading_eigenvalue(double* leading_eigenvec, double** mat, int dim){
	double val;

	val = compute_vec1_mult_mat_mult_vec2(leading_eigenvec,mat,leading_eigenvec,dim);
	val /= row_multiply_col(leading_eigenvec,leading_eigenvec, dim);

	return val;
}
void kill_mat_double(double** mat, int dim){
	double **p;
	for(p=mat;p<&mat[dim];p++){
		free(*p);
	}
	free(mat);
}
void kill_mat_int(int** mat, int dim){
	int **p;
	for(p=mat;p<&mat[dim];p++){
		free(*p);
	}
	free(mat);
}

void compute_vec_s(double* vec_s, double* eigenvec, int dim){
	double *p;
	double *q;
	for(p=vec_s, q=eigenvec; p<&vec_s[dim];p++,q++){
		if(*q>Epsilon){*p = 1.0;}
		else	{*p = 0.0;}
	}
}
void devide_according_to_s(Devision *devision, Group* group, double* vec_s){
	int *p_dev1,*p_dev2, *p_g, dim = group->size_g ;
	double *p_s;
	devision->group1.size_g=0;
	devision->group2.size_g=0;


	for(p_s=vec_s; p_s<&vec_s[dim]; p_s++){
			if(*p_s>0){
				devision->group1.size_g++;
			}
			else{
				devision->group2.size_g++;
			}
		}
	printf("\n");

	devision->group1.arr_g =(int*) calloc(devision->group1.size_g,sizeof(int));
	devision->group2.arr_g =(int*) calloc(devision->group2.size_g,sizeof(int));
	p_dev1=devision->group1.arr_g;
	p_dev2=devision->group2.arr_g;
	for(p_g=group->arr_g, p_s=vec_s; p_g<&group->arr_g[dim]; p_s++,p_g++){
		if(*p_s>Epsilon){
			*p_dev1 = *p_g;
			p_dev1++;
		}
		else{
			*p_dev2 = *p_g;
			p_dev2++;
		}
	}
}
void change_mat_to_mat_tag(double** mat, double mat_norm, int dim){
	int i;
	for(i=0;i<dim;i++){
		mat[i][i]+=mat_norm;
	}
}
void change_mat_tag_to_mat(double** mat, double mat_norm, int dim){
	int i;
		for(i=0;i<dim;i++){
			mat[i][i]-=mat_norm;
		}
}
Devision devide_group_into_two(Group* group, Graph_A* graph){
	double  *vec_s;
	double **mat_Bg_hat, *vec_eigen,val_eigen, s_Bg_s;
	double mat_norm;
	Devision devision;
	Devision *p_devision=&devision;

	int i,j;

	mat_Bg_hat = build_mat_double(group->size_g);

	vec_eigen = (double*)calloc(group->size_g,sizeof(double));

	vec_s = (double*)calloc(group->size_g,sizeof(double));
	printf("before compute_mat_Bg_hat\n");
	compute_mat_Bg_hat(mat_Bg_hat,graph,group);
	printf("after compute_mat_Bg_hat\n");
	printf("mat_Bg_hat:\n");
	for(i=0;i<group->size_g;i++){
		for(j=0;j<group->size_g;j++){
			printf("%f ",mat_Bg_hat[i][j]);
		}
		printf("\n");
	}
	mat_norm = compute_mat_norm(mat_Bg_hat,group->size_g);
	change_mat_to_mat_tag(mat_Bg_hat,mat_norm,group->size_g);
	printf("before compute_leading_eigenvec\n");
	compute_leading_eigenvec(vec_eigen,mat_Bg_hat,group->size_g);
	printf("after compute_leading_eigenvec\n");
	printf("eigen vector =( ");
	for(i=0;i<group->size_g;i++){
		printf("%f, ",vec_eigen[i]);}
	printf(")\n");
	change_mat_tag_to_mat(mat_Bg_hat,mat_norm,group->size_g);
	printf("before compute_leading_eigenvalue\n");
	val_eigen = compute_leading_eigenvalue(vec_eigen,mat_Bg_hat,group->size_g);
	printf("after compute_leading_eigenvalue\n");
	printf("eigen value = %f\n",val_eigen);
	if(val_eigen<=0){
		printf("before make_vec_of_1\n");
		make_vec_of_1s(vec_s, group->size_g);
		printf("after make_vec_of_1s\n");
	}
	else{
		printf("before compute_vec_s\n");
		compute_vec_s(vec_s,vec_eigen,group->size_g);
		printf("after compute_vec_s\n");
		s_Bg_s = compute_vec1_mult_mat_mult_vec2(vec_s,mat_Bg_hat,vec_s,group->size_g);
		if(s_Bg_s<=0){
			make_vec_of_1s(vec_s,group->size_g);
		}
	}
	printf("before devide_according_to_s\n");
	devide_according_to_s(p_devision,group,vec_s);
	printf("after devide_according_to_s\n");

	free(vec_s);
	free(vec_eigen);
	kill_mat_double(mat_Bg_hat,group->size_g);
	return devision;
}
double** build_mat_double(int dim){
	double **p, **mat;
	mat=(double**)calloc(dim,sizeof(double*));
	for(p=mat; p<&mat[dim];p++){
		*p=(double*)calloc(dim,sizeof(double));
	}

	return mat;
}
void compute_leading_eigenvec(double* eigenvec,double** mat, int dim ){
	double *vec0, *curr_vec;


	vec0 = (double*)calloc(dim,sizeof(double));
	curr_vec = (double*)calloc(dim,sizeof(double));


	generate_rand_vec0(vec0, dim);
	copy_vector(vec0, curr_vec,dim);
	printf("before generate_next_vec\n");
	generate_next_vec(eigenvec,curr_vec,mat,dim);
	printf("after generate_next_vec\n");
	while(vectors_difference_is_small(curr_vec,eigenvec, dim)==0){
		copy_vector(eigenvec,curr_vec,dim);
		generate_next_vec(eigenvec,curr_vec,mat,dim);
	}
	free(vec0);
	free(curr_vec);
}
void generate_rand_vec0(double* vec0, int dim){

	double* p;

	for (p=vec0; p<&vec0[dim]; p++){
		*p = (double) (rand()%5000)+1;
	}
}
void generate_next_vec(double* next_vec, double* curr_vec, double** mat, int dim){

	double* p_next_vec;
	double **p_mat;
	double norm = 0;
	for (p_next_vec=next_vec,p_mat=mat; p_next_vec < &next_vec[dim]; p_next_vec++,p_mat++){

		*p_next_vec = row_multiply_col(*p_mat,curr_vec,dim);
		norm += *p_next_vec * *p_next_vec;
	}
	norm = sqrt(norm);
	assert(norm!=0);
	for (p_next_vec = next_vec; p_next_vec < &next_vec[dim]; p_next_vec++){
		*p_next_vec /= norm;
	}
}
double row_multiply_col(double* row, double* col, int dim){

	double *p_row, *p_col;
	double sum = 0;
	p_col = col;

	for(p_row=row; p_row < &row[dim]; p_row++){
		sum += *p_row * *p_col;
		p_col++;
	}

	return sum;
}
void copy_vector(double* prev_bk, double* curr_vec,int dim){
	double *prev_p, *curr_p;
	curr_p = curr_vec;
	for(prev_p=prev_bk;prev_p<&prev_bk[dim];prev_p++){
		*curr_p = *prev_p;
		curr_p++;
	}
}
int vectors_difference_is_small(double* vector1, double* vector2, int dim){

	double *p1, *p2;
	p2 = vector2;

	for(p1=vector1; p1< &vector1[dim]; p1++){
		if(fabs(*p1-*p2) >= Epsilon) return 0;
		p2++;
	}

	return 1;
}
/*
int main (int argc, char* argv[]){
/*	int** mat_A = (int**)calloc(6,sizeof(int*));
	int i;
	Devision devision;
	Group group, *p_group;
	Graph_A graph, *p_graph;
	int vec_k[]={3,3,3,1,3,3};
	group.size_g = 6;
	group.arr_g =(int*) calloc(6,sizeof(int));
	group.arr_g[0] = 0;
	group.arr_g[1] = 1;
	group.arr_g[2] = 2;
	group.arr_g[3] = 3;
	group.arr_g[4] = 4;
	group.arr_g[5] = 5;



	mat_A[0] = (int*)calloc(3,sizeof(int));
	mat_A[1] = (int*)calloc(3,sizeof(int));
	mat_A[2] = (int*)calloc(3,sizeof(int));
	mat_A[3] = (int*)calloc(1,sizeof(int));
	mat_A[4] = (int*)calloc(3,sizeof(int));
	mat_A[5] = (int*)calloc(3,sizeof(int));

	mat_A[0][0]=1;
	mat_A[0][1]=2;
	mat_A[0][2]=5;

	mat_A[1][0]=0;
	mat_A[1][1]=2;
	mat_A[1][2]=4;

	mat_A[2][0]=0;
	mat_A[2][1]=1;
	mat_A[2][2]=4;

	mat_A[3][0]=5;

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

	p_graph = &graph;
	p_group = &group;
	devision = devide_group_into_two(p_group,p_graph);

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

/*	double** mat = (double**)calloc(6,sizeof(double*));
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
