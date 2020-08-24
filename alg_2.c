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
	int i,j,*pi,*pj,dim=group->size_g;
	double *vec_f_g, **mat_Bg, **p;
	vec_f_g = (double*)calloc(dim,sizeof(double));

	compute_mat_Bg(mat_Bg_hat,graph, group);
	compute_vec_f_g(vec_f_g,mat_Bg,dim);
	for(i=0;i<dim;i++){
		mat_Bg_hat[i][i] -= vec_f_g[i];
	}
	free(vec_f_g);
}
void compute_mat_Bg(double** mat_Bg, Graph_A* graph, Group* group){
	int *pi_g, *pj_g ,m,dim=group->size_g;
	double **pi_Bg, *pj_Bg;

	m = graph->m;
	build_mat_g(mat_Bg,graph,group);
	for(pi_g=group->arr_g, pi_Bg=mat_Bg; pi_g< group->arr_g + dim; pi_Bg++, pi_g++){
		for(pj_g=group->arr_g, pj_Bg=*pi_Bg; pj_g< group->arr_g + dim; pj_Bg++, pj_g++){
			*pj_Bg -=(graph->vec_k[*pi_g] * graph->vec_k[*pj_g])/m;
		}
	}
}
void build_mat_g(double** mat_g, Graph_A* graph, Group* group){
	int i;
	int *p_grp=group->arr_g, size_g=group->size_g;
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
}
Graph_A build_graph_A(FILE* file){
	int *dim,f ,**pi, *deg, *pj;
	Graph_A graph;
	f = fread(dim,sizeof(int),1,file);
	/*check f==1 */
	graph.num_nodes = *dim;
	graph.mat_A = (int**)calloc(*dim,sizeof(int*));
	graph.vec_k = (int*)calloc(*dim,sizeof(int));
	for(pi=graph.mat_A, deg=graph.vec_k; pi<&graph.mat_A[*dim]; pi++, deg++){
		f = fread(deg,sizeof(int),1,file);
		/*check f==1 */
		*pi = (int*)calloc(*deg,sizeof(int));
		f = fread(*pi,sizeof(int),*deg,file);
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
		if(*q>0){*p = 1.0;}
		else	{*p = 0.0;}
	}
}
void devide_according_to_s(Devision *devision, Group* group, double* vec_s){
	int *p_dev1,*p_dev2, *p_g, dim = group->size_g ;
	double *p_s;

	for(p_s=vec_s; p_s<&vec_s[dim]; p_s++){
			if(*p_s>0){
				devision->group1.size_g++;
			}
			else{
				devision->group2.size_g++;
			}
		}
	devision->group1.arr_g =(int*) calloc(devision->group1.size_g,sizeof(int));
	devision->group2.arr_g =(int*) calloc(devision->group2.size_g,sizeof(int));
	p_dev1=devision->group1.arr_g;
	p_dev2=devision->group2.arr_g;
	for(p_g=group->arr_g, p_s=vec_s; p_g<&group->arr_g[dim]; p_s++,p_g++){
		if(*p_s>0){
			*p_dev1 = *p_g;
			p_dev1++;
		}
		else{
			*p_dev2 = *p_g;
			p_dev2++;
		}
	}
}
Devision devide_group_into_two(Group* group, Graph_A* graph){
	double  *vec_s;
	double **mat_Bg_hat, *vec_eigen,val_eigen, **q, s_Bg_s;
	Devision devision;
	Devision *p_devision=&devision;

	mat_Bg_hat = build_mat_double(group->size_g);

	vec_eigen = (double*)calloc(group->size_g,sizeof(double));

	vec_s = (double*)calloc(group->size_g,sizeof(double));

	compute_mat_Bg_hat(mat_Bg_hat,graph,group);
	compute_leading_eigenvec(vec_eigen,mat_Bg_hat,group->size_g);
	val_eigen = compute_leading_eigenvalue(vec_eigen,mat_Bg_hat,group->size_g);
	if(val_eigen<=0){
		make_vec_of_1s(vec_s, group->size_g);
	}
	else{
		compute_vec_s(vec_s,vec_eigen,group->size_g);
		s_Bg_s = compute_vec1_mult_mat_mult_vec2(vec_s,mat_Bg_hat,vec_s,group->size_g);
		if(s_Bg_s<=0){
			make_vec_of_1s(vec_s,group->size_g);
		}
	}
	devide_according_to_s(p_devision,group,vec_s);

	free(vec_s);
	free(vec_eigen);
	kill_mat_double(mat_Bg_hat,group->size_g);
	return devision;
}
double** build_mat_double(int dim){
	double **p, *q, **mat;
	mat=(double**)calloc(dim,sizeof(double*));
	for(p=mat; p<&mat[dim];p++){
		*p=(double*)calloc(dim,sizeof(double));
	}

	return mat;
}
void compute_leading_eigenvec(double* eigenvec,double** mat, int dim ){
	double norm;
	double *vec0, *curr_vec;


	vec0 = (double*)calloc(dim,sizeof(double));
	curr_vec = (double*)calloc(dim,sizeof(double));


	generate_rand_vec0(vec0, dim);
	copy_vector(vec0, curr_vec,dim);

	generate_next_vec(eigenvec,curr_vec,mat,dim);
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
	double epsilon;
	p2 = vector2;
	epsilon = 0.00001;

	for(p1=vector1; p1< &vector1[dim]; p1++){
		if(fabs(*p1-*p2) >= epsilon) return 0;
		p2++;
	}

	return 1;
}
int main (int argc, char* argv[]){
	int** mat_A = (int**)calloc(6,sizeof(int*));
	int i;
	Devision devision;
	Group group, *p_group;
	Graph_A graph, *p_graph;
	group.size_g = 6;
	group.arr_g =(int*) calloc(6,sizeof(int));
	for(i=0;i<6;i++){
		group.arr_g[i]=i;
	}
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

	int vec_k[]={3,3,3,1,3,3};


	graph.mat_A = mat_A;
	graph.vec_k=vec_k;
	graph.num_nodes=6;
	graph.m=compute_M(vec_k,6);

	p_graph = &graph;
	p_group = &group;
	devision = devide_group_into_two(p_group,p_graph);

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
		printf("\n");

	return 0;
}
