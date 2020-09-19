#include "computations.h"
#include "algorithms.h"
#define Epsilon 0.00001

void alg3(Graph_A* graph, Partition* O)
{
	Group triv_g;
	Partition P;
	Group grp;
	Devision temp_dev;
	P.groups = (Group*)malloc(graph->num_nodes * sizeof(Group));
	P.num_of_groups = 0;

	trivial_group(&triv_g, graph);
	graph->norm = compute_graph_norm(graph, &triv_g);
	push_partition(&P, &triv_g);
	while (P.num_of_groups > 0)
	{

		pop_partition(&grp, &P);

		/*here comes algorithm 2 into part*/
		devide_group_into_two(&temp_dev, &grp, graph);


		if (temp_dev.group1.size_g == 0 || temp_dev.group2.size_g == 0)
		{
			/*add group grp into O*/

			push_partition(O, &grp);
		}
		else
		{
			if (temp_dev.group1.size_g == 1)
			{
				/*add group 1 into O*/
				printf("\n\nIn if for group 1\n\n");

				push_partition(O, &(temp_dev.group1));
			}
			else
			{
				/*add group 1 into P*/

				push_partition(&P, &(temp_dev.group1));
			}
			if (temp_dev.group2.size_g == 1)
			{
				/*add group 2 into O*/
				printf("\n\nIn if for group 2\n\n");
				push_partition(O, &(temp_dev.group2));
			}
			else
			{
				/*add group 2 into P*/

				push_partition(&P, &(temp_dev.group2));
			}
		}
		free(temp_dev.group1.arr_g);
		free(temp_dev.group2.arr_g);
		free(grp.arr_g);

	}

	free(triv_g.arr_g);
	kill_partition(&P);
}

void alg_4(double* vec_s, Graph_A* graph, Group* group)
{
	int k, i, j = 0;
	int score_max_ind, improve_max_ind;
	int* unmoved;
	double* score;
	double* improve;
	int* indices;
	double delta_Q = 0;
	int len = group->size_g;

	double rowBg_s, * vec_Bg_ii;

	unmoved = (int*)malloc(len * sizeof(int));
	indices = (int*)malloc(len * sizeof(int));
	score = (double*)malloc(len * sizeof(double));
	improve = (double*)malloc(len * sizeof(double));

	vec_Bg_ii = (double*)malloc(len * sizeof(double));

	compute_vec_Bg_ii(vec_Bg_ii, graph, group);


	do
	{
		/*step 1*/
		reset_unmoved(unmoved, len);
		improve_max_ind = 0;
		/*step 2*/
		for (i = 0; i < len; i++)
		{
			score_max_ind = -1;
			for (k = 0; k < len; k++)
			{
				if (unmoved[k] == 0)
				{
					compute_row_Bg(graph->tmp_vec, graph, group, group->arr_g[k]);
					rowBg_s = row_multiply_col(graph->tmp_vec, vec_s, len);
					score[k] = 4 * (vec_Bg_ii[k] - (vec_s[k] * rowBg_s));

					if (score_max_ind == -1) {
						score_max_ind = k;
					}
					if (score[score_max_ind] < score[k]) {
						score_max_ind = k;
					}
				}
			}

			vec_s[score_max_ind] *= -1;

			indices[i] = score_max_ind;

			if (i == 0)
			{
				improve[i] = score[score_max_ind];
			}
			else
			{
				improve[i] = improve[i - 1] + score[score_max_ind];
			}
			if (improve[improve_max_ind] < improve[i]) {
				improve_max_ind = i;
			}

			unmoved[score_max_ind] = -1;
		}

		/*step 3*/

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
	} while (delta_Q > Epsilon);

	free(unmoved);
	free(indices);
	free(score);
	free(improve);
	free(vec_Bg_ii);
}

void devide_group_into_two(Devision* devision, Group* group, Graph_A* graph) {
	double* vector, val_eigen, s_BgH_s;
	vector = (double*)malloc(group->size_g * sizeof(double));

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