#include "group.h"
#include "partition.h"

void trivial_partition(Partition* partition, int num_of_nodes) {
	int i;
	for (i = 0; i < num_of_nodes; i++) {
		partition->groups[i].arr_g = (int*)malloc(sizeof(int));
		partition->groups[i].arr_g[0] = i;
		partition->groups[i].size_g = 1;
	}
	partition->num_of_groups = num_of_nodes;
}

void push_partition(Partition* partition, Group* group) {
	Group new_group;
	copy_group(group, &new_group);
	partition->groups[(partition->num_of_groups)++] = new_group;
}

void pop_partition(Group* pop_group, Partition* partition) {
	copy_group(&(partition->groups[partition->num_of_groups - 1]), pop_group);
	free(partition->groups[--(partition->num_of_groups)].arr_g);
}

void kill_partition(Partition* partition) {
	Group* p_groups;
	for (p_groups = partition->groups; p_groups < partition->groups + partition->num_of_groups; p_groups++) {
		free(p_groups->arr_g);
	}
	free(partition->groups);
}