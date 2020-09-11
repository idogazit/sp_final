#include "group.h"

#ifndef _PARTITION_H
#define _PARTITION_H

typedef struct {
	Group group1, group2;
}Devision;

typedef struct {
	Group* groups;
	int num_of_groups;
}Partition;


void trivial_partition(Partition* partition, int num_of_nodes);
void push_partition(Partition* partition, Group* group);
void pop_partition(Group* pop_group, Partition* partition);
void kill_partition(Partition* partition);


#endif
