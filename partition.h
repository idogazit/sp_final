/*
* PARTITION Summary:
*
* A module that implements a set of groups
*
*
* trivial_partition		- Creates a the trivial partition where each node is a group
* push_partition		- Inserting a group into the partition
* pop_partition			- Removing a group from the partition
* kill_partition		- Terminates all memory allocations related to the partition
*
*/

#include "group.h"

#ifndef _PARTITION_H
#define _PARTITION_H

/*
* Type used for holding a division of the graph into two groups
*/
typedef struct {
	Group group1, group2;
}Devision;


/*
* Type used for holding a partition of the graph into a num of groups
*/
typedef struct {
	Group* groups;
	int num_of_groups;
}Partition;


/*
* This function gets an address for a partition and inserts to it the trivial partition
* where each nodes is a group of his own
*
* @param : partition - the target partition address, num_of_nodes - the amount of nodes
*/
void trivial_partition(Partition* partition, int num_of_nodes);


/*
* This function gets a partition and a group and inserts the group to the partition
* Insertion is done at the end of the array
*
* @param : partition - the destination partition address, group - the target group address
*/
void push_partition(Partition* partition, Group* group);


/*
* This function gets a partition and a group address,
* removes the group from the end of the partition and copies it to the input group addres
* Deletion is done at the end of the array
*
* @param : pop_group - the destination group address, partition - the target partition address
*/
void pop_partition(Group* pop_group, Partition* partition);


/*
* This function gets an address for a partition and frees all allocated memory of it
*
* @param : partition - the target partition address
*/
void kill_partition(Partition* partition);


#endif
