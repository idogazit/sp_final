/*
* ALGORITHMS Summary:
*
* A module that holds all three high-level and complex algorithms needed for the project
*
*
* devide_group_into_two			- Creates a division into two groups from a group (Algorithm 2)
* alg3							- Divides a network into modularity groups (Algorithm 3)
* alg_4							- Finds a modularity maximum partition (Algorithm 4)
*/


#include "computations.h"

#ifndef _ALGORITHMS_H
#define _ALGORITHMS_H

/*
* This function is the core function of Algorithm 2 - Divide A Group Into Two as described in the PDF
* The function gets an address for the result division, a group and a graph
* The function divides the group into two a saves the result division in the allocated memory address given
*
* @param: devision - the final division address, group - the current relevant group of nodes address, graph - the entire graph info address
*/
void devide_group_into_two(Devision* devision, Group* group, Graph_A* graph);


/*
* This function is the core function of Algorithm 3 - Network Division Into Modularity groups as described in the PDF
* The function gets a graph and an address to final partition
* The function starts with the trivial partition and each time tries to divide a group using Algorithm 2
* On success the function saves the final partition in the allocated memory address given
*
* @param: graph - the entire graph info address, O - final partition address
*/
void alg3(Graph_A* graph, Partition* O);


/*
* This function is the core function of Algorithm 4 - Modularity Maximization as described in the PDF
* The function gets a vecor s, a graph and a group, and each iteration moves one node from a group to the other until it finds the maximum modularity possible
*
* @param: vec_s - the target vector of +-1's address, graph - the entire graph info address, group - the current relevant group of nodes address
*/
void alg_4(double* vec_s, Graph_A* graph, Group* group);

#endif