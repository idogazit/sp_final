/*
* GROUP Summary:
*
* A module that implements a group class
*
*
* copy_group			- Copies of a group
* trivial_group			- Creates a the trivial nodes group: {0,1,...,n-1}
* output_groups			- Writes a set of groups into output file
* copy_vector_int		- Copies a vector
*
*/

#include "graph.h"

#ifndef _GROUP_H
#define _GROUP_H

/*
* Type used for containing all relevant information regarding to a group
*/
typedef struct {
	int* arr_g;
	int size_g;
}Group;


/*
* This function gets a group and copies to it to another group
*
* @param : org_group - the original group address, new_group - the destination group address
*/
void copy_group(Group* org_group, Group* new_group);


/*
* This function gets an address for a group and inserts to it the trivial group {0,1,...,n-1}
*
* @param : triv_group - the target group address, graph - the source graph address
*/
void trivial_group(Group* triv_group, Graph_A* graph);


/*
* This function gets a set of groups and writes them into the output file
* in correct way as described in the PDF.
* If an error occurs then the function exits uncleaned.
*
* @param: O - the set of groups address, num_groups - the amount of groups, out_file - output file name address
* @return: On success the function returns 0, otherwise return -1
*/
int output_groups(Group* O, int num_groups, char* out_file);


/*
* This function gets a vector and copies to it to another vector
*
* @param : org_vec - the original vector address, new_vec - the destination vec address, dim - the sizee of the vectors
*/
void copy_vector_int(int* org_vec, int* new_vec, int dim);

#endif