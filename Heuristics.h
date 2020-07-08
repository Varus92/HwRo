#include "Tsp.h"

int* swap_2opt(instance* inst, int* circuit_nodes, int node_a, int node_b, int node_c, int node_d);
double get_obj_val(instance* inst, int* nodes_cycle);
int variable_neighborhood_search(instance* inst, double time_limit);