#include "Mat_euristics.h"
#include "Tsp.h"
#include <stdlib.h>

void fix(int* edges, int dim, float percentage, instance* inst, CPXENVptr env, CPXLPptr lp)
{
	char lu = 'L';
	double one = 1.;
	double zero = 0.;
	for (int i = 0; i < dim; i++)
	{
		double drawn_n = (0. + rand()) / RAND_MAX;
		if (drawn_n > percentage)
		{
			//fix edges[i] to 1
			CPXchgbds(env, lp, 1, &edges[i], &lu, &one);	
		}
		else
		{
			CPXchgbds(env, lp, 1, &edges[i], &lu, &zero);
		}
	}
}

int* get_circuit(instance* inst, int* comp, int* succ)
{
	int (*pos)(int, int, instance*);
	pos = get_xpos(inst->model_type);

	int* circuit_nodes = get_component_array(1, succ, comp, inst->nnodes, &(inst->nnodes));

	int* circuit_pos = (int*)calloc(inst->nnodes, sizeof(int));

	for (int i = 1; i < inst->nnodes; i++)
		circuit_pos[i - 1] = pos(circuit_nodes[i - 1], circuit_nodes[i], inst);

	// chiudo il cerchio
	circuit_pos[inst->nnodes] = pos(circuit_nodes[inst->nnodes-1], circuit_nodes[0], inst);

	free(circuit_nodes);
	return circuit_pos;
}