#include "Grasp.h"

#include <stdbool.h>


bool time_limit(instance* inst);
void Grasp(instance* inst, int* succ, int* comp, int* ncomp);

void heuristicSolver(instance* inst)
{
	inst->tstart = second();
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp;
	Grasp(inst, succ, comp, &ncomp);

	//Plot soluzione

	free(comp);
}


void Grasp(instance* inst, int* succ, int* comp, int* ncomp)
{
	//Componenti iniziali
	double bestsol = 999999;
	*ncomp = 1;
	int* tmp_succ = (int*)calloc(inst->nnodes, sizeof(int));
	
	for (int i = 0; i < inst->nnodes; i++) 
	{
		tmp_succ[i] = -1;
		comp[i] = 1;
	}

	while (!time_limit(inst))
	{
		if (VERBOSE > 100) printf("Trying a new starting node...");

		//Choose starting point
		srand(time(NULL));
		int p = rand() % inst->nnodes;
		int start = p;
		bool done = false;
		int solcost = 0;

		while (!done)
		{
			//Il prossimo nodo da visitare sarà scelto tra 3 nodi possibili
			int min_dist1 = 9999999;
			int min_dist2 = 9999999;
			int min_dist3 = 9999999;
			int next_node[3] = { -1 };
			int final_node = -1;

			for (int i = 0; i < inst->nnodes; i++)
			{
				if (i != p && min_dist1 > dist(p, i, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = min_dist2;
					min_dist2 = min_dist1;
					min_dist1 = dist(p, i, inst);
					next_node[0] = i;
				}
				else if (i != p && min_dist2 > dist(p, i, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = min_dist2;
					min_dist2 = dist(p, i, inst);
					next_node[1] = i;
				}
				else if (i != p && min_dist3 > dist(p, i, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = dist(p, i, inst);
					next_node[2] = i;
				}
			}

			//Sceglamo in modo random Tra i 3 valori
			//srand(time(NULL));
			int n = rand() % 3;
			final_node = next_node[n];

			//Verifico se il circuito è chiuso
			if (final_node == -1)
			{
				tmp_succ[p] = start;
				solcost += dist(p, start, inst);
				done = true;
			}
			else //Aggiorno i valori
			{
				solcost += dist(p, final_node, inst);
				tmp_succ[p] = final_node;
				p = final_node;
			}
		}

		//Aggiorno la migliore soluzione
		if (solcost < bestsol)
		{
			for (int i = 0; i < inst->nnodes; i++)
				succ[i] = tmp_succ[i];
			bestsol = solcost;
		}

		//Reset tmp_succ
		for (int i = 0; i < inst->nnodes; i++)
			tmp_succ[i] = -1;

		inst->best_lb = solcost;
	}
	free(tmp_succ);
}


bool time_limit(instance* inst)
{
	double tspan = second() - inst->tstart;
	if (tspan > inst->timelimit)
	{
		if (VERBOSE >= 100) printf("\n\n$$$ time limit of %10.1lf sec.s expired after %10.1lf sec.s $$$\n\n", inst->timelimit, tspan);
		return true;
	}
	return false;
}
