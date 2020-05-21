#include "Grasp.h"

#include <stdbool.h>


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
	double best_sol = 9999999;
	*ncomp = 1;
	int* tmp_succ = (int*)calloc(inst->nnodes, sizeof(int));
	
	for (int i = 0; i < inst->nnodes; i++) 
	{
		tmp_succ[i] = -1;
		comp[i] = 1;
	}

	double t_span = 0; 

	while (t_span < inst->timelimit)
	{
		//if (VERBOSE > 100) printf("Trying a new starting node...");

		//Choose starting point
		srand(time(NULL));
		int n = rand() % inst->nnodes;
		int start = n;
		bool done = false;
		int sol_cost = 0;

		while (!done)
		{
			//Il prossimo nodo da visitare sarà scelto tra 3 nodi possibili
			int min_dist1 = 99999999;
			int min_dist2 = 99999999;
			int min_dist3 = 99999999;
			int next_node[3] = {-1,-1,-1};
			int final_node = -1;

			for (int i = 0; i < inst->nnodes; i++)
			{
				if (i != n && min_dist1 > dist(i, n, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = min_dist2;
					min_dist2 = min_dist1;
					min_dist1 = dist(i, n, inst);
					next_node[0] = i;
				}
				else if (i != n && min_dist2 > dist(i, n, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = min_dist2;
					min_dist2 = dist(n, i, inst);
					next_node[1] = i;
				}
				else if (i != n && min_dist3 > dist(i, n, inst) && tmp_succ[i] == -1)
				{
					min_dist3 = dist(i, n, inst);
					next_node[2] = i;
				}
			}

			//Sceglamo in modo random Tra i 3 valori
			if (next_node[0] != -1 || next_node[1] != -1 || next_node[2] != -1)
			{
				while (final_node == -1)
				{
					srand(time(NULL));
					int r = rand() % 3;
					final_node = next_node[r];
				}
			}
			
			//Verifico se il circuito è chiuso
			if (final_node == -1)
			{
				tmp_succ[n] = start;
				sol_cost += dist(n, start, inst);
				done = true;
			}
			else //Aggiorno i valori
			{
				sol_cost += dist(n, final_node, inst);
				tmp_succ[n] = final_node;
				n = final_node;
			}
		}

		//Aggiorno la migliore soluzione
		if (sol_cost < best_sol)
		{
			for (int i = 0; i < inst->nnodes; i++)
				succ[i] = tmp_succ[i];
			best_sol = sol_cost;
		}

		/*for (int i = 0; i < inst->nnodes; i++)
			tmp_succ[i] = -1;
		*/
		inst->best_lb = sol_cost;
	}
	free(tmp_succ);
}


