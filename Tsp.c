#include "Tsp.h"
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
int xpos(int i, int j, instance* inst) 
{
	if(i==j)
	{
		print_error("i==j in posx");
	}
	if(i>j)
	{
		return xpos(j, i, inst);
	}

	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	

	return pos; 
}

double dist(int i, int j, instance *inst)
{
    double dx = inst->xcoord[i] - inst->ycoord[j];
    double dy = inst->ycoord[i] - inst->ycoord[j];
    
    if( !inst->integer_costs)
        return sqrt(dx*dx+dy*dy);
    
    int distance = sqrt(dx*dx+dy*dy + 0.499999); //l'intero piu vicino

    return distance +0.0;
}

int TSPopt(instance* inst)
{
	printf("inizio initialization : \n");
	/* 1. initialization ------------------------------------------------- */

	// open cplex model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP"); //creato le strutture env & lp
	

	/* 2. build initial model  ------------------------------------------------- */
	
	printf("Inizion build initial model\n");
	build_model(inst, env, lp);
	
	//Cplex parameter settings
	
	CPXmipopt(env, lp);
	
	/* 3. final MIP run ------------------------------------------------------ */
	
	// free pools and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0;
}

void build_model(instance* inst, CPXENVptr env, CPXLPptr lp)   // basic model with asymmetric x and q var.s
/**************************************************************************************************************************/
{
	
	double zero = 0.0; // one = 1.0; 	
	char binary = 'B';
	
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i<j
	 
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i+1; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			double obj = dist(i, j, inst);			//costo
			double ub = 1.0;
			double lb = 0.0;

			if (CPXnewcols(env, lp, 1, &obj, &zero, &ub, &binary, cname)) 
			{
				print_error(" wrong CPXnewcols on x var.s");
			}

			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
			{
				print_error(" wrong position for x var.s");
			}
		}
	}
	

	for (int h = 0; h < inst->nnodes; h++)		// add the degree
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 2.0;
		char sense = 'E';
		sprintf(cname[0], "degree(%d)", h + 1);

		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
		{
			print_error(" wrong CPXnewrows [degree]");
		}
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			if (CPXchgcoef(env, lp, lastrow, xpos(h, i, inst), 1.0))
			{
				print_error(" wrong CPXchgcoef [degree]");
			}
		}
	}

	free(cname[0]);
	free(cname);

}
