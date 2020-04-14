#include "build_model_2.h"

#pragma warning(disable : 4996)

int xpos_compact(int i, int j, instance* inst);
int ypos_compact(int i, int j, instance* inst);
double dist(int i, int j, instance* inst);

void build_model2(instance* inst, CPXENVptr env, CPXLPptr lp) {
	char binary = 'B';
	char integer = 'I';
	char continuous = 'C';
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	*cname = (char*)calloc(100, sizeof(char));
	
	// aggiunta colonne
	// add binary var.s x(i,j)
	inst->x_start = CPXgetnumcols(env,lp);
	for (int i = 0; i < inst->nnodes; i++)		
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			double obj = dist(i, j, inst); // cost == distance  
			double lb = 0.0;
			double ub = (i == j) ? 0.0 : 1.0;

			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) print_error(" wrong position for x var.s "); 
		}
	}
	// add binary var.s y(i,j)
	inst->y_start = CPXgetnumcols(env, lp);
	for (int i = 0; i < inst->nnodes; i++)		
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
			double obj = 0.0; // not involved in the cost function computation
			double lb = 0.;
			double ub = (i == j || j == 0 ) ? 0.0 : inst->nnodes -1.;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) print_error(" wrong CPXnewcols on y var.s");
			if (CPXgetnumcols(env, lp) - 1 != ypos_compact(i, j, inst)) print_error(" wrong position for y var.s");
		}
	}
	
	//aggiunta righe
	// add the degree constraints
	for (int h = 0; h < inst->nnodes; h++)  // degree constraints
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = 'E';                              // 'E' for equality constraint 
		sprintf_s(cname[0], VERBOSE, "INdegree(%d)", h + 1);
		// inserisco il termine noto rhs (sommatoria di archi per ogni riga deve essere uguale a 1 in questo modello perche' usiamo grafo orientato)
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [in-degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta X(i,h)=1
			if (CPXchgcoef(env, lp, lastrow, xpos_compact(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [in-degree]");
		}		
	}
	for (int h = 0; h < inst->nnodes; h++)  // degree constraints
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = 'E';                              // 'E' for equality constraint 
		sprintf_s(cname[0], VERBOSE, "OUTdegree(%d)", h + 1);
		// inserisco il termine noto rhs (sommatoria di archi per ogni riga deve essere uguale a due nel nostro caso)
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [out-degree]");
		for (int j = 0; j < inst->nnodes; j++)
		{
			if (j == h) continue;
			// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta X(i,h)=1
			if (CPXchgcoef(env, lp, lastrow, xpos_compact(h, j, inst), 1.0)) print_error(" wrong CPXchgcoef [out-degree]");
		}
	}

	 // degree constraints sum(y_1j) = n - 1;
	int h = 0;
	int lastrow = CPXgetnumrows(env, lp);
	double rhs = inst->nnodes - 1.0;
	char sense = 'E';
	sprintf_s(cname[0], VERBOSE, "starting flow",h + 1);
	if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [starting flow]");
	for(int i=1 ; i<inst->nnodes;i++)
	{
		if (i == h) continue;
		// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta X(h,J)=1
		if (CPXchgcoef(env, lp, lastrow, ypos_compact(h,i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
	}
	
		// SUM( Y_ij) - SUM(Y_jk) = 1;

	for (int i = 1; i < inst->nnodes; i++) 
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = 'E';                              // 'E' for equality constraint 
		sprintf_s(cname[0], VERBOSE, "flow(%d)", i + 1);
		// inserisco il termine noto rhs (sommatoria di archi per ogni riga deve essere uguale a due nel nostro caso)
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree]");
		
			for (int j = 0; j < inst->nnodes; j++)  // excluding node 0 
			{
				if (i == j) continue;
				// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta X(h,J)=1
				if (CPXchgcoef(env, lp, lastrow, ypos_compact(i, j, inst), -1.0)) print_error(" wrong CPXchgcoef [degree]");
			
				// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta Y(i,j)=1
				if (CPXchgcoef(env, lp, lastrow, ypos_compact(j, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");

			}
		
	}
												//****LAZY CONTRAINTS****
	//Ex: GG formulation with directed - arc variables x_ij and x_ji-- > xpos_compact(i, j, inst)
	
	int izero = 0;
	int index[2];
	double value[2];
	
	// add lazy constraints  1.0 * n* x_ij - 1.0 * x_ij + 1 * y_ij <= 0 , for each arc (i,j) not touching node 0
	// credo che il commento sopra sia sbagliato, dovrebbe essere -(n-1)x_ij + y_ij <= 0 e va messo anche per il nodo zero
	// praticamente questo dice che ogni volta che y e' maggiore di zero va fissata la x a 1
	// si evita il nodo 0 solo quando si va a contare i flussi in entrata e uscita per permettere di chiudere il ciclo
	double n = inst->nnodes;
	rhs = 0;
	sense = 'L';
	int nnz = 2;

	for (int i = 0; i < inst->nnodes; i++) // excluding node 0 
	{
		for (int j = 0; j < inst->nnodes; j++) // excluding node 0 
		{
			if (i == j) continue;
			sprintf(cname[0], "lazyconstraints1(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);  // X_ij
			value[0] = 1. - n;
			index[1] = ypos_compact(i, j, inst);  // Y_ij
			value[1] = 1;
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}
	/*
	for (int i = 0; i < inst->nnodes; i++) // excluding node 0 
	{
		for (int j = 0; j < inst->nnodes; j++) // excluding node 0 
		{
			if (i == j) continue;
			sprintf(cname[0], "lazyconstraints1(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);  // X_ij
			value[0] = -1;
			index[1] = ypos_compact(i, j, inst);  // Y_ij
			value[1] = 1./(inst->nnodes);
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}*/

	// crea un file model.lp dove va inserire tutti valori ottimi tra i vari punti, degree di ogni punto, bounds.
	if (VERBOSE >= -100) CPXwriteprob(env, lp, "model_2.lp", NULL);

	free(cname[0]);
	free(cname);
}



int xpos_compact(int i, int j, instance* inst) 
{
	return inst->x_start + i * inst->nnodes + j;
}

int ypos_compact(int i, int j, instance* inst) 
{
	return inst->y_start + i * inst->nnodes + j;
}
