#include "build_model_1.h";
#include <math.h>

#pragma warning(disable : 4996)

//modello MTZ!
/*Attenzione, qui si usano grafi orientati, quindi  dovro' cambiare la xpos*/
void build_model_1(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	char binary = 'B';
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	*cname = (char*)calloc(100, sizeof(char));

	// aggiunta colonne: prima le x_ij
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if (i == j)
				continue;

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			double obj = dist(i, j, inst); // cost == distance  
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, inst)) print_error(" wrong position for x var.s");
		}
	}
	//aggiunta delle u_i: una variabile per ogni vertice a parte il primo
	char continuous = 'C';
	for (int i = 0; i < inst->nnodes-1; i++)
	{
		sprintf(cname[0], "u(%d)", i+2);
		double zero = 0.0; //cost = 0 per queste variabili
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, &zero, &zero, &ub, &continuous, cname)) print_error(" wrong CPXnewcols on u var.s");
		if (CPXgetnumcols(env, lp) - 1 != upos(i, inst)) print_error(" wrong position for u var.s");
	}

	//aggiunta righe: vincoli
	double rhs = 1.0;
	char sense = 'E';

	//sommatoria su i delle x_ih uguale a 1 per ogni h appartenente a V
	for (int h = 0; h < inst->nnodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);                  // 'E' for equality constraint 
		sprintf_s(cname[0], VERBOSE, "degree(%d)", h + 1);
		// inserisco il termine noto rhs (sommatoria di archi per ogni riga deve essere uguale a 1 in questo modello perche' usiamo grafo orientato)
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
		}
	}

	//sommatoria su i delle x_hi uguale a 1 per ogni h appartenente a V
	for (int h = 0; h < inst->nnodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		sprintf_s(cname[0], VERBOSE, "degree(%d)", h + 1 + inst->nnodes);
		// come sopra
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta
			if (CPXchgcoef(env, lp, lastrow, xxpos(h, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
		}
	}

	//u_j >= u_i + 1 -M(1-x_ij) lo traduco come u_j-u_i -M*x_ij >= 1 - M, quindi il termine noto diventa 1-(inst->nnodes-1) = 2 - inst->nnodes
	// tolgo la roba che ho scritto sopra e provo a mettere come equazione u_i - u_j + n x_ij <= n-1
	rhs = inst -> nnodes - 1.;
	sense = 'L';
	//vincoli sulle u
	for (int i = 1; i < inst->nnodes; i++)
		for (int j = 1; j < inst->nnodes; j++)
		{
			if (i == j) continue;


			int lastrow = CPXgetnumrows(env, lp);
			sprintf(cname[0], "vincolo u %d, %d", i, j);

			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [u]\n");

			/*if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), -inst->nnodes + 1.0)) print_error(" wrong CPXchgcoef [u1]");
			if (CPXchgcoef(env, lp, lastrow, upos(i - 1, inst), -1.0)) print_error(" wrong CPXchgcoef [u2]");
			if (CPXchgcoef(env, lp, lastrow, upos(j - 1, inst), 1.0)) print_error(" wrong CPXchgcoef [u3]");
			*/
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), inst->nnodes)) print_error(" wrong CPXchgcoef [u1]");
			if (CPXchgcoef(env, lp, lastrow, upos(i - 1, inst), 1.0)) print_error(" wrong CPXchgcoef [u2]");
			if (CPXchgcoef(env, lp, lastrow, upos(j - 1, inst), -1.0)) print_error(" wrong CPXchgcoef [u3]");

			// faccio -1 perche' u_0 non l'ho inserita
		}

	if (VERBOSE >= -100) CPXwriteprob(env, lp, "model_1.lp", NULL);

	free(cname[0]);
	free(cname);
}

int xxpos(int i, int j, instance* inst)
{
	if (i == j)
		print_error("Errore nella xxpos: i=j\n");
	int x = i * inst->nnodes + j;
	return (i < j) ? x - (i + 1) : x - i;
}

int upos(int i, instance* inst)
{
	return pow(inst->nnodes,2) - inst->nnodes + i;
}