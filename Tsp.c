#include "Tsp.h"
#include <gnuplot_c.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "build_model_1.h"
#include "build_model_2.h"
#include "Mat_euristics.h"
#include "Grasp.h"

//eps = epsilon
#define EPS 0.2

#pragma warning(disable : 4996)

void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
int build_solution(const double* xstar, instance* inst, int* succ, int* comp);
int* get_component_array(int component, int succ[], int comp[], int size, int* c_component_dim);
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* data);
static int add_secs(instance* inst, int* xstar, CPXCALLBACKCONTEXTptr context);

// required by time.h
double second();


double dist(int i, int j, instance* inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];

	if (!inst->integer_costs)
		return sqrt(dx * dx + dy * dy);

	int distance = sqrt(dx * dx + dy * dy + 0.499999); //l'intero piu vicino

	return distance + 0.0;
}


int TSPopt(instance* inst) {
	//inst->tstart = second();


	int error = 0;
	int firstLine = 0;
	inst->model_type = 3;

	CPXENVptr env = CPXopenCPLEX(&error);   //crea l'ambiente per cplex
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	switch (inst->model_type)
	{
	case 0:
		build_model(inst, env, lp);
		break;

	case 1:
		build_model_1(inst, env, lp);
		break;
	case 2:
		build_model2(inst, env, lp);
		break;
	case 3:
		heuristicSolver(inst);
	}

	printf("\nCalcolo soluzione ottima\n");

	//setta parametri per trovare la soluzione ottima

		//setto callback
	int callback = 0;
	if (callback == 1)
	{
		// cplex output
		//CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);						// MIP node display info
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);

		inst->ncols = CPXgetnumcols(env, lp);
		
		if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, my_callback, inst)) print_error("errore nel set della callback");

		// riattivo il multithreading che cplex di default disattiva quando l'utente setta una callback sua
		int ncores = 1;
		CPXgetnumcores(env, &ncores);
		CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
		
	}
	
	int use_euristic = 0; // 1 hard fixing, 0 niente, 2 local branching
	double temp_objval = CPX_INFBOUND, tstart = second();
	if (use_euristic == 1)
	{
		//launch cplex with a time limit and get incumbent

		CPXsetdblparam(env, CPX_PARAM_TILIM, 2); 							// real time
		//CPXsetdblparam(env, CPX_PARAM_DETTILIM, TICKS_PER_SECOND * timelimit);

		if (CPXmipopt(env, lp) != 0) print_error("Errore nella chiamata di mipopt");

		printf("secondi passati: %f\n", second() - tstart);
		printf("inizio %f fine %f\n", second(), tstart);

		double percentages[] = { 0.9, 0.75, 0.5, 0.15 };
		for (int i = 0; i < 4; i++)
		{
			int Ncols = CPXgetnumcols(env, lp);
			double* sol = (double*)calloc(Ncols, sizeof(double));
			CPXgetx(env, lp, sol, 0, Ncols - 1);

			int* succ = (int*)calloc(inst->nnodes, sizeof(int));
			int* comp = (int*)calloc(inst->nnodes, sizeof(int));
			int Ncomp = build_solution(sol, inst, succ, comp);

			printf("number of connected components: %d\n", Ncomp);

			printf("debug 1\n");
			int* circuit_partial = get_circuit(inst, comp, succ);

			printf("debug 2\n");
			fix(circuit_partial, inst->nnodes, percentages[i], inst, env, lp);

			printf("debug 3\n");
			CPXsetdblparam(env, CPX_PARAM_TILIM, 2);
			if (CPXmipopt(env, lp) != 0) print_error("Errore nella chiamata di mipopt");

			free(circuit_partial);
		}
	}
	else
	{
		//risolve
		if (CPXmipopt(env, lp) != 0) print_error("Errore nella chiamata di mipopt");
	}
	int Ncols = CPXgetnumcols(env, lp);
	double* sol = (double*)calloc(Ncols, sizeof(double));
	CPXgetx(env, lp, sol, 0, Ncols - 1);

	//printf("\n printing solution\n");

	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int Ncomp = build_solution(sol, inst, succ, comp);
	
	printf("soluzione trovata\n");


	h_GPC_Plot* plot;
	plot = gpc_init_xy("TSP solution",
		"X",
		"Y",
		GPC_AUTO_SCALE,
		GPC_KEY_DISABLE);

	
	printf("%d\n", Ncomp);
	if (Ncomp == 1)
		printf("Soluzione trovata, non e' necessario procedere con il LOOP\n");
	else
		if (inst->model_type == 0) {
			printf("Metodo Loop");
			Ncomp = benders_method(sol, succ, comp, Ncomp, inst, env, lp);
			if (Ncomp != 1) print_error("Errore nel benders method\n");
		}

	printf("numero componenti connesse: %d\n", Ncomp);
	
	//preparo grafico
	ComplexRect_s* points = (ComplexRect_s*)calloc(inst->nnodes, sizeof(ComplexRect_s));

	for (int i = 0; i < inst->nnodes; i++)
	{
		points[i].real = inst->xcoord[i];
		points[i].imag = inst->ycoord[i];
	}

	gpc_plot_xy(plot,
		points,
		inst->nnodes,
		"X/Y Plot",
		"points pt 2",
		"red",
		GPC_NEW);

	free(points);

	//mentre plotto conto elementi nelle varie componenti connesse per debuggare
	int count = 0;
	for (int i = 0; i < Ncomp; i++)
	{
		int dim = 0;
		int* Array = get_component_array(i + 1, succ, comp, inst->nnodes, &dim);
		count += dim;
		ComplexRect_s* coords = (ComplexRect_s*)calloc(dim + 1.0, sizeof(ComplexRect_s));
		for (int j = 0; j < dim; j++)
		{
			// *(coords + j) = inst->(xcoord+Array[j]), inst->(ycoord + Array[j]) };
			coords[j].real = inst->xcoord[Array[j]];
			coords[j].imag = inst->ycoord[Array[j]];
		}
		//inserisco arco per chiudere il ciclo
		coords[dim].real = inst->xcoord[Array[0]];
		coords[dim].imag = inst->ycoord[Array[0]];

		//plot_solution(coords, dim+1, i==0 ? GPC_NEW : GPC_ADD, plot);
		gpc_plot_xy(plot,
			coords,
			dim + 1,
			"X/Y Plot",
			"lines",
			"blue",
			GPC_ADD);

		//getchar();

		free(coords);
		free(Array);
	}

	printf("Soluzione plottata!\nnumero nodi: %d, di cui contati %d\n", inst->nnodes, count);
	getchar();
	gpc_close(plot);

	//effettua i free per evitare memory leak, alcuni sono richiesti da cplex
	free(sol);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	return error;

}

void print_error(const char* err)
{
	printf("\n debug: %s \n", err);
	fflush(NULL);
	exit(1);
}

void build_model(instance* inst, CPXENVptr env, CPXLPptr lp) {

	char binary = 'B';
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	*cname = (char*)calloc(100, sizeof(char));

	// aggiunta colonne
	// add binary var.s x(i,j) for i < j  
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);

			double obj = dist(i, j, inst); // cost == distance 

			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) print_error(" wrong position for x var.s");
		}
	}

	//aggiunta righe
	// add the degree constraints
	for (int h = 0; h < inst->nnodes; h++)  // degree constraints
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf_s(cname[0], VERBOSE, "degree(%d)", h + 1);
		// inserisco il termine noto rhs (sommatoria di archi per ogni riga deve essere uguale a due nel nostro caso)
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			// vado a dirgli quali elementi concorrono nella sommatoria, nel nostro caso tutta la riga appena aggiunta
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
		}
	}

	if (VERBOSE >= -100) CPXwriteprob(env, lp, "model.lp", NULL);

	free(cname[0]);
	free(cname);
}

// il nostro modello vuole una variabile per ogni arco, e noi le organizziamo in una matrice, come fatto in ricerca operativa 1
// la diagonale principale è vuota ovviamente perché non ci sono self loop, e ci basta tenere una sola delle due "metà" perché
// trattiamo il TSP nella versione con grafi NON orientati, prendiamo quindi la triangolare superiore (tuttavia xpos restituisce
// il numero dell'arco anche se chiediamo l'indice dell'arco nella posizione simmetrica sfruttando appunto la simmetria del problema)
int xpos(int i, int j, instance* inst) {
	if (i == j) print_error(" i == j in xpos");
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}
// ricapitolando: X_ij nel nostro modello viene codificata in cplex come X_{xpos(i,j, inst)}

void* get_xpos(int model_type)
{
	int (*pos)(int, int, instance*);
	switch (model_type)
	{
	case 0: pos = xpos;
		break;
	case 1: pos = xxpos;
		break;
	case 2: pos = xpos_compact;
		break;
	default: pos = xpos;
	}
	return pos;
}

//#define DEBUG 
int build_solution(const double* xstar, instance* inst, int* succ, int* comp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{
	//creo un puntatore a funzione per sapere quale xpos utilizzare in base al model type
	int (*pos)(int, int, instance*);
	/*switch (inst->model_type)
	{
	case 0: pos = xpos;
		break;
	case 1: pos = xxpos;
		break;
	case 2: pos = xpos_compact;
		break;
	default: pos = xpos;
	}*/
	pos = get_xpos(inst->model_type);

#ifdef DEBUG
	int* degree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			int k = pos(i, j, inst);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) print_error(" wrong xstar in build_sol()");
			if (xstar[k] > 0.5)
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		char s[100];
		sprintf(s, "wrong degree in build_sol(), error at node %d", i + 1);
		if (degree[i] != 2) print_error(s);
	}
	free(degree);
#endif

	int ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		ncomp++;
		int i = start;

		while (succ[i] != start)  // go and visit the current component
		{
			comp[i] = ncomp;
			int j;
			for (j = 0; j < inst->nnodes; j++)
			{
				if (i != j && xstar[pos(i, j, inst)] > 0.5 && comp[j] == -1)
				{
					succ[i] = j;
					i = j;
					break;
				}
			}
			if (j == inst->nnodes) // if non ho trovato niente vuol dire che sono arrivato all'inizio
			{
				succ[i] = start;  // last arc to close the cycle
			}
		}
	}

	return ncomp;
}

//funzione che data una componente connessa component restituisce un array contenenti il subtour di quella componente
int* get_component_array(int component, int succ[], int comp[], int size, int* c_component_dim) {
	int dim = 0;
	int* array = (int*)calloc(size, sizeof(int));

	//find first node
	int i;
	for (i = 0; i < size && comp[i] != component; i++);

	if (i == size)	return NULL;		//se viene scandita tutta la soluzione senza trovare la componente connessa desiderata, devo ritornare altrimenti va in loop

	int prev = -1, start = i;	//ho dovuto usare una variabile in più altrimenti i cicli di 2 nodi mi fregano
	while (succ[prev] != start && dim < size)
	{
		prev = i;
		*(array + dim) = prev;
		dim++;
		i = succ[i];
	}

	*c_component_dim = dim;
	return array;
}

int benders_method(double* xstar, int* succ, int* comp, int ncomp, instance* inst, CPXENVptr env, CPXLPptr lp)
{
	char sense = 'L';
	char** cname = (char**)malloc(sizeof(char*));
	cname[0] = (char*)calloc(100, sizeof(char));
	while (ncomp > 1)
	{
		//printf("Sono vivo! Numero componenti connesse: %d\n", ncomp);
		for (int i = 0; i < ncomp; i++)
		{
			sprintf(cname[0], "SEC %d", i + 1);
			// calcolo la dimensione del termine noto del vincolo
			double rhs = 0.;
			for (int j = 0; j < inst->nnodes; j++)
				if (comp[j] == i + 1)
				{
					rhs++;
				}
			rhs--;

			//aggiungo SEC sulla componente connessa s_i
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [u]\n");
			int lastrow = CPXgetnumrows(env, lp) - 1;

			int dim;
			int* comp_i = get_component_array(i + 1, succ, comp, inst->nnodes, &dim);
			if (comp_i == NULL)	return -1;
			// va fatto per ogni coppia appartenente alla componente connessa S_i
			for (int j = 0; j < dim; j++)
				for (int jj = j + 1; jj < dim; jj++)
				{
					if (CPXchgcoef(env, lp, lastrow, xpos(comp_i[j], comp_i[jj], inst), 1.0)) print_error(" wrong CPXchgcoef [SEC]");
				}

			free(comp_i);
		}
		//ricalcolo la soluzione ottima con i nuovi vincoli e chiamo di nuovo build solution
		CPXwriteprob(env, lp, "model.lp", NULL);

		if (CPXmipopt(env, lp) != 0) print_error("Errore nella chiamata di mipopt");
		CPXgetx(env, lp, xstar, 0, CPXgetnumcols(env, lp) - 1);
		ncomp = build_solution(xstar, inst, succ, comp);
	}
	free(cname[0]);
	free(cname);
	return ncomp;
}

int time_limit_expired(instance* inst)
{
	double tspan = second() - inst->tstart;
	if (tspan > inst->timelimit)
	{
		if (VERBOSE >= 100) printf("\n\n$$$ time limit of %10.1lf sec.s expired after %10.1lf sec.s $$$\n\n", inst->timelimit, tspan);
		//exit(0); 
		return 1;
	}
	return 0;
}

static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* data)
{
	instance* inst = (instance*)data;

	// get solution xstar
	double* xstar = (double*)calloc(inst->ncols, sizeof(double));

	// get some random information at the node (as an example)
	double objval = CPX_INFBOUND + 0.0;

	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) return 1; // xstar = current x from CPLEX-- xstar starts from position 0

	int mythread = -1;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADS, &mythread);
	double zbest;
	CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);

	//ora aggiungo i SEC su questa xstar
	int ncuts = add_secs(inst, xstar, context);

	free(xstar);
	return 0;
}


static int add_secs(instance* inst, int* xstar, CPXCALLBACKCONTEXTptr context)
{

	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = build_solution(xstar, inst, succ, comp);

	for (int i = 0; ncomp > 1 && i < ncomp; i++)
	{
		int count_nodes = 0;

		int* comp_nodes;

		// calcolo la dimensione del termine noto del vincolo
		comp_nodes = get_component_array(i + 1, succ, comp, inst->nnodes, &count_nodes);

		int zero = 0;

		int n_edges = (count_nodes * (count_nodes - 1)) / 2;
		int* indexes = (int*)calloc(n_edges, sizeof(int));
		double* values = (double*)calloc(n_edges, sizeof(double));
		int start_indexes = 0;

		double rhs = count_nodes - 1.0;
		char sense = 'L';
		int nnz = 0; //number of non-zero coefficients

		// ottengo la componente connessa S_i, costruisco il vettore di indici degli elementi in S_i e attribuisco valore 1

		for (int j = 0; j < count_nodes - 1; j++)
			for (int jj = j + 1; jj < count_nodes; jj++)
			{
				indexes[nnz] = xpos(comp_nodes[j], comp_nodes[jj], inst);
				values[nnz] = 1.0;
				nnz++;
			}

		if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &start_indexes, indexes, values)) print_error("Errore nell'aggiunta dei sec nella callback");

		free(indexes);
		free(values);
		free(comp_nodes);
	}

	free(comp);
	free(succ);

	return ncomp == 1 ? 0 : ncomp;
}