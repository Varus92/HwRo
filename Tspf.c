#include "Tspf.h"
#include "build_model_2.h"
#include <gnuplot_c.h>
#include <string.h>

//eps = epsilon
#define EPS 0.2

#pragma warning(disable : 4996)

void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
void plot_solution(ComplexRect_s Array[], int size, enum gpcNewAddGraphMode x, h_GPC_Plot*);
int build_solution(const double* xstar, instance* inst, int* succ, int* comp);
int* get_component_array(int component, int succ[], int comp[], int size, int* c_component_dim);

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
	int error = 0;
	int firstLine = 0;

	CPXENVptr env = CPXopenCPLEX(&error);   //crea l'ambiente per cplex
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
	
	printf("\n \n \n Build model= %d\n \n \n ", inst->model_type);
	
	if (inst->model_type == 0) build_model(inst, env, lp);
	if (inst->model_type == 1) build_model2(inst, env, lp);
	

	//setta parametri per trovare la soluzione ottima
	if (CPXmipopt(env, lp) != 0) print_error("Errore nella chiamata di mipopt");

	h_GPC_Plot* plot;
	plot = gpc_init_xy("TSP solution",
		"X",
		"Y",
		GPC_AUTO_SCALE,
		GPC_KEY_ENABLE);

	int Ncols = CPXgetnumcols(env, lp);
	double* sol = (double*)calloc(Ncols, sizeof(double));
	CPXgetx(env, lp, sol, 0, Ncols - 1);

	printf("\nprinting solution\n");

	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int Ncomp = build_solution(sol, inst, succ, comp);

	for (int i = 0; i < inst->nnodes; i++)
		printf("nodo %d nella componente connessa %d e successivo %d\n", i+1, *(comp+i), *(succ+i)+1);

	printf("numero componenti connesse: %d\n", Ncomp);


	//conto elementi nelle varie componenti connesse per debuggare
	int count = 0;
	for (int i = 0; i < Ncomp; i++)
	{
		int dim;
		int* Array = get_component_array(i + 1, succ, comp, inst->nnodes, &dim);
		count += dim;
		ComplexRect_s* coords = (ComplexRect_s*)calloc(dim + 1, sizeof(ComplexRect_s));
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
		plot_solution(coords, dim + 1, i == 0 ? GPC_NEW : GPC_NEW, plot);

		getchar();

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

void plot_solution(ComplexRect_s Array[], int size, enum gpcNewAddGraphMode x, h_GPC_Plot* plot)
{
	gpc_plot_xy(plot,
		Array,
		size,
		"X/Y Plot",
		"linespoints pt7 ps 2",
		"blue",
		x);
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

	// crea un file model.lp dove va inserire tutti valori ottimi tra i vari punti, degree di ogni punto, bounds.
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

//#define DEBUG 

int build_solution(const double* xstar, instance* inst, int* succ, int* comp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{

#ifdef DEBUG
	int* degree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			int k = xpos(i, j, inst);
			if (fabs(xstar[k]) > EPS&& fabs(xstar[k] - 1.0) > EPS) print_error(" wrong xstar in build_sol()");
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
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1)
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

int* get_component_array(int component, int succ[], int comp[], int size, int* c_component_dim) {
	int dim = 0;
	int* array = (int*)calloc(size, sizeof(int));

	//find first node
	int i;
	for (i = 0; i < size && comp[i] != component; i++);

	if (i == size)	return;		//se viene scandita tutta la soluzione senza trovare la componente connessa desiderata, devo ritornare altrimenti va in loop

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