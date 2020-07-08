#include "Heuristics.h"
#include "Tsp.h"
#include <time.h>

#define EPSILON_COST 0.1

double second();
int* v_duplicate(int* v, int dim);

void print_vector(int* vec, int dim)
{
	for (int i = 0; i < dim; i++)
		printf("%d ", vec[i]);
	printf("\n");
}

double calc_delta(instance* inst, int a, int b, int c, int d)
{
	return dist(a, b, inst) + dist(c, d, inst) - (dist(a, c, inst) + dist(b, d, inst));
}

void posmax(double** mat, int rows, int cols, int* max_row, int* max_col)
{
	*max_row = 0; *max_col = 1;
	for (int i = 0; i < rows; i++)
		for (int j = i+1; j < cols; j++)
			if (i != j && mat[i][j] > mat[*max_row][*max_col])
			{
				*max_row = i;
				*max_col = j;
			}
}

// assumo che gli archi correnti siano a-b e c-d, a-d e' l'arco che non posso creare per non dividere il ciclo in due subtour
int* swap_2opt(instance* inst, int* circuit_nodes, int node_a, int node_b, int node_c, int node_d)
{
	int (*pos)(int, int, instance*);
	pos = get_xpos(inst->model_type);

	int* new_circuit = (int*)calloc(inst->nnodes, sizeof(int));

	double delta = calc_delta(inst, node_a, node_b, node_c, node_d);
	if (delta > 0) //swap
	{
		int i;
		new_circuit[0] = circuit_nodes[0];
		for (i = 1; i < inst->nnodes && new_circuit[i - 1] != node_a; i++)
			new_circuit[i] = circuit_nodes[i];
		int j = i; // index for position in the new circuit

		//find position of node C
		while (circuit_nodes[i] != node_c)
			i++;

		//new circuit = {every node in order from C to B}
		while (circuit_nodes[i] != node_b)
		{
			new_circuit[j] = circuit_nodes[i];
			i--;
			j++;
		}

		//add edge b-d
		new_circuit[j++] = node_b;
		// node b will be added to the new circuit up next

		//locate D and append each node in old circuit to the new one
		for (; circuit_nodes[i] != node_d; i++);
		for (; i < inst->nnodes; i++)
		{
			new_circuit[j] = circuit_nodes[i];
			j++;
		}
	}
	else  // Non scambiare, ma bisogna comunque copiare l'array per evitare incongruenze
	{
		for (int i = 0; i < inst->nnodes; i++)
			new_circuit[i] = circuit_nodes[i];
	}
	return new_circuit;
}

int* neighborhood_search(instance* inst, int* partial_sol, int using_cplex_xstar)
{
	int convergence = 0;
	int* comp = NULL, * succ = NULL, * cycle = NULL;
	double** delta = (double**)calloc(inst->nnodes, sizeof(double*));
	comp = (int*)calloc(inst->nnodes, sizeof(int));
	succ = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		delta[i] = calloc(inst->nnodes, sizeof(double));
	}

	if (partial_sol != NULL && using_cplex_xstar == 1)
	{
		int dim;
		build_solution(partial_sol, inst, succ, comp);
		cycle = get_component_array(1, succ, comp, inst->nnodes, &dim);
		if (dim < inst->nnodes - 1)
		{
			printf("%d\n", dim);
			print_error("Attenzione, il ciclo considerato non e' hamiltoninano!\n");
		}
	}
	else
		if (partial_sol == NULL)
		{
			cycle = (int*)calloc(inst->nnodes, sizeof(int));
			// inizializzo il ciclo a qualcosa di sufficientemente stupido: 0->1->2->3->...->nnodes-1, potrei anche prendere una sequenza random in effetti
			for (int i = 0; i < inst->nnodes; i++)
				cycle[i] = i;
		}
		else
			cycle = partial_sol;

	while (convergence != 1)
	{
		//for each edge in the cycle
			//for each other edge in the cycle
				// get DELTA
		// swap the couple of edges with max delta
		for (int i = 0; i < inst->nnodes; i++)
			for (int j = i + 1; j < inst->nnodes; j++)
				delta[i][j] = calc_delta(inst, cycle[i], cycle[(i + 1) % inst->nnodes], cycle[j], cycle[(j + 1) % inst->nnodes]);

		//printf("delta matrix\n");
		//for (int i = 0; i < inst->nnodes; i++)
		//	for (int j = i+1; j < inst->nnodes; j++)
		//		printf("%f ", delta[i][j]);
		//printf("\n");
		//getchar();

		int max_i, max_j;
		posmax(delta, inst->nnodes, inst->nnodes, &max_i, &max_j);

		printf("posizioni massimo delta: %d %d, valore: %1f\nCiclo prima dello scambio:\n", max_i, max_j, delta[max_i][max_j]);
		for (int i = 0; i < inst->nnodes; i++)
			printf("%d ", cycle[i]);
		printf("\n");

		if (delta[max_i][max_j] > EPSILON_COST)
		{
			int* new_cycle = swap_2opt(inst, cycle, cycle[max_i], cycle[(max_i + 1) % inst->nnodes], cycle[max_j], cycle[(max_j + 1) % inst->nnodes]);
			
			for (int i = 0; i < inst->nnodes; i++)
				cycle[i] = new_cycle[i];
			//free(new_cycle);
		}
		else
		{
			convergence = 1;
		}
	}

	printf("Ciclo dopo la convergenza:\n");
	for (int i = 0; i < inst->nnodes; i++)
		printf("%d ", cycle[i]);
	printf("\n");

	// free di tutte le cose
	free(comp);
	free(succ);
	for (int i = 0; i < inst->nnodes; i++)
		free(delta[i]);
	free(delta);
	//if(partial_sol != NULL)
	//	free(partial_sol);
	return cycle;
}

int* v_duplicate(int* v, int dim)
{
	int* new_v = (int*)calloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++)
		new_v[i] = v[i];
	return new_v;
}

int v_search(int* v, int dim, const int value)
{
	int pos = -1;
	for (int i = 0; i < dim; i++)
		if (v[i] == value)
		{
			pos = i;
			break;
		}
	return pos;
}

int insertion_sort(int val, int* v, int dim)	// versione limitata, solo per quello che serve a me, non robusto in tutti i casi
{
	int i;
	for (i = 0; v[i] >= 0 && i < dim && v[i] < val; i++);
	if (i == dim)	//array gia' riempito
		return 0;
	if (v[i] == val)	//duplicate value
		return 1;
	for (int j = dim-1; j > i; j--)
		v[j] = v[j - 1];
	v[i] = val;
	return 0;
}

int* randperm(int dim)
{
	int* free_values = (int*)calloc(dim, sizeof(int));
	int* ret = (int*)calloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++)
		free_values[i] = i;

	for (int i = 0; i < dim; i++)
	{
		int x = rand() % (dim - i);
		ret[i] = free_values[x+i];
		//sposto il valore scelto nella posizione i cosi' non posso sceglierlo di nuovo
		int temp = free_values[i];
		free_values[i] = free_values[x + i];
		free_values[x + i] = temp;
	}

	free(free_values);
	return ret;
}

int* increment_pointer_mod_dim(int* ptr, int* start_ptr, int dim, int n)
{
	ptr = ptr + n;
	if (ptr < start_ptr)
		return start_ptr + dim - 1;
	if (ptr - start_ptr >= dim)
		return start_ptr;
	return ptr;
}

//nuovi archi sono del tipo: possible_destinations[i] <-> possible_destinations[new_dest[i]] solo se new_dest[i]!=-1
int* build_tour(int* cycle, int* new_dest, int cycle_dim, int* possible_destinations, int* subtour)
{
	int* fixed_cycle = (int*)calloc(cycle_dim, sizeof(int));
	
	int *w_node = &cycle[0], clockwise = 1;


	for (int current_node = 0; current_node < cycle_dim; current_node++)
	{
		//printf("We want to insert %d in the following vector\n", *w_node);
		//print_vector(fixed_cycle, cycle_dim);
		//getchar();
		int y = v_search(fixed_cycle, current_node, *w_node);
		if (/*current_node > 0*/ w_node != &cycle[0] && y != -1)
		{
			printf("\nSubtour identificato! ai nodi %d e %d\n", current_node, y);
			*subtour = 1;
			free(fixed_cycle);
			return NULL;
		}
		int pos = v_search(possible_destinations, 10, *w_node);
		if (pos == -1)
		{
			fixed_cycle[current_node] = *w_node;
			if (clockwise == 1)
			{
				w_node = increment_pointer_mod_dim(w_node, &cycle[0], cycle_dim, 1);
			}
			else
			{
				w_node = increment_pointer_mod_dim(w_node, &cycle[0], cycle_dim, -1);
			}
			continue;
		}

		//print_vector(fixed_cycle, cycle_dim);
		//getchar();

		// nodo critico, aggiungo lui al ciclo, aggiungo il successivo e sistemo w_node
		fixed_cycle[current_node] = *w_node;
		current_node++;

		int x;
		if ((x = v_search(fixed_cycle, current_node, possible_destinations[new_dest[pos]])) != -1)
		{
			printf("\nSubtour identificato! al nodo %d, posizione %d\n", possible_destinations[new_dest[pos]], x);
			*subtour = 1;
			return NULL;
		}

		fixed_cycle[current_node] = possible_destinations[new_dest[pos]];
		//printf("ho inserito %d\n", fixed_cycle[current_node]);
		int new_index = v_search(cycle, cycle_dim, fixed_cycle[current_node]);
		//if (new_index == -1)
		//{
		//	print_vector(cycle, cycle_dim);
		//	exit(0);
		//}
		pos = v_search(possible_destinations, 10, fixed_cycle[current_node]);
		if (pos < 5)
		{
			clockwise = 0;
			new_index = (new_index - 1) % cycle_dim;
			printf("%d\n", new_index);
		}	
		else
		{
			//printf("sto cercando %d ed e' alla posizione %d\n",fixed_cycle[current_node], pos);
			clockwise = 1;
			new_index = (new_index + 1) % cycle_dim;
		}
		w_node = &cycle[new_index];
		//print_vector(fixed_cycle, cycle_dim);
		//getchar();
	}

	//printf("Building tour:\n\n");
	//for (int i = 0; i < cycle_dim; i++)
	//	printf("%d ", fixed_cycle[i]);
	return fixed_cycle;
}

int variable_neighborhood_search(instance* inst, double time_limit)
{
	double best_obj_val = CPX_INFBOUND;
	int* partial_sol = NULL, best_sol = NULL;
	while (time_limit > 0)
	{
		double begin_time = second();
		partial_sol = neighborhood_search(inst, partial_sol, 0);

		for (int i = 0; i < inst->nnodes; i++)
			printf("%d ", partial_sol[i]);
		printf("\n");

		//questa parte e' lenta, dovrei salvare la soluzione precedente e calcolare il nuovo costo aggiungendo delta
		if (get_obj_val(inst, partial_sol) < best_obj_val)
		{
			best_obj_val = get_obj_val(inst, partial_sol);
			if(best_sol != NULL) free(best_sol);
			best_sol = v_duplicate(partial_sol, inst->nnodes);
		}

		// togli 5 archi a caso e riaggiusta a caso
		int critical_indexes[5];	//ogni arco da togliere viene identificato dal nodo iniziale e dal nodo successivo nel ciclo (cycle)
		for (int i = 0; i < 5; i++)
		{
			critical_indexes[i] = -1;	//necessario per insertion_sort()

			int done = 0;
			while (!done)
			{
				int x = rand() % inst->nnodes;
				if (x == inst->nnodes - 1 || x == 0)	//per non complicare inutilmente il codice successivo, richiedo che l'arco di ricongiungimento non venga toccato
					continue;
				if (v_search(critical_indexes, i + 1, x - 1) != -1)
					continue;
				if (v_search(critical_indexes, i + 1, x + 1) != -1)
					continue;
				done = !insertion_sort(x, critical_indexes, 5);
			}
		}
		int real_dim = 10;
		int possible_destinations[10];
		for (int i = 0; i < 5; i++)
		{
			possible_destinations[i] = partial_sol[critical_indexes[i]];
			possible_destinations[5 + i] = partial_sol[(critical_indexes[i]+1) % inst->nnodes]; //successivo di quello appena definito sopra
		}
		int dest[10];
		int subtour = 1;
		while (subtour == 1)
		{
			int* perm = randperm(10);
			subtour = 0;

			// inizializzo il vettore dest
			for (int i = 0; i < 10; i++)
				dest[i] = -1;

			for (int i = 0; i < 9; i++)
			{
				//if (dest[perm[i]] != -1)
				{
					dest[perm[i]] = perm[i + 1];
					dest[perm[i + 1]] = perm[i];
					i++;
				}
			}

			int* temp = build_tour(partial_sol, dest, inst->nnodes, possible_destinations, &subtour);

			if (subtour != 0)
			{
				printf("Subtour! %d\n", subtour);
				//getchar();
			}


			if (subtour == 0)
			{
				partial_sol = temp;
			}
			else
				if (temp != NULL)
					free(temp);
			free(perm);
		}

		int error = 0;
		printf("\n\n calcio all'ottimo locale, nuova soluzione:\n");
		for (int i = 0; i < inst->nnodes; i++)
		{
			printf("%d ", partial_sol[i]);
			if (partial_sol[i] > 47 || partial_sol[i] < 0)
				error = 1;
		}
		printf("\n\n\n\n\n");
		if (error > 0)
		{
			for (int i = 0; i < 10; i++)
				printf("%d ", dest[i]);
			printf("\n");

			for (int i = 0; i < 10; i++)
				printf("%d ", possible_destinations[i]);
			exit(1);
		}

		time_limit = time_limit - (second() - begin_time);
	}
}

double get_obj_val(instance* inst, int* nodes_cycle)
{
	int obj_val = 0;
	for (int i = 0; i < inst->nnodes; i++)
		obj_val += dist(i, (i + 1) % inst->nnodes, inst);
	return obj_val;
}

//#define debug

#ifdef debug
#define dim 30

void main()
{
	int *partial_sol = (int*)calloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++)
		partial_sol[i] = i;

	int n_scambi = 4;
	for (int k = 0; k < n_scambi; k++)
	{
		int dest[10];
		int critical_indexes[5];	//ogni arco da togliere viene identificato dal nodo iniziale e dal nodo successivo nel ciclo (cycle)
		for (int i = 0; i < 5; i++)
		{
			critical_indexes[i] = -1;	//necessario per insertion_sort()

			int done = 0;
			while (!done)
			{
				int x = rand() % dim;
				if (x == dim - 1 || x == 0)	//per non complicare inutilmente il codice successivo, richiedo che l'arco di ricongiungimento non venga toccato
					continue;
				if (v_search(critical_indexes, i + 1, x - 1) != -1)
					continue;
				if (v_search(critical_indexes, i + 1, x + 1) != -1)
					continue;
				done = !insertion_sort(x, critical_indexes, 5);
			}
		}
		int possible_destinations[10];
		for (int i = 0; i < 5; i++)
		{
			possible_destinations[i] = partial_sol[critical_indexes[i]];
			possible_destinations[5 + i] = partial_sol[(critical_indexes[i] + 1) % dim]; //successivo di quello appena definito sopra
		}

		int subtour = 1;
		while (subtour == 1)
		{
			int* perm = randperm(10);
			subtour = 0;

			for (int i = 0; i < 10; i++)
				dest[i] = -1;

			for (int i = 0; i < 9; i++)
			{
				dest[perm[i]] = perm[i + 1];
				dest[perm[i + 1]] = perm[i];
				i++;
			}



			/*************************************************/

			print_vector(dest, 10);
			print_vector(possible_destinations, 10);
			print_vector(critical_indexes, 5);

			int* temp = build_tour(partial_sol, dest, dim, possible_destinations, &subtour);
			if (temp != NULL)
				print_vector(temp, dim);

			/*************************************************/

			if (subtour != 0)
			{
				printf("Subtour! %d\n", subtour);
				//printf("Il ciclo che prendo in esame e':\n");
				//print_vector(partial_sol, dim);
				//getchar();
			}


			if (subtour == 0)
			{
				partial_sol = temp;
			}
			else
				if (temp != NULL)
					free(temp);
			free(perm);
		}
		getchar();
	}
}
#endif //debug