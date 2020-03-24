#include "Tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include<ilcplex/cplex.h>
#include<gnuplot_c.h>



void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);

void print_graph(instance* inst);

void print_error(const char *err)
{
    printf("\n debug: %s \n", err);
    fflush(NULL);
}


void free_instance(instance *inst)
{
	free(inst->demand);
	free(inst->xcoord);
	free(inst->ycoord);
	
}

int main (int argc, char **argv)
{
   
	if(argc <2)
    {
        printf("\nArgomenti passati da input errati: %s \n", argv[0]);
        exit(1);
    }

    if (VERBOSE >=2)
    {
        for(int a=0; a<argc; ++a)
        {
            printf("%s", argv[a]);
            printf("\n");
        }
    }

    instance inst;
	
    parse_command_line(argc, argv, &inst);
	read_input(&inst);
	
	if (TSPopt(&inst)) print_error(" error within TSPopt()");
		
	for (int i = 0; i < inst.nnodes; ++i) 
	{
		printf("\nnodo %2d ) ", i+1);
		printf("x=%15.7lf - ", inst.xcoord[i]);
		printf("y=%15.7lf", inst.ycoord[i]);
	}

	print_graph(&inst);

	/*
	h_GPC_Plot* plotter;
	
	//plotter = gpc_init_xy("Best solution", "X coord", "Y coord",GPC_AUTO_SCALE, GPC_KEY_DISABLE);
	plotter = gpc_init_xy("solution", (char)*inst.xcoord, (char)*inst.ycoord, GPC_AUTO_SCALE, GPC_KEY_DISABLE);

		
	//system("gnuplot > pl 'att48.tsp' using 2:3");
	gpc_close(plotter);
	*/
    
	
	free_instance(&inst);

    return 0;
}

void read_input(instance* inst) // simplified CVRP parser, not all SECTIONs detected  
{
	int active_section = 0;

	printf("file : %s", inst->input_file);

	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL) print_error(" input file not found!");

	inst->nnodes = -1;
	inst->depot = -1;
	inst->nveh = -1;

	char line[180];
	char *par_name;
	char *token1;
	char *token2;

	int do_print = (VERBOSE >= 1000);

	while (fgets(line, sizeof(line), fin) != NULL)
	{
			
		if (VERBOSE >= 2000)
		{ 
			printf("%s", line); fflush(NULL);
		}

		if (strlen(line) <= 1) continue; // skip empty lines
		
		par_name = strtok(line, " :");
		
		if (VERBOSE >= 3000) 
		{ 
			printf("parameter \"%s\" ", par_name); fflush(NULL); 
		}

		if (strncmp(par_name, "NAME", 4) == 0)
		{
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0)
		{
			active_section = 0;
			token1 = strtok(NULL, "");
			if (VERBOSE >= 10) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0)
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "CVRP", 4) != 0) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!");
			active_section = 0;
			continue;
		}


		if (strncmp(par_name, "DIMENSION", 9) == 0)
		{
			if (inst->nnodes >= 0) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if (do_print) printf(" ... nnodes %d\n", inst->nnodes);
			inst->demand = (double*)calloc(inst->nnodes, sizeof(double));
			inst->xcoord = (double*)calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double*)calloc(inst->nnodes, sizeof(double));
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "CAPACITY", 8) == 0)
		{
			token1 = strtok(NULL, " :");
			inst->capacity = atof(token1);
			if (do_print) printf(" ... vehicle capacity %lf\n", inst->capacity);
			active_section = 0;
			continue;
		}


		if (strncmp(par_name, "VEHICLES", 8) == 0)
		{
			token1 = strtok(NULL, " :");
			inst->nveh = atoi(token1);
			if (do_print) printf(" ... n. vehicles %d\n", inst->nveh);
			active_section = 0;
			continue;
		}


		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0)
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "EUC_2D", 6) != 0) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0)
		{
			if (inst->nnodes <= 0) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;
			continue;
		}


		if (strncmp(par_name, "EOF", 3) == 0)
		{
			active_section = 0;
			break;
		}

		
		if (active_section == 1) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes) print_error(" ... unknown node in NODE_COORD_SECTION section");
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if (do_print) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i + 1, inst->xcoord[i], inst->ycoord[i]);
			continue;
		}

		if (active_section == 2) // within DEMAND_SECTION
		{
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes) print_error(" ... unknown node in NODE_COORD_SECTION section");
			token1 = strtok(NULL, " :,");
			inst->demand[i] = atof(token1);
			if (do_print) printf(" ... node %4d has demand %10.5lf\n", i + 1, inst->demand[i]);
			continue;
		}

		if (active_section == 3) // within DEPOT_SECTION
		{
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes) continue;
			if (inst->depot >= 0) print_error(" ... multiple depots not supported in DEPOT_SECTION");
			inst->depot = i;
			if (do_print) printf(" ... depot node %d\n", inst->depot + 1);
			continue;
		}

		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");

	}

	fclose(fin);
}

void print_graph(instance *inst) {

	h_GPC_Plot* plotter;
	plotter = gpc_init_xy("Best solution", "X coord", "Y coord",GPC_AUTO_SCALE, GPC_KEY_DISABLE);

	if (plotter == NULL)                       // Plot creation failed - e.g is server running ?
	{
		printf("\nPlot creation failure. Please ensure gnuplot is located on your system path\n");
		exit(1);
	}

	int sample_lenght = inst->nnodes;

	typedef struct // Complex data type
	{
		double real;
		double imag;
	} ComplexRect_s;

	ComplexRect_s CArray[48];

	for (int j = 0; j < sample_lenght; j++)         // Plot a number of arrays
	{
		for (int i = 0; i < sample_lenght; i++)             // Fill the array
		{
			CArray[i].real = inst->xcoord[i];
			CArray[i].imag = inst->ycoord[i];
		}

		if (j == 0)
		{
			gpc_plot_xy(plotter,                     // Plot handle
				CArray,                    // Dataset
				sample_lenght,                 // Number of data points
				"Data graph",                // Dataset title
				"points",                  // Plot type
				"black",                 // Colour
				GPC_NEW);                  // New plot
		}
		else
		{
			gpc_plot_xy(plotter,              // Plot handle
				CArray,                 // Dataset
				sample_lenght,         // Number of data points
				"Data graph",            // Dataset title
				"points",               // Plot type
				"black",         // Colour
				GPC_ADD);              // Add plot
		}
	}


	gpc_close(plotter);


}

void parse_command_line(int argc, char** argv, instance* inst)
{
	
	if (VERBOSE >= 100)
	{
		printf(" running %s with %d parameters \n", argv[0], argc - 1);
	}
	 
	inst->model_type = 0;						// default  

	strcpy(inst->input_file, "NULL");
	
	inst->integer_costs = 0;

	inst->available_memory = 12000;   			// available memory, in MB, for Cplex execution (e.g., 12000)
	
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)        

	int help = 0; 
	
	if (argc < 1) help = 1;

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 				// input file
		
		help = 1;
	}

	if (help || (VERBOSE >= 10))		// print current parameters
	{
		printf("\n\navailable parameters (vers. 16-may-2015) --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file);
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	if (help) exit(1);
	

}





