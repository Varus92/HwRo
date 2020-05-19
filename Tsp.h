#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <stdlib.h>
#include <ilcplex/cplex.h>

#define VERBOSE 50

typedef struct
{
    //input data
    int depot;
    int nnodes;
    int nveh;
    double* demand;
    double* xcoord;
    double* ycoord;
    double capacity;

    //parameters
    int available_memory;
    int model_type;
    int max_nodes;
    int integer_costs;
    int num_threads;
    int randomseed;
    int ncols;
    double timelimit;						// overall time limit, in sec.s
    double cutoff; 							// cutoff (upper bound) for master
    char input_file[1000];
    char node_file[1000];

    //optimal data
    double tstart;
    double x_start;
    double y_start;
    double tbest;
    int best_lb;

   

} instance;


int TSPopt(instance*);
void print_error(const char*);
double dist(int, int, instance*);
void* get_xpos(int model_type);
int benders_method(double* xstar, int* succ, int* comp, int ncomp, instance* inst, CPXENVptr env, CPXLPptr lp);
int* get_component_array(int component, int succ[], int comp[], int size, int* c_component_dim);
int time_limit_expired(instance* inst);

#endif 