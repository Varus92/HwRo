#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <stdlib.h>
#include <ilcplex/cplex.h>

#define VERBOSE     50

typedef struct 
{
    //input data
    int depot;
    int nnodes;
    int nveh;
    double* demand;
    double *xcoord;
    double *ycoord;
    double capacity;

    //parameters
    int available_memory;
    int model_type;
    int max_nodes;
    int integer_costs;
    char input_file[1000];
    char node_file[1000];

    //optimal data
    double x_start;
    double y_start;
    double tbest;
    
}instance;


// double dist(int i, int j, instance* inst);
int TSPopt(instance*);
void print_error(const char*);

#endif 