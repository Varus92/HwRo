#ifndef TSP_H
#define TSP_H

//#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include<Windows.h>



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


#endif 