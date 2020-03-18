#include "Tsp.h"


double dist(int i, int j, instance *inst)
{
    double dx = inst->xcoord[i] - inst->ycoord[j];
    double dy = inst->ycoord[i] - inst->ycoord[j];
    
    if( !inst->integer_costs)
        return sqrt(dx*dx+dy*dy);
    
    int distance = sqrt(dx*dx+dy*dy + 0.499999); //l'intero piu vicino

    return distance +0.0;
}
