#ifndef Mat_euristic_h
#define Mat_euristic_h

#include "Tsp.h"

void fix(int* edges, int dim, float percentage, instance* inst, CPXENVptr env, CPXLPptr lp);
int* get_circuit(instance* inst, int* comp, int* succ);

#endif // Mat_euristic_h