#include "init.h"
#include "parameters.h"

void init_population_field(real* population) {
    unsigned shift = 0; 
    
    for(unsigned component = 0; component < parameters.discretization; ++component) {
        
        shift = parameters.num_lattices * component;
        for(unsigned i = 0; i < parameters.num_lattices; ++i) {
            population[i + shift] = weights[component];
        }
    
    }
}


