#include <fstream>
#include <string>
#include <iomanip>
#include "headers/parameters.h"
#include "headers/helper.h"

int GetIndex(int index_i, int index_j, int dim) {
    return dim * (index_i + index_j * parameters.width);
}


void WriteVectorFieldToFile(real *velocity, std::string file_name) {
    real delta_x = parameters.delta_x;
    std::string delimiter = "\t";
    real scaling_factor = 50.0;

    std::ofstream file;
    file.open(file_name);
    file << std::fixed << std::setprecision(6) << std::endl;
    file << "# X    Y   delta_x    delta_y" << std::endl;

    for (int j = 0; j < parameters.height; ++j) {
        for (int i = 0; i < parameters.width; ++i) {
            
            int index = GetIndex(i, j);
            file << real(i) << delimiter 
                 << real(j) << delimiter 
                 << scaling_factor * velocity[index] << delimiter 
                 << scaling_factor * velocity[index + parameters.num_lattices] << delimiter 
                 << std::endl;
        }
    }
    file.close();
}
