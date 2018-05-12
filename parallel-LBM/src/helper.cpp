#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <unistd.h>
#include <stdio.h>

#include "headers/gnuplot_i.h"
#include "headers/parameters.h"
#include "headers/helper.h"

int GetIndex(int index_i, int index_j, int dim) {
    return dim * (index_i + index_j * parameters.width);
}

void GnuplotCmd(gnuplot_ctrl *frame, std::string command) {
    gnuplot_cmd(frame, const_cast<char*>(command.c_str()));
}

void SetupGnuPlots(gnuplot_ctrl *velocity_frame,
                   gnuplot_ctrl *density_frame) {
    std::string command = "";

    command = std::string("set xrange [0:")
            + std::to_string(parameters.width)
            + std::string("]");
    
    GnuplotCmd(velocity_frame, command);

    command = std::string("set yrange [0:")
            + std::to_string(parameters.height)
            + std::string("]");
    
    GnuplotCmd(velocity_frame, command);

    // GnuplotCmd(velocity_frame, "set cbrange [0.9:1.3]");
    GnuplotCmd(velocity_frame, "set view map");
    GnuplotCmd(velocity_frame, "set pm3d at b map");


    // GnuplotCmd(density_frame, "set palette rbg 33,13,10");
    GnuplotCmd(density_frame, "set view map");
    GnuplotCmd(density_frame, "set pm3d at b map");
    // GnuplotCmd(density_frame, "set cbrange [0.7:1.3]");
}

void DisplayResults(real *velocity, gnuplot_ctrl *velocity_frame,
                    real *density,  gnuplot_ctrl *density_frame) {
    real delta_x = parameters.delta_x;
    int num_lattices = parameters.num_lattices;

    std::ofstream velocity_field_file;
    std::ofstream density_field_file;

    velocity_field_file.open("velocity-data.dat");
    density_field_file.open("density-data.dat");

    velocity_field_file << std::fixed << std::setprecision(6);
    density_field_file << std::fixed << std::setprecision(6);


    real scaling_factor = 40.0;
    std::string delimiter = "\t";
    for (int j = 1; j < parameters.height - 1; ++j) {
        for (int i = 1; i < parameters.width - 1; ++i) {
            int index = GetIndex(i, j);

            velocity_field_file << real(i) << delimiter
                                << real(j) << delimiter
                                << velocity[index] << delimiter
                                << std::endl;
            /*
            velocity_field_file << real(i) << delimiter
                                << real(j) << delimiter
                                << scaling_factor * velocity[index] << delimiter
                                << scaling_factor * velocity[index + num_lattices] << delimiter
                                << magnitude << delimiter
                                << std::endl;
            */
            if (density != NULL) {
                density_field_file << real(i) << delimiter
                                   << real(j) << delimiter
                                   << density[index] << delimiter
                                   << std::endl;
            }
        }
        velocity_field_file << std::endl;
        density_field_file << std::endl;
    }
    velocity_field_file.close();
    density_field_file.close();

    GnuplotCmd(velocity_frame, "splot 'velocity-data.dat' u 1:2:3");

    if (density_frame != NULL) {
        GnuplotCmd(density_frame, "splot 'density-data.dat' u 1:2:3");
    }
    //usleep(10000);
}
