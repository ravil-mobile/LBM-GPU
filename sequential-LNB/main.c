#include <stdio.h>
#include "parameters.h"
#include "stub.h"

int main() {
    char* parameter_file = "parameter.txt";
    char* boundary_file = "boundary.txt";

    ReadInputFilesStub(parameter_file,
                       boundary_file);

    return 0;
}
