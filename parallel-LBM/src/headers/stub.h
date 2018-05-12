#ifndef SEQUENTIAL_LBM_SRC_HEADERS_STUB_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_STUB_H_

void ReadInputFilesStub(struct SimulationParametes &,
                        struct BoundaryInfo &, 
                        struct Constants &,
                        char *, 
                        char *);

void InitParametersStub(SimulationParametes &,
                        Constants &);

void InitBoundaryConditionStub(BoundaryInfo &);

#endif  //  SEQUENTIAL_LBM_SRC_HEADERS_STUB_H_
