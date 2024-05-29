#pragma once
#define M_PI 3.14159265358979323846

#define USE_MPI

#define HAS_CELL_INDEX

#define WITH_SELF


//#define RUN_NRQS
//#define RUN_QS
#define RUN_OFFSET_QS

#ifdef RUN_NRQS
#define MULT_SPECIES
#endif