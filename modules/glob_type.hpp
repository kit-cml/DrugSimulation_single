#ifndef GLOB_TYPE_HPP
#define GLOB_TYPE_HPP

#include <vector>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>


// global variable for MPI.
struct mympi
{
  static char host_name[255];
  static int host_name_len;
  static int rank;
  static int size;
};

// CVode struct
// (made by Ali)
typedef struct {
    void* cvode_mem;
    N_Vector states_vec;
    SUNMatrix matrix;
    SUNLinearSolver solver;
} cvode_t;

// data structure for IC50
typedef struct row_data { double data[14]; } row_data;
typedef std::vector< row_data > drug_t;

// data structure to store
// ICaL/INaL control value
// for calculating qinward
// control means 0 concentration
// otherwise, drugs
typedef struct{
  double ical_auc_control;
  double inal_auc_control;
  double ical_auc_drug;
  double inal_auc_drug;
} qinward_t;

#endif
