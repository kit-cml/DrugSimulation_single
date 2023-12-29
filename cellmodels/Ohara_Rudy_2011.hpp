#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "cellmodel.hpp"
#include "../enums/enum_Ohara_Rudy_2011.hpp"

class Ohara_Rudy_2011 : public Cellmodel
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void initConsts ();
  void initConsts (double celltype);
  void initConsts (bool is_dutta);
  void initConsts (double celltype, bool is_dutta);
  void initConsts (double celltype, double conc, double *ic50, bool is_dutta);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
  static double set_time_step(double TIME,double time_point,double max_time_step,
    double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC);
private:
  // apply Dutta conductance scaling based on Dutta et al. 2017
  void ___applyDutta();
  // apply drug-induced equation based on drug-induced equation
  // epsilon is used to avoid NaN value data in ic50
  // (denominator is zero, app is broken)
  void ___applyDrugEffect(double conc, double *ic50, double epsilon);
  // actual initial condition function
  // that will be called by public functions
  void ___initConsts();
  // prompt the info of celltype
  // 0 is endo, 1 is epi, 2 is M cell
  void ___printCelltype(int celltype);
};


#endif

