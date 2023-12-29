#ifndef TOMEK_MODEL_ENDO_HPP
#define TOMEK_MODEL_ENDO_HPP

#include "cellmodel.hpp"
#include "../enums/enum_Tomek_model.hpp"

class Tomek_model : public Cellmodel
{
public:
  Tomek_model();
  ~Tomek_model();
  void initConsts ();
  void initConsts (double celltype);
  void initConsts (double celltype, double conc, double *ic50);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
private:
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

