#ifndef CELLMODEL_HPP
#define CELLMODEL_HPP

class Cellmodel
{
protected:
  Cellmodel(){}
public:
  bool isEctopic;
  bool isS1;
  unsigned int algebraic_size;
  unsigned int constants_size;
  unsigned int states_size;
  double *ALGEBRAIC;
  double *CONSTANTS;
  double *RATES;
  double *STATES;
  virtual ~Cellmodel() {}
  virtual void initConsts() = 0;
  virtual void initConsts(double type){}
  virtual void initConsts(bool is_dutta){}
  virtual void initConsts(double type, bool is_dutta){}
  virtual void initConsts(double type, double conc, double *hill){}
  virtual void initConsts(double type, double conc, double *hill, bool is_dutta){}
  virtual void computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC) = 0;
  virtual void solveAnalytical(double dt) {};
};

#endif
