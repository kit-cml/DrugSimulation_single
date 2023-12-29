#ifndef DRUG_SIM_HPP
#define DRUG_SIM_HPP

#include "cipa_t.hpp"
#include "param.hpp"
#include "glob_type.hpp"
#include "../cellmodels/cellmodel.hpp"

// main simulation function for drug simulation
// return TRUE if EAD happened, otherwise false
bool do_drug_sim(const double conc, row_data ic50, 
const param_t* p_param, const unsigned short sample_id, const unsigned short group_id,
Cellmodel *p_cell, cvode_t *p_cvode, qinward_t *p_qin, bool is_firsttime);

#endif
