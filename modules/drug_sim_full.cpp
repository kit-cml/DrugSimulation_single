#include "drug_sim_full.hpp"

#ifdef TOMEK_2019
#include "../cellmodels/Tomek_model.hpp"
#else
#include "../cellmodels/Ohara_Rudy_2011.hpp"
#endif
#include "commons.hpp"

#include <cmath>

bool do_drug_sim_full(const double conc, row_data ic50, 
const param_t* p_param, const unsigned short sample_id, const unsigned short group_id,
Cellmodel *p_cell, cvode_t *p_cvode, qinward_t *p_qin, bool is_firsttime)
{
  bool is_ead;
  char buffer[255];

  // looping counter
  unsigned short idx;

  // CVode variables
  double tnext, tcurr;
  int cvode_retval;
  unsigned int icount, imax;


  // files for storing results
  // time-series result
  FILE *fp_vm, *fp_ca, *fp_dvmdt, *fp_ires, *fp_inet;
  // features result
  FILE *fp_feature;
  // saving output result
  FILE *fp_output;

  // simulation parameters
#ifdef DEBUG_MODE
  bool is_print_graph = true;
  bool is_dutta = false;
  double dt = 0.5;
  double dtw = 2.0;
  const char *drug_name = "bepridil";
  const double bcl = 2000.;
  const double inet_vm_threshold = -88.0;
  const unsigned short pace_max = 1000;
  const unsigned short celltype = 0.;
  const unsigned short last_drug_check_pace = 1000;
  const unsigned int print_freq = (1./dt) * dtw;
  unsigned short pace_count = 0;
#else
  bool is_print_graph = p_param->is_print_graph;
  bool is_dutta = p_param->is_dutta;
  double dt = p_param->dt;
  double dtw = p_param->dt_write;
  const char *drug_name = p_param->drug_name;
  const double bcl = p_param->bcl;
  const double inet_vm_threshold = p_param->inet_vm_threshold;
  const unsigned short pace_max = p_param->pace_max;
  const unsigned short celltype = (unsigned short)p_param->celltype;
  const unsigned short last_drug_check_pace = 1000;
  const unsigned int print_freq = (1./dt) * dtw;
  unsigned short pace_count = 0;
#endif

  // drug features
  double inet,qinward,qnet;
  double inal_auc, ical_auc;
  double vm_repol30, vm_repol50, vm_repol90;
  double t_depol;
  double t_ca_peak, ca_amp50, ca_amp90;
  double cad50_prev, cad50_curr, cad90_prev, cad90_curr;

  // variables to store features
  // temp_result is the result of features in 1 pace,
  // will be interchanged during the simulation.
  // cipa_result is the final result of the simulation.
  cipa_t cipa_result, temp_result;

  // TRUE if vm peak more than 0 mV
  bool is_eligible_AP;

  // apply some cell initialization
#ifdef TOMEK_2019
  p_cell->initConsts( celltype, conc, ic50.data);
#else
  p_cell->initConsts( celltype, conc, ic50.data, is_dutta);
#endif
  p_cell->CONSTANTS[BCL] = bcl;

  FILE *fp_states;
  if( p_param->is_using_output > 0 ){
#ifdef TOMEK_2019
    fp_states = fopen("output_tomek.dat", "r");
#else
    fp_states = fopen("output_orudy.dat", "r");
#endif
    if( fp_states != NULL ){
      mpi_printf(0, "Using initial condition from steady state!\n");
      int idx = 0;
      while(fgets( buffer, sizeof(buffer), fp_states) != NULL){
        p_cell->STATES[idx++] = strtod(buffer, NULL);
      }
    }
    else{
      mpi_printf(0, "No initial file found! Skipped using initial file!\n");
    }
  }

  tcurr = 0.;
  tnext = dt;
  
  init_cvode(p_cvode, p_cell, tcurr, is_firsttime);

  // generate file for time-series output
  if( is_print_graph == true ){
    snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_vmcheck_smp%d.plt", 
              conc, drug_name, conc, sample_id );
    fp_vm = fopen( buffer, "w" );
    snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_dvmdt_smp%d.plt", 
              conc, drug_name, conc, sample_id);
    fp_dvmdt = fopen( buffer, "w" );
    snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_ca_i_smp%d.plt", 
              conc, drug_name, conc, sample_id );
    fp_ca = fopen( buffer, "w" );
    snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_ires_smp%d.plt", 
              conc, drug_name, conc, sample_id );
    fp_ires = fopen( buffer, "w" );
    snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_inet_smp%d.plt", 
              conc, drug_name, conc, sample_id );
    fp_inet = fopen( buffer, "w" );

    fprintf(fp_vm, "%s %s\n", "Time", "Vm");
    fprintf(fp_dvmdt, "%s %s\n", "Time", "dVm/dt");
    fprintf(fp_ca, "%s %s\n", "Time", "cai");
    fprintf(fp_ires, "%s %s %s %s %s %s %s %s\n", 
                "Time", "INa", "INaL", "ICaL", 
                "Ito", "IKr", "IKs", "IK1");
    fprintf(fp_inet, "%s %s\n", "Time", "Inet");
  }
  snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_feature_smp%d.plt", conc, p_param->drug_name, conc, sample_id );
  fp_feature = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer), "result/%.6lf/%s_%.6lf_output_smp%d.plt", conc, p_param->drug_name, conc, sample_id );
  fp_output = fopen( buffer, "w" );
  fprintf( fp_feature, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
             "Pace", "Dvm/Dt_Repol", "Max_Dvm/Dt", "Vm_Peak", "Vm_Resting","APD90", "APD50", "APDTri", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri", "Qnet", "Qinward");

  icount = 0;
  imax = (unsigned int)((pace_max * bcl)/dt);
  inet = 0.;
  is_eligible_AP = false;
  is_ead = false;

  cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
  temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );


  while(icount < imax)
  // begin simulation loop
  {
    // solve ODE
    cvode_retval = CVode(p_cvode->cvode_mem, tnext, p_cvode->states_vec, &tcurr, CV_NORMAL);
    if( cvode_retval == CV_SUCCESS ){
      icount++;
      tnext += dt;
    }
    else{
      mpi_fprintf(0, stderr, "CVode error at sample_ID %d and concentration %.6lf at rank %d\n", 
                  sample_id, conc, mympi::rank);
      break;
    }
    // get the RATES array
    p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(p_cvode->states_vec), p_cell->ALGEBRAIC);
    // calculate inet and AUC
    if(p_cell->STATES[V] > inet_vm_threshold){
      inet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*dt;
      inal_auc += p_cell->ALGEBRAIC[INaL]*dt;
      ical_auc += p_cell->ALGEBRAIC[ICaL]*dt;
    } 

    // save temporary result
    if(pace_count >= pace_max-last_drug_check_pace && icount % print_freq == 0){
      temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
      temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
      temp_result.dvmdt_data.insert( std::pair<double, double> (tcurr, p_cell->RATES[V]) );
      temp_result.inet_data.insert( std::pair<double, double> (tcurr, inet) );
      snprintf( buffer, sizeof(buffer), "%lf %lf %lf %lf %lf %lf %lf", 
              p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito], 
              p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs], p_cell->ALGEBRAIC[IK1] );
      temp_result.ires_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
      if( is_print_graph == true ){
        fprintf(fp_ca, "%lf %lf\n", tcurr, p_cell->STATES[cai]);
        fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
        fprintf(fp_dvmdt, "%lf %lf\n", tcurr, p_cell->RATES[V]);
        fprintf(fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf\n", tcurr, 
                p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito],
                p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs], p_cell->ALGEBRAIC[IK1]);
        fprintf(fp_inet, "%lf %lf\n", tcurr, inet);
      }
    }

    // execute at the beginning of new pace
    if(icount % (int)(bcl/dt) == 0 && icount > 0){

      // assuming the pace is eligible,
      // we will extract result
      if( is_eligible_AP && pace_count >= pace_max-last_drug_check_pace-1) {
        for(std::multimap<double, double>::iterator itrmap = temp_result.cai_data.begin(); 
            itrmap != temp_result.cai_data.end() ; itrmap++ ){
          // before the peak calcium
          if( itrmap->first < t_ca_peak ){
            if( itrmap->second < ca_amp50 ) cad50_prev = itrmap->first;
            if( itrmap->second < ca_amp90 ) cad90_prev = itrmap->first;
          }
          // after the peak calcium
          else{
            if( itrmap->second > ca_amp50 ) cad50_curr = itrmap->first;
            if( itrmap->second > ca_amp90 ) cad90_curr = itrmap->first;
          }
        }
        temp_result.cad50 = cad50_curr - cad50_prev;
        temp_result.cad90 = cad90_curr - cad90_prev;
        temp_result.qnet = inet/1000.0;
        temp_result.inal_auc = inal_auc;
        temp_result.ical_auc = ical_auc;
        temp_result.vm_dia = p_cell->STATES[V];
        temp_result.ca_dia = p_cell->STATES[cai];

        // replace result with steeper repolarization AP or first pace from the last 250 paces
        if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
          cipa_result = temp_result;
        }

      }// end pace eligible


      // because this is full pace simulation,
      // means that all of the result from the whole pacing will be printed
      fprintf( fp_feature, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 pace_count,
                 temp_result.dvmdt_repol,
                 temp_result.dvmdt_max,
                 temp_result.vm_peak,
                 temp_result.vm_dia,
                 temp_result.apd90,
                 temp_result.apd50,
                 temp_result.apd90-temp_result.apd50,
                 temp_result.ca_peak,
                 temp_result.ca_dia,
                 temp_result.cad90,
                 temp_result.cad50,
                 temp_result.cad90-temp_result.cad50);
      qnet = inet/1000.0;
      if( (int)ceil(conc) == 0 ) {
        p_qin->inal_auc_control = cipa_result.inal_auc;
        p_qin->ical_auc_control = cipa_result.ical_auc;
        fprintf( fp_feature, "%d %lf %lf\n", pace_count, qnet, 0.0 );
      }
      else{
        p_qin->inal_auc_drug = cipa_result.inal_auc;
        p_qin->ical_auc_drug = cipa_result.ical_auc;
        qinward =  ( (p_qin->inal_auc_drug/p_qin->inal_auc_control) + (p_qin->ical_auc_drug/p_qin->ical_auc_control) ) * 0.5;
        fprintf( fp_feature, "%d %lf %lf\n", pace_count, qnet, qinward );
      }


      // resetting inet and AUC values
      // and increase the pace count
      pace_count++;
      inet = 0.;
      inal_auc = 0.;
      ical_auc = 0.;
      temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
      t_ca_peak = tcurr;
      t_depol = (p_cell->CONSTANTS[BCL]*pace_count)+p_cell->CONSTANTS[stim_start];
      is_eligible_AP = false;


    }// end beginning pace process

    // entering the last 250 paces
    if(pace_count >= pace_max-last_drug_check_pace){
    
      // get maximum dvmdt
      if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];

      // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
      if( icount % (int)(((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+6.)) / dt) == 0 ){
        temp_result.vm_peak = p_cell->STATES[V];
        if( temp_result.vm_peak > 0. ){
          vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
          is_eligible_AP = true;
        }
        else is_eligible_AP = false;
      }
      // these operations will be executed if it's eligible AP and executed at the beginning of repolarization
      else if( is_eligible_AP && tcurr > (p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+6.) ){
        // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
        if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
          temp_result.dvmdt_repol = p_cell->RATES[V];
        }
        // get the APD90, APD50, peak calcium, 50% and 90% of amplitude of Calcium, and time of peak calcium
        if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-2 ) temp_result.apd50 = tcurr - t_depol;
        if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 ) temp_result.apd90 = tcurr - t_depol;
        if( temp_result.ca_peak < p_cell->STATES[cai] ){  
          temp_result.ca_peak = p_cell->STATES[cai];
          ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
          ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
          t_ca_peak = tcurr;
        } 
      }   
    }// end of last 250 paces operations


  }// end simulation loop

  // save the last output state
  for(idx = 0; idx < p_cell->states_size; idx++ ) fprintf( fp_output, "%lf\n", p_cell->STATES[idx] );

  // clean the memories
  fclose(fp_feature);
  if( is_print_graph == true ){
    fclose(fp_ires);
    fclose(fp_inet);
    fclose(fp_ca);
    fclose(fp_dvmdt);
    fclose(fp_vm);
  }

  return is_ead;
}
