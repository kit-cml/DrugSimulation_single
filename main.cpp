#ifdef TOMEK_2019
#include "cellmodels/Tomek_model.hpp"
#else
#include "cellmodels/Ohara_Rudy_2011.hpp"
#endif
#include "modules/cipa_t.hpp"
#include "modules/commons.hpp"
#include "modules/glob_type.hpp"
#include "modules/drug_sim.hpp"
#include "modules/drug_sim_full.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

// constants to avoid magic values
static const char *RESULT_FOLDER_PATH = "result";
static const double CONTROL_CONC = 0.;

// get the IC50 data from file
drug_t get_IC50_data_from_file(const char* file_name);
int check_IC50_content(const drug_t* ic50, const param_t* p_param);

// define MPI data structure for qinward_t to be broadcasted
MPI_Datatype create_mpi_qinward_t();

int main(int argc, char **argv)
{
  // enable real-time output in stdout
  setvbuf( stdout, NULL, _IONBF, 0 );

  // initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mympi::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &mympi::rank );
  MPI_Get_processor_name(mympi::host_name, &mympi::host_name_len);

  // buffer for writing in snprintf() function
  char buffer[255];

  // looping control and error control
  unsigned short idx;
  unsigned short sample_id;
  unsigned short group_id;
  unsigned short error_code;

  // input parameter object
  param_t *p_param;
  p_param = new param_t();
  p_param->init();
  edison_assign_params(argc,argv,p_param);
  p_param->show_val();

  // cell model and drug induced-related vars
  drug_t ic50;
  Cellmodel* p_cell;
  //std::array<double, 4> concs = {0., 33., 165., 330.};

  // Cvode vars
  cvode_t *ode_solver;
  bool cvode_firsttime;

  // variable to store AUC value of INaL and ICaL
  qinward_t *qin;

  // EAD counter;
  bool is_ead;
  unsigned short ead_counter;
  unsigned short ead_total;

  // get concentration from input
  std::vector<double> concs;
  char *token = strtok( p_param->concs, "," );
  while( token != NULL )
  { // begin data tokenizing
      concs.push_back(strtod(token, NULL));
      token = strtok(NULL, ",");
  } // end data tokenizing


  ic50 = get_IC50_data_from_file(p_param->hill_file);
  error_code = check_IC50_content(&ic50, p_param);
  if(error_code != 0) return 1;

  // make a directory for each concentration
  // and wait till all directories have been created
  if(mympi::rank == 0){
    make_directory(RESULT_FOLDER_PATH);
    snprintf( buffer, sizeof(buffer), "%s/%.6lf", RESULT_FOLDER_PATH, CONTROL_CONC);
    if(is_file_existed(RESULT_FOLDER_PATH) == 0) make_directory(buffer);
    for( idx = 0; idx < concs.size(); idx++ )
    { // begin concentration loop
      snprintf( buffer, sizeof(buffer), "%s/%.6lf", RESULT_FOLDER_PATH, concs[idx] );
      if(is_file_existed(RESULT_FOLDER_PATH) == 0) make_directory(buffer);
    } // end concentration loop
  }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef TOMEK_2019
  mpi_printf(0, "Using Tomek cell model\n");
  p_cell = new Tomek_model();
#else
  mpi_printf(0, "Using O'Hara Rudy cell model\n");
  p_cell = new Ohara_Rudy_2011();
#endif
  ode_solver = new cvode_t();
  cvode_firsttime = true;
  qin = new qinward_t();
  ead_counter = 0;
  is_ead = false;

  // execute control in order to get
  // qinward control values.
  if( mympi::rank == 0 ){
    sample_id = 0;
    group_id = 0;
    is_ead = do_drug_sim(CONTROL_CONC, ic50[0],
              p_param, sample_id, group_id,
              p_cell, ode_solver, qin, cvode_firsttime);
  }

  printf("Before INaL at rank %d: %lf\n",mympi::rank, qin->inal_auc_control);
  printf("Before ICaL at rank %d: %lf\n",mympi::rank, qin->ical_auc_control);

  // provide MPI_Datatype for broadcasting qinward_t struct.
  MPI_Datatype mpi_qinward_t = create_mpi_qinward_t();
  MPI_Bcast(qin, 1, mpi_qinward_t, 0, MPI_COMM_WORLD );

  printf("After INaL at rank %d: %lf\n",mympi::rank, qin->inal_auc_control);
  printf("After ICaL at rank %d: %lf\n",mympi::rank, qin->ical_auc_control);


  double t_begin = MPI_Wtime();
  // sample-based simulation
  if( p_param->simulation_mode == 0 )
  {
    mpi_printf(0, "Simulate sample-based simulation\n");
    for( sample_id = mympi::rank, group_id = 0;
         sample_id < ic50.size();
         sample_id += mympi::size, group_id++ )
    { // begin sample loop
      printf("Sample_ID:%hu  Count:%lu Rank:%d \nData: ", 
        sample_id, ic50.size(), mympi::rank );
    
      for( idx = 0; idx < 14; idx++ ){
        printf("%lf|", (ic50[sample_id]).data[idx]);
      }
      printf("\n");
   

      for( idx = 0; idx < concs.size(); idx++ )
      { // begin concentration loop
        // execute main simulation function
        is_ead = do_drug_sim(concs[idx], ic50[sample_id], 
                  p_param, sample_id, group_id,
                  p_cell, ode_solver, qin, cvode_firsttime);
        if(is_ead == true){
          ead_counter++;
          printf("EAD happened in drug %s sample %hu concentration %.6lf\n", p_param->drug_name, sample_id, concs[idx]);
        }
        if(cvode_firsttime == true) cvode_firsttime = false;
      } // end concentration loop
    } // end sample loop
  }
  // full-pace optimal sample simulation
  else{
#ifdef PROTOTYPE
    sample_id = 0;
    mpi_printf(0, "Simulate full-pace simulation\n");
    for( idx = 0; idx < concs.size(); idx++ )
    { // begin concentration loop
      // execute main simulation function
      is_ead = do_drug_sim_full(concs[idx], ic50[sample_id], 
                p_param, sample_id, group_id,
                p_cell, ode_solver, qin, cvode_firsttime);
      if(cvode_firsttime == true) cvode_firsttime = false;
    } // end concentration loop
#endif
    for( sample_id = mympi::rank, group_id = 0;
         sample_id < ic50.size();
         sample_id += mympi::size, group_id++ )
    { // begin sample loop
      printf("Sample_ID:%hu  Count:%lu Rank:%d \nData: ", 
        sample_id, ic50.size(), mympi::rank );
    
      for( idx = 0; idx < 14; idx++ ){
        printf("%lf|", (ic50[sample_id]).data[idx]);
      }
      printf("\n");
   

      for( idx = 0; idx < concs.size(); idx++ )
      { // begin concentration loop
        // execute main simulation function
        is_ead = do_drug_sim_full(concs[idx], ic50[sample_id], 
                  p_param, sample_id, group_id,
                  p_cell, ode_solver, qin, cvode_firsttime);
        if(is_ead == true){
          ead_counter++;
          printf("EAD happened in drug %s sample %hu concentration %.6lf\n", p_param->drug_name, sample_id, concs[idx]);
        }
        if(cvode_firsttime == true) cvode_firsttime = false;
      } // end concentration loop
    } // end sample loop

  }

  double t_end = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);

  // showing total number of EAD in each processor, and sum it all
  printf("EAD total cases in processor %d: %hu\n", mympi::rank, ead_counter);
  MPI_Reduce( &ead_counter, &ead_total, 1, MPI_UNSIGNED_SHORT, MPI_SUM, 0, MPI_COMM_WORLD );

  // post-procssing
  if( mympi::rank == 0 ){
    char buff_time[15] = {'\0'};
    FILE *fp_perf_log;
    struct tm* tm_info;
    time_t timer;

    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );
    fp_perf_log = fopen( buffer, "w" );

    // zip result folder
    snprintf( buffer, sizeof(buffer), "result_%s_%s.zip", p_param->drug_name, buff_time );
    create_zip(buffer, "result");

    // delete the content inside result folder
    remove_dir_content("result");

    // copy the zip into the result folder
    char buff_new[150] = {'\0'};
    snprintf( buffer, sizeof(buffer), "result_%s_%s.zip", p_param->drug_name, buff_time );
    snprintf( buff_new, sizeof(buff_new), "result/result_%s_%s.zip", p_param->drug_name, buff_time );
    if( rename( buffer, buff_new ) == 0 ){
      printf("Zip file %s moved succesfully to %s!!\n", buffer, buff_new);
    }
    else{
      printf("Error moving file %s to %s!!\n", buffer, buff_new);
    }

    // make a timestamp log
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );
    fp_perf_log = fopen( buffer, "w" );
    fprintf( fp_perf_log, "Iteration BCL Timestep CPU Sample Pacing_Time(Minutes) Total_EAD_Cases Celltype\n" );
    fprintf( fp_perf_log, "%lld %.6lf %.6lf %d %lu %d %lf %hu\n", (long long int)(ic50.size()*concs.size()*p_param->pace_max*p_cell->CONSTANTS[BCL]/p_param->dt), p_cell->CONSTANTS[BCL], p_param->dt, mympi::size, ic50.size(), p_param->pace_max, (t_end-t_begin)/60., ead_total, p_cell->CONSTANTS[celltype] );

    fclose(fp_perf_log);
  }


  // memory cleaning and finalize the program
  delete ode_solver;
  delete p_param;
  delete qin;
  delete p_cell;

  mpi_printf(0, "Computation Time: %lf minutes\n", (t_end-t_begin)/60.);
  mpi_printf(0, "Total EAD case: %hu\n", ead_total);

  MPI_Finalize();
}

drug_t get_IC50_data_from_file(const char* file_name)
{
  FILE *fp_drugs;
  drug_t ic50;
  char *token, buffer[255];
  row_data temp_array;
  unsigned short idx;

  if( (fp_drugs = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s in %s at rank %d\n", 
      file_name, mympi::host_name, mympi::rank);
    return ic50;
  }

  fgets(buffer, sizeof(buffer), fp_drugs); // skip header
  while( fgets(buffer, sizeof(buffer), fp_drugs) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array.data[idx++] = strtod(token, NULL);
      token = strtok(NULL, ",");
    } // end data tokenizing
    ic50.push_back(temp_array);
  } // end line reading

  fclose(fp_drugs);
  return ic50;
}

int check_IC50_content(const drug_t* ic50, const param_t* p_param)
{
  if(ic50->size() == 0){
    mpi_printf(0, "Something problem with the IC50 file!\n");
    return 1;
  }
  else if(ic50->size() > 2000){
    mpi_printf(0, "Too much input! Maximum sample data is 2000!\n");
    return 2;
  }
  else if(p_param->pace_max < 750 && p_param->pace_max > 1000){
    mpi_printf(0, "Make sure the maximum pace is around 750 to 1000!\n");
    return 3;
  }
  else if(mympi::size > ic50->size()){
    mpi_printf(0, "%s\n%s\n", 
                "Overflow of MPI Process!",
                "Make sure MPI Size is less than or equal the number of sample");
    return 4;
  }
  else{
    return 0;
  }
}

MPI_Datatype create_mpi_qinward_t()
{
  // create MPI_Datatype for qinward_t
  // Reference:
  // https://rookiehpc.github.io/mpi/docs/mpi_type_create_struct/index.html
  MPI_Datatype mpi_qinward_t;
  int length[4] = {1,1,1,1};

  // calculate displacement for each member of the struct
  MPI_Aint displacements[4];
  qinward_t dummy_qin;
  MPI_Aint base_address;

  // blame EDISON for not having MPI3.0
#if( MPI_VERSION >= 3 )
  mpi_printf(0, "This has MPI version 3 or more\n");
  MPI_Get_address(&dummy_qin, &base_address);
  MPI_Get_address(&dummy_qin.ical_auc_control, &displacements[0]);
  MPI_Get_address(&dummy_qin.inal_auc_control, &displacements[1]);
  MPI_Get_address(&dummy_qin.ical_auc_drug, &displacements[2]);
  MPI_Get_address(&dummy_qin.inal_auc_drug, &displacements[3]);
  displacements[0] = MPI_Aint_diff( displacements[0], base_address );
  displacements[1] = MPI_Aint_diff( displacements[1], base_address );
  displacements[2] = MPI_Aint_diff( displacements[2], base_address );
  displacements[3] = MPI_Aint_diff( displacements[3], base_address );

  // define the new data type mpi_qinward_t using MPI_Type_create_struct()
  MPI_Datatype block_count[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Type_create_struct(4, length, displacements, block_count, &mpi_qinward_t  );
#else
  // Alternative for MPI lower than 3.0
  // Reference:
  // https://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi
  mpi_printf(0, "This has MPI version lower than 3.0\n");
  MPI_Address(&dummy_qin, &base_address);
  MPI_Address(&dummy_qin.ical_auc_control, &displacements[0]);
  MPI_Address(&dummy_qin.inal_auc_control, &displacements[1]);
  MPI_Address(&dummy_qin.ical_auc_drug, &displacements[2]);
  MPI_Address(&dummy_qin.inal_auc_drug, &displacements[3]);

  displacements[0] = offsetof(qinward_t, ical_auc_control);
  displacements[1] = offsetof(qinward_t, inal_auc_control);
  displacements[2] = offsetof(qinward_t, ical_auc_drug);
  displacements[3] = offsetof(qinward_t, inal_auc_drug);

  // define the new data type mpi_qinward_t using MPI_Type_create_struct()
  MPI_Datatype block_count[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Type_struct(4, length, displacements, block_count, &mpi_qinward_t  );
#endif
  MPI_Type_commit(&mpi_qinward_t);

  return mpi_qinward_t;
}
