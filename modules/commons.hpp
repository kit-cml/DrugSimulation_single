#ifndef COMMONS_HPP
#define COMMONS_HPP

#include <cstdio>

#include "glob_type.hpp"
#include "param.hpp"
#include "../cellmodels/cellmodel.hpp"

// custom printf for MPI
// to avoid duplicate printing
void mpi_printf(unsigned short node_id, const char *fmt, ...);
void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...);

// CVode setup functions
void init_cvode(cvode_t* p_cvode, Cellmodel* p_cell, double tcurr, bool is_first);
void clean_cvode(cvode_t* p_cvode);

// parameter setup function
void edison_assign_params(int argc, char *argv[], param_t *p_param);

// Main interface for zipping a folder
void create_zip(const char *zip_filename, const char *target_folder);

// Main interface for extracting a zipfile
int extract_zip(const char *zip_filename, const char *target_folder);

// Zip folder recursively
// Source:
// https://github.com/kuba--/zip
void zip_walk(struct zip_t *zip, const char *path);

// delete folder that have contents
// Source:
// https://stackoverflow.com/a/54956690/981481
void remove_dir_content(const char *path);

// create a directory.
// supporting different OS.
int make_directory(const char* dirname );

// checking file availability
// supporting different OS.
int is_file_existed(const char* pathname);

#endif
