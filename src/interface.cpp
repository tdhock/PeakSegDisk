/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegFPOPLog.h"
#include <R.h>//for error
#include <R_ext/Rdynload.h>//for registering routines.
#include <Rinternals.h>//for STRSXP

static int PeakSegFPOP_nargs = 3;
static R_NativePrimitiveArgType PeakSegFPOP_types[] = {STRSXP, STRSXP, STRSXP};
void PeakSegFPOP_interface
(char **file_vec, char **pen_vec, char **temp_vec){
  char *bedGraph = file_vec[0];
  char *penalty = pen_vec[0];
  char *db = temp_vec[0];
  int status = PeakSegFPOP_disk(bedGraph, penalty, db);
  if(status==ERROR_PENALTY_NOT_FINITE){
    Rf_error("penalty=%s but must be finite", penalty);
  }
  if(status==ERROR_PENALTY_NEGATIVE){
    Rf_error("penalty=%s must be non-negative", penalty);
  }
  if(status==ERROR_UNABLE_TO_OPEN_BEDGRAPH){
    Rf_error("unable to open input file for reading %s", bedGraph);
  }
  if(status==ERROR_NOT_ENOUGH_COLUMNS){
    Rf_error("each line of input data file %s should have exactly four columns", bedGraph);
  }
  if(status==ERROR_NON_INTEGER_DATA){
    Rf_error("fourth column of input data file %s should be integer", bedGraph);
  }
  if(status==ERROR_INCONSISTENT_CHROMSTART_CHROMEND){
    Rf_error("there should be no gaps (columns 2-3) in input data file %s", bedGraph);
  }
  if(status==ERROR_WRITING_COST_FUNCTIONS){
    Rf_error("unable to write to cost function database file %s", db);
  }
  if(status==ERROR_WRITING_LOSS_OUTPUT){
    Rf_error("unable to write to loss output file %s_penalty=%s_loss.tsv",
	  bedGraph, penalty);
  }
  if(status==ERROR_WRITING_SEGMENTS_OUTPUT){
    Rf_error("unable to write to segments output file %s_penalty=%s_segments.bed",
	  bedGraph, penalty);
  }
  if(status==ERROR_NO_DATA){
    Rf_error("input file %s contains no data", bedGraph);
  }
  if(status==ERROR_PENALTY_NOT_NUMERIC){
    Rf_error
      ("penalty string '%s' is not numeric; it should be convertible to double",
       penalty);
  }
  if(status != 0){
    Rf_error("error code %d", status);
  }
}
  
R_CMethodDef cMethods[] = {
  {"PeakSegFPOP_interface",
   (DL_FUNC) &PeakSegFPOP_interface, PeakSegFPOP_nargs, PeakSegFPOP_types
  },
  {NULL, NULL, 0, NULL}
};

extern "C" {
  void R_init_PeakSegDisk(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }
}
