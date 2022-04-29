#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_HI(DllInfo *info)
{
  /* Register routines */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);

  /* Disable symbol search */
  R_useDynamicSymbols(info, TRUE);
}
