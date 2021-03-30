#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mudepth(void *, void *, void *, void *, void *, void *, void *);
extern void KDEloc2(void *, void *, void *, void *, void *, void *, void *);
extern void KDElocscale(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void KDErepeated(void *, void *, void *, void *, void *, void *, void *);
extern void KDErepeatedbw(void *, void *, void *, void *, void *, void *, void *);
extern void KDEsymloc(void *, void *, void *, void *, void *, void *, void *);
extern void KDEsymloc1comp(void *, void *, void *, void *, void *, void *, void *);
extern void KDEsymloc2(void *, void *, void *, void *, void *, void *, void *);
extern void multinompost(void *, void *, void *, void *, void *);
extern void mvwkde_adaptbw(void *, void *, void *, void *, void *, void *, void *, void *);
extern void mvwkde_samebw(void *, void *, void *, void *, void *, void *, void *, void *);
extern void newz(void *, void *, void *, void *, void *);
extern void normpost(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void npMSL_Estep(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void npMSL_Estep_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void npMSL_Mstep(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void npMSL_Mstep_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"mudepth",      (DL_FUNC) &mudepth,       7},
  {"KDEloc2",        (DL_FUNC) &KDEloc2,         7},
  {"KDElocscale",    (DL_FUNC) &KDElocscale,    10},
  {"KDErepeated",    (DL_FUNC) &KDErepeated,     7},
  {"KDErepeatedbw",  (DL_FUNC) &KDErepeatedbw,   7},
  {"KDEsymloc",      (DL_FUNC) &KDEsymloc,       7},
  {"KDEsymloc1comp", (DL_FUNC) &KDEsymloc1comp,  7},
  {"KDEsymloc2",     (DL_FUNC) &KDEsymloc2,      7},
  {"multinompost",   (DL_FUNC) &multinompost,    5},
  {"mvwkde_adaptbw", (DL_FUNC) &mvwkde_adaptbw,  8},
  {"mvwkde_samebw",  (DL_FUNC) &mvwkde_samebw,   8},
  {"newz",           (DL_FUNC) &newz,            5},
  {"normpost",       (DL_FUNC) &normpost,       10},
  {"npMSL_Estep",   (DL_FUNC) &npMSL_Estep,     15},
  {"npMSL_Estep_bw",(DL_FUNC) &npMSL_Estep_bw,  15},
  {"npMSL_Mstep",   (DL_FUNC) &npMSL_Mstep,     13},
  {"npMSL_Mstep_bw",(DL_FUNC) &npMSL_Mstep_bw,  13},
  {NULL, NULL, 0}
};

void R_init_mixtools(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}