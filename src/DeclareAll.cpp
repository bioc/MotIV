#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP motifMatch (SEXP  cc, SEXP align, SEXP top, SEXP go, SEXP ge, SEXP inputPWM, SEXP inputDB, SEXP inputScores);

R_CallMethodDef callMethods[]  = {
	{"motifMatch", (DL_FUNC)&motifMatch, 8},
	{NULL, NULL, 0}
};
