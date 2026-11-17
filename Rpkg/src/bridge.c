#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "pmlr.h"

#ifdef _WIN32
#include <windows.h>
static HMODULE pmlr_library = NULL;
static HMODULE julia_runtime = NULL;
#define PMLR_SEP "\\"
static void *lookup_symbol(const char *name) {
    return (void *)GetProcAddress(pmlr_library, name);
}
#else
#include <dlfcn.h>
static void *pmlr_library = NULL;
static void *julia_runtime = NULL;
#define PMLR_SEP "/"
static void *lookup_symbol(const char *name) {
    return dlsym(pmlr_library, name);
}
#endif

typedef void (*init_julia_fn)(int, char **);
typedef int32_t (*fit_fn)(const int64_t *, const double *,
                         int32_t, int32_t, int32_t, int32_t, int32_t, int32_t,
                         double, int32_t, int32_t, int32_t, double, double,
                         const double *, int32_t);
typedef int32_t (*length_fn)(void);
typedef int32_t (*copy_fn)(double *, int32_t);
typedef int32_t (*error_fn)(unsigned char *, int32_t);

static fit_fn call_fit = NULL;
static length_fn result_length = NULL;
static copy_fn copy_result = NULL;
static error_fn last_error = NULL;
static int initialized = 0;

static void libpath(char *dest, size_t size, const char *root, const char *relative) {
    int written = snprintf(dest, size, "%s" PMLR_SEP "%s", root, relative);
    if (written < 0 || (size_t)written >= size)
        Rf_error("Julia bundle path is too long");
}

SEXP C_pmlr_init(SEXP rootS) {
    if (initialized)
        return ScalarLogical(1);
    if (TYPEOF(rootS) != STRSXP || LENGTH(rootS) != 1)
        Rf_error("Expected the compiled Julia bundle directory");

    const char *root = CHAR(STRING_ELT(rootS, 0));
    char path[8192];
    init_julia_fn init_julia_ptr;

#ifdef _WIN32
    libpath(path, sizeof(path), root, "bin\\pmlr.dll");
    char bindir[8192];
    libpath(bindir, sizeof(bindir), root, "bin");
    SetDllDirectoryA(bindir);
    pmlr_library = LoadLibraryExA(path, NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
    if (pmlr_library == NULL)
        Rf_error("Cannot load %s (Windows error %lu)", path, (unsigned long)GetLastError());
#else
#ifdef __APPLE__
    libpath(path, sizeof(path), root, "lib/libjulia.dylib");
#else
    libpath(path, sizeof(path), root, "lib/libjulia.so");
#endif
    julia_runtime = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
    if (!julia_runtime)
        Rf_error("Cannot load Julia runtime: %s", dlerror());
#ifdef __APPLE__
    libpath(path, sizeof(path), root, "lib/libpmlr.dylib");
#else
    libpath(path, sizeof(path), root, "lib/libpmlr.so");
#endif
    pmlr_library = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
    if (!pmlr_library)
        Rf_error("Cannot load compiled model: %s", dlerror());
#endif

    init_julia_ptr = (init_julia_fn)lookup_symbol("init_julia");
    call_fit = (fit_fn)lookup_symbol("pmlr_fit");
    result_length = (length_fn)lookup_symbol("pmlr_result_length");
    copy_result = (copy_fn)lookup_symbol("pmlr_copy_result");
    last_error = (error_fn)lookup_symbol("pmlr_last_error");
    if (!init_julia_ptr || !call_fit || !result_length || !copy_result || !last_error)
        Rf_error("The compiled library is missing one or more required exported symbols");

    char arg0[] = "R";
    char *argv[] = {arg0, NULL};
    init_julia_ptr(1, argv);
    initialized = 1;
    return ScalarLogical(1);
}

static int32_t intarg(SEXP arg, const char *label) {
    if (TYPEOF(arg) != INTSXP || LENGTH(arg) != 1 || INTEGER(arg)[0] == NA_INTEGER)
        Rf_error("Invalid integer argument: %s", label);
    return (int32_t)INTEGER(arg)[0];
}

static double dblarg(SEXP arg, const char *label) {
    if (TYPEOF(arg) != REALSXP || LENGTH(arg) != 1 || !R_FINITE(REAL(arg)[0]))
        Rf_error("Invalid numeric argument: %s", label);
    return REAL(arg)[0];
}

SEXP C_pmlr_fit(SEXP xS, SEXP covS, SEXP qS, SEXP kS, SEXP cS,
                SEXP tolS, SEXP maxitS, SEXP maxrestS, SEXP replixS,
                SEXP fixS, SEXP loadingS, SEXP latentS) {
    if (!initialized)
        Rf_error("Call Julia runtime initialization first");
    if (TYPEOF(xS) != INTSXP && TYPEOF(xS) != REALSXP)
        Rf_error("X must be numeric or integer");
    if (TYPEOF(covS) != REALSXP || !Rf_isMatrix(covS))
        Rf_error("ics must be a numeric matrix");
    if (!Rf_isMatrix(xS))
        Rf_error("X must be a matrix");

    SEXP dimx = Rf_getAttrib(xS, R_DimSymbol);
    SEXP dimcov = Rf_getAttrib(covS, R_DimSymbol);
    int32_t n = (int32_t)INTEGER(dimx)[0];
    int32_t p = (int32_t)INTEGER(dimx)[1];
    int32_t ncov = (int32_t)INTEGER(dimcov)[1];
    if (INTEGER(dimcov)[0] != n)
        Rf_error("X and ics must have the same number of rows");

    int32_t q = intarg(qS, "q");
    int32_t k = intarg(kS, "K");
    int32_t c = intarg(cS, "c");
    int32_t maxit = intarg(maxitS, "maxit");
    int32_t maxrest = intarg(maxrestS, "maxNumRestart");
    int32_t replix = intarg(replixS, "replix");
    double tol = dblarg(tolS, "tol");
    double fix = dblarg(fixS, "fixthrehsolds");
    double latent = dblarg(latentS, "latentmeans");

    R_xlen_t cells = XLENGTH(xS);
    int64_t *x = (int64_t *)R_alloc((size_t)cells, sizeof(int64_t));
    for (R_xlen_t i = 0; i < cells; ++i) {
        double val = (TYPEOF(xS) == REALSXP) ? REAL(xS)[i] : (double)INTEGER(xS)[i];
        if (!R_FINITE(val) || val != floor(val) || val < 1.0 || val > c)
            Rf_error("X must contain integer categories between 1 and c, without NA");
        x[i] = (int64_t)val;
    }

    int32_t has_loading = 0;
    const double *loading = NULL;
    if (TYPEOF(loadingS) == REALSXP && Rf_isMatrix(loadingS)) {
        SEXP dl = Rf_getAttrib(loadingS, R_DimSymbol);
        if (INTEGER(dl)[0] != p || INTEGER(dl)[1] != q)
            Rf_error("L_str must be a p x q matrix");
        has_loading = 1;
        loading = REAL(loadingS);
    }

    int32_t status = call_fit(x, REAL(covS), n, p, q, k, c, ncov,
                              tol, maxit, maxrest, replix, fix, latent,
                              loading, has_loading);
    if (status != 0) {
        unsigned char message[4096] = {0};
        last_error(message, (int32_t)sizeof(message));
        Rf_error("Julia model error: %s", (char *)message);
    }
    int32_t nout = result_length();
    if (nout <= 0)
        Rf_error("The Julia model returned an empty result");
    SEXP result = PROTECT(Rf_allocVector(REALSXP, nout));
    if (copy_result(REAL(result), nout) != 0) {
        UNPROTECT(1);
        Rf_error("Could not copy the model result from Julia to R");
    }
    UNPROTECT(1);
    return result;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_pmlr_init", (DL_FUNC)&C_pmlr_init, 1},
    {"C_pmlr_fit", (DL_FUNC)&C_pmlr_fit, 12},
    {NULL, NULL, 0}
};

void R_init_PMLCFAEFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
