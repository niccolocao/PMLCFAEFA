#include <dlfcn.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

typedef void (*init_julia_fn)(int, char **);
typedef int32_t (*fit_fn)(const int64_t *, const double *, int32_t, int32_t, int32_t,
                          int32_t, int32_t, int32_t, double, int32_t, int32_t,
                          int32_t, double, double, const double *, int32_t, int32_t);
typedef int32_t (*length_fn)(void);
typedef int32_t (*copy_fn)(double *, int32_t);
typedef int32_t (*error_fn)(unsigned char *, int32_t);

static void *require_symbol(void *lib, const char *symbol) {
    void *ptr = dlsym(lib, symbol);
    if (!ptr) {
        fprintf(stderr, "Missing Julia symbol %s: %s\n", symbol, dlerror());
        exit(2);
    }
    return ptr;
}
static void *load_library(const char *root, const char *relative) {
    char path[PATH_MAX];
    int n = snprintf(path, sizeof path, "%s/%s", root, relative);
    if (n < 0 || (size_t)n >= sizeof path) {
        fprintf(stderr, "Julia library path too long\n");
        exit(2);
    }
    void *lib = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
    if (!lib) {
        fprintf(stderr, "Cannot open %s: %s\n", path, dlerror());
        exit(2);
    }
    return lib;
}
int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: pmlr-runner BUNDLE_DIR INPUT_FILE OUTPUT_FILE\n");
        return 2;
    }
    FILE *in = fopen(argv[2], "r");
    if (!in) { perror("Cannot read input file"); return 2; }
    int32_t n, p, q, k, c, ncov, maxit, maxrest, replix, has_l;
    double tol, fix, latent;
    if (fscanf(in, " %" SCNd32 " %" SCNd32 " %" SCNd32 " %" SCNd32 " %" SCNd32
               " %" SCNd32 " %" SCNd32 " %" SCNd32 " %" SCNd32 " %" SCNd32
               " %lf %lf %lf", &n, &p, &q, &k, &c, &ncov, &maxit, &maxrest,
               &replix, &has_l, &tol, &fix, &latent) != 13 || n < 2 || p < 2 ||
               ncov < 1 || q < 1 || k < 1 || c < 2 || n > 1000000 || p > 10000 ||
               ncov > 10000 || maxit < 1 || maxrest < 1 || replix < 1 ||
               has_l < 0 || has_l > 1) {
        fprintf(stderr, "Malformed or unsupported input header\n");
        fclose(in);
        return 2;
    }
    size_t xlen = (size_t)n * (size_t)p;
    size_t covlen = (size_t)n * (size_t)ncov;
    int64_t *x = malloc(xlen * sizeof *x);
    double *cov = malloc(covlen * sizeof *cov);
    double *loading = has_l ? malloc((size_t)p * (size_t)q * sizeof *loading) : NULL;
    if (!x || !cov || (has_l && !loading)) {
        fprintf(stderr, "Could not allocate input buffers\n");
        fclose(in);
        return 2;
    }
    for (size_t i = 0; i < xlen; i++) {
        if (fscanf(in, " %" SCNd64, &x[i]) != 1) {
            fprintf(stderr, "Malformed integer data at index %zu\n", i);
            fclose(in);
            return 2;
        }
    }
    for (size_t i = 0; i < covlen; i++) {
        if (fscanf(in, " %lf", &cov[i]) != 1) {
            fprintf(stderr, "Malformed covariates at index %zu\n", i);
            fclose(in);
            return 2;
        }
    }
    if (has_l) {
        for (size_t i = 0; i < (size_t)p * (size_t)q; i++) {
            if (fscanf(in, " %lf", &loading[i]) != 1) {
                fprintf(stderr, "Malformed loading structure at index %zu\n", i);
                fclose(in);
                return 2;
            }
        }
    }
    fclose(in);
    fprintf(stderr, "[pmlr] parsed inputs n=%d p=%d q=%d K=%d\n", n, p, q, k); fflush(stderr);
    fprintf(stderr, "[pmlr] loading Julia libraries\n"); fflush(stderr);
    void *runtime = load_library(argv[1], "lib/libjulia.so");
    void *model = load_library(argv[1], "lib/libpmlr.so");
    init_julia_fn init_julia = (init_julia_fn)require_symbol(model, "init_julia");
    fit_fn fit = (fit_fn)require_symbol(model, "pmlr_fit");
    length_fn length = (length_fn)require_symbol(model, "pmlr_result_length");
    copy_fn copy = (copy_fn)require_symbol(model, "pmlr_copy_result");
    error_fn last_error = (error_fn)require_symbol(model, "pmlr_last_error");
    char julia_name[] = "pmlr-runner";
    char *julia_argv[] = {julia_name, NULL};
    fprintf(stderr, "[pmlr] initializing Julia\n"); fflush(stderr);
    init_julia(1, julia_argv);
    fprintf(stderr, "[pmlr] Julia initialization returned\n"); fflush(stderr);
    typedef void *(*julia_eval_fn)(const char *);
    typedef void *(*julia_exception_fn)(void);
    julia_eval_fn eval_julia = (julia_eval_fn)dlsym(runtime, "jl_eval_string");
    julia_exception_fn julia_exception = (julia_exception_fn)dlsym(runtime, "jl_exception_occurred");
    if (!eval_julia || !julia_exception) {
        fprintf(stderr, "Julia C API is missing evaluation functions\n");
        return 3;
    }
    eval_julia("try; let mods = values(Base.loaded_modules); opt = first(m for m in mods if nameof(m) === :Optim); nls = first(m for m in mods if nameof(m) === :NLSolversBase); if !isdefined(opt, :only_fg!); f = getfield(nls, :only_fg!); Core.eval(opt, :(const only_fg! = $f)); end; end; catch e; Base.showerror(stderr, e, catch_backtrace()); println(stderr); end");
    eval_julia("any(m -> nameof(m) === :Optim && isdefined(m, :only_fg!), values(Base.loaded_modules)) || error(\"Optim.only_fg! compatibility shim did not register\")");
    if (julia_exception() != NULL) {
        fprintf(stderr, "Cannot activate Optim/NLSolversBase compatibility alias\n");
        return 3;
    }
    fprintf(stderr, "[pmlr] Optim compatibility alias ready\n"); fflush(stderr);
    fprintf(stderr, "[pmlr] invoking Julia pmlr_fit\n"); fflush(stderr);
    int32_t code = fit(x, cov, n, p, q, k, c, ncov, tol, maxit, maxrest,
                       replix, fix, latent, loading, has_l, 0);
    if (code != 0) {
        unsigned char message[16384] = {0};
        last_error(message, (int32_t)sizeof message);
        fprintf(stderr, "Julia model returned an error: %s\n", message);
        return 3;
    }
    int32_t len = length();
    if (len <= 0 || len > 200000000) {
        fprintf(stderr, "Unexpected result length: %" PRId32 "\n", len);
        return 3;
    }
    double *result = malloc((size_t)len * sizeof *result);
    if (!result || copy(result, len) != 0) {
        fprintf(stderr, "Cannot copy Julia result\n");
        return 3;
    }
    FILE *out = fopen(argv[3], "w");
    if (!out) { perror("Cannot create output file"); return 3; }
    fprintf(out, "%" PRId32 "\n", len);
    for (int32_t i = 0; i < len; i++) fprintf(out, "%.17g\n", result[i]);
    fclose(out);
    free(result);
    free(loading);
    free(cov);
    free(x);
    return 0;
}
