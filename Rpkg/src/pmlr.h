#ifndef PMLR_H
#define PMLR_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
int32_t pmlr_fit(const int64_t *x, const double *covariates,
    int32_t n, int32_t p, int32_t q, int32_t k, int32_t c, int32_t ncov,
    double tol, int32_t maxit, int32_t maxrest, int32_t replix,
    double fixthresholds, double latentmeans, const double *loading, int32_t has_loading);
int32_t pmlr_result_length(void);
int32_t pmlr_copy_result(double *destination, int32_t capacity);
int32_t pmlr_last_error(unsigned char *destination, int32_t capacity);
#ifdef __cplusplus
}
#endif
#endif
