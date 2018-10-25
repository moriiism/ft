#ifndef MORIIISM_FT_LASSO_SUB_FFT_CV_H_
#define MORIIISM_FT_LASSO_SUB_FFT_CV_H_

#include "fftw3.h"
#include "mib_blas.h"
#include "mir_hist1d_nerr.h"
#include "mir_hist1d_ope.h"

void GetLassoFft(const double* const h_arr, long nrow,
                 const long* const index_arr,
                 long ncol,
                 double delta_freq,
                 long nfreq_fft,
                 double freq_lo_fft,
                 double freq_up_fft,
                 double lambda,
                 double* const out_arr);

#endif // MORIIISM_FT_LASSO_SUB_FFT_CV_H_
