#ifndef MORIIISM_FT_LASSO_SUB_FFT_H_
#define MORIIISM_FT_LASSO_SUB_FFT_H_

#include "fftw3.h"
#include "mib_blas.h"
#include "mir_hist1d_nerr.h"
#include "mir_hist1d_ope.h"

void GetProxMapFft(const double* const h_arr, long nrow,
                   const long* const index_arr, // nrow
                   const double* const y_arr, long ncol,
                   double delta_freq,
                   double lconst,
                   double lambda,
                   double* const pL_arr);  // ncol

double QFuncFft(const double* const h_arr, long nrow,
                const long* const index_arr, // nrow
                const double* const x_arr, long ncol,
                const double* const y_arr, // ncol
                double delta_freq,
                double lconst,
                double lambda);

void DiffFFuncFft(const double* const h_arr, long nrow,
                  const long* const index_arr, // nrow
                  const double* const x_arr, long ncol,
                  double delta_freq,
                  double* const out_arr); // ncol

double FFuncFft(const double* const h_arr, long nrow,
                const long* const index_arr, // nrow
                const double* const x_arr, long ncol,
                double delta_freq);

double GFuncFft(const double* const x_arr, long ncol, double lambda);

double FGFuncFft(const double* const h_arr, long nrow,
                 const long* const index_arr, // nrow
                 const double* const x_arr, long ncol,
                 double delta_freq,
                 double lambda);

double GetLconstFft(const double* const h_arr, long nrow,
                    const long* const index_arr, // nrow
                    const double* const x_arr, long ncol,
                    double delta_freq,
                    double lconst,
                    double eta,
                    double lambda);

void AxFft(const double* const x_arr, long ncol, 
           const long* const index_arr, long nrow,
           double delta_freq,
           double* const out_arr); // nrow

void ATyFft(const double* const y_arr, long nrow, 
            const long* const index_arr, // nrow
            double delta_freq,
            long ncol,
            double* const out_arr); // ncol

void GenFft(const double* const real_arr, long nbin,
            const double* const image_arr,
            double** const out_real_ptr,
            double** const out_image_ptr,
            double** const power_ptr);

#endif // MORIIISM_FT_LASSO_SUB_FFT_H_
