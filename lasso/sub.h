#ifndef MORIIISM_FT_LASSO_SUB_H_
#define MORIIISM_FT_LASSO_SUB_H_

#include "fftw3.h"
#include "mib_blas.h"
#include "mir_hist1d_nerr.h"
#include "mir_hist1d_ope.h"

void LoadRootLc(string infile, string tag, double time_lo, double time_up,
                long* const nbin_lc_ptr,
                double* const xlo_lc_ptr,
                double* const xup_lc_ptr,
                double* const delta_lc_ptr,
                HistDataNerr1d* const hd1d_out);

void Gd2dNormAndStd(const GraphData2d* const gd2d,
                    GraphDataNerr2d* const gd2d_out);

void GenMatLasso(const double* const time_arr, long ntime,
                 double freq_lo, double delta_freq, long nfreq,
                 double** const mat_arr_ptr);

void GetProxMap(const double* const A_mat_arr, long nrow, long ncol,
                const double* const h_arr,
                const double* const y_arr,
                double lconst,
                double lambda,
                double* const pL_arr);

double GetSoftThres_OverX(double val, double lambda);

double QFunc(const double* const A_mat_arr, long nrow, long ncol,
             const double* const h_arr,
             const double* const x_arr,
             const double* const y_arr,
             double lconst,
             double lambda);

void DiffFFunc(const double* const A_mat_arr, long nrow, long ncol,
               const double* const h_arr,
               const double* const x_arr,
               double* const out_arr);

double FFunc(const double* const A_mat_arr, long nrow, long ncol,
             const double* const h_arr,
             const double* const x_arr);

double GFunc(const double* const x_arr, long ncol, double lambda);

double FGFunc(const double* const A_mat_arr, long nrow, long ncol,
              const double* const h_arr,
              const double* const x_arr,
              double lambda);

double GetLconst(const double* const A_mat_arr, long nrow, long ncol,
                 const double* const h_arr,
                 const double* const x_arr,
                 double lconst,
                 double eta,
                 double lambda);

void GenFFT(long nbin, const double* const val_arr,
            double** const out_real_ptr,
            double** const out_image_ptr,
            double** const power_ptr);

#endif // MORIIISM_FT_LASSO_SUB_H_
