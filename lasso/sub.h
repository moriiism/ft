#ifndef MORIIISM_FT_LASSO_SUB_H_
#define MORIIISM_FT_LASSO_SUB_H_

//#include "mi_rand.h"
//#include "mi_sort.h"
//#include "mi_time.h"
//#include "mif_fits.h"
//#include "mir_math.h"
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

void GenMatLasso(const double* const time_arr, int ntime,
                 double freq_lo, double delta_freq, int nfreq,
                 double** const mat_arr_ptr);

void GetProxMap(const double* const A_mat_arr, int nrow, int ncol,
                const double* const h_arr,
                const double* const y_arr,
                double lconst,
                double lambda,
                double* const pL_arr);

double GetSoftThres_OverX(double val, double lambda);

double QFunc(const double* const A_mat_arr, int nrow, int ncol,
             const double* const h_arr,
             const double* const x_arr,
             const double* const y_arr,
             double lconst,
             double lambda);

void DiffFFunc(const double* const A_mat_arr, int nrow, int ncol,
               const double* const h_arr,
               const double* const x_arr,
               double* const out_arr);

double FFunc(const double* const A_mat_arr, int nrow, int ncol,
             const double* const h_arr,
             const double* const x_arr);

double GFunc(const double* const x_arr, int ncol, double lambda);

double FGFunc(const double* const A_mat_arr, int nrow, int ncol,
              const double* const h_arr,
              const double* const x_arr,
              double lambda);

double GetLconst(const double* const A_mat_arr, int nrow, int ncol,
                 const double* const h_arr,
                 const double* const y_arr,
                 double lconst,
                 double eta,
                 double lambda);



#endif // MORIIISM_FT_LASSO_SUB_H_
