#include "sub_fft.h"
#include "sub.h"

void GetProxMapFft(const double* const h_arr, long nrow,
                   const int* const index_arr, // nrow
                   const double* const y_arr, long ncol,
                   double delta_freq,
                   double lconst,
                   double lambda,
                   double* const pL_arr)  // ncol
{
    double* tmp_arr = new double[ncol];
    DiffFFuncFft(h_arr, nrow,
                 index_arr,
                 y_arr, ncol,
                 delta_freq,
                 tmp_arr);
    double* z_arr = new double[ncol];
    for(long icol = 0; icol < ncol; icol ++){
        z_arr[icol] = y_arr[icol] - tmp_arr[icol] / lconst;
    }
    for(long icol = 0; icol < ncol/2; icol ++){
        double xval = lconst * sqrt( pow(z_arr[icol], 2) + pow(z_arr[icol + ncol/2], 2) );
        double factor = GetSoftThres_OverX(xval, lambda);
        pL_arr[icol]          = factor * z_arr[icol];
        pL_arr[icol + ncol/2] = factor * z_arr[icol + ncol/2];
    }
    delete [] tmp_arr;
    delete [] z_arr;
}

double QFuncFft(const double* const h_arr, long nrow,
                const int* const index_arr, // nrow
                const double* const x_arr, long ncol,
                const double* const y_arr, // ncol
                double delta_freq,
                double lconst,
                double lambda)
{
    double term1 = FFuncFft(h_arr, nrow,
                            index_arr,
                            y_arr, ncol,
                            delta_freq);
    double* x_y_arr = new double [ncol];
    for(long icol = 0; icol < ncol; icol ++){
        x_y_arr[icol] = x_arr[icol];
    }
    daxpy_(ncol, -1.0, const_cast<double*>(y_arr), 1, x_y_arr, 1);
    double* diff_f_arr = new double [ncol];
    DiffFFuncFft(h_arr, nrow,
                 index_arr,
                 y_arr, ncol,
                 delta_freq,
                 diff_f_arr);
    double term2 = ddot_(ncol, x_y_arr, 1, diff_f_arr, 1);
    double term3 = lconst / 2.0 * ddot_(ncol, x_y_arr, 1, x_y_arr, 1);
    double term4 = GFuncFft(x_arr, ncol, lambda);
    double ans = term1 + term2 + term3 + term4;
    delete [] x_y_arr;
    delete [] diff_f_arr;
    return(ans);
}

void DiffFFuncFft(const double* const h_arr, long nrow,
                  const int* const index_arr, // nrow
                  const double* const x_arr, long ncol,
                  double delta_freq,
                  double* const out_arr) // ncol
{
    double* tmp_arr = new double[nrow];
    for(long irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] = h_arr[irow];
    }
    AxFft(x_arr, ncol, 
          index_arr, nrow,
          delta_freq,
          tmp_arr); // nrow
    for(long irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] -= h_arr[irow];
    }
    ATyFft(tmp_arr, nrow, 
           index_arr,
           delta_freq,
           ncol,
           out_arr);
    for(long icol = 0; icol < ncol; icol ++){
        out_arr[icol] *= 2;
    }
    
    //char* trans1 = new char [1];
    //strcpy(trans1, "N");
    //dgemv_(trans1, nrow, ncol, 1.0, const_cast<double*>(A_mat_arr), nrow,
    //       const_cast<double*>(x_arr), 1,
    //       -1.0, tmp_arr, 1);
    //char* trans2 = new char [1];
    //strcpy(trans2, "T");
    //dgemv_(trans2, nrow, ncol, 2.0, const_cast<double*>(A_mat_arr), nrow,
    //       tmp_arr, 1,
    //       1.0, out_arr, 1);
    // delete [] tmp_arr;

}

double FFuncFft(const double* const h_arr, long nrow,
                const int* const index_arr, // nrow
                const double* const x_arr, long ncol,
                double delta_freq)
                
{
    // tmp.vec = A.mat %*% x.vec - h.vec
    //char* trans = new char [1];
    //strcpy(trans, "N");
    //dgemv_(trans, nrow, ncol, 1.0, const_cast<double*>(A_mat_arr), nrow,
    //       const_cast<double*>(x_arr), 1,
    //       -1.0, tmp_arr, 1);

    double* tmp_arr = new double[nrow];
    AxFft(x_arr, ncol, 
          index_arr, nrow,
          delta_freq,
          tmp_arr);
    for(long irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] -= h_arr[irow];
    }
    double ans = ddot_(nrow, tmp_arr, 1, tmp_arr, 1);
    delete [] tmp_arr;
    return(ans);
}

double GFuncFft(const double* const x_arr, long ncol, double lambda)
{
    double ans = 0.0;
    for(long icol = 0; icol < ncol/2; icol ++){
        ans += sqrt( pow(x_arr[icol], 2) + pow(x_arr[icol + ncol/2], 2) );
    }
    ans *= lambda;
    return(ans);
}

double FGFuncFft(const double* const h_arr, long nrow,
                 const int* const index_arr, // nrow
                 const double* const x_arr, long ncol,
                 double delta_freq,
                 double lambda)
{
    double ans = 0.0;
    ans = FFuncFft(h_arr, nrow,
                   index_arr,
                   x_arr, ncol,
                   delta_freq)
        + GFuncFft(x_arr, ncol, lambda);
    return(ans);
}

double GetLconstFft(const double* const h_arr, long nrow,
                    const int* const index_arr, // nrow
                    const double* const x_arr, long ncol,
                    double delta_freq,
                    double lconst,
                    double eta,
                    double lambda)
{
    double lconst_new = lconst;
    long ik_max = 1000;
    double* pL_arr = new double [ncol];
    for(long icol = 0; icol < ncol; icol ++){
        pL_arr[icol] = 0.0;
    }
    for(long ik = 0; ik < ik_max; ik ++){
        lconst_new = pow(eta, ik) * lconst;
        GetProxMapFft(h_arr, nrow,
                      index_arr,
                      x_arr, ncol,
                      delta_freq,
                      lconst_new,
                      lambda,
                      pL_arr);
        double FGFunc_val = FGFuncFft(h_arr, nrow,
                                      index_arr,
                                      pL_arr, ncol,
                                      delta_freq,
                                      lambda);
        double QFunc_val  = QFuncFft(h_arr, nrow,
                                     index_arr,
                                     pL_arr, ncol,
                                     x_arr,
                                     delta_freq,
                                     lconst_new,
                                     lambda);
        if(FGFunc_val <= QFunc_val){
            break;
        }
    }
    delete [] pL_arr;
    return(lconst_new);
}


void AxFft(const double* const x_arr, long ncol, 
           const int* const index_arr, long nrow,
           double delta_freq,
           double* const out_arr) // nrow

{
    double* real_arr = new double[ncol];
    double* image_arr = new double[ncol];
    for(long ibin = 0; ibin < ncol; ibin ++){
        real_arr[ibin] = 0.0;
        image_arr[ibin] = 0.0;
    }
    for(long ibin = 0; ibin < ncol / 2; ibin ++){
        real_arr[ibin] = x_arr[ibin];
        image_arr[ibin] = -1 * x_arr[ibin + ncol/2];
    }
    double* out_real_arr = NULL;
    double* out_image_arr = NULL;
    double* power_arr = NULL;
    GenFft(real_arr, ncol,
           image_arr,
           &out_real_arr,
           &out_image_arr,
           &power_arr);
    for(long ibin = 0; ibin < nrow; ibin ++){
        out_arr[ibin] = 2 * delta_freq * out_real_arr[index_arr[ibin]];
    }
    delete [] real_arr;
    delete [] image_arr;
    delete [] out_real_arr;
    delete [] out_image_arr;
    delete [] power_arr;
}

void ATyFft(const double* const y_arr, long nrow, 
            const int* const index_arr, // nrow
            double delta_freq,
            long ncol,
            double* const out_arr) // ncol
{
    double* real_arr = new double[ncol];
    double* image_arr = new double[ncol];
    for(long ibin = 0; ibin < ncol; ibin ++){
        real_arr[ibin] = 0.0;
        image_arr[ibin] = 0.0;
    }
    for(long ibin = 0; ibin < nrow; ibin ++){
        real_arr[index_arr[ibin]] = y_arr[ibin];
    }
    double* out_real_arr = NULL;
    double* out_image_arr = NULL;
    double* power_arr = NULL;
    GenFft(real_arr, ncol,
           image_arr,
           &out_real_arr,
           &out_image_arr,
           &power_arr);

    for(long ibin = 0; ibin < ncol / 2; ibin ++){
        out_arr[ibin] = 2 * delta_freq * out_real_arr[ibin];
        out_arr[ibin + ncol/2] = 2 * delta_freq * out_image_arr[ibin];
    }
    delete [] real_arr;
    delete [] image_arr;
    delete [] out_real_arr;
    delete [] out_image_arr;
    delete [] power_arr;
}


void GenFft(const double* const real_arr, long nbin,
            const double* const image_arr,
            double** const out_real_ptr,
            double** const out_image_ptr,
            double** const power_ptr)
{
    long nbin_half = nbin / 2;
    
    fftw_complex* in = NULL;
    fftw_complex* out = NULL;
    fftw_plan plan = NULL;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);

    plan = fftw_plan_dft_1d(nbin, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(long ibin = 0; ibin < nbin; ibin ++){
        in[ibin][0] = 0.0;
        in[ibin][1] = 0.0;
        out[ibin][0] = 0.0;
        out[ibin][1] = 0.0;
    }

    for(long ibin = 0; ibin < nbin; ibin ++){
        in[ibin][0] = real_arr[ibin];
        in[ibin][1] = image_arr[ibin];
    }
  
    fftw_execute(plan);
    double* out_real  = new double [nbin];
    double* out_image = new double [nbin];
    for(long ibin = 0; ibin < nbin; ibin++ ){
        out_real[ibin]  = out[ibin][0];
        out_image[ibin] = out[ibin][1];
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    // calculate Power Spectrum
    double* freq  = new double [nbin_half + 1];
    double* power = new double [nbin_half + 1];
    double nbin_pow2 = pow(nbin, 2);
    freq[0] = 0.0;
    power[0] = (pow(out_real[0], 2) + pow(out_image[0], 2)) / nbin_pow2;
    for(long ibin = 1; ibin < nbin_half; ibin++){
        freq[ibin] = ibin / nbin;
        power[ibin] = 2.0 * (pow(out_real[ibin], 2) + pow(out_image[ibin], 2)) / nbin_pow2;
    }
    freq[nbin_half] = 1.0 / 2.0;
    power[nbin_half] = (pow(out_real[nbin_half], 2) + pow(out_image[nbin_half], 2)) / nbin_pow2;

    delete [] freq;

    *out_real_ptr = out_real;
    *out_image_ptr = out_image;
    *power_ptr = power;
}



