#include "sub_fft.h"

void GetProxMap(const double* const A_mat_arr, int nrow, int ncol,
                const double* const h_arr,  // nrow
                const double* const y_arr,  // ncol
                double lconst,
                double lambda,
                double* const pL_arr)  // ncol
{
    double* tmp_arr = new double[ncol];
    DiffFFunc(A_mat_arr, nrow, ncol,
              h_arr, y_arr, tmp_arr);
    double* z_arr = new double[ncol];
    for(int icol = 0; icol < ncol; icol ++){
        z_arr[icol] = y_arr[icol] - tmp_arr[icol] / lconst;
    }
    for(int icol = 0; icol < ncol/2; icol ++){
        double xval = lconst * sqrt( pow(z_arr[icol], 2) + pow(z_arr[icol + ncol/2], 2) );
        double factor = GetSoftThres_OverX(xval, lambda);
        pL_arr[icol]          = factor * z_arr[icol];
        pL_arr[icol + ncol/2] = factor * z_arr[icol + ncol/2];
    }
    delete [] tmp_arr;
    delete [] z_arr;
}

// OK
double GetSoftThres_OverX(double val, double lambda)
{
    double ans = 0.0;
    if(lambda < val){
        ans = 1.0 - lambda / val;
    } else if(-1 * lambda <= val && val <= lambda){
        ans = 0.0;
    } else if(val < -1 * lambda){
        ans = 1.0 + lambda / val;
    }
    return(ans);
}

// OK
double QFunc(const double* const A_mat_arr, int nrow, int ncol,
             const double* const h_arr,  // nrow
             const double* const x_arr,  // ncol
             const double* const y_arr,  // ncol
             double lconst,
             double lambda)
{
    double term1 = FFunc(A_mat_arr, nrow, ncol, h_arr, y_arr);
    double* x_y_arr = new double [ncol];
    for(int icol = 0; icol < ncol; icol ++){
        x_y_arr[icol] = x_arr[icol];
    }
    daxpy_(ncol, -1.0, const_cast<double*>(y_arr), 1, x_y_arr, 1);
    double* diff_f_arr = new double [ncol];
    DiffFFunc(A_mat_arr, nrow, ncol,
              h_arr,
              y_arr,
              diff_f_arr);
    double term2 = ddot_(ncol, x_y_arr, 1, diff_f_arr, 1);
    double term3 = lconst / 2.0 * ddot_(ncol, x_y_arr, 1, x_y_arr, 1);
    double term4 = GFunc(x_arr, ncol, lambda);
    double ans = term1 + term2 + term3 + term4;
    delete [] x_y_arr;
    delete [] diff_f_arr;
    return(ans);
}

// OK
void DiffFFunc(const double* const A_mat_arr, int nrow, int ncol,
               const double* const h_arr,  // nrow
               const double* const x_arr,  // ncol
               double* const out_arr)      // ncol
{
    double* tmp_arr = new double[nrow];
    for(int irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] = h_arr[irow];
    }
    for(int icol = 0; icol < ncol; icol ++){
        out_arr[icol] = 0.0;
    }
    char* trans1 = new char [1];
    strcpy(trans1, "N");
    dgemv_(trans1, nrow, ncol, 1.0, const_cast<double*>(A_mat_arr), nrow,
           const_cast<double*>(x_arr), 1,
           -1.0, tmp_arr, 1);
    char* trans2 = new char [1];
    strcpy(trans2, "T");
    dgemv_(trans2, nrow, ncol, 2.0, const_cast<double*>(A_mat_arr), nrow,
           tmp_arr, 1,
           1.0, out_arr, 1);
    delete [] tmp_arr;
    delete [] trans1;
    delete [] trans2;
}

// OK
double FFunc(const double* const A_mat_arr, int nrow, int ncol,
             const double* const h_arr,  // nrow
             const double* const x_arr)  // ncol
{
    double* tmp_arr = new double[nrow];
    for(int irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] = h_arr[irow];
    }
    // tmp.vec = A.mat %*% x.vec - h.vec
    char* trans = new char [1];
    strcpy(trans, "N");
    dgemv_(trans, nrow, ncol, 1.0, const_cast<double*>(A_mat_arr), nrow,
           const_cast<double*>(x_arr), 1,
           -1.0, tmp_arr, 1);
    double ans = ddot_(nrow, tmp_arr, 1, tmp_arr, 1);
    delete [] tmp_arr;
    delete [] trans;
    return(ans);
}

// OK
double GFunc(const double* const x_arr, int ncol, double lambda)
{
    double ans = 0.0;
    for(int icol = 0; icol < ncol/2; icol ++){
        ans += sqrt( pow(x_arr[icol], 2) + pow(x_arr[icol + ncol/2], 2) );
    }
    ans *= lambda;
    return(ans);
}

// OK
double FGFunc(const double* const A_mat_arr, int nrow, int ncol,
              const double* const h_arr, // nrow
              const double* const x_arr, // ncol
              double lambda)
{
    double ans = 0.0;
    ans = FFunc(A_mat_arr, nrow, ncol, h_arr, x_arr)
        + GFunc(x_arr, ncol, lambda);
    return(ans);
}

double GetLconstFFT(int nrow, int ncol,
                    const double* const h_arr,  // nrow
                    const double* const x_arr,  // ncol
                    double lconst,
                    double eta,
                    double lambda)
{
    double lconst_new = lconst;
    int ik_max = 1000;
    double* pL_arr = new double [ncol];
    for(int icol = 0; icol < ncol; icol ++){
        pL_arr[icol] = 0.0;
    }
    for(int ik = 0; ik < ik_max; ik ++){
        lconst_new = pow(eta, ik) * lconst;
        GetProxMap(A_mat_arr, nrow, ncol,
                   h_arr, x_arr, lconst_new, lambda,
                   pL_arr);
        double FGFunc_val = FGFunc(A_mat_arr, nrow, ncol,
                                   h_arr, pL_arr, lambda);
        double QFunc_val  = QFunc(A_mat_arr, nrow, ncol,
                                  h_arr, pL_arr, x_arr, lconst_new, lambda);
        if(FGFunc_val <= QFunc_val){
            break;
        }
    }
    delete [] pL_arr;
    return(lconst_new);
}






