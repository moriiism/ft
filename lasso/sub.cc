#include "sub.h"

void LoadRootLc(string infile, string tag, double time_lo, double time_up,
                long* const nbin_lc_ptr,
                double* const xlo_lc_ptr,
                double* const xup_lc_ptr,
                double* const delta_lc_ptr,
                HistDataNerr1d* const hd1d_out)
{
    TFile* tfile = new TFile(infile.c_str());
    TH1D* th1d = static_cast<TH1D*>(tfile->Get(tag.c_str()));
    long nbin_lc  = th1d->GetNbinsX();
    double xlo_lc = th1d->GetXaxis()->GetBinLowEdge(1);
    double xup_lc = th1d->GetXaxis()->GetBinUpEdge(th1d->GetNbinsX());
    double delta_lc = (xup_lc - xlo_lc) / nbin_lc;

    double xlo_hd1d = time_lo;
    double xup_hd1d = time_up;
    long nbin_hd1d  = (long) floor( ( (xup_hd1d - xlo_hd1d) + delta_lc * 0.1 ) / delta_lc );
    hd1d_out->Init(nbin_hd1d, xlo_hd1d, xup_hd1d);
    for(long ibin = 0; ibin < nbin_hd1d; ibin ++){
        double xval = xlo_hd1d + delta_lc * (ibin + 0.5);
        long ibin_lc = (long) floor( (xval - xlo_lc) / delta_lc );
        hd1d_out->SetOvalElm(ibin, th1d->GetBinContent(ibin_lc + 1));
    }
    // output light curve
    //char qdplc[kLineSize];
    //sprintf(qdplc, "%s/%s_lc.qdp",
    // argval->GetOutdir().c_str(),
    //        argval->GetOutfileHead().c_str());
    //MirQdpTool::MkQdp(hd1d, qdplc, "x,y");
    
    *nbin_lc_ptr = nbin_lc;
    *xlo_lc_ptr = xlo_lc;
    *xup_lc_ptr = xup_lc;
    *delta_lc_ptr = delta_lc;
}

// normalize and standardize for Gd2d
void Gd2dNormAndStd(const GraphData2d* const gd2d,
                    GraphDataNerr2d* const gd2d_out)
{
    // normalize and standardize
    //    mean = mean(h.vec)
    //    sd   = sd(h.vec)
    //    h.vec = h.vec - mean
    //    h.vec = h.vec / sd
    double mean   = MirMath::GetAMean(gd2d->GetNdata(), gd2d->GetOvalArr()->GetVal());
    double stddev = MirMath::GetStddev(gd2d->GetNdata(), gd2d->GetOvalArr()->GetVal());

    DataArrayNerr1d* da1d_tmp1 = new DataArrayNerr1d;
    DataArrayNerr1d* da1d_tmp2 = new DataArrayNerr1d;
    DataArray1dOpe::GetScale(gd2d->GetOvalArr(), 1.0, -1 * mean, da1d_tmp1);
    DataArray1dOpe::GetScale(da1d_tmp1, 1.0 / stddev, 0.0, da1d_tmp2);

    gd2d_out->Init(gd2d->GetNdata());
    gd2d_out->SetXvalArr(gd2d->GetXvalArr());
    gd2d_out->SetOvalArr(da1d_tmp2);
    delete da1d_tmp1;
    delete da1d_tmp2;
}

void GenMatLasso(const double* const time_arr, long ntime,
                 double freq_lo, double delta_freq, long nfreq,
                 double** const mat_arr_ptr)
{
    long nrow = ntime;
    long ncol = 2 * nfreq;
    
    // lasso
    // matrix
    // (nrow, ncol)
    double* mat_arr = new double [nrow * ncol];
    for(long irow = 0; irow < nrow; irow ++){
        for(long icol = 0; icol < ncol; icol ++){
            long ibin = icol * nrow + irow;
            mat_arr[ibin] = 0.0;
        }
    }
    for(long irow = 0; irow < nrow; irow ++){
        double time = time_arr[irow];
        for(long icol = 0; icol < ncol/2; icol ++){
            long ibin = icol * nrow + irow;
            long ibin2 = (icol + ncol/2) * nrow + irow;
            double freq = freq_lo + delta_freq * (icol + 0.5);
            mat_arr[ibin]  = 2 * delta_freq * cos(2 * M_PI * freq * time);
            mat_arr[ibin2] = 2 * delta_freq * sin(2 * M_PI * freq * time);
        }
    }
    *mat_arr_ptr = mat_arr;
}

// OK
void GetProxMap(const double* const A_mat_arr, long nrow, long ncol,
                const double* const h_arr,  // nrow
                const double* const y_arr,  // ncol
                double lconst,
                double lambda,
                double* const pL_arr)  // ncol
{
    for(long icol = 0; icol < ncol; icol ++){
        pL_arr[icol] = 0.0;
    }
    double* tmp_arr = new double[ncol];
    DiffFFunc(A_mat_arr, nrow, ncol,
              h_arr, y_arr, tmp_arr);
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
double QFunc(const double* const A_mat_arr, long nrow, long ncol,
             const double* const h_arr,  // nrow
             const double* const x_arr,  // ncol
             const double* const y_arr,  // ncol
             double lconst,
             double lambda)
{
    double term1 = FFunc(A_mat_arr, nrow, ncol, h_arr, y_arr);
    double* x_y_arr = new double [ncol];
    for(long icol = 0; icol < ncol; icol ++){
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
void DiffFFunc(const double* const A_mat_arr, long nrow, long ncol,
               const double* const h_arr,  // nrow
               const double* const x_arr,  // ncol
               double* const out_arr)      // ncol
{
    double* tmp_arr = new double[nrow];
    for(long irow = 0; irow < nrow; irow ++){
        tmp_arr[irow] = h_arr[irow];
    }
    for(long icol = 0; icol < ncol; icol ++){
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
           0.0, out_arr, 1);
    delete [] tmp_arr;
    delete [] trans1;
    delete [] trans2;
}

// OK
double FFunc(const double* const A_mat_arr, long nrow, long ncol,
             const double* const h_arr,  // nrow
             const double* const x_arr)  // ncol
{
    double* tmp_arr = new double[nrow];
    for(long irow = 0; irow < nrow; irow ++){
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
double GFunc(const double* const x_arr, long ncol, double lambda)
{
    double ans = 0.0;
    for(long icol = 0; icol < ncol/2; icol ++){
        ans += sqrt( pow(x_arr[icol], 2) + pow(x_arr[icol + ncol/2], 2) );
    }
    ans *= lambda;
    return(ans);
}

// OK
double FGFunc(const double* const A_mat_arr, long nrow, long ncol,
              const double* const h_arr, // nrow
              const double* const x_arr, // ncol
              double lambda)
{
    double ans = 0.0;
    ans = FFunc(A_mat_arr, nrow, ncol, h_arr, x_arr)
        + GFunc(x_arr, ncol, lambda);
    return(ans);
}

// OK
double GetLconst(const double* const A_mat_arr, long nrow, long ncol,
                 const double* const h_arr,  // nrow
                 const double* const x_arr,  // ncol
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


// ------------------------------------------

void GenFFT(long nbin, const double* const val_arr,
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

    // plan = fftw_plan_dft_1d(nbin, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan = fftw_plan_dft_1d(nbin, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(long ibin = 0; ibin < nbin; ibin ++){
        in[ibin][0] = 0.0;
        in[ibin][1] = 0.0;
        out[ibin][0] = 0.0;
        out[ibin][1] = 0.0;
    }

    for(long ibin = 0; ibin < nbin; ibin ++){
        in[ibin][0] = val_arr[ibin];
        in[ibin][1] = 0.0;
    }
  
    fftw_execute(plan);
    double* output    = new double [2 * (nbin_half + 1)];
    double* out_real  = new double [nbin_half + 1];
    double* out_image = new double [nbin_half + 1];
  
    for(long ibin = 0; ibin < nbin_half + 1; ibin++ ){
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
        // printf("%ld %e\n", ibin, power[ibin]);
    }
    freq[nbin_half] = 1.0 / 2.0;
    power[nbin_half] = (pow(out_real[nbin_half], 2) + pow(out_image[nbin_half], 2)) / nbin_pow2;

    delete [] output;
    delete [] out_real;
    delete [] out_image;
    delete [] freq;

    *out_real_ptr = out_real;
    *out_image_ptr = out_image;
    *power_ptr = power;
}



