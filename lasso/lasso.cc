#include "fftw3.h"
#include "mi_iolib.h"
#include "mir_hist1d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_lasso.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

//
// unit: sec, Hz
//

// delta_time = 0.1 msec = 1e-4 sec
// time_lo
// time_up

// --- LASSO case ----
// ntime      = N ( GTI = 1)
// freq_lo
// freq_up
// nfreq      = M / 2
// delta_freq = (freq_up - freq_lo) / nfreq
// A_mat_arr : nrow x ncol =  ntime x (2 * nfreq) = N x M
// light curve in time domain: h_arr : nsize = ntime = N
// light curve in freq domain: x_arr : nsize = 2 * nfreq = M

// --- FFT case ----
// ntime_fft = (time_up - time_lo) / delta_time
// delta_freq_fft = 1 / (ntime_fft * delta_time) = 1 / (time_up - time_lo)
// freq_lo_fft = 0
// freq_up_fft = delta_freq_fft * ntime_fft / 2 = 1 / (2 * delta_time)
//             = freq_nyquist
// nfreq_fft   = (freq_up_fft - freq_lo_fft) / delta_freq_fft = ntime_fft / 2
// A_mat_fft_arr : nrow_fft x ncol_fft = ntime_fft x (2 * nfreq_fft) = ntime_fft x ntime_fft
// light curve in time domain: h_fft_arr : nsize = ntime_fft
// light curve in freq domain: x_fft_arr : nsize = 2 * nfreq_fft = ntime_fft


int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValLasso* argval = new ArgValLasso;
    argval->Init(argc, argv);

    // make output directory
    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str(),
            argval->GetProgname().c_str());
    FILE* fp_log = fopen(logfile, "w");
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    argval->Print(fp_log);
    argval->Print(stdout);


    // load light curve 
    long nbin_lc_org  = 0;
    double xlo_lc_org = 0.0;
    double xup_lc_org = 0.0;
    double delta_lc_org = 0.0;
    HistDataNerr1d* hd1d_lc = new HistDataNerr1d;
    LoadRootLc(argval->GetInfile(), "lcurve",
               argval->GetTimeLo(),
               argval->GetTimeUp(),
               &nbin_lc_org, &xlo_lc_org, &xup_lc_org, &delta_lc_org,
               hd1d_lc);

    // load gti light curve 
    long nbin_gti_org  = 0;
    double xlo_gti_org = 0.0;
    double xup_gti_org = 0.0;
    double delta_gti_org = 0.0;
    HistDataNerr1d* hd1d_gti = new HistDataNerr1d;
    LoadRootLc(argval->GetGtifile(), "gti",
               argval->GetTimeLo(),
               argval->GetTimeUp(),
               &nbin_gti_org, &xlo_gti_org, &xup_gti_org, &delta_gti_org,
               hd1d_gti);
    MirQdpTool::MkQdp(hd1d_gti, "temp_gti.qdp", "x,y");

    
    double delta_time = hd1d_lc->GetXvalBinWidth();
    // LASSO
    int ntime = MirMath::GetSum(hd1d_gti->GetOvalArr()->GetNdata(),
                                hd1d_gti->GetOvalArr()->GetVal());
    double delta_freq = (argval->GetFreqUp() - argval->GetFreqLo()) / argval->GetNfreq();
    int nrow = ntime;
    int ncol = 2 * argval->GetNfreq();
    // FFT
    int ntime_fft = hd1d_lc->GetNbinX();
    double delta_freq_fft = 1.0 / (hd1d_lc->GetXvalUp() - hd1d_lc->GetXvalLo());
    double freq_lo_fft = 0.0;
    double freq_up_fft = 1.0 / (2 * delta_time);
    int nfreq_fft = ntime_fft / 2;
    int nrow_fft = ntime_fft;
    int ncol_fft = ntime_fft;

    printf("delta_time = %e\n", delta_time);
    printf("ntime = %d\n", ntime);
    printf("delta_freq = %e\n", delta_freq);
    printf("nrow = %d\n", nrow);
    printf("ncol = %d\n", ncol);
    printf("ntime_fft = %d\n", ntime_fft);
    printf("delta_freq_fft = %e\n", delta_freq_fft);
    printf("freq_lo_fft = %e\n", freq_lo_fft);
    printf("freq_up_fft = %e\n", freq_up_fft);
    printf("nfreq_fft = %d\n", nfreq_fft);
    printf("nrow_fft = %d\n", nrow_fft);
    printf("ncol_fft = %d\n", ncol_fft);
    

    // light curve of gti = 1
    GraphDataNerr2d* gd2d_lc = new GraphDataNerr2d;
    int itime = 0; 
    gd2d_lc->Init(ntime);
    for(int ibin = 0; ibin < hd1d_lc->GetNbinX(); ibin ++){
        if(1 == hd1d_gti->GetOvalElm(ibin)){
            gd2d_lc->SetPoint(itime,
                              hd1d_lc->GetHi1d()->GetBinCenter(ibin),
                              hd1d_lc->GetOvalElm(ibin));
            itime ++;
        }
    }

    printf("nbin ntime = %d %ld\n",  hd1d_lc->GetNbinX(), ntime);
    
    GraphDataNerr2d* gd2d_norm_lc = new GraphDataNerr2d;
    Gd2dNormAndStd(gd2d_lc, gd2d_norm_lc);
    MirQdpTool::MkQdp(gd2d_lc, "lc.qdp", "x,y");
    MirQdpTool::MkQdp(gd2d_norm_lc, "lc_norm.qdp", "x,y");

//    abort();
    

//    // -------------------------------------
//
//    long nbin = hd1d_lc->GetNbinX();
//    long nbin_half = nbin / 2;
//    
//    fftw_complex* in = NULL;
//    fftw_complex* out = NULL;
//    fftw_plan plan = NULL;
//    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);
//    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);
//
//    plan = fftw_plan_dft_1d(nbin, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//
//    for(long ibin = 0; ibin < nbin; ibin ++){
//        in[ibin][0] = 0.0;
//        in[ibin][1] = 0.0;
//        out[ibin][0] = 0.0;
//        out[ibin][1] = 0.0;
//    }
//
//    for(long ibin = 0; ibin < nbin; ibin ++){
//        in[ibin][0] = hd1d_lc->GetOvalElm(ibin);
//        in[ibin][1] = 0.0;
//    }
//  
//    fftw_execute(plan);
//    double* output    = new double [2 * (nbin_half + 1)];
//    double* out_real  = new double [nbin_half + 1];
//    double* out_image = new double [nbin_half + 1];
//  
//    for(long ibin = 0; ibin < nbin_half + 1; ibin++ ){
//        out_real[ibin]  = out[ibin][0];
//        out_image[ibin] = out[ibin][1];
//    }
//
//    fftw_destroy_plan(plan);
//    fftw_free(in);
//    fftw_free(out);
//
//
//    // calculate Power Spectrum
//    double* freq  = new double [nbin_half + 1];
//    double* power = new double [nbin_half + 1];
//    double nbin_pow2 = pow(nbin, 2);
//    freq[0] = 0.0;
//    power[0] = (pow(out_real[0], 2) + pow(out_image[0], 2)) / nbin_pow2;
//    for(long ibin = 1; ibin < nbin_half; ibin++){
//        freq[ibin] = ibin / nbin;
//        power[ibin] = 2.0 * (pow(out_real[ibin], 2) + pow(out_image[ibin], 2)) / nbin_pow2;
//
//        printf("%ld %e\n", ibin, power[ibin]);
//    }
//    freq[nbin_half] = 1.0 / 2.0;
//    power[nbin_half] = (pow(out_real[nbin_half], 2) + pow(out_image[nbin_half], 2)) / nbin_pow2;
//
//
//    // -----------------------------------

    double* A_mat_arr = NULL;
    GenMatLasso(gd2d_norm_lc->GetXvalArr()->GetVal(), ntime,
                argval->GetFreqLo(), delta_freq, argval->GetNfreq(),
                &A_mat_arr);
    double* h_arr = new double [nrow];
    for(int irow = 0; irow < nrow; irow++){
        h_arr[irow] = gd2d_norm_lc->GetOvalElm(irow);
    }
    double* x_arr = new double [ncol];
    double* y_arr = new double [ncol];
    for(int icol = 0; icol < ncol; icol++){
        x_arr[icol] = 0.0;
        y_arr[icol] = 0.0;
    }
    
    double tolerance = 1.0e-10;
    double eta = 1.2;
    double lconst = 1.0e-3;
    double lconst_pre = lconst;
    double kiter_max = 500;
    double cost = 0.0;
    double cost_pre = cost;
    for(int kiter = 0; kiter < kiter_max; kiter ++){
        printf("kiter = %d, cost = %e, lconst = %e\n", kiter, cost, lconst);
        lconst = GetLconst(A_mat_arr, nrow, ncol,
                           h_arr, y_arr,
                           lconst, eta, argval->GetLambda());
        GetProxMap(A_mat_arr, nrow, ncol,
                   h_arr,
                   y_arr,
                   lconst,
                   argval->GetLambda(),
                   x_arr);
        cost = FGFunc(A_mat_arr, nrow, ncol,
                      h_arr, x_arr, argval->GetLambda());

        printf("ratio = %e\n", (cost_pre - cost) / cost);
        
        if(kiter > 1 && fabs((cost_pre - cost) / cost) < tolerance){
            printf("kiter = %d, cost = %e\n", kiter, cost);
            break;
        }
        dcopy_(ncol, x_arr, 1, y_arr, 1);
        cost_pre = cost;
    }

    HistDataNerr1d* hd1d_lc_in_freq = new HistDataNerr1d;
    hd1d_lc_in_freq->Init(2 * argval->GetNfreq(), argval->GetFreqLo(), argval->GetFreqUp());
    hd1d_lc_in_freq->SetOvalArr(2 * argval->GetNfreq(), x_arr);
    MirQdpTool::MkQdp(hd1d_lc_in_freq, "lc_in_freq.qdp", "x,y");
    
    delete argval;
    
    return status;
}

