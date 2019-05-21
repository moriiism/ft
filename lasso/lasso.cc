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
    MirQdpTool::MkQdp(hd1d_lc, "hd1d_lc.qdp", "x,y");
    
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
    MirQdpTool::MkQdp(hd1d_gti, "hd1d_gti.qdp", "x,y");

    
    double delta_time = hd1d_lc->GetXvalBinWidth();
    // LASSO
    long ntime = MirMath::GetSum(hd1d_gti->GetOvalArr()->GetNdata(),
                                 hd1d_gti->GetOvalArr()->GetVal());
    double delta_freq = (argval->GetFreqUp() - argval->GetFreqLo()) / argval->GetNfreq();
    long nrow = ntime;
    long ncol = 2 * argval->GetNfreq();
    // FFT
    long ntime_fft = hd1d_lc->GetNbinX();
    double delta_freq_fft = 1.0 / (hd1d_lc->GetXvalUp() - hd1d_lc->GetXvalLo());
    double freq_lo_fft = 0.0;
    double freq_up_fft = 1.0 / (2 * delta_time);
    long nfreq_fft = ntime_fft / 2;
    long nrow_fft = ntime_fft;
    long ncol_fft = ntime_fft;

    printf("delta_time = %e\n", delta_time);
    printf("ntime = %ld\n", ntime);
    printf("delta_freq = %e\n", delta_freq);
    printf("nrow = %ld\n", nrow);
    printf("ncol = %ld\n", ncol);
    printf("ntime_fft = %ld\n", ntime_fft);
    printf("delta_freq_fft = %e\n", delta_freq_fft);
    printf("freq_lo_fft = %e\n", freq_lo_fft);
    printf("freq_up_fft = %e\n", freq_up_fft);
    printf("nfreq_fft = %ld\n", nfreq_fft);
    printf("nrow_fft = %ld\n", nrow_fft);
    printf("ncol_fft = %ld\n", ncol_fft);
    

    // light curve of gti = 1
    GraphDataNerr2d* gd2d_lc = new GraphDataNerr2d;
    long itime = 0; 
    gd2d_lc->Init(ntime);
    for(long ibin = 0; ibin < hd1d_lc->GetNbinX(); ibin ++){
        if(1 == hd1d_gti->GetOvalElm(ibin)){
            gd2d_lc->SetPoint(itime,
                              hd1d_lc->GetHi1d()->GetBinCenter(ibin),
                              hd1d_lc->GetOvalElm(ibin));
            itime ++;
        }
    }
    printf("nbin ntime = %ld %ld\n",  hd1d_lc->GetNbinX(), ntime);
    
    GraphDataNerr2d* gd2d_norm_lc = new GraphDataNerr2d;
    Gd2dNormAndStd(gd2d_lc, gd2d_norm_lc);
    MirQdpTool::MkQdp(gd2d_lc, "lc.qdp", "x,y");
    MirQdpTool::MkQdp(gd2d_norm_lc, "lc_norm.qdp", "x,y");

    // -------------------------------

    double* power_arr = NULL;
    double* fft_real_arr = NULL;
    double* fft_image_arr = NULL;
    GenFFT(hd1d_lc->GetNbinX(), hd1d_lc->GetOvalArr()->GetVal(),
           &fft_real_arr, &fft_image_arr, &power_arr);

    FILE* fp = fopen("fft_power.qdp", "w");
    for(long ipow = 0; ipow < hd1d_lc->GetNbinX() / 2; ipow ++){
        fprintf(fp, "%ld  %e\n", ipow, power_arr[ipow]);
    }
    fclose(fp);

    fp = fopen("fft_real_image.qdp", "w");
    for(long ipow = 0; ipow < hd1d_lc->GetNbinX() / 2; ipow ++){
        fprintf(fp, "%ld  %e\n", ipow, fft_real_arr[ipow]);
    }
    for(long ipow = 0; ipow < hd1d_lc->GetNbinX() / 2; ipow ++){
        fprintf(fp, "%ld  %e\n", ipow + hd1d_lc->GetNbinX() / 2, fft_image_arr[ipow]);
    }    
    fclose(fp);

    printf("out fft\n");
    
    double* A_mat_arr = NULL;
    //GenMatLasso(gd2d_norm_lc->GetXvalArr()->GetVal(), ntime,
    //            argval->GetFreqLo(), delta_freq, argval->GetNfreq(),
    //            &A_mat_arr);
    GenMatLasso(gd2d_lc->GetXvalArr()->GetVal(), ntime,
                argval->GetFreqLo(), delta_freq, (long) argval->GetNfreq(),
                &A_mat_arr);
    
    double* h_arr = new double [nrow];
    for(long irow = 0; irow < nrow; irow++){
        // h_arr[irow] = gd2d_norm_lc->GetOvalElm(irow);
        h_arr[irow] = gd2d_lc->GetOvalElm(irow);
    }
    double* x_arr = new double [ncol];
    double* y_arr = new double [ncol];
    for(long icol = 0; icol < ncol; icol++){
        x_arr[icol] = 0.0;
        y_arr[icol] = 0.0;
    }
    
    double tolerance = 1.0e-10;
    double eta = 1.2;
    double lconst = 1.0e-3;
    double kiter_max = 10000;
    double cost = 0.0;
    double cost_pre = cost;
    for(long kiter = 0; kiter < kiter_max; kiter ++){
        //for(int icol = 0; icol < ncol; icol++){
        //    x_arr[icol] = 0.0;            
        //}
        double lconst_new = GetLconst(A_mat_arr, nrow, ncol,
                                      h_arr, y_arr,
                                      lconst, eta, argval->GetLambda());
        GetProxMap(A_mat_arr, nrow, ncol,
                   h_arr,
                   y_arr,
                   lconst_new,
                   argval->GetLambda(),
                   x_arr);
        cost = FGFunc(A_mat_arr, nrow, ncol,
                      h_arr, x_arr, argval->GetLambda());

        double ratio_cost_improve = (cost_pre - cost) / cost;
        printf("kiter = %ld, cost = %e, lconst_new = %e, ratio_cost_improve = %e\n",
               kiter, cost, lconst_new, ratio_cost_improve);
        
        if(kiter > 1 && fabs(ratio_cost_improve) < tolerance){
            printf("kiter = %ld, cost = %e\n", kiter, cost);
            break;
        }
        dcopy_(ncol, x_arr, 1, y_arr, 1);
        cost_pre = cost;
        lconst = lconst_new;
    }

    HistDataNerr1d* hd1d_lc_in_freq = new HistDataNerr1d;
    hd1d_lc_in_freq->Init(2 * argval->GetNfreq(), argval->GetFreqLo(), argval->GetFreqUp());
    hd1d_lc_in_freq->SetOvalArr(2 * argval->GetNfreq(), x_arr);
    // MirQdpTool::MkQdp(hd1d_lc_in_freq, "lc_in_freq.qdp", "x,y");

    fp = fopen("lasso_power.qdp", "w");
    for(long icol = 0; icol < argval->GetNfreq(); icol ++){
        fprintf(fp, "%ld  %e\n",
                icol, pow(hd1d_lc_in_freq->GetOvalElm(icol), 2)
                + pow(hd1d_lc_in_freq->GetOvalElm(icol + argval->GetNfreq()), 2) );
    }
    fclose(fp);

    fp = fopen("lasso_real_image.qdp", "w");
    for(long icol = 0; icol < argval->GetNfreq(); icol ++){
        fprintf(fp, "%ld  %e\n", icol, hd1d_lc_in_freq->GetOvalElm(icol));
    }
    for(long icol = 0; icol < argval->GetNfreq(); icol ++){
        fprintf(fp, "%ld  %e\n", icol + argval->GetNfreq(),
                hd1d_lc_in_freq->GetOvalElm(icol + argval->GetNfreq()) );
    }
    fclose(fp);

    // check light curve
    // rec_arr = A.mat %*% x_arr
    double* rec_arr = new double[nrow];
    for(long irow = 0; irow < nrow; irow ++){
        rec_arr[irow] = 0.0;
    }
    char* trans = new char [1];
    strcpy(trans, "N");
    dgemv_(trans, nrow, ncol, 1.0, const_cast<double*>(A_mat_arr), nrow,
           const_cast<double*>(x_arr), 1,
           0.0, rec_arr, 1);
    delete [] trans;

    FILE* fp_compare = fopen("compare_orglc_rec.qdp", "w");
    fprintf(fp_compare, "skip sing\n");
    for(long idata = 0; idata < gd2d_lc->GetNdata(); idata ++){
        fprintf(fp_compare, "%e  %e  %e \n",
                gd2d_lc->GetXvalElm(idata), gd2d_lc->GetOvalElm(idata), rec_arr[idata]);
    }        
    fclose(fp_compare);
    
    delete argval;
    
    return status;
}

