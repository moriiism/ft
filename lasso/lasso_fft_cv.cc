#include "fftw3.h"
#include "mi_iolib.h"
#include "mi_time.h"
#include "mir_hist1d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_lasso_fft_cv.h"
#include "sub_fft_cv.h"
#include "sub_fft.h"
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

    double time_st = MiTime::GetTimeSec();
    
    ArgValLassoFftCv* argval = new ArgValLassoFftCv;
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
    // N
    long ntime = MirMath::GetSum(hd1d_gti->GetOvalArr()->GetNdata(),
                                 hd1d_gti->GetOvalArr()->GetVal());
    // FFT
    long ntime_fft = hd1d_lc->GetNbinX();
    double delta_freq_fft = 1.0 / (hd1d_lc->GetXvalUp() - hd1d_lc->GetXvalLo());
    double freq_lo_fft = 0.0;
    double freq_up_fft = 1.0 / (2 * delta_time);
    long nfreq_fft = ntime_fft / 2;
    long nrow_fft = ntime_fft;
    long ncol_fft = ntime_fft;

    long nrow = ntime;
    long ncol = ncol_fft;
    double delta_freq = delta_freq_fft;

    printf("delta_time = %e\n", delta_time);
    printf("ntime = %ld\n", ntime);
    printf("ntime_fft = %ld\n", ntime_fft);
    printf("delta_freq_fft = %e\n", delta_freq_fft);
    printf("freq_lo_fft = %e\n", freq_lo_fft);
    printf("freq_up_fft = %e\n", freq_up_fft);
    printf("nfreq_fft = %ld\n", nfreq_fft);
    printf("nrow_fft = %ld\n", nrow_fft);
    printf("ncol_fft = %ld\n", ncol_fft);

    // light curve of gti = 1
    GraphDataNerr2d* gd2d_lc = new GraphDataNerr2d;
    // index of gti = 1
    long* index_arr = new long [ntime];
    long itime = 0; 
    gd2d_lc->Init(ntime);
    for(long ibin = 0; ibin < hd1d_lc->GetNbinX(); ibin ++){
        if(1 == hd1d_gti->GetOvalElm(ibin)){
            gd2d_lc->SetPoint(itime,
                              hd1d_lc->GetHi1d()->GetBinCenter(ibin),
                              hd1d_lc->GetOvalElm(ibin));
            index_arr[itime] = ibin;
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
        double freq = freq_lo_fft + delta_freq * ipow;
        fprintf(fp, "%e  %e\n", freq, power_arr[ipow]);
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

    double* lambda_arr = new double [argval->GetNlambda()];
    double* cv_variance_arr = new double [argval->GetNlambda()];
    for(int ilambda = 0; ilambda < argval->GetNlambda(); ilambda ++){
        lambda_arr[ilambda] = 0.0;
        cv_variance_arr[ilambda] = 0.0;
    }
    
    double loglambda_lo = log10(argval->GetLambdaLo());
    double loglambda_up = log10(argval->GetLambdaUp());
    int nlambda = argval->GetNlambda();
    double delta_loglambda = (loglambda_up - loglambda_lo) / nlambda;
    for(int ilambda = 0; ilambda < argval->GetNlambda(); ilambda ++){
        double loglambda = loglambda_lo + delta_loglambda * ilambda;
        double lambda = pow(10, loglambda);
        lambda_arr[ilambda] = lambda;

//        GetLassoFft(gd2d_lc->GetOvalArr()->GetVal(), nrow,
//                    index_arr,
//                    ncol,
//                     delta_freq,
//                    nfreq_fft,
//                    freq_lo_fft,
//                    freq_up_fft,
//                    lambda);

        int nfold = argval->GetNfold();
        
        int rand_seed = 1;
        TRandom3* trand = new TRandom3(rand_seed);
        double* rand_arr = new double[nrow];
        for(long irow = 0; irow < nrow; irow++){
            double rand = trand->Rndm();
            rand_arr[irow] = rand;
        }
        delete trand;
        long* index_rand_arr = new long[nrow];
        MiSort::Sort(nrow, rand_arr, index_rand_arr, 1);
        delete [] rand_arr;

        int* nbin_fold_arr = new int[nfold];
        for(int ifold = 0; ifold < nfold; ifold ++){
            nbin_fold_arr[ifold] = 0;
        }
        for(long irow = 0; irow < nrow; irow ++){
            int ifold_modulo = irow % nfold;
            nbin_fold_arr[ifold_modulo] ++;
        }

        //for(int ifold = 0; ifold < nfold; ifold ++){
        //    printf("nbin_fold_arr[%d] = %d\n", ifold, nbin_fold_arr[ifold]);
        //}

        double* variance_arr = new double[nfold];
        for(int ifold_cv = 0; ifold_cv < nfold; ifold_cv ++){
            long* index_vl_arr = new long [nbin_fold_arr[ifold_cv]];
            long* index_tr_arr = new long [nrow - nbin_fold_arr[ifold_cv]];
            long index_vl_st = 0;
            for(int ifold = 0; ifold < ifold_cv; ifold ++){
                index_vl_st += nbin_fold_arr[ifold];
            }
            long index_vl_ed = index_vl_st + nbin_fold_arr[ifold_cv] - 1;

            long ivl = 0;
            long itr = 0;
            for(long irow = 0; irow < nrow; irow ++){
                if(index_vl_st <= irow && irow <= index_vl_ed){
                    index_vl_arr[ivl] = index_rand_arr[irow];
                    ivl ++;
                } else {
                    index_tr_arr[itr] = index_rand_arr[irow];
                    itr ++;
                }
            }
            //printf("ivl = %ld\n", ivl);
            //printf("itr = %ld\n", itr);

            std::sort(index_vl_arr, index_vl_arr + nbin_fold_arr[ifold_cv]);
            std::sort(index_tr_arr, index_tr_arr + nrow - nbin_fold_arr[ifold_cv]);
            //for(long ivl = 0; ivl < nbin_fold_arr[ifold_cv]; ivl ++){
                // printf("index_vl_arr[%ld] = %ld\n", ivl, index_vl_arr[ivl]);
            //}
            long* index_fold_vl_arr = new long [nbin_fold_arr[ifold_cv]];
            long* index_fold_tr_arr = new long [nrow - nbin_fold_arr[ifold_cv]];
            for(long itr = 0; itr < nrow - nbin_fold_arr[ifold_cv]; itr ++){
                index_fold_tr_arr[itr] = index_arr[ index_tr_arr[itr] ];
            }
            for(long ivl = 0; ivl < nbin_fold_arr[ifold_cv]; ivl ++){
                index_fold_vl_arr[ivl] = index_arr[ index_vl_arr[ivl] ];
            }
            delete [] index_vl_arr;
            delete [] index_tr_arr;

            // light curve of gti = 1 and index_fold_tr_arr
            GraphDataNerr2d* gd2d_tr_lc = new GraphDataNerr2d;
            gd2d_tr_lc->Init(nrow - nbin_fold_arr[ifold_cv]);
            for(long ibin = 0; ibin < gd2d_tr_lc->GetNdata(); ibin ++){
                gd2d_tr_lc->SetPoint(ibin,
                                     hd1d_lc->GetHi1d()->GetBinCenter( index_fold_tr_arr[ibin] ),
                                     hd1d_lc->GetOvalElm( index_fold_tr_arr[ibin] ));
            }

            double* x_arr = new double [ncol];
            GetLassoFft(gd2d_tr_lc->GetOvalArr()->GetVal(),
                        nrow - nbin_fold_arr[ifold_cv],
                        index_fold_tr_arr,
                        ncol,
                        delta_freq,
                        nfreq_fft,
                        freq_lo_fft,
                        freq_up_fft,
                        lambda,
                        x_arr);
            
            // --------------------

            double* reconst_tr_arr = new double [nrow - nbin_fold_arr[ifold_cv]];
            AxFft(x_arr, ncol,
                  index_fold_tr_arr,
                  nrow - nbin_fold_arr[ifold_cv],
                  delta_freq,
                  reconst_tr_arr);

            double* reconst_arr = new double [nrow];
            AxFft(x_arr, ncol,
                  index_arr,
                  nrow,
                  delta_freq,
                  reconst_arr);

            double* reconst_vl_arr = new double [nbin_fold_arr[ifold_cv]];
            AxFft(x_arr, ncol,
                  index_fold_vl_arr,
                  nbin_fold_arr[ifold_cv],
                  delta_freq,
                  reconst_vl_arr);
            
            FILE* fp_out = fopen("compare.qdp", "w");
            fprintf(fp_out, "skip sing\n");
            fprintf(fp_out, "\n");
            for(long ibin = 0; ibin < nrow - nbin_fold_arr[ifold_cv]; ibin ++){
                fprintf(fp_out, "%e  %e\n",
                        hd1d_lc->GetHi1d()->GetBinCenter( index_fold_tr_arr[ibin] ),
                        reconst_tr_arr[ibin]);
            }
            fprintf(fp_out, "\n");
            fprintf(fp_out, "no\n");
            fprintf(fp_out, "\n");
            for(long ibin = 0; ibin < nbin_fold_arr[ifold_cv]; ibin ++){
                fprintf(fp_out, "%e  %e\n",
                        hd1d_lc->GetHi1d()->GetBinCenter( index_fold_vl_arr[ibin] ),
                        reconst_vl_arr[ibin]);
            }
            fprintf(fp_out, "\n");
            fprintf(fp_out, "no\n");
            fprintf(fp_out, "\n");
            for(long ibin = 0; ibin < nrow; ibin ++){
                fprintf(fp_out, "%e  %e\n",
                        hd1d_lc->GetHi1d()->GetBinCenter( index_arr[ibin] ),
                        reconst_arr[ibin]);
            }
            fprintf(fp_out, "\n");
            fprintf(fp_out, "no\n");
            fprintf(fp_out, "\n");
            for(long ibin = 0; ibin < ntime; ibin ++){
                fprintf(fp_out, "%e  %e\n",
                        gd2d_lc->GetXvalElm(ibin),
                        gd2d_lc->GetOvalElm(ibin));
            }
            fclose(fp_out);

            // validation error
            double variance = 0.0;
            for(long ibin = 0; ibin < nbin_fold_arr[ifold_cv]; ibin ++){
                variance += pow( reconst_vl_arr[ibin] - hd1d_lc->GetOvalElm( index_fold_vl_arr[ibin] ), 2);
            }
            variance_arr[ifold_cv] = variance;

            delete [] x_arr;
            delete [] reconst_tr_arr;
            delete [] reconst_vl_arr;            
            delete [] reconst_arr;

            delete gd2d_tr_lc;
            delete [] index_fold_vl_arr;
            delete [] index_fold_tr_arr;
        }

        double ave_variance = 0.0;
        for(int ifold = 0; ifold < nfold; ifold ++){
            // printf("%d %e\n", ifold, variance_arr[ifold]);
            ave_variance += variance_arr[ifold];
        }
        ave_variance /= nfold;
        printf("lambda, ave_variance = %e %e\n", lambda, ave_variance);

        cv_variance_arr[ilambda] = ave_variance;

        delete [] variance_arr;
        delete [] nbin_fold_arr;        
        delete [] index_rand_arr;
    }

    FILE* fp_cv_out = fopen("cv.qdp", "w");
    for(int ilambda = 0; ilambda < argval->GetNlambda(); ilambda ++){
        fprintf(fp_cv_out, "%e %e\n", lambda_arr[ilambda], cv_variance_arr[ilambda]);
    }
    fclose(fp_cv_out);

    int cv_best_index = 0;
    double cv_variance_min = 1e10;
    for(int ilambda = 0; ilambda < argval->GetNlambda(); ilambda ++){
        if(cv_variance_arr[ilambda] < cv_variance_min){
            cv_variance_min = cv_variance_arr[ilambda];
            cv_best_index = ilambda;
        }
    }

    MiIolib::Printf2(fp_log, "best lambda: %e\n", lambda_arr[cv_best_index]);
    double lambda_best = lambda_arr[cv_best_index];

    double* x_best_arr = new double [ncol];
    GetLassoFft(gd2d_lc->GetOvalArr()->GetVal(), nrow,
                index_arr,
                ncol,
                delta_freq,
                nfreq_fft,
                freq_lo_fft,
                freq_up_fft,
                lambda_best,
                x_best_arr);
    

    double time_ed = MiTime::GetTimeSec();
    double calc_time = time_ed - time_st;
    MiIolib::Printf2(fp_log, "calc_time = %e sec, %e min, %e hour\n",
                     calc_time, calc_time/60., calc_time/3600.);
    
    delete argval;
    
    return status;
}

