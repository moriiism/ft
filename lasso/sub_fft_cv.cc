#include "sub_fft_cv.h"
#include "sub_fft.h"
#include "sub.h"

void GetLassoFft(const double* const h_arr, long nrow,
                 const long* const index_arr,
                 long ncol,
                 double delta_freq,
                 long nfreq_fft,
                 double freq_lo_fft,
                 double freq_up_fft,
                 double lambda,
                 double* const out_arr)
{
    //double* h_arr = new double [nrow];
    //for(long irow = 0; irow < nrow; irow++){
    //    // h_arr[irow] = gd2d_norm_lc->GetOvalElm(irow);
    //    h_arr[irow] = gd2d_lc->GetOvalElm(irow);
    // }

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
        double lconst_new = GetLconstFft(h_arr, nrow,
                                         index_arr,
                                         y_arr, ncol,
                                         delta_freq,
                                         lconst,
                                         eta,
                                         lambda);
        GetProxMapFft(h_arr, nrow,
                      index_arr,
                      y_arr, ncol,
                      delta_freq,
                      lconst_new,
                      lambda,
                      x_arr);
        cost = FGFuncFft(h_arr, nrow,
                         index_arr,
                         x_arr, ncol,
                         delta_freq,
                         lambda);
        double ratio_cost_improve = (cost_pre - cost) / cost;
        //printf("kiter = %ld, cost = %e, lconst_new = %e, ratio_cost_improve = %e\n",
        //       kiter, cost, lconst_new, ratio_cost_improve);
        
        if(kiter > 1 && fabs(ratio_cost_improve) < tolerance){
            printf("kiter = %ld, cost = %e\n", kiter, cost);
            break;
        }
        dcopy_(ncol, x_arr, 1, y_arr, 1);
        cost_pre = cost;
        lconst = lconst_new;
    }

    HistDataNerr1d* hd1d_lc_in_freq = new HistDataNerr1d;
    hd1d_lc_in_freq->Init(2 * nfreq_fft, freq_lo_fft, freq_up_fft);
    hd1d_lc_in_freq->SetOvalArr(2 * nfreq_fft, x_arr);
    // MirQdpTool::MkQdp(hd1d_lc_in_freq, "lc_in_freq.qdp", "x,y");


    char outqdp[kLineSize];
    sprintf(outqdp, "lasso_power_%1.1e.qdp", lambda);
    FILE* fp = fopen(outqdp, "w");
    for(long icol = 0; icol < nfreq_fft; icol ++){
        double freq = freq_lo_fft + delta_freq * icol;
        fprintf(fp, "%e  %e\n",
                freq, pow(hd1d_lc_in_freq->GetOvalElm(icol), 2)
                + pow(hd1d_lc_in_freq->GetOvalElm(icol + nfreq_fft), 2) );
    }
    fclose(fp);

    fp = fopen("lasso_real_image.qdp", "w");
    for(long icol = 0; icol < nfreq_fft; icol ++){
        fprintf(fp, "%ld  %e\n", icol, hd1d_lc_in_freq->GetOvalElm(icol));
    }
    for(long icol = 0; icol < nfreq_fft; icol ++){
        fprintf(fp, "%ld  %e\n", icol + nfreq_fft,
                hd1d_lc_in_freq->GetOvalElm(icol + nfreq_fft) );
    }
    fclose(fp);

    for(long icol = 0; icol < ncol; icol++){
        out_arr[icol] = x_arr[icol];
    }

    delete hd1d_lc_in_freq;
    delete [] x_arr;
    delete [] y_arr;
}

