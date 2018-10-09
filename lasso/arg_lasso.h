#ifndef MORIIISM_FT_LASSO_ARG_LASSO_H_
#define MORIIISM_FT_LASSO_ARG_LASSO_H_

#include "mi_base.h"

class ArgValLasso : public MiArgBase{
public:
    ArgValLasso() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        gtifile_(""),
        time_lo_(0.0),
        time_up_(0.0),
        freq_lo_(0.0),
        freq_up_(0.0),
        nfreq_(0),
        lambda_(0.0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValLasso(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile()   const {return infile_;};
    string GetGtifile()   const {return gtifile_;};
    double GetTimeLo()   const {return time_lo_;};
    double GetTimeUp()   const {return time_up_;};
    double GetFreqLo()   const {return freq_lo_;};
    double GetFreqUp()   const {return freq_up_;};
    int    GetNfreq()    const {return nfreq_;};
    double GetLambda()   const {return lambda_;};
    string GetOutdir()   const {return outdir_;};
    string GetOutfileHead()   const {return outfile_head_;};

private:
    string progname_;
    string infile_;
    string gtifile_;
    double time_lo_;
    double time_up_;
    double freq_lo_;
    double freq_up_;
    int    nfreq_;
    double lambda_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_FT_LASSO_ARG_LASSO_H_
