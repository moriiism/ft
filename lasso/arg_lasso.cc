#include "arg_lasso.h"

// public

void ArgValLasso::Init(int argc, char* argv[])
{
    progname_ = "lasso";

    option long_options[] = {
        {"debug",      required_argument, NULL, 'd'},
        {"help",       required_argument, NULL, 'h'},
        {"verbose",    required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default
    SetOption(argc, argv, long_options);

    if(0 < g_flag_verbose){
        printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    }
    int narg = 10;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    infile_       = argv[iarg]; iarg++;
    gtifile_      = argv[iarg]; iarg++;
    time_lo_      = atof(argv[iarg]); iarg++;
    time_up_      = atof(argv[iarg]); iarg++;
    freq_lo_      = atof(argv[iarg]); iarg++;
    freq_up_      = atof(argv[iarg]); iarg++;
    nfreq_        = atoi(argv[iarg]); iarg++;
    lambda_       = atof(argv[iarg]); iarg++;
    outdir_       = argv[iarg]; iarg++;
    outfile_head_ = argv[iarg]; iarg++;    
}

void ArgValLasso::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);
    
    fprintf(fp, "%s: progname_      : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: infile_        : %s\n", __func__, infile_.c_str());
    fprintf(fp, "%s: gtifile_       : %s\n", __func__, gtifile_.c_str());    
    fprintf(fp, "%s: time_lo_       : %e\n", __func__, time_lo_);
    fprintf(fp, "%s: time_up_       : %e\n", __func__, time_up_);
    fprintf(fp, "%s: freq_lo_       : %e\n", __func__, freq_lo_);
    fprintf(fp, "%s: freq_up_       : %e\n", __func__, freq_up_);
    fprintf(fp, "%s: nfreq_         : %d\n", __func__, nfreq_);
    fprintf(fp, "%s: lambda_        : %e\n", __func__, lambda_);
    fprintf(fp, "%s: outdir_        : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_  : %s\n", __func__, outfile_head_.c_str());    
}

void ArgValLasso::Null()
{
    progname_     = "";
    infile_       = "";
    gtifile_      = "";    
    time_lo_      = 0.0;
    time_up_      = 0.0;    
    freq_lo_      = 0.0;
    freq_up_      = 0.0;
    nfreq_        = 0;
    lambda_       = 0.0;
    outdir_       = "";
    outfile_head_ = "";    
}

void ArgValLasso::SetOption(int argc, char* argv[], option* long_options)
{
    if(0 < g_flag_verbose){
        MPrintInfo("start...");
    }
    // option default
    g_flag_debug   = 0;
    g_flag_help    = 0;
    g_flag_verbose = 0;
    while (1) {
        int option_index = 0;
        int retopt = getopt_long(argc, argv, "dhv",
                                 long_options, &option_index);
        if(-1 == retopt)
            break;
        switch (retopt) {
        case 0:
            // long option
            break;
        case 'd':
            g_flag_debug = atoi(optarg);
            printf("%s: g_flag_debug = %d\n", __func__, g_flag_debug);
            break;
        case 'h':
            g_flag_help = atoi(optarg);
            printf("%s: g_flag_help = %d\n", __func__, g_flag_help);
            if(0 != g_flag_help){
                Usage(stdout);
            }
            break;
        case 'v':
            g_flag_verbose = atoi(optarg);
            printf("%s: g_flag_verbose = %d\n", __func__, g_flag_verbose);
            break;
        case '?':
            printf("%s: retopt (= %c) is invalid flag.\n",
                   __func__, retopt);
            Usage(stdout);
            break;

        default:
            printf("%s: error: getopt returned character code 0%o ??\n",
                   __func__, retopt);
            abort();
        }
    }
    if(0 < g_flag_verbose){
        MPrintInfo("done.");
    }
}

void ArgValLasso::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "infile  gtifile  time_lo  time_up  freq_lo  freq_up  nfreq  lambda  "
            "outdir  outfile_head \n",
            progname_.c_str());
    exit(1);
}
