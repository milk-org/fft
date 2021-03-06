#ifndef _FFT_H
#define _FFT_H


void __attribute__((constructor)) libinit_fft();



errno_t import_wisdom();

int fft_setoffsets(long o1, long o2);

errno_t init_fftw_plans(int mode);

errno_t init_fftw_plans0();

errno_t export_wisdom();

int permut(const char *ID_name);

//void permutfliphv(const char *ID_name);

imageID do1dfft(const char *in_name, const char *out_name);

imageID do1drfft(const char *in_name, const char *out_name);

imageID do1dffti(const char *in_name, const char *out_name);

imageID do2dfft(const char *in_name, const char *out_name);

imageID do2dffti(const char *in_name, const char *out_name);

int pupfft(const char *ID_name_ampl, const char *ID_name_pha,
           const char *ID_name_ampl_out, const char *ID_name_pha_out, const char *options);

imageID do2drfft(const char *in_name, const char *out_name);

imageID do2drffti(const char *in_name, const char *out_name);

imageID fft_correlation(const char *ID_name1, const char *ID_name2,
                     const char *ID_nameout);

int fftzoom(const char *ID_name, const char *ID_out, long factor);

int fftczoom(const char *ID_name, const char *ID_out, long factor);

int test_fftspeed(int nmax);


imageID fft_DFT(
    const char *IDin_name,
    const char *IDinmask_name,
    const char *IDout_name,
    const char *IDoutmask_name,
    double      Zfactor,
    int         dir,
    long        kin
);


imageID fft_DFTinsertFPM(
    const char *pupin_name,
    const char *fpmz_name,
    double      zfactor,
    const char *pupout_name
);

long fft_DFTinsertFPM_re(const char *pupin_name, const char *fpmz_name,
                         double zfactor, const char *pupout_name);

int fft_image_translate(const char *ID_name, const char *ID_out, double xtransl,
                        double ytransl);

#endif
