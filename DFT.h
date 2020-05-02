/** @file DFT.h
 */

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
