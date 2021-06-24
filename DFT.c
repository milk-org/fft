/** @file DFT.c
 */


#include <math.h>

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"


# ifdef _OPENMP
# include <omp.h>
# endif

/* ----------------- CUSTOM DFT ------------- */

//
// Zfactor is zoom factor
// dir = -1 for FT, 1 for inverse FT
// kin in selects slice in IDin_name if this is a cube
//
imageID fft_DFT(
    const char *IDin_name,
    const char *IDinmask_name,
    const char *IDout_name,
    const char *IDoutmask_name,
    double      Zfactor,
    int         dir,
    long        kin
)
{
    imageID IDin;
    imageID IDout;
    imageID IDinmask;
    imageID IDoutmask;

    uint32_t NBptsin;
    uint32_t NBptsout;

    uint32_t xsize, ysize;
    uint32_t ii, jj, k, kout;
    double val;
    double re, im;
    float pha;

    uint_fast16_t *iiinarray;
    uint_fast16_t *jjinarray;
    double *xinarray;
    double *yinarray;
    double *valinamp;
    double *valinpha;
    float *cosvalinpha;
    float *sinvalinpha;

    uint_fast16_t *iioutarray;
    uint_fast16_t *jjoutarray;
    double *xoutarray;
    double *youtarray;

    float cospha, sinpha;

    long IDcosXX, IDcosYY, IDsinXX, IDsinYY;

    // list of active coordinates
    uint_fast16_t *iiinarrayActive;
    uint_fast16_t *jjinarrayActive;
    uint_fast16_t *iioutarrayActive;
    uint_fast16_t *jjoutarrayActive;
    uint_fast16_t pixiiin, pixiiout, pixjjin, pixjjout;

    uint_fast8_t pixact;
    uint_fast16_t NBpixact_iiin, NBpixact_jjin;
    uint_fast16_t NBpixact_iiout, NBpixact_jjout;

    float *XinarrayActive;
    float *YinarrayActive;
    float *XoutarrayActive;
    float *YoutarrayActive;
    uint_fast16_t iiin, jjin, iiout, jjout;

    float cosXX, sinXX, cosYY, sinYY, cosXY, sinXY;





    IDin = image_ID(IDin_name);

    IDinmask = image_ID(IDinmask_name);
    xsize = data.image[IDinmask].md[0].size[0];
    ysize = data.image[IDinmask].md[0].size[1];
    iiinarrayActive = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * xsize);
    jjinarrayActive = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * ysize);
    iioutarrayActive = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * xsize);
    jjoutarrayActive = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * ysize);





    NBptsin = 0;
    NBpixact_iiin = 0;
    for(ii = 0; ii < xsize; ii++)
    {
        pixact = 0;
        for(jj = 0; jj < ysize; jj++)
        {
            val = data.image[IDinmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                pixact = 1;
                NBptsin ++;
            }
        }
        if(pixact == 1)
        {
            iiinarrayActive[NBpixact_iiin] = ii;
            NBpixact_iiin++;
        }
    }

    NBpixact_jjin = 0;
    for(jj = 0; jj < ysize; jj++)
    {
        pixact = 0;
        for(ii = 0; ii < xsize; ii++)
        {
            val = data.image[IDinmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                pixact = 1;
            }
        }
        if(pixact == 1)
        {
            jjinarrayActive[NBpixact_jjin] = jj;
            NBpixact_jjin++;
        }
    }

    XinarrayActive = (float *) malloc(sizeof(float) * NBpixact_iiin);
    YinarrayActive = (float *) malloc(sizeof(float) * NBpixact_jjin);

    for(pixiiin = 0; pixiiin < NBpixact_iiin; pixiiin++)
    {
        iiin = iiinarrayActive[pixiiin];
        XinarrayActive[pixiiin] = (1.0 * iiin / xsize - 0.5);
    }
    for(pixjjin = 0; pixjjin < NBpixact_jjin; pixjjin++)
    {
        jjin = jjinarrayActive[pixjjin];
        YinarrayActive[pixjjin] = (1.0 * jjin / ysize - 0.5);
    }

    printf("DFT (factor %f, slice %ld):  %u input points (%ld %ld)-> ", Zfactor,
           kin, NBptsin, NBpixact_iiin, NBpixact_jjin);

    iiinarray = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * NBptsin);
    jjinarray = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * NBptsin);
    xinarray = (double *) malloc(sizeof(double) * NBptsin);
    yinarray = (double *) malloc(sizeof(double) * NBptsin);
    valinamp = (double *) malloc(sizeof(double) * NBptsin);
    valinpha = (double *) malloc(sizeof(double) * NBptsin);
    cosvalinpha = (float *) malloc(sizeof(float) * NBptsin);
    sinvalinpha = (float *) malloc(sizeof(float) * NBptsin);
    k = 0;



    for(ii = 0; ii < xsize; ii++)
        for(jj = 0; jj < ysize; jj++)
        {
            val = data.image[IDinmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                iiinarray[k] = ii;
                jjinarray[k] = jj;
                xinarray[k] = 1.0 * ii / xsize - 0.5;
                yinarray[k] = 1.0 * jj / xsize - 0.5;
                re = data.image[IDin].array.CF[kin * xsize * ysize + jj * xsize + ii].re;
                im = data.image[IDin].array.CF[kin * xsize * ysize + jj * xsize + ii].im;
                valinamp[k] = sqrt(re * re + im * im);
                valinpha[k] = atan2(im, re);
                cosvalinpha[k] = cosf(valinpha[k]);
                sinvalinpha[k] = sinf(valinpha[k]);
                k++;
            }
        }




    IDoutmask = image_ID(IDoutmask_name);

    NBptsout = 0;
    NBpixact_iiout = 0;
    for(ii = 0; ii < xsize; ii++)
    {
        pixact = 0;
        for(jj = 0; jj < ysize; jj++)
        {
            val = data.image[IDoutmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                pixact = 1;
                NBptsout ++;
            }
        }
        if(pixact == 1)
        {
            iioutarrayActive[NBpixact_iiout] = ii;
            NBpixact_iiout++;
        }
    }

    NBpixact_jjout = 0;
    for(jj = 0; jj < ysize; jj++)
    {
        pixact = 0;
        for(ii = 0; ii < xsize; ii++)
        {
            val = data.image[IDoutmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                pixact = 1;
            }
        }
        if(pixact == 1)
        {
            jjoutarrayActive[NBpixact_jjout] = jj;
            NBpixact_jjout++;
        }
    }
    XoutarrayActive = (float *) malloc(sizeof(float) * NBpixact_iiout);
    YoutarrayActive = (float *) malloc(sizeof(float) * NBpixact_jjout);

    for(pixiiout = 0; pixiiout < NBpixact_iiout; pixiiout++)
    {
        iiout = iioutarrayActive[pixiiout];
        XoutarrayActive[pixiiout] = (1.0 / Zfactor) * (1.0 * iiout / xsize - 0.5) *
                                    xsize;
    }

    for(pixjjout = 0; pixjjout < NBpixact_jjout; pixjjout++)
    {
        jjout = jjoutarrayActive[pixjjout];
        YoutarrayActive[pixjjout] = (1.0 / Zfactor) * (1.0 * jjout / ysize - 0.5) *
                                    ysize;
    }

    printf("%u output points (%ld %ld) \n", NBptsout, NBpixact_iiout,
           NBpixact_jjout);



    iioutarray = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * NBptsout);
    jjoutarray = (uint_fast16_t *) malloc(sizeof(uint_fast16_t) * NBptsout);
    xoutarray = (double *) malloc(sizeof(double) * NBptsout);
    youtarray = (double *) malloc(sizeof(double) * NBptsout);
    kout = 0;
    for(ii = 0; ii < xsize; ii++)
        for(jj = 0; jj < ysize; jj++)
        {
            val = data.image[IDoutmask].array.F[jj * xsize + ii];
            if(val > 0.5)
            {
                iioutarray[kout] = ii;
                jjoutarray[kout] = jj;
                xoutarray[kout] = (1.0 / Zfactor) * (1.0 * ii / xsize - 0.5) * xsize;
                youtarray[kout] = (1.0 / Zfactor) * (1.0 * jj / ysize - 0.5) * ysize;
                kout++;
            }
        }

    IDout = create_2DCimage_ID(IDout_name, xsize, ysize);




    IDcosXX = create_2Dimage_ID("_cosXX", xsize, xsize);
    IDsinXX = create_2Dimage_ID("_sinXX", xsize, xsize);
    IDcosYY = create_2Dimage_ID("_cosYY", ysize, ysize);
    IDsinYY = create_2Dimage_ID("_sinYY", ysize, ysize);



    printf(" <");
    fflush(stdout);

#ifdef _OPENMP
    printf("Using openMP %d", omp_get_max_threads());
    #pragma omp parallel default(shared) private(pixiiout, pixiiin, iiout, iiin, pha, cospha, sinpha)
    {
        #pragma omp for
#endif

        for(pixiiout = 0; pixiiout < NBpixact_iiout; pixiiout++)
        {
            iiout = iioutarrayActive[pixiiout];
            for(pixiiin = 0; pixiiin < NBpixact_iiin; pixiiin++)
            {
                iiin = iiinarrayActive[pixiiin];
                pha = 2.0 * dir * M_PI * (XinarrayActive[pixiiin] * XoutarrayActive[pixiiout]);
                cospha = cosf(pha);
                sinpha = sinf(pha);

                data.image[IDcosXX].array.F[iiout * xsize + iiin] = cospha;
                data.image[IDsinXX].array.F[iiout * xsize + iiin] = sinpha;

            }
        }
# ifdef _OPENMP
    }
# endif

    printf("> ");
    fflush(stdout);




    printf(" <");
    fflush(stdout);

# ifdef _OPENMP
    #pragma omp parallel default(shared) private(pixjjout, pixjjin, jjout, jjin, pha, cospha, sinpha)
    {
        #pragma omp for
# endif
        for(pixjjout = 0; pixjjout < NBpixact_jjout; pixjjout++)
        {
            jjout = jjoutarrayActive[pixjjout];
            for(pixjjin = 0; pixjjin < NBpixact_jjin; pixjjin++)
            {
                jjin = jjinarrayActive[pixjjin];
                pha = 2.0 * dir * M_PI * (YinarrayActive[pixjjin] * YoutarrayActive[pixjjout]);
                cospha = cosf(pha);
                sinpha = sinf(pha);

                data.image[IDcosYY].array.F[jjout * ysize + jjin] = cospha;
                data.image[IDsinYY].array.F[jjout * ysize + jjin] = sinpha;

            }
        }
# ifdef _OPENMP
    }
# endif
    printf("> ");
    fflush(stdout);


    // DFT


    printf(" <");
    fflush(stdout);


#ifdef _OPENMP
    printf(" -omp- %d ", omp_get_max_threads());
    fflush(stdout);
    #pragma omp parallel default(shared) private(kout, k, pha, re, im, cospha, sinpha, iiin, jjin, iiout, jjout, cosXX, cosYY, sinXX, sinYY, cosXY, sinXY)
    {
        #pragma omp master
        {
            printf(" [%d thread(s)] ", omp_get_num_threads());
            fflush(stdout);
        }

        #pragma omp for
#endif

        for(kout = 0; kout < NBptsout; kout++)
        {
            iiout = iioutarray[kout];
            jjout = jjoutarray[kout];

            re = 0.0;
            im = 0.0;
            for(k = 0; k < NBptsin; k++)
            {
                iiin = iiinarray[k];
                jjin = jjinarray[k];

                cosXX = data.image[IDcosXX].array.F[iiout * xsize + iiin];
                cosYY = data.image[IDcosYY].array.F[jjout * ysize + jjin];

                sinXX = data.image[IDsinXX].array.F[iiout * xsize + iiin];
                sinYY = data.image[IDsinYY].array.F[jjout * ysize + jjin];

                cosXY = cosXX * cosYY - sinXX * sinYY;
                sinXY = sinXX * cosYY + cosXX * sinYY;

                cospha = cosvalinpha[k] * cosXY - sinvalinpha[k] * sinXY;
                sinpha = sinvalinpha[k] * cosXY + cosvalinpha[k] * sinXY;


                re += valinamp[k] * cospha;
                im += valinamp[k] * sinpha;
            }
            data.image[IDout].array.CF[jjoutarray[kout]*xsize + iioutarray[kout]].re = re /
                    Zfactor;
            data.image[IDout].array.CF[jjoutarray[kout]*xsize + iioutarray[kout]].im = im /
                    Zfactor;
        }

#ifdef _OPENMP
    }
#endif
    printf("> ");
    fflush(stdout);

    free(cosvalinpha);
    free(sinvalinpha);

    delete_image_ID("_cosXX", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("_sinXX", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("_cosYY", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("_sinYY", DELETE_IMAGE_ERRMODE_WARNING);

    free(XinarrayActive);
    free(YinarrayActive);
    free(XoutarrayActive);
    free(YoutarrayActive);

    free(iiinarrayActive);
    free(jjinarrayActive);
    free(iioutarrayActive);
    free(jjoutarrayActive);

    free(iiinarray);
    free(jjinarray);
    free(xinarray);
    free(yinarray);
    free(valinamp);
    free(valinpha);


    free(iioutarray);
    free(jjoutarray);
    free(xoutarray);
    free(youtarray);

    return IDout;
}





/**
 * @brief Use DFT to insert Focal Plane Mask
 *
 *  Pupil convolution by complex focal plane mask of limited support
 *  typically used with fpmz = zoomed copy of 1-fpm
 *
 * High resolution focal plane mask using DFT
 *
 * Forces computation over pixels >0.5 in _DFTmask00 if it exists
 *
 */
imageID fft_DFTinsertFPM(
    const char *pupin_name,
    const char *fpmz_name,
    double      zfactor,
    const char *pupout_name
)
{
    double eps = 1.0e-16;
    imageID ID;
    imageID IDpupin_mask;
    imageID IDfpmz;
    imageID IDfpmz_mask;
    long xsize, ysize, zsize;
    imageID IDin, IDout;
    long ii, jj, k;
    double re, im, rein, imin, amp, pha, ampin, phain, amp2;
    double x, y;
    double total = 0;
    imageID IDout2D;
    int FORCE_IMZERO = 0;
    //double imresidual = 0.0;
    double tx, ty, tcx, tcy;
    long size2;

    long ID_DFTmask00;


    if(variable_ID("_FORCE_IMZERO") != -1)
    {
        FORCE_IMZERO = 1;
        printf("---------------FORCING IMAGINARY PART TO ZERO-------------\n");
    }


    ID_DFTmask00 = image_ID("_DFTmask00");


    printf("zfactor = %f\n", zfactor);

    IDin = image_ID(pupin_name);
    xsize = data.image[IDin].md[0].size[0];
    ysize = data.image[IDin].md[0].size[1];
    if(data.image[IDin].md[0].naxis > 2)
    {
        zsize = data.image[IDin].md[0].size[2];
    }
    else
    {
        zsize = 1;
    }
    printf("zsize = %ld\n", zsize);
    size2 = xsize * ysize;

    IDout = create_3DCimage_ID(pupout_name, xsize, ysize, zsize);


    for(k = 0; k < zsize; k++) // increment slice (= wavelength)
    {
        //
        // Create default input mask for DFT
        // if amplitude above threshold value, turn pixel "on"
        //
        IDpupin_mask = create_2Dimage_ID("_DFTpupmask", xsize, ysize);
        for(ii = 0; ii < xsize * ysize; ii++)
        {
            re = data.image[IDin].array.CF[k * size2 + ii].re;
            im = data.image[IDin].array.CF[k * size2 + ii].im;
            amp2 = re * re + im * im;
            if(amp2 > eps)
            {
                data.image[IDpupin_mask].array.F[ii] = 1.0;
            }
            else
            {
                data.image[IDpupin_mask].array.F[ii] = 0.0;
            }
        }
        //
        // If _DFTmask00 exists, make corresponding pixels = 1
        //
        if(ID_DFTmask00 != -1)
            for(ii = 0; ii < xsize * ysize; ii++)
            {
                if(data.image[ID_DFTmask00].array.F[ii] > 0.5)
                {
                    data.image[IDpupin_mask].array.F[ii] = 1.0;
                }
            }

        //
        // Construct focal plane mask for DFT
        // If amplitude >eps, turn pixel ON, save result in _fpmzmask
        //
        IDfpmz = image_ID(fpmz_name);
        IDfpmz_mask = create_2Dimage_ID("_fpmzmask", xsize, ysize);
        for(ii = 0; ii < xsize * ysize; ii++)
        {
            re = data.image[IDfpmz].array.CF[k * size2 + ii].re;
            im = data.image[IDfpmz].array.CF[k * size2 + ii].im;
            amp2 = re * re + im * im;
            if(amp2 > eps)
            {
                data.image[IDfpmz_mask].array.F[ii] = 1.0;
            }
            else
            {
                data.image[IDfpmz_mask].array.F[ii] = 0.0;
            }
        }

        //	save_fits("_DFTpupmask", "!_DFTpupmask.fits");



        fft_DFT(pupin_name, "_DFTpupmask", "_foc0", "_fpmzmask", zfactor, -1, k);

        ID = image_ID("_foc0");
        total = 0.0;
        tx = 0.0;
        ty = 0.0;
        tcx = 0.0;
        tcy = 0.0;
        for(ii = 0; ii < xsize; ii++)
            for(jj = 0; jj < ysize; jj++)
            {
                x = 1.0 * ii - 0.5 * xsize;
                y = 1.0 * jj - 0.5 * ysize;
                re = data.image[IDfpmz].array.CF[k * size2 + jj * xsize + ii].re;
                im = data.image[IDfpmz].array.CF[k * size2 + jj * xsize + ii].im;
                amp = sqrt(re * re + im * im);
                pha = atan2(im, re);

                rein = data.image[ID].array.CF[jj * xsize + ii].re;
                imin = data.image[ID].array.CF[jj * xsize + ii].im;
                ampin = sqrt(rein * rein + imin * imin);
                phain = atan2(imin, rein);

                ampin *= amp;
                total += ampin * ampin;
                phain += pha;

                data.image[ID].array.CF[jj * xsize + ii].re = ampin * cos(phain);
                data.image[ID].array.CF[jj * xsize + ii].im = ampin * sin(phain);

                tx += x * ampin * sin(phain) * ampin;
                ty += y * ampin * sin(phain) * ampin;
                tcx += x * x * ampin * ampin;
                tcy += y * y * ampin * ampin;
            }
        printf("TX TY = %.18lf %.18lf", tx / tcx, ty / tcy);
        if(FORCE_IMZERO == 1)   // Remove tip-tilt in focal plane mask imaginary part
        {
            tx = 0.0;
            ty = 0.0;
            for(ii = 0; ii < xsize; ii++)
                for(jj = 0; jj < ysize; jj++)
                {
                    x = 1.0 * ii - 0.5 * xsize;
                    y = 1.0 * jj - 0.5 * ysize;

                    re = data.image[ID].array.CF[jj * xsize + ii].re;
                    im = data.image[ID].array.CF[jj * xsize + ii].im;
                    amp = sqrt(re * re + im * im);

                    data.image[ID].array.CF[jj * xsize + ii].im -= amp * (x * tx / tcx + y * ty /
                            tcy);
                    tx += x * data.image[ID].array.CF[jj * xsize + ii].im * amp;
                    ty += y * data.image[ID].array.CF[jj * xsize + ii].im * amp;
                }
            printf("  ->   %.18lf %.18lf", tx / tcx, ty / tcy);

            mk_amph_from_complex("_foc0", "_foc0_amp", "_foc0_pha", 0);
            save_fl_fits("_foc0_amp", "!_foc_amp.fits");
            save_fl_fits("_foc0_pha", "!_foc_pha.fits");
            delete_image_ID("_foc0_amp", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("_foc0_pha", DELETE_IMAGE_ERRMODE_WARNING);
        }
        printf("\n");


        data.FLOATARRAY[0] = (float) total;

        /*  if(FORCE_IMZERO==1) // Remove tip-tilt in focal plane mask imaginary part
        {
        imresidual = 0.0;
        ID = image_ID("_foc0");
        ID1 = create_2Dimage_ID("imresidual", xsize, ysize);
        for(ii=0; ii<xsize*ysize; ii++)
        {
        data.image[ID1].array.F[ii] = data.image[ID].array.CF[ii].im;
        imresidual += data.image[ID].array.CF[ii].im*data.image[ID].array.CF[ii].im;
        data.image[ID].array.CF[ii].im = 0.0;
        }
        printf("IM RESIDUAL = %lf\n", imresidual);
        save_fl_fits("imresidual", "!imresidual.fits");
        delete_image_ID("imresidual");
        }
        */


        if(0) // TEST
        {
            /// @warning This internal test could crash the process as multiple write operations to the same filename may occurr: leave option OFF for production

            mk_amph_from_complex("_foc0", "tmp_foc0_a", "tmp_foc0_p", 0);
            save_fl_fits("tmp_foc0_a", "!_DFT_foca.fits");
            save_fl_fits("tmp_foc0_p", "!_DFT_focp.fits");
            delete_image_ID("tmp_foc0_a", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("tmp_foc0_p", DELETE_IMAGE_ERRMODE_WARNING);
        }


        /* for(ii=0; ii<xsize; ii++)
        for(jj=0; jj<ysize; jj++)
         {
         x = 1.0*ii-xsize/2;
         y = 1.0*jj-ysize/2;
         r = sqrt(x*x+y*y);
         if(r<150.0)
         data.image[IDpupin_mask].array.F[jj*xsize+ii] = 1.0;
         }*/

        fft_DFT("_foc0", "_fpmzmask", "_pupout2D", "_DFTpupmask", zfactor, 1, 0);

        //	save_fits("_DFTpupmask", "!test_DFTpupmask.fits");//TEST

        IDout2D = image_ID("_pupout2D");
        for(ii = 0; ii < xsize * ysize; ii++)
        {
            data.image[IDout].array.CF[k * xsize * ysize + ii].re =
                data.image[IDout2D].array.CF[ii].re / (xsize * ysize);
            data.image[IDout].array.CF[k * xsize * ysize + ii].im =
                data.image[IDout2D].array.CF[ii].im / (xsize * ysize);
        }
        delete_image_ID("_pupout2D", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("_foc0", DELETE_IMAGE_ERRMODE_WARNING);

        delete_image_ID("_DFTpupmask", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("_fpmzmask", DELETE_IMAGE_ERRMODE_WARNING);
    }

    return IDout;
}




//
// pupil convolution by real focal plane mask of limited support
// typically used with fpmz = zoomed copy of 1-fpm
// high resolution focal plane mask using DFT
// zoom factor
//
//
//
imageID fft_DFTinsertFPM_re(
    const char *pupin_name,
    const char *fpmz_name,
    double      zfactor,
    const char *pupout_name
)
{
    double eps = 1.0e-10;
    imageID ID;
    imageID IDpupin_mask;
    imageID IDfpmz;
    imageID IDfpmz_mask;
    long xsize, ysize;
    imageID IDin, IDout;
    long ii;
    double re, im, rein, imin, amp, ampin, phain, amp2;
    double total = 0;
    imageID ID_DFTmask00;

    IDin = image_ID(pupin_name);
    xsize = data.image[IDin].md[0].size[0];
    ysize = data.image[IDin].md[0].size[1];


    ID_DFTmask00 = image_ID("_DFTmask00");

    printf("zfactor = %f\n", zfactor);

    IDpupin_mask = create_2Dimage_ID("_DFTpupmask", xsize, ysize);
    for(ii = 0; ii < xsize * ysize; ii++)
    {
        re = data.image[IDin].array.CF[ii].re;
        im = data.image[IDin].array.CF[ii].im;
        amp2 = re * re + im * im;
        if(amp2 > eps)
        {
            data.image[IDpupin_mask].array.F[ii] = 1.0;
        }
        else
        {
            data.image[IDpupin_mask].array.F[ii] = 0.0;
        }
    }
    //  save_fl_fits("_DFTpupmask", "!_DFTpupmask.fits");

    if(ID_DFTmask00 != -1)
        for(ii = 0; ii < xsize * ysize; ii++)
        {
            if(data.image[ID_DFTmask00].array.F[ii] > 0.5)
            {
                data.image[IDpupin_mask].array.F[ii] = 1.0;
            }
        }



    IDfpmz = image_ID(fpmz_name);
    IDfpmz_mask = create_2Dimage_ID("_fpmzmask", xsize, ysize);
    for(ii = 0; ii < xsize * ysize; ii++)
    {
        amp = fabs(data.image[IDfpmz].array.F[ii]);
        if(amp > eps)
        {
            data.image[IDfpmz_mask].array.F[ii] = 1.0;
        }
        else
        {
            data.image[IDfpmz_mask].array.F[ii] = 0.0;
        }
    }

    fft_DFT(pupin_name, "_DFTpupmask", "_foc0", "_fpmzmask", zfactor, -1, 0);

    ID = image_ID("_foc0");
    total = 0.0;
    for(ii = 0; ii < xsize * ysize; ii++)
    {
        amp = data.image[IDfpmz].array.F[ii];

        rein = data.image[ID].array.CF[ii].re;
        imin = data.image[ID].array.CF[ii].im;
        ampin = sqrt(rein * rein + imin * imin);
        phain = atan2(imin, rein);

        ampin *= amp;
        total += ampin * ampin;

        data.image[ID].array.CF[ii].re = ampin * cos(phain);
        data.image[ID].array.CF[ii].im = ampin * sin(phain);
    }

    data.FLOATARRAY[0] = (float) total;


    if(1) // TEST
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        mk_amph_from_complex("_foc0", "tmp_foc0_a", "tmp_foc0_p", 0);
        WRITE_FULLFILENAME(fname, "%s/_DFT_foca", data.SAVEDIR);
        save_fl_fits("tmp_foc0_a", fname);
        WRITE_FULLFILENAME(fname, "%s/_DFT_focp", data.SAVEDIR);
        save_fl_fits("tmp_foc0_p", fname);
        delete_image_ID("tmp_foc0_a", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("tmp_foc0_p", DELETE_IMAGE_ERRMODE_WARNING);
    }

    /* for(ii=0; ii<xsize; ii++)
      for(jj=0; jj<ysize; jj++)
        {
    x = 1.0*ii-xsize/2;
    y = 1.0*jj-ysize/2;
    r = sqrt(x*x+y*y);
    if(r<150.0)
      data.image[IDpupin_mask].array.F[jj*xsize+ii] = 1.0;
      }*/

    fft_DFT("_foc0", "_fpmzmask", pupout_name, "_DFTpupmask", zfactor, 1, 0);

    IDout = image_ID(pupout_name);
    for(ii = 0; ii < xsize * ysize; ii++)
    {
        data.image[IDout].array.CF[ii].re /= xsize * ysize;
        data.image[IDout].array.CF[ii].im /= xsize * ysize;
    }

    delete_image_ID("_foc0", DELETE_IMAGE_ERRMODE_WARNING);

    delete_image_ID("_DFTpupmask", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("_fpmzmask", DELETE_IMAGE_ERRMODE_WARNING);

    return IDout;
}



