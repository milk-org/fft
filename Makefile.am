lib_LTLIBRARIES = libfft.la
libfft_la_SOURCES = fft.c fft.h

#libfft_la_CFLAGS = -fopenacc -foffload=nvptx-none

#AM_CFLAGS = -fopenacc -foffload=nvptx-none

AM_CPPFLAGS = -I@abs_top_srcdir@/src -DFFTCONFIGDIR=\"@abs_top_srcdir@/src/fft/config\" -lfftw3_omp
