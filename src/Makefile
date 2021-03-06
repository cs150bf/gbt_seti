FC = gfortran
FFLAGS = -g -fPIC
OS = Linux
BINDIR = $(SETI_GBT)/bin
CUDA = /usr/local/cuda
CUDAC = $(CUDA)/bin/nvcc
# Linux is the first choice
ifeq ($(OS),Linux)
        LIBSUFFIX = .so
        LIBCMD = -shared
# else assume Darwin (i.e. OSX)
else
        LIBSUFFIX = .dylib
        LIBCMD = -dynamiclib
endif

CFLAGS = -O3 -Wall -DFOLD_USE_INTRINSICS -I$(OPT64)/include -I/usr/include/cfitsio -I/usr/include/mysql -I$(CUDA)/include 
PROGS = test_psrfits_read \
	raw_quant raw_read rawdopplersearch raw2wu \
	treedop raw_head two2eight rawfilterbank rawfilterbank_singleband_xy raw_ascii2 
OBJS  = guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
	guppi_params.o guppi_time.o median.o setimysql.o\
	write_psrfits.o read_psrfits.o misc_utils.o \
	fold.o polyco.o hget.o hput.o sla.o downsample.o barycenter.o chkio.o
LIBS = -L$(OPT64)/lib -I$(OPT64)/include -L$(CUDA)/lib64 -lfftw3f_threads -lcfitsio -L$(SETI_GBT)/lib -lsla -lm -lgsl -lgslcblas -lpthread -I/usr/include/cfitsio -lfftw3f -L/usr/lib64/mysql -lmysqlclient -lcufft -lcudart
INSTALL_DIR = ../bin

all: slalib $(PROGS) gpu
clean: 
	rm -f $(PROGS) *~ *.o
	exit 0
install:$(PROGS)
	mkdir -p $(INSTALL_DIR) && \
	cp -f $(PROGS) $(INSTALL_DIR)
slalib: libsla$(LIBSUFFIX)
	cd slalib ; $(FC) -o sla_test sla_test.f -fno-second-underscore -L$(SETI_GBT)/lib -lsla
	slalib/sla_test
libsla$(LIBSUFFIX):
	@if [ -z "${SETI_GBT}" ]; \
	then echo "Error: environment variable SETI_GBT is not set!"; exit 2; \
	else echo "SETI_GBT="${SETI_GBT}; fi
	cd slalib ; $(FC) $(FFLAGS) -fno-second-underscore -c -I. *.f *.F
	rm slalib/sla_test.o
	cd slalib ; $(FC) $(LIBCMD) -o $(SETI_GBT)/lib/libsla$(LIBSUFFIX) -fno-second-underscore *.o
gpu: 
	$(CUDAC) -c guppi2spectra_gpu.cu -o guppi2spectra_gpu.o
	$(CC) $(CFLAGS) -fno-strict-aliasing guppi2spectra.c -o $(BINDIR)/gpuspec $(OBJS) guppi2spectra_gpu.o $(LIBS) $(THREAD_LIBS)
	$(CC) $(CFLAGS) -fno-strict-aliasing guppi2wu.c -o $(BINDIR)/guppi2wu $(OBJS) guppi2spectra_gpu.o $(LIBS) $(THREAD_LIBS)

.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $(BINDIR)/$@ $(OBJS) $(LIBS) $(THREAD_LIBS)
