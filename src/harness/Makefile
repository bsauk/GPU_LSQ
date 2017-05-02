MAGMADIR     =  /home/bsauk/programs/magma-2.2.0
CUDADIR      = /usr/local/cuda-7.5
OPENBLASDIR  = /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -g -Wall
LDFLAGS       = 

LIBBLAS     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

MAGMA_CFLAGS = -DADD_ -g -G -I$(MAGMADIR)/include -I$(CUDADIR)/include 
MAGMA_LIBS   = -L$(MAGMADIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
		$(MAGMADIR)/lib/libmagma.so -lcublas -lcudart $(LIBBLAS) -lcusolver 

%.o: %.cu
	$(CUDADIR)/bin/nvcc $(MAGMA_CFLAGS) $(MAGMA_LIBS)  -c $<

subset:	magma_sub.o dn_sub.o
	/usr/local/cuda-7.5/bin/nvcc magma_sub.o dn_sub.o find_sub.cpp lsq.cpp selection.cpp benchmark.cpp -o sub.out -lineinfo $(MAGMA_LIBS) 

clean:
	-rm -f *.out
	-rm -f *.o

all:
	g++ -g lsq_demo.cpp lsq.cpp -o sequential_gold

#subset:
#	g++ -g find_sub.cpp lsq.cpp selection.cpp benchmark.cpp -o sub


subset2:
	g++ -g find_sub.cpp lsq.cpp selection.cpp benchmark.cpp -o test.out

gpusub: gpuLSQ.o
	$(CUDADIR)/bin/nvcc -g find_sub.cpp lsq.cpp subsetGPU.cpp benchmark.cpp gpuLSQ.o -o gpusub.out

