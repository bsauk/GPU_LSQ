################################################################################################
# The following commented out code is meant for later testing when I download other programs
#MAGMADIR     =  /home/bsauk/programs/magma-2.2.0
#CUDADIR      = /usr/local/cuda-7.5
#OPENBLASDIR  = /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -g -Wall
LDFLAGS       = 
NVCCFLAGS     = -O0 -g -G -m64 --gpu-architecture compute_35
LDFLAGS       = -L/usr/local/cuda/lib64/ -lcudart

#LIBBLAS     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#MAGMA_CFLAGS = -DADD_ -g -G -I$(MAGMADIR)/include -I$(CUDADIR)/include 
#MAGMA_LIBS   = -L$(MAGMADIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
		$(MAGMADIR)/lib/libmagma.so -lcublas -lcudart $(LIBBLAS) -lcusolver 
#################################################################################################
#%.o: %.cu
#	$(CUDADIR)/bin/nvcc $(MAGMA_CFLAGS) $(MAGMA_LIBS)  -c $<

#subset:	magma_sub.o dn_sub.o
#	/usr/local/cuda-7.5/bin/nvcc magma_sub.o dn_sub.o find_sub.cpp lsq.cpp selection.cpp benchmark.cpp -o sub.out -lineinfo $(MAGMA_LIBS) 
EXECUTABLE := gpusub
OBJDIR=objs
OBJS=$(OBJDIR)/subsetGPU.o $(OBJDIR)/gpuLSQ.o

default: $(EXECUTABLE)

clean:
	-rm -f *.out
	rm -rf $(OBJDIR) $(EXECUTABLE)

dirs: 
	mkdir -p $(OBJDIR)/

$(EXECUTABLE): dirs $(OBJS)
	g++ -g -o $@ $(OBJS) $(LDFLAGS)

$(OBJDIR)/%.o: %.cpp
	g++ $< -g -c -o $@

$(OBJDIR)/%.o: %.cu
	nvcc $< $(NVCCFLAGS) -c -o $@

gold: 
	g++ -g lsq.cpp find_sub.cpp selection.cpp benchmark.cpp -o gold
