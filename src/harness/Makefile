################################################################################################
# The following commented out code is meant for later testing when I download other programs
# The NVCCFLAGS and NVCC lines may need to be edited based on computer architecture being used. 

CC            = gcc
LD            = gcc
#CFLAGS        = -g -Wall
CFLAGS        = -O3 -Wall
LDFLAGS       = 
#NVCCFLAGS     = -O0 -g -G -m64 --gpu-architecture compute_35
NVCCFLAGS     = -O3 -m64 --gpu-architecture compute_35
LDFLAGS       = -L/usr/local/cuda/lib64/ -lcudart
NVCC          = /usr/local/cuda-7.5/bin/nvcc #This is used for testing on own computer
#NVCC         = nvcc  # This is used for gates and latedays
#LIBBLAS     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#################################################################################################
#%.o: %.cu
#	$(CUDADIR)/bin/nvcc $(MAGMA_CFLAGS) $(MAGMA_LIBS)  -c $<

EXECUTABLE := gpusub
OBJDIR=objs
OBJS=$(OBJDIR)/subsetGPU.o $(OBJDIR)/gpuLSQ.o $(OBJDIR)/lsq.o $(OBJDIR)/find_sub.o

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
	$(NVCC) $< $(NVCCFLAGS) -c -o $@

gold: 
	g++ -g lsq.cpp find_sub.cpp selection.cpp benchmark.cpp -o gold
