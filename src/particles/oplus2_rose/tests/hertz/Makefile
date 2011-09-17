# Added by Nao for local setting
HERTZ_COMMONDIR :=../../../hertz_framework/common

COMMONDIR := ../common
OP2ROSE := ../../translator/op2rose

NVCC := nvcc
#Modified by Nao
#INCLUDEFLAGS := -I$(COMMONDIR) -I$(HERTZ_COMMONDIR)
INCLUDEFLAGS := -I$(COMMONDIR) -I$(HERTZ_COMMONDIR) -I./

#Modified by Nao
#NVCCFLAGS := -arch sm_13 $(INCLUDEFLAGS) -D OP_DIAGS=0 #-g -O2
LIBS := -L$(HERTZ_COMMONDIR)
NVCCFLAGS := -arch sm_13 $(INCLUDEFLAGS) $(LIBS) -D OP_DIAGS=0 #-g -O2

VPATH = $(COMMONDIR)

#SOURCES = rose_hertz.cpp op_par.cu cutil.cpp \
#	cmd_arg_reader.cpp stopwatch.cpp bank_checker.cpp stopwatch_linux.cpp $(HERTZ_COMMONDIR)/libcommon.so \

SOURCES = rose_hertz.cpp op_par.cu cutil.cpp \
	cmd_arg_reader.cpp stopwatch.cpp bank_checker.cpp stopwatch_linux.cpp

all: src2src build

build: $(SOURCES)
	$(NVCC) $(NVCCFLAGS) -DBSIZE=512 $^ -o op2_cpu_timer -lcommon
	cp rose_hertz.cpp rose_hertz.cu
	$(NVCC) $(NVCCFLAGS) -D GPU_TIMER \
    op_par.cu ../common/cutil.cpp \
	  ../common/cmd_arg_reader.cpp ../common/stopwatch.cpp ../common/bank_checker.cpp ../common/stopwatch_linux.cpp \
    $(HERTZ_COMMONDIR)/libcommon.so \
    rose_hertz.cu -o op2_gpu_timer
	$(NVCC) $(NVCCFLAGS) -lrt $^ -o op2_posix_timer -lcommon

#    scp op2_cpu_timer ni10@login-0.cx1.hpc.ic.ac.uk:/work/ni10/particles/oplus2_rose/tests/hertz
#    scp op2_gpu_timer ni10@login-0.cx1.hpc.ic.ac.uk:/work/ni10/particles/oplus2_rose/tests/hertz
#    scp op2_posix_timer ni10@login-0.cx1.hpc.ic.ac.uk:/work/ni10/particles/oplus2_rose/tests/hertz


src2src: clean
	@if [ -z "${HERTZ_COMMONDIR}" ]; then echo "ERROR: set HERTZ_COMMONDIR env variable first."; exit 1; fi
	$(OP2ROSE) $(INCLUDEFLAGS) hertz.cpp

clean: 
	rm -f op_par.cu kernels.h res_kernel.cu update_kernel.cu rose_hertz.cpp \
	op2_cpu_timer op2_gpu_timer op2_posix_timer rose_hertz.cu blank.cpp
