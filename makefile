CXX=/usr/local/cuda/bin/nvcc
CPPFLAGS=-lm -arch=compute_70 -code=sm_70,sm_72 -lineinfo -lnvToolsExt

srun_flags=-p gpu -o log/out.log -e log/error.log
jobid_find=s/.*job \([[:digit:]]*\) queued.*/\1/p

REAL=data/real/data_100k_arcmin.txt
FAKE=data/real/flat_100k_arcmin.txt

all : darkmatter

darkmatter: src/darkmatter.cu src/*.h
	$(CXX) $(CPPFLAGS) -o $@ src/darkmatter.cu

.PHONY: run
run: darkmatter
	@mkdir -p log
	@mkdir -p out
	$(eval PROG_ID=$(shell /bin/bash -c "srun $(srun_flags) $^ $(REAL) $(FAKE) 2>&1 | tee /dev/tty |sed -n '$(jobid_find)'"))
	@./report_run.sh $(PROG_ID)

.PHONY: debug
debug: darkmatter
	@mkdir -p log
	srun -p gpu --pty cuda-gdb --args $^ $(REAL) $(FAKE)

.PHONY: time
time: darkmatter
	srun -p gpu --pty time ./$^ $(REAL) $(FAKE)

.PHONY: prof
prof: darkmatter
	srun -p gpu --pty nvprof ./$^ $(REAL) $(FAKE)
	
.PHONY: vprof
vprof: darkmatter
	srun -p gpu --pty nvprof -f -o darkmatter.prof --cpu-profiling on ./$^ $(REAL) $(FAKE)
	#srun -p gpu --pty nvprof -f -o darkmatter.prof --analysis-metrics ./$^ $(REAL) $(FAKE)
.PHONY : clean
clean: 
	rm -rf build
	rm -f src/**/*.o
	rm -rf  log
