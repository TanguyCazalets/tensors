CPP_FLAGS=-O3 -march=native -mtune=native -std=c++11 -fopenmp -mavx

all: create-ttmat create-ttvec compare-ttvec ttmatvec ttmatvec-opti ttmatvec-omp ttmatvec-task

ttmatvec: ttmatvec.cpp ttmat.cpp ttmat.h ttvec.cpp ttvec.h
	g++ ${CPP_FLAGS} ttmatvec.cpp ttmat.cpp ttvec.cpp -lm -o ttmatvec	

ttmatvec-opti: ttmatvec.cpp ttmat_opti_max.cpp ttmat.h ttvec.cpp ttvec.h
	g++ ${CPP_FLAGS} ttmatvec.cpp ttmat_opti.cpp ttvec.cpp -lm -o ttmatvec_opti

ttmatvec-omp: ttmatvec.cpp ttmat_omp.cpp ttmat.h ttvec.cpp ttvec.h
	g++ ${CPP_FLAGS} ttmatvec.cpp ttmat_omp.cpp ttvec.cpp -lm -o ttmatvec_omp

ttmatvec-task: ttmatvec.cpp ttmat_task.cpp ttmat.h ttvec.cpp ttvec.h
	g++ ${CPP_FLAGS} ttmatvec.cpp ttmat_task.cpp ttvec.cpp -lm -o ttmatvec_task

create-ttmat: create-ttmat.cpp
	g++ ${CPP_FLAGS} create-ttmat.cpp -o create-ttmat

create-ttvec: create-ttvec.cpp
	g++ ${CPP_FLAGS} create-ttvec.cpp -o create-ttvec

compare-ttvec: compare-ttvec.cpp ttvec.cpp
	g++ ${CPP_FLAGS} compare-ttvec.cpp ttvec.cpp -lm -o compare-ttvec

clean:
	rm -f ttmatvec ttmatvec-opti ttmatvec-omp ttmatvec-task create-ttmat create-ttvec compare-ttvec
