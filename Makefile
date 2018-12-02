all:
	g++ create-ttmat.cpp -o create-ttmat
	g++ create-ttvec.cpp -o create-ttvec
	g++ compare-ttvec.cpp ttvec.cpp -lm -o compare-ttvec
	g++ ttmatvec.cpp ttmat.cpp ttvec.cpp -lm -o ttmatvec	
