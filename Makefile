all:
	g++ create-ttmat.cpp -o create-ttmat
	g++ create-ttvec.cpp -o create-ttvec
	g++ ttmatvec.cpp ttmat.cpp ttvec.cpp -o ttmatvec	
