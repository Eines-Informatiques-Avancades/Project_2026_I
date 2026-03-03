main : constants.o nonbonded.o initconfig.o bonded.o energy.o main.o
	gfortran -o main main.o constants.o nonbonded.o initconfig.o bonded.o energy.o

main.o : main.f90 constants.o nonbonded.o initconfig.o bonded.o energy.o
	gfortran -c main.f90

constants.o : constants.f90
	gfortran -c constants.f90

nonbonded.o : nonbonded.f90 constants.o
	gfortran -c nonbonded.f90

bonded.o : bonded.f90 constants.o
	gfortran -c bonded.f90

energy.o : constants.o nonbonded.o bonded.o
	gfortran -c energy.f90

initconfig.o : initconfig.f90 constants.o
	gfortran -c initconfig.f90