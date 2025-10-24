FFLAGS =  -O  #-fopenmp 
FF = -c
FC = gfortran

LIBBLAS = -lblas
LIBLAPACK = -llapack
LIBFFT = # -lfftw3

LIBS = $(LIBFFT) $(LIBLAPACK) $(LIBBLAS)

OBJ= declarations.o declfft.o OpenFiles.o timer_gnu.o algebra.o io.o spin_parameters.o KrondProd.o \
     H_QD.o QME_F.o Transport_F.o SpinTools.o Matrix_Coeff.o output.o Floquet_ESR.o 

Floquet_ESR.out: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o Floquet_ESR_v7.4.0.out

declarations.o: declarations.f90
	$(FC) $(FFLAGS) $(FF)  declarations.f90

declfft.o: declfft.f90
	$(FC) $(FFLAGS) $(FF) declfft.f90 

algebra.o: algebra.f90
	$(FC) $(FFLAGS) $(FF) algebra.f90

io.o: io.f90
	$(FC) $(FFLAGS) $(FF) io.f90

OpenFiles.o: OpenFiles.f90
	$(FC) $(FFLAGS) $(FF) OpenFiles.f90

timer_gnu.o: timer_gnu.f90
	$(FC) $(FFLAGS) $(FF) timer_gnu.f90

spin_parameters.o: spin_parameters.f90
	$(FC) $(FFLAGS) $(FF) spin_parameters.f90

KrondProd.o: KrondProd.f90
	$(FC) $(FFLAGS) $(FF) KrondProd.f90

H_QD.o: H_QD.f90
	$(FC) $(FFLAGS) $(FF) H_QD.f90

QME_F.o: QME_F.f90
	$(FC) $(FFLAGS) $(FF) QME_F.f90

Matrix_Coeff.o: Matrix_Coeff.f90
	$(FC) $(FFLAGS) $(FF) Matrix_Coeff.f90

Transport_F.o: Transport_F.f90
	$(FC) $(FFLAGS) $(FF) Transport_F.f90

SpinTools.o: SpinTools.f90
	$(FC) $(FFLAGS) $(FF) SpinTools.f90

output.o: output.f90
	$(FC) $(FFLAGS) $(FF) output.f90

Floquet_ESR.o: Floquet_ESR.f90
	$(FC) $(FFLAGS) $(FF) Floquet_ESR.f90

clean:
	@echo "CLEANING up!"
	rm -f *.o *.mod *.out

