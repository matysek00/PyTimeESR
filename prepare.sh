cd PyTimeESR/Calculations/fortran/

python -m numpy.f2py -c -m exec --dep lapack --dep blas declarations.f90 OpenFiles.f90 algebra.f90 spin_parameters.f90 KrondProd.f90 H_QD.f90 QME_F.f90 Matrix_Coeff.f90
mv exec.cpython-313-x86_64-linux-gnu.so exec.so

#
cd ../../..