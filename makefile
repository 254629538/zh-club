ifort  -O3 -qopenmp  -check bounds Tmatrix.f90 definemu.f90 gaussint.F90 gsldint.f90 solven_nonlin.F90 interpolation.f90 -I/opt/intel/compilers_and_libraries_2017.3.191/linux/mkl/include/intel64/lp64 -L/opt/intel/compilers_and_libraries_2017.3.191/linux/mkl/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -liomp5  ~/桌面/try/Cuba-4.1/libcuba.a  -g -traceback
