#pragma once
#include "D:\Intel\oneAPI\mkl\latest\include\mkl.h"

extern "C"  _declspec(dllexport)
void SplineInterpolation(MKL_INT nx, MKL_INT ny, double* x, double* y, double* scoeff, MKL_INT nsite, double* site, MKL_INT ndorder, MKL_INT * dorder, double* approximation, int maxiter, int maxiter_step, int& stopCriteria, double& resFinal, int& ndoneIter, double rs, double* result, double& ret);