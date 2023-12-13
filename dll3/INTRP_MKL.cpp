#include "pch.h"
#include <mkl.h>
#include <mkl_df_functions.h>


extern "C" __declspec(dllexport)
void SplineInterpolation(MKL_INT nx, MKL_INT ny, double* x, double* y, double* scoeff, MKL_INT nsite, double* site, MKL_INT ndorder, MKL_INT * dorder, double* result, int& ret)
{
	try
	{
		int status;
		DFTaskPtr task;
		status = dfdNewTask1D(&task, nx, x, DF_UNIFORM_PARTITION, ny, y, DF_MATRIX_STORAGE_ROWS);
		if (status != DF_STATUS_OK) { ret = -1; return; }
		status = dfdEditPPSpline1D(task, DF_PP_CUBIC, DF_PP_NATURAL, DF_BC_FREE_END, NULL, DF_NO_IC, NULL, scoeff, DF_NO_HINT);
		if (status != DF_STATUS_OK) { ret = -1; return; }
		status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
		if (status != DF_STATUS_OK) { ret = -1; return; }
		status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP, nsite, site, DF_UNIFORM_PARTITION, ndorder, dorder, NULL, result, DF_MATRIX_STORAGE_ROWS, NULL);
		if (status != DF_STATUS_OK) { ret = -1; return; }
		status = dfDeleteTask(&task);
		if (status != DF_STATUS_OK) { ret = -1; return; }

		ret = 0;
	}
	catch (...)
	{
		ret = -1;
	}
}