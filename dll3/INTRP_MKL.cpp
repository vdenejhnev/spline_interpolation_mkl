#include "pch.h"
#include <mkl.h>
#include <mkl_df_functions.h>
#include "get_values.h"

extern "C" __declspec(dllexport)

void SplineInterpolation(MKL_INT nx, MKL_INT ny, double* x, double* y, double* scoeff, MKL_INT nsite, double* site, MKL_INT ndorder, MKL_INT * dorder, double* approximation, int maxiter, int maxiter_step, int& stopCriteria, double &resFinal, int& ndoneIter, double rs, double* result, double& ret)
{

	int status;
	DFTaskPtr task;
	_TRNSP_HANDLE_t handle = NULL;

	try 
	{
		const double eps[6] = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12 };
		MKL_INT RCI_Request = 0;
		MKL_INT checkInfo[4];
		MKL_INT stopCriteria = NULL;
		double resInitial = NULL;
		double* fjac = new double[nx * ny];
		double jac_eps = 1.0E-8;

		MKL_INT status = dtrnlsp_init(&handle, &nx, &ny, approximation, eps, &maxiter, &maxiter_step, &rs);
		if (status != TR_SUCCESS) { ret = -1; return; }
		status = dtrnlsp_check(&handle, &nx, &ny, x, y, eps, checkInfo);
		if (status != TR_SUCCESS) { ret = -1; return; }


		while (true) 
		{
			try
			{
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

			ret = dtrnlsp_solve(&handle, y, x, &RCI_Request);
			if (status != TR_SUCCESS) { ret = -1; return; }

			if (RCI_Request == 0) continue;
			else if (RCI_Request == 1) {
				getValues(&ny, &nx, approximation, y);
			}
			else if (RCI_Request == 2)
			{
				ret = djacobi(getValues, &nx, &ny, fjac, approximation, &jac_eps);
				if (ret != TR_SUCCESS) { ret = -1; return; }
			}
			else if (RCI_Request >= -6 && RCI_Request <= -1) break;
			else { ret = -1; return; }
		}

			 
		ret = dtrnlsp_get(&handle, &ndoneIter, &stopCriteria ,&resInitial, &resFinal);
		{ ret = -1; return; }

		ret = dtrnlsp_delete(&handle);
		{ ret = -1; return; }

		if (fjac != NULL) delete[] fjac;
		ret = -0;
	}
	catch (...) 
	{
		ret = -1;
	}
}