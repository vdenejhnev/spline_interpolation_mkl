#pragma once
void getValues(int* m, int* n, double* x, double* y)
{
	for (int i = 0; i <= *n; i++) {
		y[i] = x[i] * x[i];
	}

	return;
}