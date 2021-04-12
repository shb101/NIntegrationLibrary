#include <stdio.h>
#include <stdlib.h>

// Author: Shobhit Bhatnagar
// Date: 12/04/2021

// 1D Numerical Integration with the Trapezoidal Rule
__declspec(dllexport) double NIntTrapz(double x0, double x1, double* F, int N) {
	int i;
	double I, dx = (x1 - x0) / N;
	I = 0;
	for (i = 0; i < N; i++)
		I += dx * (F[i] + F[i + 1]) / 2;
	return I;
}

// 1D Numerical Integration with the Simpson's Rule
__declspec(dllexport) double NIntSimpson(double x0, double x1, double* F, int N) {
	int i;
	double I, dx = (x1 - x0) / N;
	I = 0;
	for (i = 0; i < N / 2; i++)
		I += dx * (F[2 * i] + 4 * F[2 * i + 1] + F[2 * i + 2]) / 3;
	if (N % 2 != 0) {
		I += dx * (F[N - 1] + F[N]) / 2;
	}
	return I;
}

// 2D Numerical Integration with the Trapezoidal Rule
__declspec(dllexport) double NIntTrapz2D(double x0, double x1, double* y0, double* y1, double** F, int N, int M) {
	int i, j;
	double I, dx = (x1 - x0) / N, dy, dA;
	I = 0;
	for (i = 0; i < N; i++) {
		dy = (y1[i] - y0[i]) / M;
		dA = dy * dx;
		for (j = 0; j < M; j++) {
			I += dA * (F[i][j] + F[i + 1][j] + F[i][j+1] + F[i + 1][j+1]) / 4;
		}
	}
	return I;
}

// 3D Numerical Integration with the Trapezoidal Rule
__declspec(dllexport) double NIntTrapz3D(double x0, double x1, double* y0, double* y1, double** z0, double** z1, double*** F, int N, int M, int L) {
	int i, j, k;
	double I, dx = (x1 - x0) / N, dy, dz, dV;
	I = 0;
	for (i = 0; i < N; i++) {
		dy = (y1[i] - y0[i]) / M;
		for (j = 0; j < M; j++) {
			dz = (z1[i][j] - z0[i][j]) / L;
			dV = dx * dy * dz;
			for (k = 0; k < L; k++)
				I += dV * (F[i][j][k] + F[i + 1][j][k] + F[i][j + 1][k] + F[i + 1][j + 1][k]
					+ F[i][j][k + 1] + F[i + 1][j][k + 1] + F[i][j + 1][k + 1] + F[i + 1][j + 1][k + 1]) / 8;
		}
	}
	return I;
}

// 2D Numerical Integration with the Simpson's Rule
__declspec(dllexport) double NIntSimpson2D(double x0, double x1, double* y0, double* y1, double** F, int N, int M) {
	int i;
	double I, * I1;
	I1 = (double*)malloc((N + 1) * sizeof(double));
	for (i = 0; i <= N; i++) {
		I1[i] = NIntSimpson(y0[i], y1[i], F[i], M);
	}
	I = NIntSimpson(x0, x1, I1, N);
	free(I1);
	return I;
}

// 3D Numerical Integration with Simpson's Rule
__declspec(dllexport) double NIntSimpson3D(double x0, double x1, double* y0, double* y1, double** z0, double** z1, double*** F, int N, int M, int L) {
	int i, j;
	double I, * I2, * I1;
	I2 = (double*)malloc((N + 1) * sizeof(double));
	I1 = (double*)malloc((M + 1) * sizeof(double));
	for (i = 0; i <= N; i++) {
		for (j = 0; j <= M; j++) {
			I1[j] = NIntSimpson(z0[i][j], z1[i][j], F[i][j], L);
		}
		I2[i] = NIntSimpson(y0[i], y1[i], I1, M);
	}
	I = NIntSimpson(x0, x1, I2, N);
	free(I2);
	free(I1);
	return I;
}