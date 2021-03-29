#include <iostream>

// 1D Numerical Integration with the Trapezoidal Rule
double NIntTrapz(double x0, double x1, double* F, int N) {
	int i;
	double I, dx = (x1 - x0) / N;
	I = 0;
	for (i = 0; i < N; i++)
		I += dx * (F[i] + F[i+1]) / 2;
	return I;
}

// 1D Numerical Integration with the Simpson's Rule
double NIntSimpson(double x0, double x1, double* F, int N) {
	int i;
	double I, dx = (x1 - x0) / N;
	I = 0;
	for (i = 0; i < N/2; i++)
		I += dx * (F[2*i] + 4 * F[2*i + 1] + F[2*i + 2])/3;
	if (N % 2 != 0) {
		I += dx * (F[N - 1] + F[N]) / 2;
	}
	return I;
}

// 2D Numerical Integration with Trapezoidal Rule
double NIntTrapz2D(double x0, double x1, double y0, double y1, double** F, int N, int M) {
	int i, j;
	double dx, dy, I;
	dx = (x1 - x0) / N;
	dy = (y1 - y0) / M;
	I = 0;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++)
			I += dx * dy * (F[i][j] + F[i + 1][j] + F[i][j+1]+F[i+1][j+1]) / 4;
	}
	return I;
}

void Read1DArray(int N, double* F) {
	int i;
	for (i = 0; i <= N; i++)
		std::cin >> F[i];
}

int main() {
	int algo, N, i;
	double x0, x1, * F, I;
	std::cin >> algo >> N;
	F = new double[N + 1];

	// Read the upper and lower bounds
	std::cin >> x0 >> x1;
	
	// Read the function
	Read1DArray(N, F);
	// If algo = 0 then algorithm selected is Trapezoidal rule; if algo = 1, then it is Simpson's rule
	if (algo == 0)
		I = NIntTrapz(x0, x1, F, N);
	else if (algo == 1)
		I = NIntSimpson(x0, x1, F, N);
	// Print the result
	std::cout << I << std::endl;
	// Free the memory
	delete[] F;
	return 0;
}