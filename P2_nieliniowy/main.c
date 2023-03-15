//This is the program created on the Numerical Methods lab on Warsaw University on Technology in March2023
//This is Problem 2 from lab1. You can find the description here http://ccfd.github.io/courses/metnum_lab1.html

#define _CRT_SECURE_NO_WARNINGS

#include "gauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void initial_guess(double* x) {						//Initial guess of x's
	x[0] = 2.72;
	x[1] = 4.2;
	x[2] = 3;
	x[3] = 6;
}

void print_vector(int N, double* b) {
	//printf("Vect:	");
	for (int i = 0; i < N; i++)
		printf("%.17lf ", b[i]);			//The precision is set
	printf("\n");
	
}

void print_matrix(int N, double** A) {
	printf("Matrix:\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%.2lf\t", A[i][j]);
		printf("\n");
	}
}

void Constraints(double* x, double* F) {			//Calculate the values of the constraints functions
		double alpha = 1;
		F[0] = x[0] - 5 * cos(alpha);
		F[1] = x[1] - 5 * sin(alpha);
		F[2] = x[2] * x[2] - 2 * x[0] * x[2] + x[0] * x[0] + x[3] * x[3] - 2 * x[1] * x[3] + x[1] * x[1] - 4;
		F[3] = -6 * x[2] + x[2] * x[2] + x[3] * x[3] - 27;
}

void JacobiMatrix(double** J, double* x) {			//Calculates Jacobi matrix analitically for a given x-vector
	J[0][0] = 1; J[0][1] = 0; J[0][2] = 0; J[0][3] = 0;
	J[1][0] = 0; J[1][1] = 1; J[1][2] = 0; J[1][3] = 0;
	J[2][0] = -2*x[2] + 2*x[0]; J[2][1] = -2 * x[3] + 2 * x[1]; J[2][2] = -2 * x[0] + 2 * x[2]; J[2][3] = -2 * x[1] + 2 * x[3];
	J[3][0] = 0; J[3][1] = 0; J[3][2] = -6 + 2*x[2]; J[3][3] = 2*x[3];
}


void NewtonRaphson(double* x, int N, int n_iter) {				//Whole procedure of the Newton Raphson method
	
	double** J, * J_row;							//ALLOCation of the Jacobi matrix
	J = (double**)malloc(N * sizeof(double*));
	J_row = (double*)malloc(N * N * sizeof(double));
	for (int i = 0; i < N; i++) {
		J[i] = &J_row[i * N];
	}
	double* F;										//ALLOCation of the Constraints vector
	F = (double*)malloc(N * sizeof(double));
	double* h;										//ALLOCation of the h vector
	h = (double*)malloc(N * sizeof(double));

	for (int i = 0; i < n_iter; i++) {
		JacobiMatrix(J, x);								//CALCulation of the first Jacobi matrix 
		Constraints(x, F);								//CALCulation of the values of first Constraints vector
		for (int k = 0; k < N; k++)		F[k] = -F[k];

		gauss(N, J, h, F);
		for (int k = 0; k < N; k++)		x[k] = x[k] + h[k];

		//print_vector(N, x);
		//print_matrix(N, J);
		//print_vector(N, F);
		//print_vector(N, h);
		
	}
	printf("\nThe calculated x vector equals to:	"); print_vector(N, x);
	printf("\nThe final residua vector equals to:	");  print_vector(N, F);

	free(F);
	free(h);
	free(J);
	free(J_row);
}




int main() {
	int n_iter = 8;									//Number of iterations
	int N = 4;										//Size of the system of equations

	double* x;										//ALLOCation of x vector
	x = (double*)malloc(N * sizeof(double));
	initial_guess(x);								//INITialization of the first approximation

	NewtonRaphson(x, N, n_iter);					//Termination of Newton Raphson iterations
	free(x);

	return 0;
}