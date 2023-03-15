//This is the program created on the Numerical Methods lab on Warsaw University on Technology in March2023
//This is Problem 2 from lab1. You can find the description here http://ccfd.github.io/courses/metnum_lab1.html

#define _CRT_SECURE_NO_WARNINGS

#include "gauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
	J[2][0] = -2*x[2] + 2*x[0]; J[2][1] = -2 * x[3] + 2 * x[1]; J[2][2] = -2 * x[0] + 2 * x[2]; J[0][3] = -2 * x[1] + 2 * x[3];
	J[3][0] = 0; J[3][1] = 0; J[3][2] = 0; J[3][3] = 0;
}

void NewtonRaphson(double* x) {


	double* F;										//ALLOCation of the Constraints vector
	F = (double*)malloc(N * sizeof(double));
	Constraints(x, F);								//CALCulation of the values of first Constraints vector

	gauss(N, J, x, F);


	/*for (int i = 0; i < 5; i++) {

		gauss(N, J, x, F);
		JacobiMatrix(J, x);
		Constraints(x, F);

	}*/
}
void initial_guess(double* x) {						//Initial guess of x's
	x[0] = 4.2;
	x[1] = 2.72;
	x[2] = 6;
	x[3] = 3.7;
}

void print_vector(int N, double* b) {				//Prints a vector
	for (int i = 0; i < N; i++)
		printf("%lf \n", b[i]);
}



int main() {
	int N = 4;

	double* x;										//ALLOCation of old and new x's
	x = (double*)malloc(N * sizeof(double));
	initial_guess(x);								//INITialization of the first approximation

	double** J, * J_row;							//ALLOCation of the Jacobi matrix
	J = (double**)malloc(N * sizeof(double*));
	J_row = (double*)malloc(N * N * sizeof(double));
	for (int i = 0; i < N; i++) {
		J[i] = &J_row[i * N];
	}
	JacobiMatrix(J, x);								//CALCulation of the first Jacobi matrix 

	NewtonRaphson(x);								//Termination of Newton Raphson iterations
	
	print_vector(N, F); 
	print_vector(N, x); 

	return 0;
}