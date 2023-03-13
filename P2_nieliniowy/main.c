//This is the program created on the Numerical Methods lab on Warsaw University on Technology in March2023
//This is Problem 2 from lab1. You can find the description here http://ccfd.github.io/courses/metnum_lab1.html

#define _CRT_SECURE_NO_WARNINGS

#include "gauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Constraints(double* x, double* F);		//To do

void JacobiMatrix(double** J, double* x) {
	J[0][0] = 1; J[0][1] = 0; J[0][2] = 0; J[0][3] = 0;
	J[1][0] = 0; J[1][1] = 1; J[1][2] = 0; J[1][3] = 0;
	J[2][0] = -2*x[2] + 2*x[0]; J[2][1] = -2 * x[3] + 2 * x[1]; J[2][2] = -2 * x[0] + 2 * x[2]; J[0][3] = -2 * x[1] + 2 * x[3];
	J[3][0] = 0; J[3][1] = 0; J[3][2] = 0; J[3][3] = 0;
}

void NewtonRaphson(double* x);		//To do

void initial_guess(double* x) {		     //Initial guess of x's
	x[0] = 2;
	x[1] = 3.5;
	x[2] = 2.7;
	x[3] = 5.5;
}

void print_vector(int N, double* b) {
	for (int i = 0; i < N; i++)
		printf("%lf \n", b[i]);
}



int main() {
	int N = 4;

	double* x_old; double* x_new;					//allocation of old and new x's
	x_old = (double*)malloc(N * sizeof(double));
	x_new = (double*)malloc(N * sizeof(double));
	initial_guess(x_old);							//initialization of the first approximation

	double** J, * J_row;							//allocation of the Jacobi matrix
	J = (double**)malloc(N * sizeof(double*));
	J_row = (double*)malloc(N * N * sizeof(double));
	for (int i = 0; i < N; i++) {
		J[i] = &J_row[i * N];
	}
	JacobiMatrix(J, x_old);							//initiation of the first Jacobi matrix





	

	print_vector(N, x_new);

	return 0;
}