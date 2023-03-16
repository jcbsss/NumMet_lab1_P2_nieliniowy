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

void JacobiMatrix(double** J, double* x) {		//METHOD1//Calculates Jacobi matrix analitically for a given x-vector
	J[0][0] = 1; J[0][1] = 0; J[0][2] = 0; J[0][3] = 0;
	J[1][0] = 0; J[1][1] = 1; J[1][2] = 0; J[1][3] = 0;
	J[2][0] = -2*x[2] + 2*x[0]; J[2][1] = -2 * x[3] + 2 * x[1]; J[2][2] = -2 * x[0] + 2 * x[2]; J[2][3] = -2 * x[1] + 2 * x[3];
	J[3][0] = 0; J[3][1] = 0; J[3][2] = -6 + 2*x[2]; J[3][3] = 2*x[3];
}
void JacobiMatrixFD(double** J, double* x, int N) {  //METHOD2//Calculates appoximated Jacobi matrix, using Central Finite Difference method
	double epsilon = 1e-8;
	double* x_l;									//ALLOCation of x_LEFT vector
	x_l = (double*)malloc(N * sizeof(double));
	double* x_r;									//ALLOCation of x_RIGHT vector
	x_r = (double*)malloc(N * sizeof(double));
	double* F_l;									//ALLOCation of F_LEFT vector
	F_l = (double*)malloc(N * sizeof(double));
	double* F_r;									//ALLOCation of F_RIGHT vector
	F_r = (double*)malloc(N * sizeof(double));
	
	for (int j = 0; j < N; j++) {
		for (int k = 0; k < N; k++) {
			if (k != j) { 
				x_l[k] = x[k]; 
				x_r[k] = x[k];}	//copied vector of x, except of one variable, because we are taking only partial derivative
			else {
				x_l[j] = x[j] - epsilon; 
				x_r[j] = x[j] + epsilon;} //left and right value of the j-th variable, needed for finite difference
			}
		Constraints(x_l, F_l);
		Constraints(x_r, F_r);
		for (int i = 0; i < N; i++) {
			J[i][j] = (F_r[i] - F_l[i]) / (2*epsilon);
		}
	}
	
	free(x_l);
	free(x_r);
	free(F_l);
	free(F_r);
}

void NewtonRaphson(double* x, int N, int n_iter, int jacobian_meth) {	//Whole procedure of the Newton Raphson method
	
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
		switch (jacobian_meth) {
		case 1: JacobiMatrix(J, x);	break;		//CALCulation of the first Jacobi matrix analytically
		case 2: JacobiMatrixFD(J, x, N); break;	//CALCulation of the first Jacobi matrix using finite difference method
		}

		Constraints(x, F);							//CALCulation of the values of first Constraints vector
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
	int n_iter = 5;									//Number of iterations
	int N = 4;										//Size of the system of equations

	double* x;										//ALLOCation of x vector
	x = (double*)malloc(N * sizeof(double));
	initial_guess(x);								//INITialization of the first approximation

	NewtonRaphson(x, N, n_iter, 2);					//Termination of Newton Raphson iterations
	free(x);

	return 0;
}