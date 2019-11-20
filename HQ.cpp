/***********************************************************************
!SUBMITTED WITH THE MANUSCRIPT:                                        *
!IMPROVING COMPUTATIONAL EFFICIENCY IN LARGE LINEAR INVERSE PROBLEMS   *
!AN EXAMPLE FROM CARBON DIOXIDE FLUX ESTIMATION                        *
!VINEET YADAV[1] AND ANNA M. MICHALAK[1]                               *
![1]{DEPARTMENT OF GLOBAL ECOLOGY, CARNEGIE INSTITUTION FOR SCIENCE,   *
!STANFORD, CALIFORNIA, USA 94305}                                      *
!***********************************************************************
!WRITTEN BY : VINEET YADAV                                             *
!DATE : 2/28/2013                                                      *
!PURPOSE : COMPUTE HQ FROM THE DIRECT AN INDIRECT METHOD WHERE Q CAN   *
!DECOMPOSED AS A KRONECKER PRODUCT OF TWO MATRICES                     *
!NOTE: THIS IS A COMPLETELY SERIAL FUNCTION                            *
!***********************************************************************
!COMPILATION PROCEDURE:						       					   *
!TESTED WITH COMPILERS: INTEL (ICPC), GCC AND MICROSOFT COMPILERS      *
!THERE ARE NO MEMORY LEAKS IN THIS PROGRAM			                   *
!*****************NOTES*************************************************
!                                                                      *
!NOTE 1: THIS CODE DEMONSTRATES MATRIX MULTIPLICATION BETWEEN AN       *
!ARBITRARY MATRIX H AND AN ARBITRARY SQUARE SYMMETRIC MATRIX Q THAT CAN*
!BE EXPRESSED AS A KRONECKER PRODUCT OF TWO MATRICES.                  *
!                                                                      *
!NOTE 2: THIS CODE ASSUMES THAT D AND E MATRIX ARE SYMMETRIC WHICH IS  *
!TRUW IF Q IS A COVARIANCE MATRIX                                      *
!                                                                      *
!NOTE 3: THIS IS A PARTLY SERIAL AND PARTLY PARALLEL PROGRAM . EXCEPT  *
!MATRIX MULTIPLICATION AND TRANSPOSE FUNCTIONS WHICH UTILIZES OPENMP   *
!OTHERS FUCNTIONS ARE SERIAL FUNCTIONS.COMPUTATIONAL EFFICIENCY OF THE *
!INDIRECT METHOD IN COMPARISON TO THE DIRECT METHOD IS TESTED USING    *
!RANDOM NUMBERS. TO MAKE THIS PROGRAM COMPLETELY SERIAL DO NOT COMPILE *
!WITH OPENMP DIRECTIVES	                                               *
!                                                                      *
!NOTE 4: ALL VARIABLES ARE AS DEFINED IN THE MANUSCRIPT                *
!                                                                      *
!NOTE 5: RATHER THAN THIS CODE MATLAB PROTOTYPE CODE SHOULD BE USED TO *
!UNDERSTAND THE EQUATIONS IN THE MANUSCRIPT AND THERE IMPLEMENTATION.  *
!                                                                      *
!NOTE 6: SINCE WE DO NOT HAVE ANY CONTROL OVER HOW MATLAB MANAGES      *
!MEMORY AND THE ITS EFFICIENCY IN IMPLEMENTING FOR LOOPS IN COMPARISON *
!TO VECTOR OPERATIONS; WE HAVE PROVIDED THE C++ AND FORTRAN CODE TO    *
!PRACTICALLY SHOW THE COMPUTATIONAL EFFICIENCY OF THE PROPOSED         *
!MATRIX MULTIPLICATION METHOD                                          *
!                                                                      *
!NOTE 7: THE MATRIX MULTIPLICATION IN THIS PROGRAM CAN BE SPED UP BY   *
!USING ATLAS/INTEL MKL LIBRARY SUBROUTINE DGEMM. HOWEVER WE HAVE CODED *
!A CRUDE MATRIX MULTIPLICATION METHOD JUST FOR PERFORMANCE COMPARISON  *
!                                                                      *
!NOTE 8: ALL ARRAYS ARE INITIALIZED AS SINGLE DIMENSION ARRAYS AS THEY *
!ARE LAID OUT IN MEMORY AND ARE USED TO PERFORM ALL THESE CALCULATIONS.*
!THE RESULTS ARE HOWEVER PRINTED AS TWO DIMENSION ARRAYS. THIS SAVES   *
!MEMORY AS POINTERS TO POINTERS ARE NOT REQUIRED FOR INITILIZING ARRAYS*
!MOREOVER THIS ALSO ALLOWS TO BETTER MAINTAIN CACHE COHERENCY.         *
!                                                                      *
!NOTE 9: DUE TO ROW WISE STORAGE IN C++ WE ACTUALLY PERFORM KRON(D,E)  *
!*H'AND TRANSFORM THE RESULTS FOR THE INDIRECT METHOD TO BETTER	       *
!OPTIMIZE FOR CACHE COHERENCE                                          *
!                                                                      *
!***********CODE BEGINS:TYPSET AND DEFINITIONS**************************/

#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string>
#include <iostream>

using namespace std;

/* A FUNCTION TO PRINT DOUBLE ARRAYS*/
void print_array_1d(double * __restrict Array, int rows, int cols)
	/* Inputs: double *Array, a pointer to a double array
	int rows: Number of rows in the array
	int cols: Number of cols in the array*/
{
	for (int i = 0; i <rows; ++i){
		for (int j = 0; j<cols; ++j){
			cout <<fixed<<setw(2)<<setprecision(2)<<Array[i*cols+j]<< " " ;
		}
		cout << "\n";
	}
}

/* A FUNCTION TO INITIALIZE DOUBLE ARRAYS WITH CONSTANTS*/
void inital_1cons_double(double * __restrict Array, int elements ,double initial)
	/* Inputs: double *Array, a pointer to a double array
	int elements: Number of elements in a vector or a 2d array stored as a vector
	double initial: The number by which to initialize the array*/
{
	for (int i = 0; i < elements; ++i){
		Array[i] = initial;
	}
}

/* A FUNCTION TO PERFORM OUT OF PLACE TRANSPOSE OF ARRAYS
THIS IS NOT AN OPTIMIZED FUNCTION BUT DOES ITS JOB*/
void transpose(double * __restrict Array, double * __restrict trans, int rows, int cols)
	/* Inputs: double *Array, a pointer to a double array that would be transposed
	double *Array: A pointer to an array that will store the result after the transpose operation
	int rows: Total rows in the array that would be transposed
	int cols: Total cols in the array that would be transposed*/
{
#pragma omp parallel for shared (Array,trans)
	for (int i = 0; i <rows; ++i){
		for (int j = 0; j <cols; ++j)
			trans[i+j*rows] = Array[i*cols+j] ;

	}
}

/* A FUNCTION TO INITIALIZE ARRAYS WITH RANDOM NUMBERS*/
void inital_1rand(double * __restrict Array, int elements, double low, double high)
	/* Inputs: double *Array, a pointer to a double array
	int elements: Number of elements in a vector or a 2d array stored as a vector
	double low: The lower limit for generating random numbers
	double high: The upper limit for generating random numbers*/
{
	double range = (high-low);
	for (int i = 0; i < elements; ++i){
		Array[i] = rand() * range / (RAND_MAX + low) ;
	}
}

/* A FUNCTION TO COMPUTE KRONECKER PRODUCT OF TWO ARRAYS*/
void kron_prod_single(double * __restrict D, double * __restrict E, double * __restrict Result, int drows, int dcols, int erows, int ecols)
	/* EX: A*kron(D,E)=Result
	NOTE: In this function 2d arrays are stored as 1D array in the memory
	Inputs: double *D, a pointer to a double array (D)
	Inputs: double *E, a pointer to a double array (E)
	Inputs: double *Result, a pointer to a double array (Result)
	int drows: total rows in D matrix
	int dcols: total cols in D matrix
	int erows: total rows in E matrix
	int ecols: total cols in E matrix*/

{
	/* The logic is complicated in this function as we are dealing with two dimension
	arrays stored in row major order and we treat them simply as vectors and also tries to
	maintain cache coherence*/
	int m ;
	int n = 0;
	int x = 0 ;
	for (int i = 0; i <= drows-1; i++){
		for (int j = (i*erows); j <= (i*erows)+(erows-1); j++){
			m=i*dcols;
			for (int k = (j*dcols); k <= (j*dcols)+(dcols-1); k++){
				n=(j*ecols)-(i*x);
				for (int l = (k*ecols); l <= (k*ecols)+(ecols-1); l++){
					Result[l] = D[m]*E[n];
					n = n+1;
				}
				m = m+1;
			}
		}
		x=n;
	}
}

/* A PARALLEL FUNCTION TO COMPUTE MATRIX MULTIPLICATION OF TWO ARRAYS*/
void mat_mul(double* __restrict A, double* __restrict B, double* __restrict Result, int arows, int acols, int bcols)
	/* EX: A*B=Result
	NOTE: In this function 2d arrays are stored as 1D array in the memory
	Inputs: double *A, a pointer to a double array (A)
	Inputs: double *B, a pointer to a double array (B)
	Inputs: double *Result, a pointer to a double array (Result)
	int arows: total rows in A matrix
	int acols: total cols in A matrix
	int bcols: total cols in B matrix */
{
	int idx = 0;
	double sum = 0.0;
    #pragma omp parallel for shared (A,B,Result)
	/* this loop computes all rows all columns*/
	for (int i = 0; i <arows; ++i){
		idx = (i*acols);
		/* this loop computes one row all columns*/
		for (int j = 0; j<bcols; ++j){
			sum = 0.0;
			/* this loop computes one row one column*/
			for (int k = 0; k < acols; ++k){
				sum += A[idx+k]*B[(k*bcols)+j];
			}
			Result[i*bcols+j]=sum;
		}
	}
}

/* A FUNCTION TO COMPUTE HQ INDIRECTLY*/
void hq (double *h, double *d, double *e, double *Result, int n, int p, int r)
	/* EX: h*kron(d,e)=Result. This function performs matrix multiplication of h and q
	indirectly
	NOTE: In this function 2d arrays are stored as 1D array in the memory
	Inputs: double *h, a pointer to a double array (h)
	Inputs: double *d, a pointer to a double array (d) which is temporal covariance matrix
	Inputs: double *e, a pointer to a double array (e) which is spatial covariance matrix
	Inputs: double *Result, a pointer to a double array (Result) which would have
	oputput from matrix multiplication of h*kron(d,e)
	int n: total observation or first dimension of h matrix
	int p: total rows=cols of the square symmetric matrix D
	int r: total rows=cols of the square symmetric matrix E */
{
	/* initilize an array that would contain results*/
	double *he;
	double *hd;
	int m = 0;
	int t = -1;
	int cons = n*r;
	/* allocate memory*/
	hd=new double[n*r];// Each column block of h multiplied by a d value
	he=new double[n*r];
	inital_1cons_double(hd, n*r ,0.0);

	/* this loop repeats this process for all columns in d matrix
	Since memory is dynamically allocated inside this functions a
	try and catch block is included to avoid memory leaks */

	try{for (int i = 0; i<=(p-1); i++){
		/* this loop iterates over all the blocks*/
		for (int j = 0; j<=(p-1); j++) {
			m=0;
			t=t+1;
			/* if an h[k] within a column block it multiplies it with
			the a constant d value*/
			for (int k = (j*cons); k <=(j*cons)+(cons-1); k++){
				/* if condition to avoid multiplication with zeros
				and only do addition if the d matrix has ones*/
				if ( d[t] != 0 && d[t] !=1 ){
					hd[m]+= h[k]*d[t];
				}
				else if (d[t] == 1) {
					hd[m]+=h[k];
				}
				m=m+1;
			}
		}
		/* matrix multiplication of e matrix by hd*/
		mat_mul(e,hd,he,r,r,n);
		m=0;
		/*place results at correct locations*/
		for (int l = (i*cons); l <=(i*cons)+(cons-1); l++){
			Result[l]=he[m];
			m=m+1;;
		}
		/* again initilize hd with zeros to for next round*/
		inital_1cons_double(hd, n*r ,0.0);
	}
	/* recover memory and set pointers to null in try block*/
	delete [] hd; hd = NULL;
	delete [] he; he = NULL;
	}
	/* recover memory and set pointers to null in catch block*/
	catch(...){
		delete [] hd; hd = NULL;
		delete [] he; he = NULL;
		throw;
	}

}
/* MAIN FUNCTION*/
int main(void)
{
	/* Main function to compute hq for direct and indirect method*/

	/*pointers to arrays initialized as null*/
	double *D = NULL ; // temporal covariance matrix
	double *Temp = NULL; // temporary array 1 for storing transpose of a matrix and matrix multiplication
	double *Temp2 = NULL ; // temporary array 2 for storing transpose of a matrix and matrix multiplication
	double *E = NULL ; // spatial covariance matrix
	double *H = NULL ; // H matrix or forward operator in inverse problems
	double *Q = NULL ; // full covariance matrix that is kron(D,E)
	double *HQ = NULL ; // HQ matrix from the direct method
	double *HQI = NULL; // HQI HQ matrix from the indirect method

	/*initialize integers that will hold dimensions of arrays*/
	int n = 0;
	int p = 0;
	int r = 0;

	string PrintCon ;

	/*take user input for dimensions*/
	cout<<"Enter First Dimension of H Matrix"<<'\n';
	cin>>n;

	cout<<"Enter First Dimension of Square Symmetric D or Temporal Covariance Matrix "<<'\n';
	cin>>p;

	cout<<"Enter First Dimension of Square Symmetric E or Spatial Covariance Matrix "<<'\n';
	cin>>r;

	cout<<"Warning: Option of Printing Matrices is to Check Results for Very Small Matrices"<<'\n';
	cout<<"Do Not Use This Option to Print Large Matrices"<<'\n';
	cout<<"Do You Want to Print Matrix D, E, H , HQ_Direct, HQ_Indirect: small case yes or no "<<'\n';
	cin.ignore();
	getline(cin,PrintCon,'\n');;

	cout<<"*******************************************************"<<'\n';
	cout<<"The Parameters You Entered"<<'\n';
	cout<<"Size of H:"<<"["<<n<<","<<p*r<<"]"<<'\n';
	cout<<"Temporal Covariance Size:"<<"["<<p<<","<<p<<"]"<<'\n';
	cout<<"Spatial Covariance Size:"<<"["<<r<<","<<r<<"]"<<'\n';
	cout<<"*******************************************************"<<'\n';


	/*initialize arrays with random numbers*/
	srand ( time(NULL) ); // Initilize a random seed

	/*Compute and make a symmetric temporal covariance matrix D */
	Temp = new double[p*p];
	inital_1rand(Temp, p*p, 1.0, 2.0);
	Temp2 = new double[p*p];
	transpose(Temp,Temp2,p,p);
	D = new double[p*p];
	mat_mul(Temp,Temp2,D,p,p,p);
	delete [] Temp; Temp = NULL;
	delete [] Temp2; Temp2 = NULL;

	/*Compute and make a symmetric spatial covariance matrix E */
	Temp = new double[r*r];
	inital_1rand(Temp, r*r, 0.0, 1.0);
	Temp2 = new double[r*r];
	transpose(Temp,Temp2,r,r);
	E = new double[r*r];
	mat_mul(Temp,Temp2,E,r,r,r);
	delete [] Temp; Temp = NULL;
	delete [] Temp2; Temp2 = NULL;

	/* Compute Q and HQ from the direct method*/
	H = new double[n*p*r];
	inital_1rand(H, n*p*r, 0.0, 1.0);
	Q = new double[p*r*p*r];
	kron_prod_single(D,E,Q,p,p,r,r);
	HQ = new double[n*p*r];
	clock_t tStart = clock();
	mat_mul(H,Q,HQ,n,p*r,p*r);
	printf("Time taken to do HQ Directly: %.8fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	/*Compute HQ from the indirect method*/
	Temp = new double[n*p*r];
	transpose(H,Temp,n,p*r);
	HQI = new double[n*p*r];
	clock_t mStart = clock();
	hq(Temp,D,E,HQI,n,p,r);
	printf("Time taken to do HQ Indirectly: %.8fs\n", (double)(clock() - mStart)/CLOCKS_PER_SEC);
	delete [] Temp; Temp = NULL;


	/*print matrices*/
	if (PrintCon == "yes"){
		cout<<'\n'<<"Now Printing All the Matrices"<<'\n';
		cout<<'\n'<<"Printing H:" <<'\n';
		print_array_1d(H, n, p*r);
		cout<<'\n'<<"Printing D:" <<'\n';
		print_array_1d(D, p, p);
		cout<<'\n'<<"Printing E:" <<'\n';
		print_array_1d(E, r, r);
		cout<<'\n'<<"Printing Q:" <<'\n';
		print_array_1d(Q, p*r, p*r);
		cout<<'\n'<<"Printing HQ From Direct Method:" <<'\n';
		print_array_1d(HQ,n,p*r);//
		cout<<'\n'<<"Printing HQ From Indirect Method:" <<'\n';
		Temp = new double[n*p*r];
		transpose(HQI,Temp,p*r,n);
		print_array_1d(Temp,n,p*r);
		delete [] Temp; Temp=NULL;
	}

	/*avoid memory leaks and dangling pointers*/
	delete [] D; D=NULL;
	delete [] E; E=NULL;
	delete [] H; H=NULL;
	delete [] Q; Q=NULL;
	delete [] HQ; HQ=NULL;
	delete [] HQI; HQI=NULL;

	/*print total time taken by the program*/
	printf("Total Time Taken By The Complete Program: %.8fs\n", (double)(clock() - mStart)/CLOCKS_PER_SEC);
	return 0;

}



