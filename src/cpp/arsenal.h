
#ifndef arsenal_h
#define arsenal_h

#include "stdlib.h"
#include <vector>
#include <string>

using namespace std;

void printline();
void print_processor();

// used in ParameterReader
double stringToDouble(string);
vector<double> stringToDoubles(string);
string toLower(string str);
string trim(string str);

vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);


// matrix operations for modified momentum transformation
void matrix_multiplication(double ** A, const double x[], double y[], int n, int m);	// compute A.x, store in y
void vector_copy(const double a[], double c[], int n);									// store a in c
void vector_addition(const double a[], const double b[], double c[], int n);			// compute a + b, store in c
void vector_subtraction(const double a[], const double b[], double c[], int n);			// compute a - b, store in c

// free arrays
void free_2D(double ** M, int n);
void free_3D(double *** M, int n, int m);

#endif
