
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cpuid.h>
#include <iomanip>
#include <cstdarg>
#include <stdio.h>
#include "Macros.h"
#include "Arsenal.h"

using namespace std;

void printline()
{
  cout << "\n--------------------------------------------------\n" << endl;
}

void print_processor()
{
#ifdef PRINT_PROCESSOR
  char CPUBrandString[0x40];              // copied this routine online
  unsigned int CPUInfo[4] = {0,0,0,0};

  __cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
  unsigned int nExIds = CPUInfo[0];

  memset(CPUBrandString, 0, sizeof(CPUBrandString));

  for (unsigned int i = 0x80000000; i <= nExIds; ++i)
  {
      __cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

      if (i == 0x80000002)
          memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
      else if (i == 0x80000003)
          memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
      else if (i == 0x80000004)
          memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
  }

  printf("\n\nProcessor = %s\n\n", CPUBrandString);
#endif
}


//**********************************************************************
vector<double> stringToDoubles(string str)
// Return a vector of doubles from the string "str". "str" should
// be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  vector<double> valueList;
  double val;
  sst >> val;
  while (sst.eof()==false)
  {
    valueList.push_back(val);
    sst >> val;
  }
  return valueList;
}

//**********************************************************************
double stringToDouble(string str)
// Return the 1st doubles number read from the string "str". "str" should be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  double val;
  sst >> val;
  return val;
}


//**********************************************************************
vector< vector<double>* >* readBlockData(istream &stream_in)
// Return a nested vector of vector<double>* object. Each column of data
// is stored in a vector<double> array and the collection is the returned
// object. Data are read from the input stream "stream_in". Each line
// of data is processed by the stringToDoubles function. Note that the
// data block is dynamicall allocated and is not release within the
// function.
// Note that all "vectors" are "new" so don't forget to delete them.
// Warning that also check if the last line is read correctly. Some files
// are not endded properly and the last line is not read.
{
  vector< vector<double>* >* data;
  vector<double> valuesInEachLine;
  long lineSize;
  long i; // temp variable
  char buffer[99999]; // each line should be shorter than this

  // first line:
  stream_in.getline(buffer,99999);
  valuesInEachLine = stringToDoubles(buffer);
  // see if it is empty:
  lineSize = valuesInEachLine.size();
  if (lineSize==0)
  {
    // empty:
    cout << "readBlockData warning: input stream has empty first row; no data read" << endl;
    return NULL;
  }
  else
  {
    // not empty; allocate memory:
    data = new vector< vector<double>* >(lineSize);
    for (i=0; i<lineSize; i++) (*data)[i] = new vector<double>;
  }

  // rest of the lines:
  while (stream_in.eof()==false)
  {
    // set values:
    for (i=0; i<lineSize; i++) (*(*data)[i]).push_back(valuesInEachLine[i]);
    // next line:
    stream_in.getline(buffer,99999);
    valuesInEachLine = stringToDoubles(buffer);
  }

  return data;
}


//**********************************************************************
void releaseBlockData(vector< vector<double>* >* data)
// Use to delete the data block allocated by readBlockData function.
{
  if (data)
  {
    for (unsigned long i=0; i<data->size(); i++) delete (*data)[i];
    delete data;
  }
}


//**********************************************************************
string toLower(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  for (string::iterator it=tmp.begin(); it<=tmp.end(); it++) *it = tolower(*it);
  return tmp;
}

//**********************************************************************
string trim(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  long number_of_char = 0;
  for (size_t ii=0; ii<str.size(); ii++)
    if (str[ii]!=' ' && str[ii]!='\t')
    {
      tmp[number_of_char]=str[ii];
      number_of_char++;
    }
  tmp.resize(number_of_char);
  return tmp;
}



void matrix_multiplication(double ** A, const double x[], double y[], int n, int m)
{
  // multiplication y_i = A_ij . x_j (stores result in y)

  for(int i = 0; i < n; i++)    // rows
  {
    y[i] = 0;                   // initialize to zero

    for(int j = 0; j < m; j++)  // columns
    {
      y[i] += A[i][j] * x[j];
    }
  }
}

void vector_copy(const double a[], double c[], int n)
{
  // copy a to c
  for(int i = 0; i < n; i++)
  {
    c[i] = a[i];
  }
}

void vector_addition(const double a[], const double b[], double c[], int n)
{
  // addition c_i = a_i + b_i (stores result in c)
  for(int i = 0; i < n; i++)
  {
    c[i] = a[i] + b[i];
  }
}

void vector_subtraction(const double a[], const double b[], double c[], int n)
{
  // subtraction c_i = a_i - b_i (stores result in c)
  for(int i = 0; i < n; i++)
  {
    c[i] = a[i] - b[i];
  }
}


void free_2D(double ** M, int n)
{
  for(int i = 0; i < n; i++)
  {
    free(M[i]);
  }
  free(M);
}


void free_3D(double *** M, int n, int m)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m; j++)
    {
      free(M[i][j]);
    }
    free(M[i]);
  }
  free(M);
}



