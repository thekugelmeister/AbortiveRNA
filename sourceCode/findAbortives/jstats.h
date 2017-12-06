/* jstats.h
   statistics functions made available for use in searching for abortive binding
 */

#ifndef JSTATS
#define JSTATS

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

double jMean(int a[], int len);
double jSSD(int a[], int len);
double jSTT(int a[], int len, int val);

#endif
