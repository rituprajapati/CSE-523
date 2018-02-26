#include<math.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <iomanip>

using namespace std;
#define ldouble  long double

vector< ldouble > phi_norms;
const double pi = 3.1415926535897;

#include "tensor.cpp"
#include "helperFunction.cpp"
#include "quadrature.cpp"
#include "testFunction.cpp"
#include "twoscalecoeffs.cpp"

int autorefine = 0;	