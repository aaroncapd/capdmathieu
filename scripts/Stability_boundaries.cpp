/////////////////////////////////////////////////////////////////////////////
//
/// @file Stability_boundaries.cpp
/// @author Aaron Fernandez
//
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <stdio.h>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
#define LEN 256

IVector functioneval(interval a, interval b, interval eps)
{
    // This is vector field for the Mathieu system
    IMap mathieu("time:t;par:a,b,eps;var:x1,x2,x3,x4;fun:x2,-x1*(a+b*(cos(t)+eps*cos(2*t))),x4,-x3*(a+b*(cos(t)+eps*cos(2*t)));");
    mathieu.setParameter("a",a);
    mathieu.setParameter("b",b);
    mathieu.setParameter("eps",eps);

    // Integrator
    int order = 25;
    IOdeSolver solver(mathieu,order);
    ITimeMap tm(solver);

    // Initial conditions
    IVector y0(4);
    y0[0] = interval(1.);
    y0[1] = interval(0.);
    y0[2] = interval(0.);
    y0[3] = interval(1.);

    C0HOTripletonSet u( y0 );
    IVector y = tm(2.*interval::pi(),u);

    // Resulting values for the period map
    IVector result(3);
    result[0] = abs(y[0]+y[3])-2.;
    result[1] = y[1];
    result[2] = y[2];

    return result;
}

// Bisection method to discard 2-cubes
void bisection(interval a, interval b, interval eps, double tol, FILE* fp)
{
    double bdiam = diam(b).sup();
    double adiam = diam(a).sup();
    IVector val = functioneval(a,b,eps);

    if (val[1].contains(0.) || val[2].contains(0.)) {

    if (adiam < tol && bdiam < tol) {
        fprintf(fp, "%f %f %f %f\n",a.inf(),a.sup(),b.inf(),b.sup());
    } else if (bdiam < tol) {
        bisection(interval(a.inf(),mid(a).sup()),b,eps,tol,fp);
        bisection(interval(mid(a).inf(),a.sup()),b,eps,tol,fp);
    } else if (adiam < tol) {
        bisection(a,interval(b.inf(),mid(b).sup()),eps,tol,fp);
        bisection(a,interval(mid(b).inf(),b.sup()),eps,tol,fp);
    } else {
        bisection(interval(a.inf(),mid(a).sup()),interval(b.inf(),mid(b).sup()),eps,tol,fp);
        bisection(interval(mid(a).inf(),a.sup()),interval(b.inf(),mid(b).sup()),eps,tol,fp);
        bisection(interval(a.inf(),mid(a).sup()),interval(mid(b).inf(),b.sup()),eps,tol,fp);
        bisection(interval(mid(a).inf(),a.sup()),interval(mid(b).inf(),b.sup()),eps,tol,fp);
    }
    }
}

int main()
{
  cout.precision(16);
  try{
    // Initial parameters
    double tol = 1e-2;
    interval a = interval(-1.,5.);
    interval b = interval(0.,6.);
    interval eps = interval(1.);

    // Open the file for writing
    FILE* fp;
    fp = fopen("/home/aaron/Descargas/capd-capdDynSys-5.1.2/projectStarter/test.txt","w");

    // Bisection method starts
    bisection(a,b,eps,tol,fp);

    // Close the file
    fclose(fp);

  }catch(exception& e)
  {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
