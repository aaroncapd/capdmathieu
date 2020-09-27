/////////////////////////////////////////////////////////////////////////////
//
/// @file Uniqueness_boundaries.cpp
/// @author Aaron Fernandez
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
#define LEN 256

IVector functioneval_a(interval a, interval b, interval eps)
{
    // This is vector field for the deriatives respect "a" of the Mathieu system
    IMap mathieua("time:t;par:a,b,eps;var:x1,x2,x3,x4,y1,y2,y3,y4;fun:x2,-x1*(a+b*(cos(t)+eps*cos(2*t))),x4,-x3*(a+b*(cos(t)+eps*cos(2*t))),y3,y4,-x1-y1*(a+b*(cos(t)+eps*cos(2*t))),-x3-y2*(a+b*(cos(t)+eps*cos(2*t)));");
    mathieua.setParameter("a",a);
    mathieua.setParameter("b",b);
    mathieua.setParameter("eps",eps);

    // Integrator
    int order = 25;
    IOdeSolver solver(mathieua,order);
    ITimeMap tm(solver);

    // Initial conditions
    IVector y0(8);
    y0[0] = interval(1.);
    y0[1] = interval(0.);
    y0[2] = interval(0.);
    y0[3] = interval(1.);
    y0[4] = interval(0.);
    y0[5] = interval(0.);
    y0[6] = interval(0.);
    y0[7] = interval(0.);

    C0HOTripletonSet u( y0 );
    IVector y = tm(2.*interval::pi(),u);

    IVector result(4);
    result[0] = y[1];
    result[1] = y[2];
    result[2] = y[6];
    result[3] = y[5];

    return result;
}

// Newton interval method
void newton(interval a, interval b, interval eps, double tol)
{
    double adiam = diam(a).sup();
    cout << "----------------------------------------------------" << endl;
    cout << "Interval a = " << a << "\nDiameter of a = " << adiam << endl;
    IVector mideval = functioneval_a(mid(a),b,eps);
    IVector deriv = functioneval_a(a,b,eps);

    if (adiam < tol) {
        cout << "Newton's method has reached convergence." << endl;
    } else {

    // If the enclosure of the derivative contains a 0
    if (deriv[2].contains(0.)) {
        cout << "Derivative enclosure: " << deriv[2] << endl;
        cout << "Derivative contains 0. Dividing into 2 intervals." << endl;

        // First interval
        interval newa;
        interval newt = mid(a)-mideval[0]/(interval(deriv[2].inf(),-tol));
        cout << "Interval N(a) (right) = " << newt << endl;
        if (intersection(a,newt,newa)) {
            if (!newt.contains(a)) {
                newton(newa,b,eps,tol);
            } else { cout << "Newton ends. Same new interval." << endl;}
        } else { cout << "Newton ends. No intersection." << endl;}

        // Second interval
        newt = mid(a)-mideval[0]/(interval(tol,deriv[2].sup()));
        cout << "Interval N(a) (left) = " << newt << endl;
        if (intersection(a,newt,newa)) {
            if (!newt.contains(a)) {
                newton(newa,b,eps,tol);
            } else { cout << "Newton ends. Same new interval." << endl;}
        } else { cout << "Newton ends. No intersection." << endl;}

    // Normal case withour the derivative containing 0
    } else {
        interval newa;
        interval newt = mid(a)-mideval[0]/(deriv[2]);
        cout << "Interval N(a) = " << newt << endl;
        if (intersection(a,newt,newa)) {
            if (!newt.contains(a)) {
                newton(newa,b,eps,tol);
            } else { cout << "Newton ends. Same new interval." << endl;}
        } else { cout << "Newton ends. No intersection." << endl;}
    }
    }
}

int main()
{
  cout.precision(15);
  try{
    // Initial parameters
    double tol = 1e-12;
    interval a = interval(99.,100.)/100;
    interval b = interval(0.2);
    interval eps = interval(0.2);

    // Newton method starts
    newton(a,b,eps,tol);

  }catch(exception& e)
  {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
