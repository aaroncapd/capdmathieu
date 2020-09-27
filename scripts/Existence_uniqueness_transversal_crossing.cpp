/////////////////////////////////////////////////////////////////////////////
//
/// @file Existence_uniqueness_transversal_crossing.cpp
/// @author Aaron Fernandez
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
#define LEN 256

IVector functioneval_ab(interval a, interval b, interval eps)
{
    // This is vector field for the deriatives respect "a" and "b" of the Mathieu system
    IMap mathieuab("time:t;par:a,b,eps;var:x1,x2,x3,x4,y1,y2,y3,y4,y5,y6,y7,y8;fun:x2,-x1*(a+b*(cos(t)+eps*cos(2*t))),x4,-x3*(a+b*(cos(t)+eps*cos(2*t))),y3,y4,-x1-y1*(a+b*(cos(t)+eps*cos(2*t))),-x3-y2*(a+b*(cos(t)+eps*cos(2*t))),y7,y8,-x1*(cos(t)+eps*cos(2*t))-y5*(a+b*(cos(t)+eps*cos(2*t))),-x3*(cos(t)+eps*cos(2*t))-y6*(a+b*(cos(t)+eps*cos(2*t)));");
    mathieuab.setParameter("a",a);
    mathieuab.setParameter("b",b);
    mathieuab.setParameter("eps",eps);

    // Integrator
    int order = 25;
    IOdeSolver solver(mathieuab,order);
    ITimeMap tm(solver);

    // Initial conditions
    IVector y0(12);
    y0[0] = interval(1.);
    y0[1] = interval(0.);
    y0[2] = interval(0.);
    y0[3] = interval(1.);
    y0[4] = interval(0.);
    y0[5] = interval(0.);
    y0[6] = interval(0.);
    y0[7] = interval(0.);
    y0[8] = interval(0.);
    y0[9] = interval(0.);
    y0[10] = interval(0.);
    y0[11] = interval(0.);

    C0HOTripletonSet u( y0 );
    IVector y = tm(2.*interval::pi(),u);

    IVector result(6);
    result[0] = y[2];
    result[1] = y[1];
    result[2] = y[5];
    result[3] = y[9];
    result[4] = y[6];
    result[5] = y[10];

    return result;
}

// Krawczyk method
void krawczyk(interval a, interval b, interval eps, double tol)
{
    double adiam = diam(a).sup();
    double bdiam = diam(b).sup();
    cout << "-------------------------------------" << endl;
    cout << "2-cube a x b: " << a << "x" << b << endl;
    cout << "Widths (a,b): (" << adiam << "," << bdiam << ")" << endl;
    IVector mideval = functioneval_ab(mid(a),mid(b),eps);
    IVector deriv = functioneval_ab(a,b,eps);

    if (adiam < tol && bdiam < tol) {
        cout << "Krawczyk method has reached convergence." << endl;
    } else {
        IVector X(2);
        X[0] = a;
        X[1] = b;

        IVector midp(2);
        midp[0] = mid(a);
        midp[1] = mid(b);

        IVector midfp(2);
        midfp[0] = mid(mideval[0]);
        midfp[1] = mid(mideval[1]);

        IMatrix C(2,2);
        C[0][0] = mid(deriv[2]);
        C[0][1] = mid(deriv[3]);
        C[1][0] = mid(deriv[4]);
        C[1][1] = mid(deriv[5]);

        IMatrix dF_X(2,2);
        dF_X[0][0] = deriv[2];
        dF_X[0][1] = deriv[3];
        dF_X[1][0] = deriv[4];
        dF_X[1][1] = deriv[5];

        IMatrix I = IMatrix::Identity(2);

        IVector DP1(2);
        IMatrix DP2(2,2);
        capd::matrixAlgorithms::gauss(C,midfp,DP1);
        capd::matrixAlgorithms::gauss(C,dF_X,DP2);

        IVector newab;
        IVector kraw = midp - DP1 +(I - DP2)*(X-midp);
        cout << "2-cube K(a,b): " << kraw[0] << "x" << kraw[1] << endl;

        if (intersection(a,kraw[0],newab[0]) && intersection(b,kraw[1],newab[1])) {
            if (!kraw[0].contains(a) || !kraw[1].contains(b)) {
                krawczyk(newab[0],newab[1],eps,tol);
            } else { cout << "Krawczyk ends. Same new interval." << endl;}
        } else { cout << "Krawczyk ends. No intersection." << endl;}
    }
}

int main()
{
  cout.precision(16);
  try{
    // Initial parameters
    double tol = 1e-14;
    interval a = interval(1024.,1028.)/1000.;
    interval b = interval(391.,408.)/1000.;
    interval eps = interval(0.2);

    // Krawczyk method starts
    krawczyk(a,b,eps,tol);

  }catch(exception& e)
  {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
