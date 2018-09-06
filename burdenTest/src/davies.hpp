#ifndef davies_hpp
#define davies_hpp

#include <stdio.h>

//Carry out integration with nterm terms, at stepsize interv. If not main then multiply integrand by 1.0 - exp(-.5 * tausq * u^2)


    double qf(double lb[], double nc[], int n[], int r, double sigma, double c, int lim, double acc, double trace[], int ifault);
    double ln1(double x, bool first);
    void order();
    double errbd(double u, double cx);
    double ctff(double accx, double upn);
    double truncation(double u, double tausq);
    void findu(double utx, double accx);
    void integrate(int nterm, double interv, double tausq, bool main);
    double cfe(double x);




#endif /* davies_hpp */
