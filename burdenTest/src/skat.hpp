//
//  skat.hpp
//  burdenTest
//
//  Created by Corin Thummel on 7/30/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#ifndef skat_hpp
#define skat_hpp

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

using namespace std;

class skat
{
    
public:
    skat(gsl_matrix* geno, gsl_vector* maf, gsl_matrix* covariates, gsl_vector* pheno);
    
    
    
private:
    void logisticRegression();
    void linearRegression();
    
    int X_Count;
    int subjectCount;
    int variantCount;
    double testStatistic;
    double regression;
    bool isBinary;
    
    void setWeights(gsl_vector *maf);
    void makeKernel(string kernel_type);
    void setTestStatistic();
    gsl_matrix* matrixInverse(gsl_matrix* m);

    
    gsl_matrix *genoMatrix;
    gsl_matrix *weightMatrix;
    gsl_matrix *kernel;
    gsl_matrix *X;
    gsl_vector *pheno;
    gsl_vector *coeff;
    gsl_vector *res;
};





#endif /* skat_hpp */
