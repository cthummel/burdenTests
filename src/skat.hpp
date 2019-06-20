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
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

using namespace std;

class skat
{
    
public:
    skat(gsl_matrix* geno, gsl_vector* maf, gsl_matrix* covariates, gsl_vector* pheno);
    double getPvalue(){return pvalue;}
    double getQ(){return testStatistic;}
    double getRvQ(){return rvStat;}
    
private:
    int X_Count;
    int subjectCount;
    int variantCount;
    int caseCount;
    double rvStat;
    double testStatistic;
    double regression;
    double pvalue;
    bool isBinary;
    bool rvtestVersion = true;
    
    void setWeights(gsl_vector *maf);
    void makeKernel(string kernel_type);
    void setTestStatistic();
    void qDistribution();
    void setPvalue();
    gsl_vector* logisticRegression();
    double linearRegression();
    gsl_matrix* matrixInverse(gsl_matrix* m);
    int sqrtMatrix(gsl_matrix *input);
    double dotProductCheck(gsl_vector* left, gsl_vector *right);
    bool checkVariants(gsl_matrix *genotype, gsl_vector *phenotype);

    
    gsl_matrix *genoMatrix = nullptr;
    gsl_matrix *weightMatrix = nullptr;
    gsl_matrix *kernel = nullptr;
    gsl_matrix *X = nullptr;
    gsl_vector *pheno = nullptr;
    gsl_vector *coeff = nullptr;
    gsl_vector *eigenvalues = nullptr;
};





#endif /* skat_hpp */
