//
//  skat.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/30/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "skat.hpp"
#include "davies.cpp"

using namespace std;

skat::skat(gsl_matrix* geno, gsl_vector* maf, gsl_vector* covariates)
{
    //Initialize variables.
    genoMatrix = geno;
    X = covariates;
    weightMatrix = gsl_matrix_calloc(maf->size, maf->size);
    kernel = gsl_matrix_alloc(geno->size2, geno->size2);
    
    string kernel_type = "linear";
    
    //Run Test.
    setWeights(maf);
    makeKernel(kernel_type);
    setTestStatistic();
    
    
    //Cleanup after test.
    gsl_matrix_free(weightMatrix);
    gsl_matrix_free(kernel);
}

//Creates the mxm matrix of weights.
void skat::setWeights(gsl_vector *maf)
{
    for(int i = 0; i < maf->size; i++)
    {
        double tempWeight = gsl_ran_beta_pdf(gsl_vector_get(maf, i), 1, 25);
        gsl_matrix_set(weightMatrix, i, i, tempWeight * tempWeight);
    }
}

//Remember matrix multiplication goes right to left. So G'WG needs to calculate G'(WG) in a sense. (v = rows, s = columns)
//sxv vxv vxs -> sxv vxs -> sxs is kernal's final dimensions.
void skat::makeKernel(string kernel_type)
{
    if(kernel_type == "linear")
    {
        gsl_matrix *tempkernel = gsl_matrix_alloc(genoMatrix->size1, genoMatrix->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, weightMatrix, genoMatrix, 0, tempkernel);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, genoMatrix, tempkernel, 0.0, kernel);
    }
    else if(kernel_type == "quad")
    {
        
    }
    else if (kernel_type == "IBS")
    {
        
    }
    
}

void skat::setTestStatistic()
{
    //Q = (y - u) * K * (y - u)'
    // 1xs sxs sx1 -> 1xs sx1 -> 1x1 = Q
    gsl_vector *pheno = gsl_vector_alloc(genoMatrix->size2);
    gsl_vector *tempstat = gsl_vector_alloc(genoMatrix->size2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, weightMatrix, pheno, 0.0, tempstat);
    gsl_blas_ddot(pheno, tempstat, &testStatistic);
    
}












