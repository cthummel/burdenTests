//
//  skat.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/30/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "skat.hpp"


skat::skat(gsl_matrix* geno, gsl_vector* maf, gsl_vector* covariates)
{
    weightMatrix = gsl_matrix_calloc(maf->size, maf->size);
    
    
}

//Creates the nxn matrix of weights.
void skat::setWeights(gsl_vector *maf)
{
    for(int i = 0; i < maf->size; i++)
    {
        gsl_matrix_set(weightMatrix, i, i, gsl_ran_beta_pdf(gsl_vector_get(maf, i), 1, 25));
    }
}

void skat::makeKernel(gsl_matrix *geno)
{
    gsl_matrix *tempkernel = gsl_matrix_alloc(geno->size1, geno->size1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, geno, weightMatrix, 0, tempkernel);
}

