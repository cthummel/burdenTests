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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

class skat
{
    
public:
    skat(gsl_matrix* geno, gsl_vector* maf, gsl_vector* covariates);
    
    
    
private:
    void setWeights(gsl_vector *maf);
    void makeKernel(gsl_matrix *geno);
    
    
    gsl_matrix *weightMatrix;
    gsl_matrix *kernel;
    
};





#endif /* skat_hpp */
