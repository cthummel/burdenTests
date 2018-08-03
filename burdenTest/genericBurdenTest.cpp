//
//  This burden test only works for continuous phenotype data.
//
//
//  Created by Corin Thummel on 7/13/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include <iostream>
#include "genericBurdenTest.hpp"


genericBurdenTest::genericBurdenTest(gsl_matrix* G, gsl_vector* inputMaf, gsl_vector* ptype)
{
    pvalue = 0;
    expectedPhenotype = 0;
    weights = gsl_vector_alloc(inputMaf->size);
    scores = gsl_vector_alloc(inputMaf->size);
    
    for(int i = 0; i < ptype->size; i++)
    {
        expectedPhenotype += gsl_vector_get(ptype, i);
    }
    expectedPhenotype = expectedPhenotype / ptype->size;
    
    setWeights(inputMaf);
    setScores(G, ptype);
    testStatistic();
}

//Set weights for variants based off Beta(maf; 1, 25). Can make it customizable later.
void genericBurdenTest::setWeights(gsl_vector* maf)
{
    for(int i = 0; i < weights->size; ++i)
    {
        gsl_vector_set(weights, i, gsl_ran_beta_pdf(gsl_vector_get(maf, i), 1, 25));
        //weights[i] = gsl_ran_beta_pdf(maf[i],1,25);
    }
}


void genericBurdenTest::setScores(gsl_matrix* genotypeMatrix, gsl_vector* phenotype)
{
    for(int i = 0; i < genotypeMatrix->size1; i++)
    {
        for(int j = 0; j < genotypeMatrix->size2; j++)
        {
            double tempscore = gsl_vector_get(scores, j) + (gsl_matrix_get(genotypeMatrix, i, j) * (gsl_vector_get(phenotype, i) - expectedPhenotype));
            gsl_vector_set(scores, j, tempscore);
        }
    }
}

void genericBurdenTest::testStatistic()
{
    double tempStat = 0;
    for(int i = 0; i < scores->size; i++)
    {
        tempStat += gsl_vector_get(scores, i) * gsl_vector_get(weights, i);
    }
    testStat = tempStat * tempStat;
    pvalue = gsl_cdf_chisq_Q(testStat, 1);
}
