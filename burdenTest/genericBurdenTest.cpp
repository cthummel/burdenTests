//
//  This burden test only works for continuous phenotype data.
//
//
//  Created by Corin Thummel on 7/13/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include <iostream>
#include "genericBurdenTest.hpp"


genericBurdenTest::genericBurdenTest(std::vector<std::vector<int> > G, std::vector<double> inputMaf, std::vector<double> ptype)
{
    pvalue = 0;
    expectedPhenotype = 0;
    weights = std::vector<double>(inputMaf.size());
    scores = std::vector<double>(inputMaf.size());
    
    for(int i = 0; i < ptype.size(); ++i)
    {
        expectedPhenotype += ptype[i];
    }
    expectedPhenotype = expectedPhenotype / ptype.size();
    
    setWeights(inputMaf);
    setScores(G, ptype);
    testStatistic();
}

//Set weights for variants based off Beta(maf; 1, 25). Can make it customizable later.
void genericBurdenTest::setWeights(std::vector<double> maf)
{
    for(int i = 0; i < weights.size(); ++i)
    {
        //std::cout << "AF= " << maf[i] << std::endl;
        weights[i] = gsl_ran_beta_pdf(maf[i],1,25);
    }
}


void genericBurdenTest::setScores(std::vector<std::vector<int> > genotypeMatrix, std::vector<double> phenotype)
{
    for(int i = 0; i < genotypeMatrix.size(); ++i)
    {
        for(int j = 0; j < genotypeMatrix[i].size(); ++j)
        {
            scores[j] += genotypeMatrix[i][j] * (phenotype[i] - expectedPhenotype);
        }
    }
}

void genericBurdenTest::testStatistic()
{
    double tempStat = 0;
    for(int i = 0; i < scores.size(); i++)
    {
        tempStat += scores[i] * weights[i];
    }
    testStat = tempStat * tempStat;
    pvalue = gsl_cdf_chisq_Q(testStat, 1);
}
