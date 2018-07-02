//
//  wsbt.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "wsbt.hpp"


wsbt::wsbt()
{
    setWeights();
    setScores();
    testStatistic();
}


void wsbt::setWeights()
{
    for(int i = 0; i < weights.size(); ++i)
    {
        //weights[i] = gsl_ran_beta_pdf(maf[i],1,25);
    }
}

void wsbt::setScores()
{
    for(int j = 0; j < scores.size(); ++j)
    {
        for(int i = 0; i < genotypeMatrix[j].size(); ++i)
        {
            scores[j] += genotypeMatrix[j][i] * (phenotype[i] - expectedPhenotype);
        }
    }
}

void wsbt::testStatistic()
{
    double tempStat = 0;
    for(int i = 0; i < scores.size(); i++)
    {
        tempStat += scores[i] * weights[i];
    }
    testStat = tempStat * tempStat;
    //pvalue = gsl_cdf_chisq_P(1.0, testStat);
    
}
