//
//  genericBurdenTest.hpp
//  burdenTest
//
//  Created by Corin Thummel on 7/13/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#ifndef genericBurdenTest_hpp
#define genericBurdenTest_hpp

#include <stdio.h>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


class genericBurdenTest
{
public:
    genericBurdenTest(gsl_matrix* G, gsl_vector* maf, gsl_vector* ptype);
    
    double getPvalue(){return pvalue;}
    gsl_vector* getScores(){return scores;}
    gsl_vector* getWeights(){return weights;}
    
private:
    void setWeights(gsl_vector* maf);
    void setScores(gsl_matrix* G, gsl_vector* ptype);
    void testStatistic();
    
    double pvalue;
    double expectedPhenotype;
    double testStat;
    
    gsl_vector* scores;
    gsl_vector* weights;
};

#endif /* genericBurdenTest_hpp */
