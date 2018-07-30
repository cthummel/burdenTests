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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


class genericBurdenTest
{
public:
    genericBurdenTest(std::vector<std::vector<int> > G, gsl_vector* maf, std::vector<double> ptype);
    
    double getPvalue(){return pvalue;}
    gsl_vector* getScores(){return scores;}
    gsl_vector* getWeights(){return weights;}
    
private:
    void setWeights(gsl_vector* maf);
    void setScores(std::vector<std::vector<int> > G, std::vector<double> ptype);
    void testStatistic();
    
    double pvalue;
    double expectedPhenotype;
    double testStat;
    
    gsl_vector* scores;
    gsl_vector* weights;
};

#endif /* genericBurdenTest_hpp */
