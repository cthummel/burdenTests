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
    genericBurdenTest(std::vector<std::vector<int> > G, std::vector<double> maf, std::vector<double> ptype);
    
    double getPvalue(){return pvalue;}
    std::vector<double> getScores(){return scores;}
    std::vector<double> getWeights(){return weights;}
    
private:
    void setWeights(std::vector<double> maf);
    void setScores(std::vector<std::vector<int> > G, std::vector<double> ptype);
    void testStatistic();
    
    double pvalue;
    double expectedPhenotype;
    double testStat;
    
    std::vector<double> scores;
    std::vector<double> weights;
};

#endif /* genericBurdenTest_hpp */
