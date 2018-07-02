//
//  wsbt.hpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#ifndef wsbt_hpp
#define wsbt_hpp

#include <stdio.h>
#include <vector>
#include "third/gsl-2.4/gsl/gsl_randist.h"


class wsbt
{
public:
    wsbt();
    ~wsbt();
    
    double getPvalue(){return pvalue;}
    std::vector<double> getScores(){return scores;}
    
    
    
private:
    void setWeights();
    void setScores();
    void testStatistic();
    
    double pvalue;
    double expectedPhenotype;
    double testStat;
    
    std::vector<std::vector<int> > genotypeMatrix;
    
    std::vector<double> phenotype;
    std::vector<double> scores;
    std::vector<double> weights;
    std::vector<double> maf;
};

#endif /* wsbt_hpp */
