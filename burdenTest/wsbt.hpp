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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


class wsbt
{
public:
    wsbt(std::vector<std::vector<int> > G, std::vector<double> maf, std::vector<double> ptype);
    //~wsbt();
    
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
    
    //std::vector<std::vector<int> > genotypeMatrix;
    
    //std::vector<double> phenotype;
    std::vector<double> scores;
    std::vector<double> weights;
    //std::vector<double> maf;
};

#endif /* wsbt_hpp */
