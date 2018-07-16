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
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

using namespace std;

class wsbt
{
public:
    wsbt(vector<vector<int> > unaffectedGenotype, vector<vector<int> > affectedGenotype);
    
    double getPvalue(){return pvalue;}
    std::vector<double> getWeights(){return weights;}
    
private:

    void setWeights();
    void setScores();
    double testStatistic();
    
    double pvalue;
    double testStat;
    
    
    
    vector<double> weights;
    vector<double> scores;
    vector<int> combinedIndividualCount;
    vector<int> unaffectedAlleleCount;
    vector<int> unaffectedIndvidualCount;
    vector<double> testStatistics;
    
    vector<vector<int> > unaffectedGenotype;
    vector<vector<int> > affectedGenotype;
};

#endif /* wsbt_hpp */
