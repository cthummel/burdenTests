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
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_matrix.h>
#include <gsl/gsl_permutation.h>

using namespace std;

class wsbt
{
public:
    wsbt(vector<vector<int> > unaffectedGenotype, vector<vector<int> > affectedGenotype);
    wsbt(gsl_matrix* totalGenotype, int affectedCount);
    
    double getPvalue(){return pvalue;}
    std::vector<double> getWeights(){return weights;}
    
private:
    void setWeights();
    void setScores();
    double testStatistic();
    
    double pvalue;
    double testStat;
    int affectedCount;
    
    vector<double> weights;
    //vector<double> scores;
    vector<double> testStatistics;
    
    vector<vector<int> > unaffectedGenotype;
    vector<vector<int> > affectedGenotype;
    gsl_matrix* totalGenotype;
    gsl_vector* scores;
};

#endif /* wsbt_hpp */
