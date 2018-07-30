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
#include <fstream>
#include <iostream>
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
    wsbt(gsl_matrix* totalGenotype, int affectedCount, gsl_vector* maf);
    
    double getPvalue(){return pvalue;}
    gsl_vector* getWeights(){return weights;}
    
private:
    void setWeights();
    void setScores();
    double testStatistic();
    
    double pvalue;
    int affectedCount;
    
    gsl_matrix* totalGenotype;
    gsl_vector* scores;
    gsl_vector* testStatistics;
    gsl_vector* weights;
};

#endif /* wsbt_hpp */
