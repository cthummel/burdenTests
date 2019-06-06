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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_movstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>
using namespace std;

class wsbt
{
public:
    wsbt(gsl_matrix* totalGenotype, int affectedCount, string gene, bool);
    ~wsbt();
    
    double getPvalue(){return normpvalue;}
    double getScore(){return gsl_vector_get(scores, 0);}
    double getTestStat(){return testStat;}
    gsl_vector* getWeights(){return initialWeights;}
    
private:
    void recalculate();
    double testStatistic();
    double exactTestStatistic();
    void shuffleMatrix();

    double normpvalue;
    double testStat;
    int affectedCount;
    int totalSubjects;
    bool verbose = false;
    string geneName;

    gsl_rng *r = nullptr;
    gsl_permutation *subjectPerm = nullptr;
    
    gsl_matrix* totalGenotype = nullptr;
    gsl_matrix* changedGenotype = nullptr;
    gsl_vector* scores = nullptr;
    gsl_vector* testStatistics = nullptr;
    gsl_vector* weights = nullptr;
    gsl_vector* initialWeights = nullptr;

};

#endif /* wsbt_hpp */
