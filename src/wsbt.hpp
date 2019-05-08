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
    wsbt(gsl_matrix* totalGenotype, int affectedCount, string gene);
    ~wsbt();
    
    double getPvalue(){return normpvalue;}
    double getPermPvalue(){return permpvalue;}
    gsl_vector* getWeights(){return initialWeights;}
    
private:
    void recalculate();
    double testStatistic();
    void shuffleMatrix();

    double normpvalue;
    double permpvalue;
    int affectedCount;
    int totalSubjects;
    bool verbose = false;

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