//
//  wsbt.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//
#include <iostream>
#include "wsbt.hpp"


using namespace std;

wsbt::wsbt(vector<vector<int> > GtypeU, vector<vector<int> > GtypeA)
{
    const int permutationCount = 1000;
    testStat = 0;
    pvalue = 0;
    unaffectedGenotype = GtypeU;
    affectedGenotype = GtypeA;
    weights.resize(unaffectedGenotype[0].size());
    scores.resize(unaffectedGenotype.size() + affectedGenotype.size());
    testStatistics.resize(permutationCount);
    
    for(int k = 0; k < permutationCount; k++)
    {
        setWeights();
        setScores();
        testStatistics[k] = testStatistic();
        
    }
    
    double testStatMean = gsl_stats_mean(testStatistics.data(), 1, testStatistics.size());
    double testStatSigma = gsl_stats_sd(testStatistics.data(), 1, testStatistics.size());
    
}

//Set weights for variants based off the proportion of unaffected individuals with the variant.
void wsbt::setWeights()
{
    for(int j = 0; j < unaffectedGenotype[0].size(); j++)
    {
        int mutantAllelesU = 0;
        int indivudualsU = 0;
        int totalVariant = 0;
        
        for(int i = 0; i < unaffectedGenotype.size(); i++)
        {
            mutantAllelesU += unaffectedGenotype[i][j];
            if(unaffectedGenotype[i][j] != 0)
            {
                indivudualsU++;
            }
        }
        for(int i = 0; i < affectedGenotype.size(); i++)
        {
            if(affectedGenotype[i][j] != 0)
            {
                totalVariant++;
            }
        }
        
        totalVariant += indivudualsU;
        int unaffectedRatio = (mutantAllelesU + 1) / (2 * indivudualsU + 2);
        
        weights[j] = totalVariant * unaffectedRatio * (1 - unaffectedRatio);
        
    }
    
}

//Sets the score for each indivudal by checking their genotype for each variant and weighting it with the variant weight.
void wsbt::setScores()
{
    for(int i = 0; i < unaffectedGenotype.size(); i++)
    {
        int tempscore = 0;
        for(int j = 0; j < unaffectedGenotype[0].size(); j++)
        {
            tempscore += unaffectedGenotype[i][j] / weights[j];
        }
        scores[i] = tempscore;
    }
    for(int i = 0; i < affectedGenotype.size(); i++)
    {
        int tempscore = 0;
        for(int j = 0; j < affectedGenotype[0].size(); j++)
        {
            tempscore += affectedGenotype[i][j] / weights[j];
        }
        scores[i + unaffectedGenotype.size()] = tempscore;
    }
}

//Ranks the individual scores.
double wsbt::testStatistic()
{
    unsigned long int n = scores.size();
    testStat = 0;
    
    gsl_vector *v = gsl_vector_alloc(n);
    for(int i = 0; i < scores.size(); i++)
    {
        gsl_vector_set(v, 1, scores[i]);
    }
    
    gsl_permutation * perm = gsl_permutation_alloc(n);
    gsl_permutation * rank = gsl_permutation_alloc(n);
    
    //It sorts and then stores in perm the index of the score that would have been placed there.
    //In other words, perm index is rank and perm[index] is the index of the score that was matched there.
    gsl_sort_vector_index (perm, v);
    
    //Scores was initially set up that all the unaffected genotype scores are in the front of the vector and
    //the affected genotype scores are at the end.
    for(int i = 0; i < scores.size(); i++)
    {
        if(gsl_permutation_get(perm, i) >= unaffectedGenotype.size())
        {
            testStat += i;
        }
    }
    
    //maybe easier to get rank from this.
    //gsl_permutation_inverse (rank, perm);
    
    gsl_vector_free(v);
    gsl_permutation_free(perm);
    gsl_permutation_free(rank);
    
    return testStat;
}
