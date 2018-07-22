//
//  wsbt.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//
#include <iostream>
#include <fstream>
#include "wsbt.hpp"


using namespace std;
/*
wsbt::wsbt(vector<vector<int> > GtypeU, vector<vector<int> > GtypeA)
{
    const int permutationCount = 1000;
    unsigned long int totalSubjects = unaffectedGenotype.size() + affectedGenotype.size();
    
    testStat = 0;
    pvalue = 0;
    unaffectedGenotype = GtypeU;
    affectedGenotype = GtypeA;
    
    weights.resize(unaffectedGenotype[0].size());
    //scores.resize(unaffectedGenotype.size() + affectedGenotype.size());
    testStatistics.resize(permutationCount);
    
    vector<vector<int> > total = GtypeU;
    total.insert(total.end(), GtypeA.begin(), GtypeA.end());
    
 
    gsl_vector *subjects = gsl_vector_alloc(totalSubjects);
    for(int i = 0; i < totalSubjects; i++)
    {
        gsl_vector_set(subjects, 1, i);
    }
    gsl_permutation *subjectPerm = gsl_permutation_alloc(totalSubjects);
    gsl_sort_vector_index (subjectPerm, subjects);
 
    
    
    for(int k = 0; k < permutationCount; k++)
    {
        setWeights();
        setScores();
        testStatistics[k] = testStatistic();
        
        //Need to shuffle unaffected/affected status
        vector<vector<int> > tempGtypeU(unaffectedGenotype.size());
        vector<vector<int> > tempGtypeA(affectedGenotype.size());
        
        
        
        
        //unaffectedGenotype = shuffledunaffected;
        //affectedGenotype = shuffledaffected;
        
        
    }
    
    //double testStatMean = gsl_stats_mean(testStatistics.data(), 1, testStatistics.size());
    //double testStatSigma = gsl_stats_sd(testStatistics.data(), 1, testStatistics.size());
    
}
*/

wsbt::wsbt(gsl_matrix* totalGtype, int aCount)
{
    const int permutationCount = 1000;
    unsigned long int totalSubjects = totalGtype->size2;
    
    ofstream outfile;
    affectedCount = aCount;
    testStat = 0;
    pvalue = 0;
    totalGenotype = totalGtype;
    testStatistics.resize(permutationCount);
    weights.resize(totalGtype->size1);
    scores = gsl_vector_alloc(totalGtype->size2);
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_permutation *subjectPerm = gsl_permutation_calloc(totalSubjects);
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lastTime = startTime;
    for(int k = 0; k < permutationCount; k++)
    {
        setWeights();
        setScores();
        /*
        cout << "First 10 of subjectPerm: ";
        for(int i = 0; i < 5; i++)
        {
            cout << subjectPerm->data[i] << " ";
        }
        cout << endl;
        cout << "First 10 of weights: ";
        for(int i = 0; i < 20; i++)
        {
            cout << weights[i] << " ";
        }
        cout << endl;
        cout << "First 10 of scores: ";
        for(int i = 0; i < 20; i++)
        {
            cout << scores->data[i] << " ";
        }
        cout << endl;
         */
        
        if(k == 0)
        {
            outfile.open("output.txt");
            outfile << "Weights: ";
            for(int i = 0; i < totalGtype->size1; i++)
            {
                outfile << weights[i] << " ";
            }
            outfile << endl;
            outfile << "Scores: ";
            for(int i = 0; i < totalGtype->size2; i++)
            {
                outfile << gsl_vector_get(scores, i) << " ";
            }
            outfile << endl;
            outfile.close();
        }
        testStatistics[k] = testStatistic();
        currentTime = chrono::high_resolution_clock::now();
        cout << "Permutation "<< k <<" took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lastTime).count() / 1000.0 << " Seconds."<< endl;
        lastTime = currentTime;
        cout << "TestStatistic for permutation " << k << " is " << testStatistics[k] << endl;
        double testStatMean = gsl_stats_mean(testStatistics.data(), 1, k+1);
        cout << "TestStatisticMean for permutation " << k << " is " << testStatMean << endl;
        double testStatSigma = gsl_stats_sd(testStatistics.data(), 1, k+1);
        cout << "TestStatisticSigma for permutation " << k << " is " << testStatSigma << endl;
        double zscore = (testStatistics[0] - testStatMean) / testStatSigma;
        cout << "Zscore for permutation " << k << " is " << zscore << endl;
        pvalue = gsl_ran_ugaussian_pdf(zscore);
        cout << "Pvalue for permutation " << k << " is " << pvalue << endl;
        
        
        //cout << "P-value at permutation " << k << " is " << pvalue << endl;
        //Permutes the data between affected and unaffected status.
        //subjectPerm = gsl_permutation_alloc(totalSubjects);
        gsl_ran_shuffle(r, subjectPerm->data, totalSubjects, sizeof(size_t));
        gsl_permute_matrix(subjectPerm, totalGenotype);
        
        /*
        for(int i = 0; i < 10; i++)
        {
            for(int j = 0; j < 10; j++)
            {
                cout << gsl_matrix_get(totalGenotype, i, j) << " ";
            }
            cout << endl;
        }
         */
        
    }
    outfile.open("statoutput.txt");
    outfile << "test statistics: ";
    for(int k = 0; k < permutationCount; k++)
    {
        outfile << testStatistics[k] << " ";
    }
    outfile << endl;
    outfile.close();
    
    cout << endl;
    gsl_rng_free(r);
    gsl_permutation_free(subjectPerm);
    gsl_vector_free(scores);
    gsl_matrix_free(totalGenotype);
}

void wsbt::setWeights()
{
    //Variants are in row(i) with subject in column(j).
    for(int i = 0; i < totalGenotype->size1; i++)
    {
        double mutantAllelesU = 0;
        double indivudualsU = 0;
        double totalVariant = 0;
        
        for(int j = 0; j < totalGenotype->size2; j++)
        {
            if(j < affectedCount && gsl_matrix_get(totalGenotype, i, j) != -1)
            {
                totalVariant++;
            }
            if(j >= affectedCount && gsl_matrix_get(totalGenotype, i, j) != -1)
            {
                mutantAllelesU += gsl_matrix_get(totalGenotype, i, j);
                indivudualsU++;
            }
        }
        totalVariant += indivudualsU;
        //cout << "Mutant Alleles: " << mutantAllelesU << endl;
        //cout << "Unaffected People: " << indivudualsU << endl;
        double upper = mutantAllelesU + 1.0;
        double lower = (2.0 * indivudualsU) + 2.0;
        double unaffectedRatio = upper / lower;
        /*
        if (unaffectedRatio < 0 || unaffectedRatio > 1)
        {
            cout << "The ratio is not between 0 and 1. It is "  << unaffectedRatio << endl;
            cout << "In this case the mutantAlleles = " << mutantAllelesU << endl;
            cout << "And the unaffected people = " << indivudualsU << endl;
            for(int j = 0; j < totalGenotype->size2; j++)
            {
                if (gsl_matrix_get(totalGenotype, i, j) > 2)
                {
                    cout << "Genotype = " << gsl_matrix_get(totalGenotype, i, j) << " at (" << i << "," << j << ")" << endl;
                }
            }
        }
         */
        //cout << "unaffectedRatio: " << unaffectedRatio << endl;
        double nancheck = sqrt(totalVariant * unaffectedRatio * (1.0 - unaffectedRatio));
        if (isnan(nancheck))
        {
            cout << "Weight for variant " << i <<  " is trying to sqrt a negative" << endl;
        }
        if (nancheck == 0)
        {
            cout << "Weight is somehow 0" << endl;
        }
        weights[i] = nancheck;
        //cout << "Weight: " << weights[i] << endl;
    }
}

void wsbt::setScores()
{
    for(int j = 0; j < totalGenotype->size2; j++)
    {
        double tempscore = 0.0;
        for(int i = 0; i < totalGenotype->size1; i++)
        {
            double genoData = gsl_matrix_get(totalGenotype, i, j);
            if(genoData == -1)
            {
                //This means the genotype in the file was ./.
            }
            else
            {
                tempscore = tempscore + (genoData / weights[i]);
            }
        }
        gsl_vector_set(scores, j, tempscore);
    }
}

double wsbt::testStatistic()
{
    //We have many subjects with the same score (namely 0). We need to assign them the average of their rank and then boost the following back to normal.
    testStat = 0;
    gsl_permutation * perm = gsl_permutation_calloc(totalGenotype->tda);
    gsl_vector * rank = gsl_vector_alloc(totalGenotype->size2);
    gsl_sort_vector_index(perm, scores);
    for(int j = 0; j < totalGenotype->size2; j++)
    {
        gsl_vector_set(rank, j, gsl_vector_get(scores, j));
    }
    gsl_sort_vector(rank);
    double currentRank = 1;
    //double currentScore = 0;
    double totalTiedSubjects = 0;
    
    
    for(int j = 0; j < scores->size; j++)
    {
        //Tie ranker
        if (j + 1 < scores->size)
        {
            if (gsl_vector_get(scores, perm->data[j]) == gsl_vector_get(scores, perm->data[j+1]))
            {
                for(int i = j; gsl_vector_get(scores, perm->data[i]) == gsl_vector_get(scores, perm->data[i+1]); i++)
                {
                    totalTiedSubjects++;
                    //This checks if the next tie is also the last entry in the vector.
                    if (i + 1 == scores->size - 1)
                    {
                        break;
                    }
                }
                //Calculate Average
                double averageRank = 0;
                for(int k = 0; k < totalTiedSubjects; k++)
                {
                    averageRank += currentRank + k;
                }
                averageRank = averageRank / (totalTiedSubjects);
                for(int i = 0; i < totalTiedSubjects; i++)
                {
                    gsl_vector_set(rank, i + j, averageRank);
                }
                
                //-1 to place us at the last tie.
                j += totalTiedSubjects - 1;
                currentRank += totalTiedSubjects;
                totalTiedSubjects = 1;
            }
            else
            {
                gsl_vector_set(rank, j, currentRank);
                currentRank++;
            }
        }
        //If the last element is unique.
        else
        {
            gsl_vector_set(rank, j, currentRank);
        }
    }
    
    cout << "Affected subjects have individual ranks: ";
    for(int j = 0; j < totalGenotype->tda; j++)
    {
        //Checks if the element at rank(j) was an affected subject and if so adds its rank to the test stat.
        if(gsl_permutation_get(perm, j) < affectedCount)
        {
            cout << gsl_vector_get(rank, j) << " ";
            testStat += gsl_vector_get(rank, j);
        }
        
    }
    cout << endl;
    cout << endl;
    gsl_permutation_free(perm);
    gsl_vector_free(rank);
    return testStat;
}
