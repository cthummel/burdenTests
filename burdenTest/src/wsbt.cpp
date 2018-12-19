//
//  wsbt.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//
#include "wsbt.hpp"
#include <chrono>

using namespace std;

wsbt::wsbt(gsl_matrix* totalGtype, int aCount, gsl_vector *inputMaf)
{
    const int permutationCount = 100;
    unsigned long int totalSubjects = totalGtype->size2;
    bool verbose = false;
    
    ofstream outfile;
    affectedCount = aCount;
    normpvalue = 0;
    totalGenotype = totalGtype;
    
    testStatistics = gsl_vector_alloc(permutationCount);
    scores = gsl_vector_alloc(totalGenotype->size2);
    weights = gsl_vector_alloc(totalGenotype->size1);
    initialWeights = gsl_vector_alloc(totalGenotype->size1);
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_permutation *subjectPerm = gsl_permutation_calloc(totalSubjects);
    
    //Timing
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lastTime = startTime;

    //First run through data.
    double testStat;
    setWeights();
    setScores();
    testStat = testStatistic();
    gsl_vector_memcpy(initialWeights, weights);
    gsl_ran_shuffle(r, subjectPerm->data, totalSubjects, sizeof(size_t));
    for (int i = 0; i < affectedCount; i++)
    {
        gsl_matrix_swap_columns(totalGenotype, i, subjectPerm->data[i]);
    }

    //Permutations to get the mean and standard deviation of the test statistic.
    for(int k = 0; k < permutationCount; k++)
    {
        setWeights();
        setScores();
        gsl_vector_set(testStatistics, k, testStatistic());
        
        double testStatMean = gsl_stats_mean(testStatistics->data, 1, k+1);
        double testStatSigma = gsl_stats_sd(testStatistics->data, 1, k+1);
        double zscore = (testStat - testStatMean) / testStatSigma;
        normpvalue = gsl_cdf_ugaussian_P(zscore);
        if (normpvalue > .5)
        {
            normpvalue = 1 - normpvalue;
        }
        normpvalue = normpvalue * 2;
        
        //Reporting output. 
        currentTime = chrono::high_resolution_clock::now();
        if(verbose || k + 1 == permutationCount)
        {
            cout << "Permutation " << k + 1 << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count() / 1000.0 << " seconds." << endl;
            cout << "TestStatistic for permutation " << k + 1 << " is " << gsl_vector_get(testStatistics, k) << endl;
            cout << "TestStatisticMean for permutation " << k + 1 << " is " << testStatMean << endl;
            cout << "TestStatisticSigma for permutation " << k + 1 << " is " << testStatSigma << endl;
            cout << "Zscore for permutation " << k + 1 << " is " << zscore << endl;
            cout << "Pvalue for permutation " << k + 1 << " is " << normpvalue << endl << endl;
        }
        lastTime = currentTime;
        
        
        //Permutes the minimum required columns of the genotype matrix for the next permutation.
        gsl_ran_shuffle(r, subjectPerm->data, totalSubjects, sizeof(size_t));
        if(affectedCount < totalSubjects / 2)
        {
            //Shuffles affected status.
            for (int i = 0; i < affectedCount; i++)
            {
                gsl_matrix_swap_columns(totalGenotype, i, subjectPerm->data[i]);
            }
        }
        else
        {
            //Shuffles unaffected status.
            for (int i = affectedCount; i < totalSubjects; i++)
            {
                gsl_matrix_swap_columns(totalGenotype, i, subjectPerm->data[i]);
            }
        }
    }
    int extremeCount = 1;
    for(int i = 1; i < permutationCount; i++)
    {
        if(testStat > 0)
        {
            if(testStat <= gsl_vector_get(testStatistics, i))
            {
                extremeCount++;
            }
        }
        else
        {
            if(testStat >= gsl_vector_get(testStatistics, i))
            {
                extremeCount++;
            }
        }
    }
    permpvalue = (1.0 * extremeCount) / (permutationCount + 1);

    gsl_rng_free(r);
    gsl_permutation_free(subjectPerm);
    gsl_vector_free(scores);
    gsl_vector_free(weights);
    gsl_vector_free(testStatistics);
}

void wsbt::setWeights()
{
    //Variants are in row(i) with subject in column(j).
    for (int i = 0; i < totalGenotype->size1; i++)
    {
        double mutantAllelesU = 0;
        double indivudualsU = 0;
        double totalVariant = 0;

        for (int j = 0; j < totalGenotype->size2; j++)
        {
            if (j < affectedCount && gsl_matrix_get(totalGenotype, i, j) != -1)
            {
                totalVariant++;
            }
            if (j >= affectedCount && gsl_matrix_get(totalGenotype, i, j) != -1)
            {
                mutantAllelesU += gsl_matrix_get(totalGenotype, i, j);
                indivudualsU++;
            }
        }
        totalVariant += indivudualsU;
        double upper = mutantAllelesU + 1.0;
        double lower = (2.0 * indivudualsU) + 2.0;
        double unaffectedRatio = upper / lower;
        double nancheck = sqrt(totalVariant * unaffectedRatio * (1.0 - unaffectedRatio));

        gsl_vector_set(weights, i, nancheck);
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
            else if (genoData > 0)
            {
                tempscore = tempscore + (genoData / gsl_vector_get(weights, i));
            }
        }
        gsl_vector_set(scores, j, tempscore);
    }
}

double wsbt::testStatistic()
{
    //We have many subjects with the same score (namely 0). We need to assign them the average of their rank and then boost the following back to normal.
    double testStat = 0;
    double currentRank = 1;
    double totalTiedSubjects = 0;
    gsl_permutation * perm = gsl_permutation_calloc(totalGenotype->tda);
    gsl_vector * rank = gsl_vector_alloc(totalGenotype->size2);
    gsl_sort_vector_index(perm, scores);

    for(int j = 0; j < scores->size; j++)
    {
        //Ranks ties in score.
        if (j + 1 < scores->size)
        {
            if (gsl_vector_get(scores, perm->data[j]) == gsl_vector_get(scores, perm->data[j+1]))
            {
                //Count the number of places that tie in rank.
                for(int i = j; gsl_vector_get(scores, perm->data[i]) == gsl_vector_get(scores, perm->data[i+1]); i++)
                {
                    totalTiedSubjects++;
                    //This checks if the next tie is also the last entry in the vector.
                    if (i + 1 == scores->size - 1)
                    {
                        break;
                    }
                }
                //Calculate average rank.
                double averageRank = 0;
                for(int k = 0; k < totalTiedSubjects; k++)
                {
                    averageRank += currentRank + k;
                }
                averageRank = averageRank / (totalTiedSubjects);
                
                //Save calculated ranks in vector.
                for(int i = 0; i < totalTiedSubjects; i++)
                {
                    gsl_vector_set(rank, i + j, averageRank);
                }
                
                //Reset.
                j += totalTiedSubjects - 1;  //I used -1 since we are about to j++ in the for loop.
                currentRank += totalTiedSubjects;
                totalTiedSubjects = 1;
            }
            //Unique score gets unique rank.
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
    if(verbose)
    {
        cout << "Affected subjects have individual scores: ";
        for (int j = 0; j < affectedCount; j++)
        {
            cout << scores->data[j] << " ";
        }
        cout << endl;
    }
    
    if (verbose)
    {
        cout << "Affected subjects have individual ranks: ";
    }
    
    for(int j = 0; j < totalGenotype->tda; j++)
    {
        //Checks if the element at rank(j) was an affected subject and if so adds its rank to the test stat.
        if(gsl_permutation_get(perm, j) < affectedCount)
        {
            if(verbose)
            {
                cout << gsl_vector_get(rank, j) << " ";
            }
            testStat += gsl_vector_get(rank, j);
        }
        
    }
    if(verbose)
    {
        cout << endl;
    }
    gsl_permutation_free(perm);
    gsl_vector_free(rank);
    return testStat;
}
