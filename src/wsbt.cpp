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

wsbt::wsbt(gsl_matrix *totalGtype, int aCount, string gene, bool exactPvalueCalculation)
{
    const int permutationCount = 1000;
    bool verbose = false;
    bool exact = exactPvalueCalculation;

    ofstream outfile;
    affectedCount = aCount;
    normpvalue = 0;
    totalGenotype = totalGtype;
    geneName = gene;

    testStatistics = gsl_vector_calloc(permutationCount);
    scores = gsl_vector_calloc(totalGenotype->size2);
    weights = gsl_vector_calloc(totalGenotype->size1);
    initialWeights = gsl_vector_calloc(totalGenotype->size1);
    changedGenotype = gsl_matrix_calloc(totalGenotype->size1, affectedCount);
    gsl_matrix_set_all(changedGenotype, 1);

    r = gsl_rng_alloc(gsl_rng_default);
    subjectPerm = gsl_permutation_calloc(totalGenotype->size2);

    //Timing
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lastTime = startTime;
    
    
    if (exact)
    {
        recalculate();
        testStat = exactTestStatistic();
        normpvalue = gsl_cdf_ugaussian_P(testStat);
        //Low p-values correspond to when the variant is expressed more in sample than background.

        backMean = gsl_stats_mean(scores->data, 1, scores->size) * scores->size;
        for(int i = 0; i < affectedScores.size(); i++)
        {
            backMean -= affectedScores[i];
        }
        backMean = backMean / (scores->size - affectedCount);
    }
    else
    {
        //First run through data.
        recalculate();
        testStat = testStatistic();
        shuffleMatrix();
        //gsl_vector_memcpy(initialWeights, weights);

        //Permutations to get the mean and standard deviation of the test statistic.
        for (int k = 0; k < permutationCount; k++)
        {
            recalculate();
            gsl_vector_set(testStatistics, k, testStatistic());

            double testStatMean = gsl_stats_mean(testStatistics->data, 1, k + 1);
            double testStatSigma = gsl_stats_sd(testStatistics->data, 1, k + 1);
            double zscore = (testStat - testStatMean) / testStatSigma;
            normpvalue = gsl_cdf_ugaussian_P(zscore);
            //if (normpvalue > .5)
            //{
                normpvalue = 1 - normpvalue;
            //}
            //normpvalue = normpvalue * 2;

            //Reporting output.
            currentTime = chrono::high_resolution_clock::now();

            if (verbose || k + 1 == permutationCount)
            {
                cout << "Gene name: " << gene << endl;
                cout << "Subject count: " << totalGtype->size2 << endl;
                cout << "Variant count: " << totalGtype->size1 << endl;
                cout << k + 1 << " permutations took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
                cout << "TestStatistic is " << testStat << endl;
                cout << "TestStatisticMean for permutation " << k + 1 << " is " << testStatMean << endl;
                cout << "TestStatisticSigma for permutation " << k + 1 << " is " << testStatSigma << endl;
                cout << "Zscore for permutation " << k + 1 << " is " << zscore << endl;
                cout << "Pvalue for permutation " << k + 1 << " is " << normpvalue << endl
                     << endl;
            }
            lastTime = currentTime;

            //Shuffle for next permutation.
            shuffleMatrix();
        }
    }
}

//Only deletes gsl objects directly related to this instance of wsbt.
wsbt::~wsbt()
{
    gsl_rng_free(r);
    gsl_permutation_free(subjectPerm);
    gsl_vector_free(scores);
    gsl_vector_free(weights);
    gsl_vector_free(initialWeights);
    gsl_vector_free(testStatistics);
    gsl_matrix_free(changedGenotype);
}


double wsbt::testStatistic()
{
    //We have many subjects with the same score (namely 0). We need to assign them the average of their rank and then boost the following back to normal.
    double testStat = 0;
    double currentRank = 1;
    double totalTiedSubjects = 1;
    gsl_permutation *perm = gsl_permutation_calloc(totalGenotype->size2);
    gsl_vector *rank = gsl_vector_calloc(totalGenotype->size2);
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
    
    for(int j = 0; j < totalGenotype->size2; j++)
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

double wsbt::exactTestStatistic()
{
    //We have many subjects with the same score (namely 0). We need to assign them the average of their rank and then boost the following back to normal.
    double currentRank = 1;
    double totalTiedSubjects = 1;
    vector<int> tieTracker;

    for(int i = 0; i < affectedCount; i++)
    {
        affectedScores.push_back(gsl_vector_get(scores, i));
    }

    gsl_permutation *perm = gsl_permutation_calloc(totalGenotype->size2);
    gsl_vector *rank = gsl_vector_calloc(totalGenotype->size2);
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

                //Keep track of number of tied ranks for standard deviation calculation.
                tieTracker.push_back(totalTiedSubjects);
                
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

        cout << "Affected subjects have individual ranks: ";
        for (int j = 0; j < totalGenotype->size2; j++)
        {
            //Checks if the element at rank(j) was an affected subject and if so adds its rank to the test stat.
            if (gsl_permutation_get(perm, j) < affectedCount)
            {
                cout << gsl_vector_get(rank, j) << " ";
            }
        }
        cout << endl;
    }

    //cout << "Finished calculating ranks." << endl;
    //cout << "Total of " << tieTracker.size() << " tied ranks." << endl;
    
    for(int j = 0; j < totalGenotype->size2; j++)
    {
        //Checks if the element at rank(j) was an affected subject and if so adds its rank to U1.
        if(gsl_permutation_get(perm, j) < affectedCount)
        {
            U1 += gsl_vector_get(rank, j);
        }
        else
        {
            U2 += gsl_vector_get(rank, j);
        }
    }
    
    double n1 = affectedCount;
    double n2 = totalGenotype->size2 - affectedCount;
    double result = 0;

    U1 -= n1 * (n1 + 1) / 2.0;
    U2 -= n2 * (n2 + 1) / 2.0;
    //cout << "(U1, U2): " << "(" << U1 << ", " << U2 << ")" << endl;
    
    /*
    if(U1 < U2)
    {
        //cout << "U1 < U2: " << U1 << "<" << U2 << endl;
        result = U1;
    }
    else
    {
        //cout << "U2 < U1: " << U2 << "<" << U1 << endl;
        result = U2;
    }
    */
    
    
    result = U2;

    //The data set will likely always have more than 20 total subjects but we should implement the other distribution calculation later anyway.
    if(totalGenotype->size2 > 20)
    {
        //The following calculation does the tie correction for the sd. In the case there arent many ties, this is just slightly slower to run.
        //cout << result << " - " << n1 * n2 / 2.0 << endl;
        result -= n1 * n2 / 2.0;

        double sd = 0;
        double tieSummation = 0;
        for(int i = 0; i < tieTracker.size(); i++)
        {
            //cout << tieTracker[i] << ", ";
            tieSummation += (pow(tieTracker[i], 3) - tieTracker[i]);
        }

        tieSummation = tieSummation / (totalGenotype->size2 * 1.0 * (totalGenotype->size2 - 1));
        sd = (n1 * n2 / 12.0) *  (totalGenotype->size2 + 1 - tieSummation);
        sd = sqrt(sd);
        //cout << "Sd: " << sd << endl;
        result = result / sd;
    }
    
    gsl_permutation_free(perm);
    gsl_vector_free(rank);
    return result;
}


void wsbt::recalculate()
{
    for(int i = 0; i < changedGenotype->size1; i++)
    {
        for(int j = 0; j < affectedCount; j++)
        {
            //Found a weight to recalculate.
            if(gsl_matrix_get(changedGenotype, i, j) == 1)
            {
                //First we remove the old weight's effect on all of the scores.
                if(gsl_vector_get(weights, i) != 0)   //Just want to skip the following loop in the first weight calculation.
                {
                    for (int k = 0; k < totalGenotype->size2; k++)
                    {
                        if (gsl_matrix_get(totalGenotype, i, k) > 0)
                        {
                            double oldWeightedScoreComponent = gsl_matrix_get(totalGenotype, i, k) / gsl_vector_get(weights, i);
                            gsl_vector_set(scores, k, gsl_vector_get(scores, k) - oldWeightedScoreComponent);
                        }
                    }
                }

                //Next we recalculate the weight for the variant.
                int totalVariant = 0;
                int mutantAllelesU = 0;
                int indivudualsU = 0;
                for (int k = 0; k < totalGenotype->size2; k++)
                {
                    if (k < affectedCount)
                    {
                        if (gsl_matrix_get(totalGenotype, i, k) > -1)
                        {
                            totalVariant++;
                        }
                    }
                    else
                    {
                        //We look for variants that have more than 0 alleles because we need individuals genotyped for the variant.
                        if(gsl_matrix_get(totalGenotype, i, k) > -1)
                        {
                            mutantAllelesU += gsl_matrix_get(totalGenotype, i, k);
                            indivudualsU++;
                        }
                    }
                }
                totalVariant += indivudualsU;
                double upper = mutantAllelesU + 1.0;
                double lower = (2.0 * indivudualsU) + 2.0;
                double unaffectedRatio = upper / lower;
                if(unaffectedRatio >= 1)
                {
                    cout << "Uh oh. Qi is " << unaffectedRatio << endl;
                }
                double nancheck = sqrt(totalVariant * unaffectedRatio * (1.0 - unaffectedRatio));
                gsl_vector_set(weights, i, nancheck);

                //Finally we update all scores using the new weight.
                for(int k = 0; k < scores->size; k++)
                {
                    if(gsl_matrix_get(totalGenotype, i, k) > 0)
                    {
                        double newWeightedScoreComponent = gsl_matrix_get(totalGenotype, i, k) / gsl_vector_get(weights, i);
                        gsl_vector_set(scores, k, gsl_vector_get(scores, k) + newWeightedScoreComponent);
                    }
                }
                //Only need to recalculate once per variant.
                break;
            }
        }
    }
}

void wsbt::shuffleMatrix()
{
    //Permutes the minimum required columns of the genotype matrix for the next permutation.
    gsl_ran_shuffle(r, subjectPerm->data, totalGenotype->size2, sizeof(size_t));
    //Check what needs to be updated for the new round since if we dont need to update much we can save a lot of time.
    for (int i = 0; i < totalGenotype->size1; i++)
    {
        for (int j = 0; j < affectedCount; j++)
        {
            //Check if current affected is the same as the new affected for variant k.
            if (gsl_matrix_get(totalGenotype, i, j) == gsl_matrix_get(totalGenotype, i, subjectPerm->data[j]))
            {
                gsl_matrix_set(changedGenotype, i, j, 0);
            }
            else
            {
                gsl_matrix_set(changedGenotype, i, j, 1);
            }
        }
    }
    if (affectedCount < totalGenotype->size2 / 2)
    {
        //Shuffles affected status.
        for (int j = 0; j < affectedCount; j++)
        {
            gsl_matrix_swap_columns(totalGenotype, j, subjectPerm->data[j]);
            gsl_vector_swap_elements(scores, j, subjectPerm->data[j]);
        }
    }
    else
    {
        //Shuffles unaffected status.
        for (int j = affectedCount; j < totalGenotype->size2; j++)
        {
            gsl_matrix_swap_columns(totalGenotype, j, subjectPerm->data[j]);
            gsl_vector_swap_elements(scores, j, subjectPerm->data[j]);
        }
    }
}


void wsbt::driverOutput()
{
    cout << "Gene name: " << geneName << endl;
    cout << "Variant Count: " << totalGenotype->size1 << endl;
    if (affectedCount == 1)
    {
        cout << "Score: " << affectedScores[0] << endl;
    }
    else
    {
        cout << "Scores: ";
        for (int i = 0; i < affectedScores.size(); i++)
        {
            cout << affectedScores[i] << " ";
        }
        cout << endl;
    }
    cout << "Mean Background Score: " << backMean << endl;
    cout << "U_Case: " << U1 << endl;
    cout << "U_Control: " << U2 << endl;
    cout << "Test Statistic: " << testStat << endl;
    cout << "Exact P-value: " << normpvalue << endl;
    cout << endl;
}


