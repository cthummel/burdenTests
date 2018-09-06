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
    
    ofstream outfile;
    affectedCount = aCount;
    pvalue = 0;
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
    
    //Permutations to get the mean and standard deviation of the test statistic.
    for(int k = 0; k < permutationCount; k++)
    {
        setWeights();
        setScores();
        
        //Just to save the initial weights and scores for the data set.
        if(k == 0)
        {
            //writeOutput out("test5.vcf", weights);
            outfile.open("output.txt");
            for(int i = 0; i < totalGtype->size1; i++)
            {
                outfile << "Weight " << i << ": " << gsl_vector_get(weights, i) << " " << gsl_ran_beta_pdf(gsl_vector_get(inputMaf, i),1,25) << endl;
            }
            outfile << endl;
            outfile << "Scores: ";
            for(int i = 0; i < totalGtype->size2; i++)
            {
                outfile << gsl_vector_get(scores, i) << " ";
            }
            outfile << endl;
            outfile.close();
            gsl_vector_memcpy(initialWeights, weights);
        }
        gsl_vector_set(testStatistics, k, testStatistic());
        currentTime = chrono::high_resolution_clock::now();
        cout << "Permutation "<< k <<" took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lastTime).count() / 1000.0 << " seconds."<< endl;
        lastTime = currentTime;
        cout << "TestStatistic for permutation " << k << " is " << gsl_vector_get(testStatistics, k) << endl;
        double testStatMean = gsl_stats_mean(testStatistics->data, 1, k+1);
        cout << "TestStatisticMean for permutation " << k << " is " << testStatMean << endl;
        double testStatSigma = gsl_stats_sd(testStatistics->data, 1, k+1);
        cout << "TestStatisticSigma for permutation " << k << " is " << testStatSigma << endl;
        double zscore = (gsl_vector_get(testStatistics, 0) - testStatMean) / testStatSigma;
        cout << "Zscore for permutation " << k << " is " << zscore << endl;
        //Two-sided P-value calculation.
        pvalue = 2 * gsl_cdf_ugaussian_P(zscore);
        cout << "Pvalue for permutation " << k << " is " << pvalue << endl;
        cout << endl;
        
        //Permutes the columns of the genotype matrix for the next permutation.
        gsl_ran_shuffle(r, subjectPerm->data, totalSubjects, sizeof(size_t));
        for(int i = 0; i < affectedCount; i++)
        {
            gsl_matrix_swap_columns(totalGenotype, i, subjectPerm->data[i]);
        }
        //gsl_permute_matrix(subjectPerm, totalGenotype);
    }
    
    cout << endl;
    gsl_rng_free(r);
    gsl_permutation_free(subjectPerm);
    gsl_vector_free(scores);
    gsl_vector_free(weights);
    gsl_vector_free(testStatistics);
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
        double upper = mutantAllelesU + 1.0;
        double lower = (2.0 * indivudualsU) + 2.0;
        double unaffectedRatio = upper / lower;
        
	double nancheck = sqrt(totalVariant * unaffectedRatio * (1.0 - unaffectedRatio));
        /*
	if (isnan(nancheck))
        {
            cout << "Weight for variant " << i <<  " is trying to sqrt a negative" << endl;
        }
        if (nancheck == 0)
        {
            cout << "Weight is somehow 0" << endl;
        }
	*/
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
    cout << "Affected subjects have individual scores: ";
    for(int j = 0; j < affectedCount; j++)
    {
        cout << scores->data[j] << " ";
    }
    cout << endl;
    
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
    gsl_permutation_free(perm);
    gsl_vector_free(rank);
    return testStat;
}
