#include <gtest/gtest.h>
#include <vector>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "genericBurdenTest.cpp"
#include "input.cpp"

using namespace std;


vector<double> testRank(gsl_vector* inputScores)
{
    gsl_permutation* perm = gsl_permutation_calloc(inputScores->size);
    gsl_sort_vector_index(perm, inputScores);
    
    vector<double> ranks(inputScores->size);
    double currentRank = 1;
    double totalTiedSubjects = 1;
    
    /*
    cout << "Input Scores:";
    for(int j = 0; j < inputScores->size; j++)
    {
        cout << " " << gsl_vector_get(inputScores, j);
    }
    cout << endl;
    cout << "perm:";
    for(int j = 0; j < inputScores->size; j++)
    {
        cout << " " << gsl_permutation_get(perm, j);
    }
    cout << endl;
    */
    
    for(int j = 0; j < inputScores->size; j++)
    {
        //Tie ranker
        if (j + 1 < inputScores->size)
        {
            if (gsl_vector_get(inputScores, perm->data[j]) == gsl_vector_get(inputScores, perm->data[j+1]))
            {
                for(int i = j; gsl_vector_get(inputScores, perm->data[i]) == gsl_vector_get(inputScores, perm->data[i+1]); i++)
                {
                    totalTiedSubjects++;
                    //This checks if the next tie is also the last entry in the vector.
                    if (i + 1 == inputScores->size - 1)
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
                    ranks[i + j] = averageRank;
                }
                
                //-1 to place us at the last tie.
                j += totalTiedSubjects - 1;
                currentRank += totalTiedSubjects;
                totalTiedSubjects = 1;
            }
            else
            {
                ranks[j] = currentRank;
                currentRank++;
            }
        }
        //If the last element is unique.
        else
        {
            ranks[j] = currentRank;
        }
    }
    return ranks;
}

genericBurdenTest setupGenericBurdenTest(int subjectCount, int variantCount)
{
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    
    //Initializes the matrix of size subjectCount x variantCount with entry 1.
    vector<vector<int> > G(subjectCount, vector<int>(variantCount, 1));
    gsl_vector *maf = gsl_vector_alloc(variantCount);
    vector<double> ptype(subjectCount);
    
    gsl_vector_set(maf, 0, .1);
    ptype[0] = 0;
    for(int j = 1; j < variantCount; j++)
    {
        //maf[j] = maf[j-1] * .5;
        gsl_vector_set(maf, j, gsl_vector_get(maf, j-1) * .9);
        //maf[j] = maf[j-1] * .9;
        ptype[j] = j;
    }

    for(int i = 0; i < subjectCount; i++)
    {
        ptype[i] = gsl_ran_gaussian(r, 10) + 175;
    }
    return genericBurdenTest(G, maf, ptype);
}

TEST(Mathtest, allDifferentRank)
{
    ASSERT_EQ(4, 2+2);
    gsl_vector* scores = gsl_vector_alloc(10);
    for(int i = 0; i < 10; i++)
    {
        gsl_vector_set(scores, i, i+1);
    }
    
    vector<double> result = testRank(scores);
    
    for(int i = 0; i < 10; i++)
    {
        EXPECT_EQ(i+1, result[i]);
    }
    
    
    
}

TEST(Mathtest, allSameRank)
{
    gsl_vector* scores = gsl_vector_alloc(10);
    for(int i = 0; i < 10; i++)
    {
        gsl_vector_set(scores, i, 1);
    }
    
    vector<double> result = testRank(scores);
    
    for(int i = 0; i < 10; i++)
    {
        EXPECT_EQ((11.0/2.0), result[i]);
    }
    
    
}

TEST(Mathtest, mixedRank)
{
    gsl_vector* scores = gsl_vector_alloc(10);
    for(int i = 0; i < 10; i++)
    {
        if (i < 5)
        {
            gsl_vector_set(scores, i, 1);
        }
        else
        {
            gsl_vector_set(scores, i, 2);
        }
       
    }
    
    vector<double> result = testRank(scores);
    
    for(int i = 0; i < 10; i++)
    {
        if (i < 5)
        {
            EXPECT_EQ((3), result[i]);
        }
        else
        {
            EXPECT_EQ((8), result[i]);
        }
        
    }
}

TEST(Mathtest, mixedRank2)
{
    gsl_vector* scores = gsl_vector_alloc(11);
    for(int i = 0; i < 10; i++)
    {
        if (i < 5)
        {
            gsl_vector_set(scores, i, i);
        }
        else
        {
            gsl_vector_set(scores, i, 8);
        }
        
    }
    gsl_vector_set(scores, 10, 9);
    
    vector<double> result = testRank(scores);
    
    for(int i = 0; i < 10; i++)
    {
        if (i < 5)
        {
            EXPECT_EQ((i+1), result[i]);
        }
        else
        {
            EXPECT_EQ((8), result[i]);
        }
        
    }
    ASSERT_EQ(11, result[10]);
}

/*
TEST(setWeights,GenericBurdenTest)
{
    genericBurdenTest test = setupGenericBurdenTest(10, 10);
    vector<double> tempWeights = test.getWeights();
    ASSERT_EQ(tempWeights[0], gsl_ran_beta_pdf(.1, 1, 25));
    ASSERT_NE(tempWeights[0], gsl_ran_beta_pdf(.11, 1, 25));
    
    double temp = .1;
    for(int i = 1; i < 10; i++)
    {
        temp = temp * .9;
        ASSERT_EQ(tempWeights[i], gsl_ran_beta_pdf(temp, 1, 25));
    }
}

TEST(setScores, GenericBurdenTest)
{
    genericBurdenTest test = setupGenericBurdenTest(10, 10);
    vector<double> tempscore = test.getScores();
    vector<double> ptype;
    ptype.resize(10);
    for(int j = 0; j < 10; j++)
    {
        ASSERT_EQ(ptype[j], tempscore[j]);
    }
    
}

TEST(mafInput, fileInput)
{
    string vcfFile = "test3.vcf";
    string phenoFile = "test.pheno";
    readInput result("burden", "-vcf", vcfFile, vcfFile, phenoFile);
    vector<double> maf = result.getMaf();
    ASSERT_EQ(0.000199681, maf[0]);
    ASSERT_EQ(0.000798722, maf[1]);
    ASSERT_EQ(0.000199681, maf[2]);
}

TEST(generalRun, GenericBurdenTest)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    string vcfFile = "test3.vcf";
    string phenoFile = "test.pheno";
    readInput result("burden", "-vcf", vcfFile, vcfFile, phenoFile);
    vector<double> ptype(2505);

    for(int i = 0; i < 2505; i++)
    {
        ptype[i] = gsl_ran_gaussian(r, 10) + 175;
    }
    genericBurdenTest test(result.getGenotype(), result.getMaf(), ptype);
    cout << "Test P-value: " << test.getPvalue() << endl;
    cout << "Test Scores: ";
    for(int i = 0; i < test.getScores().size(); i++)
    {
        cout << test.getScores()[1] << " ";
    }
    cout << endl;
    cout << "Test Weights: ";
    for(int i = 0; i < test.getWeights().size(); i++)
    {
        cout << test.getWeights()[1] << " ";
    }
    cout << endl;
    gsl_rng_free(r);
    
}
*/


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
