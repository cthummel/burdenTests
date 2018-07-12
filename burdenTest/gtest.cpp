#include <gtest/gtest.h>
#include <vector>
#include "wsbt.cpp"
#include "input.cpp"

using namespace std;

TEST(Mathtest, basicmath)
{
  ASSERT_EQ(4, 2+2);
}

wsbt setupWSBT(int sCount, int vCount)
{
    //gsl_rng * r;
    
    //Initializes the matrix of size sCount x vCount with entry 1.
    vector<vector<int> > G(sCount, vector<int>(vCount, 1));
    vector<double> maf(vCount);
    vector<double> ptype(sCount);
    
    maf[0] = .1;
    ptype[0] = 0;
    for(int j = 1; j < vCount; j++)
    {
        //maf[j] = maf[j-1] * .5;
        maf[j] = maf[j-1] * .9;
        ptype[j] = j;
    }

    for(int i = 0; i < sCount; i++)
    {
        
    }
    
    return wsbt(G, maf, ptype);
    
    
}

TEST(setWeights,WeightedSumsBurdenTest)
{
    wsbt test = setupWSBT(10, 10);
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

TEST(setScores, WeightedSumsBurdenTest)
{
    wsbt test = setupWSBT(10, 10);
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
    readInput result = readInput(vcfFile, phenoFile);
    vector<double> maf = result.getMaf();
    ASSERT_EQ(0.000199681, maf[0]);
    ASSERT_EQ(0.000798722, maf[1]);
    ASSERT_EQ(0.000199681, maf[2]);
}



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
