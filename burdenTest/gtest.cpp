#include <gtest/gtest.h>
#include <vector>
#include "wsbt.cpp"
#include "input.cpp"

using namespace std;

TEST(Mathtest, basicmath)
{
  ASSERT_EQ(4, 2+2);
}

wsbt setupWSBT(int subjectCount, int variantCount)
{
    gsl_rng * r;
    
    //Initializes the matrix of size subjectCount x variantCount with entry 1.
    vector<vector<int> > G(subjectCount, vector<int>(variantCount, 1));
    vector<double> maf(variantCount);
    vector<double> ptype(subjectCount);
    
    maf[0] = .1;
    ptype[0] = 0;
    for(int j = 1; j < variantCount; j++)
    {
        //maf[j] = maf[j-1] * .5;
        maf[j] = maf[j-1] * .9;
        ptype[j] = j;
    }

    for(int i = 0; i < subjectCount; i++)
    {
        //ptype[i] = gsl_ran_gaussian(r, 10) + 175;
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
    readInput result = readInput(vcfFile, "-vcf", phenoFile);
    vector<double> maf = result.getMaf();
    ASSERT_EQ(0.000199681, maf[0]);
    ASSERT_EQ(0.000798722, maf[1]);
    ASSERT_EQ(0.000199681, maf[2]);
}

TEST(generalRun, WeightedSumsBurdenTest)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    string vcfFile = "test3.vcf";
    string phenoFile = "test.pheno";
    readInput result(vcfFile, "-vcf", phenoFile);
    vector<double> ptype(2505);

    for(int i = 0; i < 2505; i++)
    {
        ptype[i] = gsl_ran_gaussian(r, 10) + 175;
    }
    wsbt test(result.getGenotype(), result.getMaf(), ptype);
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



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
