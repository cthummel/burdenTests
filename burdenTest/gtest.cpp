#include <gtest/gtest.h>
#include <vector>
#include "wsbt.cpp"
//#include "input.cpp"

TEST(Mathtest, basicmath)
{
  ASSERT_EQ(4, 2+2);
}

void setupWSBT()
{
    
}

TEST(WeightedSumsBurdenTest, setWeights)
{
    std::vector<std::vector<int> > G;
    std::vector<double> maf = ;
    std::vector<double> ptype = ;
    
    
    
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
