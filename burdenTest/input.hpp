//
//  input.hpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#ifndef input_hpp
#define input_hpp

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>

using namespace std;

class readInput
{
    
public:
    readInput(string vcfFile, string phenoFile);
    
    vector<vector<int> > getGenotype(){return genotypeMatrix;}
    vector<double> getMaf(){return maf;}

private:
    int variantCount;
    int subjectCount;
    
    vector<vector<int> > genotypeMatrix;
    vector<double> maf;
};

#endif /* input_hpp */
