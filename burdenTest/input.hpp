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
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class readInput
{
    
public:
    readInput(string inputFile);
    
    
    
private:
    vector<vector<int> > genotypeMatrix;
    vector<double> maf;
    
    
    
    
};





#endif /* input_hpp */
