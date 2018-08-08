//
//  output.hpp
//  burdenTest
//
//  Created by Corin Thummel on 7/30/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#ifndef output_hpp
#define output_hpp

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <regex>
#include <gsl/gsl_vector.h>

using namespace std;

class writeOutput
{
    
public:
    writeOutput(string filename, string test_type, gsl_vector* weights);
    
    
private:
    ofstream outfile;
    ifstream infile;
    string bcftools_location = "externals/bcftools/src/bcftools_project/bcftools";
    
    
    
};







#endif /* output_hpp */
