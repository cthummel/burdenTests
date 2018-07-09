//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"


readInput::readInput(string filename)
{
    
    std::string genoCommand, mafCommand;
    genoCommand = "./third/vcflib/bin/vcfkeepgeno " + filename + " GT DP";
    mafCommand = "./third/vcflib/bin/vcfkeepgeno genotempfile.vcf AF";
    
    
    ofstream tempfile;
    tempfile.open("genotempfile.vcf");
    //tempfile << system(genoCommand);
    tempfile.close();
    
    tempfile.open("finaltempfile.vcf");
    //tempfile << system(mafCommand);
    tempfile.close();
};
