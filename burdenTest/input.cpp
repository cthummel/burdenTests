//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"
using namespace std;

readInput::readInput(string vcfFile, string vcfType, string phenoFile)
{
    variantCount = 0;
    subjectCount = 0;
    
    string line;
    string genoCommand, mafCommand;
    smatch match;
    ifstream tempfile;
    
    regex gMatch("([0-9])(\\/|\\|)([0-9])");
    regex mafMatch(";AF=(0\\.\\d*)");
    regex subjectCountMatch(";NS=(\\d*)");
    regex headerMatch("^#");
    
    genoCommand = "vcftools -" + vcfType + " "  + vcfFile + " --extract-FORMAT-info GT";
    mafCommand = "vcftools -" + vcfType + " "  + vcfFile + " --get-INFO AF";
    
    //Runs vcftools command to separate genotype data.
    system(genoCommand.c_str());
    
    //Obtain variant counts and subject counts.
    tempfile.open(vcfFile);
    for(int j = 0; getline(tempfile, line); j++)
    {
        if(regex_search(line, match, subjectCountMatch))
        {
            subjectCount = stoi(match[1]);
        }
        
        if(!regex_search(line, match, headerMatch))
        {
            variantCount++;
        }
    }
    //cout << "variantCount: " << variantCount << endl;
    //cout << "subjectCount: " << subjectCount << endl;
    tempfile.close();
    
    //Parsing genotype data
    tempfile.open("out.GT.FORMAT");
    genotypeMatrix = vector<vector<int> >(subjectCount, vector<int>(variantCount, 1));
    
    for(int j = 0; getline(tempfile, line); j++)
    {
        for(int i = 0; regex_search(line, match, gMatch); i++)
        {
            //Since the first line of the file is just headers, we input into j-1
            genotypeMatrix[i][j-1] = stoi(match[1]) + stoi(match[3]);
            line = match.suffix();
        }
    }
    tempfile.close();
    
    //system(mafCommand.c_str());
    //tempfile.open("out.INFO");
    
    //Parse the maf
    tempfile.open(vcfFile);
    maf = vector<double>(variantCount);
    
    int mafPos = 0;
    for(int j = 0; getline(tempfile, line); j++)
    {
        if(!regex_search(line, match, headerMatch))
        {
            if(regex_search(line, match, mafMatch))
            {
                maf[mafPos] = stod(match[1]);
                mafPos++;
            }
        }
    }
    tempfile.close();
    
    
};
