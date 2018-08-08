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
#include <gsl/gsl_matrix.h>

using namespace std;

class readInput
{
    
public:
    readInput(string testType, string vcfType, string vcfFile1, string vcfFile2, string phenoFile);
    
    void readVcfInitialInfo(string filename);
    void readGenotype(string filename, gsl_matrix *inputMatrix);
    void readMaf(string filename);
    void makePositionFile(string filename);
    
    gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
    gsl_vector* getMaf(){return maf;}
    gsl_vector* getPheno(){return pheno;}
    int getCaseCount(){return caseCount;}

private:
    int variantCount;
    int subjectCount;
    int caseCount;
    string vcfType;
    string bcftools_location = "externals/bcftools/src/bcftools_project/bcftools";
    regex gMatch;
    regex altAlleleCountMatch;
    regex mafMatch;
    regex subjectCountMatch;
    regex variantCountMatch;
    regex headerMatch;
    regex posMatch;

    ifstream inputFile;
    
    
    gsl_matrix* genotypeGslMatrix;
    gsl_vector* maf;
    gsl_vector* pheno;
};

#endif /* input_hpp */
