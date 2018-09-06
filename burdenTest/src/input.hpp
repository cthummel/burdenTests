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
    readInput(string testType, string vcfType, string vcfFile1, string region, string phenoFile, string covFile);
    
    void readVcfInitialInfo(string filename);
    void readGenotype(string filename, gsl_matrix *inputMatrix, int subjectOffset);
    void readPhenotype(string phenoFile);
    void readMaf(string filename);
    void makePositionFile(string filename);
    void mergeData(string filename);
    void bcfInput(string filename);
    
    int getCaseCount(){return caseCount;}
    gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
    gsl_matrix* getCovariates(){return covariates;}
    gsl_vector* getMaf(){return maf;}
    gsl_vector* getPheno(){return pheno;}
    

private:
    int variantCount;
    int subjectCount;
    int caseCount;
    string vcfType;
    string bcftools_loc = "../externals/bin/bcftools";
    string bgzip_loc = "../externals/bin/bgzip";
    string externals_loc = "../externals/bin/";
    regex gMatch;
    regex altAlleleCountMatch;
    regex mafMatch;
    regex subjectCountMatch;
    regex variantCountMatch;
    regex headerMatch;
    regex posMatch;

    ifstream inputFile;
    
    
    gsl_matrix* genotypeGslMatrix;
    gsl_matrix* covariates;
    gsl_vector* maf;
    gsl_vector* pheno;
};

#endif /* input_hpp */
