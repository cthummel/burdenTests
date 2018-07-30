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
    void readGenotype(string filename);
    void readGenotype(string filename, gsl_matrix *inputMatrix);
    void readMaf(string filename);
    void makePositionFile(string filename);
    
    vector<vector<int> > getGenotype(){return genotypeMatrix;}
    gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
    gsl_vector* getMaf(){return maf;}
    int getCaseCount(){return caseCount;}

private:
    int variantCount;
    int subjectCount;
    unsigned long int caseCount;
    string vcfType;
    //const string backgroundVcf = "gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz";
    regex gMatch;
    regex altAlleleCountMatch;
    regex mafMatch;
    regex subjectCountMatch;
    regex variantCountMatch;
    regex headerMatch;
    regex posMatch;

    ifstream inputFile;
    
    
    gsl_matrix *genotypeGslMatrix;
    vector<vector<int> > genotypeMatrix;
    gsl_vector *maf;
};

#endif /* input_hpp */
