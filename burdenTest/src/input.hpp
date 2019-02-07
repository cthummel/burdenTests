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
#include <map>
#include <gsl/gsl_matrix.h>
#include "geneInfo.cpp"

using namespace std;

class readInput
{
  public:
    readInput(string dir, string testType, string vcfType, string vcfFile1, string region, string phenoFile, string covFile);
    ~readInput();
    
    int getCaseCount(){return caseCount;}
    gsl_matrix* getCovariates(){return covariates;}
    gsl_vector* getPheno(){return pheno;}
    map<string, string> getRegions(){return regions;}
    vector<geneId> * getGenes(){return &info;}
    

private:
    void readVcfInitialInfo(string filename, string region, int thread_ID);
    void readVcfInitialInfo(string filename);
    void readCaseCount(string filename);
    void readPhenotype(string phenoFile);
    void readCovariates(string covFile);
    void makePositionFile(string filename);
    void readGenes(string filename);
    void buildGeneInfo(string filename);
    void buildPosMap(string filename);
    void matchGenes();

    int variantCount = 0;
    int subjectCount = 0;
    int caseCount = 0;
    bool preMerged;
    string userFile;
    string backFile;
    string variantRegion;
    string testType;
    string vcfType;
    string externals_loc = "../../externals/bin/";
    string testDir;
    regex subjectCountMatch;
    regex variantCountMatch;
    
    vector<geneId> info;
    map<string, vector<int> > posMap;
    map<string, string> regions;
    map<string, pair<int,int>> genePosMap;

    gsl_matrix* covariates = nullptr;
    gsl_vector* pheno = nullptr;
};

#endif /* input_hpp */
