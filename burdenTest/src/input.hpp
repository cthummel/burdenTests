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

using namespace std;

class readInput
{
    struct geneId
    {
        string geneName;        //Column 1
        string transcriptName;  //Column 2
        string geneChrom;       //Column 3
        bool strand;            //Column 4
        int txStartPos;         //Column 5
        int txEndPos;           //Column 6
        int codingStartPos;     //Column 7
        int codingEndPos;       //Column 8
        int exonCount;          //Column 9
        gsl_vector_int *exonStarts; //Column 10
        gsl_vector_int *exonEnds;   //Column 11
        gsl_matrix *genoData;
        gsl_vector *maf;
    };

  public:
    readInput(string dir, string testType, string vcfType, string vcfFile1, string region, string phenoFile, string covFile);
    readInput(bool mergeData, string userVcf, string backFilename, string region, int cases, int thread);
    
    void readVcfInitialInfo(string filename, string region, int thread_ID);
    void readVcfInitialInfo(string filename);
    void readCaseCount(string filename);
    void readPhenotype(string phenoFile);
    void readMaf(string filename);
    void readMaf(string filename, string region, string outfile);
    void makePositionFile(string filename);
    void mergeData(string user, string background, string region, string outfile);
    void bcfInput(string filename);
    void bcfInput(string filename, string back, string region, string outfile);
    void readGenes(string filename);
    
    
    int getCaseCount(){return caseCount;}
    gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
    gsl_matrix* getCovariates(){return covariates;}
    gsl_vector* getMaf(){return maf;}
    gsl_vector* getMaf(string geneName){return geneMaf[geneName];}
    gsl_vector* getPheno(){return pheno;}
    map<string, gsl_matrix*> getGeneSubsets(){return genes;}
    map<string, string> getRegions(){return region;}
    

private:
    void parseGenes(string chromosome, vector<string> *geneName);
    void buildGeneInfo(string filename);
    void loadGene(geneId gene, vector<pair<int, gsl_vector*> > userData);
    void buildPosMap(string filename);
    void matchGenes();
    void variantMatchGene();
    bool testReadFromStream(string filename, string region);
    void test(string filename);
    string exec(const char* cmd);
    int variantCount;
    int subjectCount;
    int caseCount;
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
    map<string, gsl_matrix *> genes;
    map<string, gsl_vector *> geneMaf;
    map<string, string> region;
    map<string, pair<int,int>> genePosMap;

    gsl_matrix* genotypeGslMatrix;
    gsl_matrix* covariates;
    gsl_vector* maf;
    gsl_vector* pheno;
};

#endif /* input_hpp */
