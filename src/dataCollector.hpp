#ifndef dataCollector_hpp
#define dataCollector_hpp

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <gsl/gsl_matrix.h>
#include "geneInfo.cpp"

using namespace std;

class dataCollector
{
    public:
      dataCollector(bool mergeData, string userVcf, string backFilename, string region, string test_type, int caseCount, int thread);
      ~dataCollector();

      gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
      gsl_vector* getMaf(){return maf;}
      int getCaseUniqueVariantCount(){return caseUniqueVariantCount;}
      int getBackgroundUniqueVariantCount(){return backgroundUniqueVariantCount;}

    private:
      void readVcfInitialInfo(string filename, string region, string outfile);
      void readMaf(string filename, string region, string outfile);
      void bcfInput(string filename, string back, string region, string outfile);

    int variantCount = 0;
    int subjectCount = 0;
    int caseCount = 0;
    int caseUniqueVariantCount = 0;
    int backgroundUniqueVariantCount = 0;
    bool preMerged = false;
    string userFile = "";
    string backFile = "";
    string variantRegion = "";
    string testType = "";
    string externals_loc = "../../externals/bin/";
    regex subjectCountMatch;
    regex variantCountMatch;

    gsl_matrix* genotypeGslMatrix = nullptr;
    gsl_vector* maf = nullptr;
};





#endif