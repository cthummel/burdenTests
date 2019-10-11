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
      dataCollector(bool mergeData, bool CADD, string userVcf, string backFilename, string region, string test_type, int caseCount, int thread);
      ~dataCollector();

      gsl_matrix_short* getShortGslGenotype(){return shortGenotypeGslMatrix;}
      gsl_matrix* getGslGenotype(){return genotypeGslMatrix;}
      gsl_vector* getMaf(){return maf;}
      gsl_vector* getCADDWeights(){return CADDWeights;}
      int getCaseUniqueVariantCount(){return caseUniqueVariantCount;}
      int getBackgroundUniqueVariantCount(){return backgroundUniqueVariantCount;}
      double getAverageImpact(){return averageImpact;}

    private:
      void readVcfInitialInfo(string filename, string region, string outfile);
      void readMaf(string filename, string region, string outfile);
      void readCADD(string filename, string back, string region, string rawfile, string annoFile, string outfile);
      void shortBcfInput(string filename, string back, string region, string outfile);
      void doubleBcfInput(string filename, string back, string region, string outfile);
      void annotationParser(string filename, string back, string region, string outfile);
      void weightImport(string region, string outfile);
      void setAverageImpact(vector<int> userUniqueVariants, int uniqueCount, string region, string outfile);

    int variantCount = 0;
    int subjectCount = 0;
    int caseCount = 0;
    int caseUniqueVariantCount = 0;
    int backgroundUniqueVariantCount = 0;
    double averageImpact = 0;
    bool preMerged = false;
    string userFile = "";
    string backFile = "";
    string variantRegion = "";
    string testType = "";
    string externals_loc = "../../externals/bin/";
    regex subjectCountMatch;
    regex variantCountMatch;
    vector<int> userUniques;

    gsl_matrix_short* shortGenotypeGslMatrix = nullptr;
    gsl_matrix* genotypeGslMatrix = nullptr;
    gsl_vector* CADDWeights = nullptr;
    gsl_vector* maf = nullptr;
};





#endif