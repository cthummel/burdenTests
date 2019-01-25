#include <vector>
#include <string>
#include <gsl/gsl_vector.h>

#ifndef geneInfo_cpp
#define geneInfo_cpp

struct geneId
{
    std::string geneName;       //Column 1
    std::string transcriptName; //Column 2
    std::string geneChrom;      //Column 3
    bool strand;                //Column 4
    int txStartPos;             //Column 5
    int txEndPos;               //Column 6
    int codingStartPos;         //Column 7
    int codingEndPos;           //Column 8
    int exonCount;              //Column 9
    gsl_vector_int *exonStarts; //Column 10
    gsl_vector_int *exonEnds;   //Column 11
    std::string region;
};

#endif


