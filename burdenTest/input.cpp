//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"
using namespace std;

readInput::readInput(string testType, string inputVcfType, string vcfFile1, string vcfFile2, string phenoFile)
{
    variantCount = 0;
    subjectCount = 0;
    vcfType = inputVcfType;
    
    //string line;
    string genoCommand = "vcftools -" + vcfType + " "  + vcfFile1 + " --extract-FORMAT-info GT";
    string mafCommand = "vcftools -" + vcfType + " "  + vcfFile1 + " --get-INFO AF";
    string subsetCommand = "vcftools -" + vcfType + " " + vcfFile1 + " --positions " + vcfFile2 + " --recode";
    //smatch match;
    
    gMatch = regex("([0-9])(\\/|\\|)([0-9])");
    mafMatch = regex(";AF=(0\\.\\d*)");
    subjectCountMatch = regex(";NS=(\\d*)");
    headerMatch = regex("^#");
    posMatch = regex("^\\d{1,2}\\s\\d*");
    
    
    
    /*
    //Obtain variant counts and subject counts.
    readVcfInitialInfo(vcfFile1);
    
    cout << "variant count: " << variantCount <<endl;
    cout << "subject count: " << subjectCount <<endl;
    
    //Initilizing genotype data structures.
    genotypeMatrix = vector<vector<int> >(subjectCount, vector<int>(variantCount, 1));
    genotypeGslMatrix = gsl_matrix_calloc(variantCount, subjectCount);
    cout << "gsl_matrix rows: " << genotypeGslMatrix->size1 << endl;
    cout << "gsl_matrix columns: " << genotypeGslMatrix->size2 << endl;
    cout << "gsl_matrix max columns: " << genotypeGslMatrix->tda << endl;
    
    //Initilizing maf data structure.
    maf = vector<double>(variantCount);
    */
    
    //Switching on test type to get proper input.
    if(testType == "wsbt")
    {
        //Read in user data.
        readVcfInitialInfo(vcfFile1);
        gsl_matrix *affectedGenotype = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Subset background data to match user data and initialize background matrix.
        makePositionFile(vcfFile1);
        readVcfInitialInfo("out.recode.vcf");
        gsl_matrix *unaffectedGenotype = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Initialize combined matrix and read in data.
        genotypeGslMatrix = gsl_matrix_calloc(variantCount, affectedGenotype->tda + unaffectedGenotype->tda);
        readGenotype(vcfFile1, affectedGenotype);
        readGenotype("out.recode.vcf", unaffectedGenotype);
        
        //Combining the two matricies with affected in front and unaffected after.
        for(int i = 0; i < affectedGenotype->size1; i++)
        {
            for(int j = 0; j < affectedGenotype->size2; j++)
            {
                gsl_matrix_set(genotypeGslMatrix, i, j, gsl_matrix_get(affectedGenotype, i, j));
            }
        }
        
        unsigned long int affectedSubjectCount = affectedGenotype->size2;
        
        for(int i = 0; i < unaffectedGenotype->size1; i++)
        {
            for(int j = 0; j < unaffectedGenotype->size2; j++)
            {
                gsl_matrix_set(genotypeGslMatrix, i, affectedSubjectCount + j, gsl_matrix_get(unaffectedGenotype, i, j));
            }
        }
    }
    else if(testType == "burden")
    {
        //Get basic info from VCF
        readVcfInitialInfo(vcfFile1);
        
        cout << "variant count: " << variantCount <<endl;
        cout << "subject count: " << subjectCount <<endl;
        
        //Initilizing genotype data structures.
        genotypeMatrix = vector<vector<int> >(subjectCount, vector<int>(variantCount, 1));
        genotypeGslMatrix = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Initilizing maf data structure.
        maf = vector<double>(variantCount);
        
        //Runs vcftools command to separate genotype data.
        //system(genoCommand.c_str());
        readGenotype(vcfFile1);
        
        //Parse the maf
        readMaf(vcfFile1);
    }
    else if (testType == "cast")
    {
        
    }
    else if (testType == "skat")
    {
        
    }
    else if (testType == "skato")
    {
        
    }
    
}

void readInput::readVcfInitialInfo(string filename)
{
    string line;
    smatch match;
    //regex subjectCountMatch(";NS=(\\d*)");
    //regex headerMatch("^#");
    
    if(!inputFile.is_open())
    {
        inputFile.open(filename);
        for(int j = 0; getline(inputFile, line); j++)
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
        inputFile.close();
    }
}


void readInput::readGenotype(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;

        string genoCommand = "vcftools -" + vcfType + " "  + filename + " --extract-FORMAT-info GT";
        system(genoCommand.c_str());
        inputFile.open("out.GT.FORMAT");

        for(int j = 0; getline(inputFile, line); j++)
        {
            for(int i = 0; regex_search(line, match, gMatch); i++)
            {
                //Since the first line of the file is just headers, we input into j-1
                gsl_matrix_set(genotypeGslMatrix, j-1, i, stoi(match[1]) + stoi(match[3]));
                genotypeMatrix[i][j-1] = stoi(match[1]) + stoi(match[3]);
                line = match.suffix();
            }
        }
        inputFile.close();
    }
}

void readInput::readGenotype(string filename, gsl_matrix *inputMatrix)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;
        
        string genoCommand = "vcftools -" + vcfType + " "  + filename + " --extract-FORMAT-info GT";
        system(genoCommand.c_str());
        inputFile.open("out.GT.FORMAT");
        
        for(int j = 0; getline(inputFile, line); j++)
        {
            for(int i = 0; regex_search(line, match, gMatch); i++)
            {
                //Since the first line of the file is just headers, we input into j-1
                gsl_matrix_set(inputMatrix, j-1, i, stoi(match[1]) + stoi(match[3]));
                line = match.suffix();
            }
        }
        inputFile.close();
    }
}

void readInput::readMaf(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;
        //regex mafMatch(";AF=(0\\.\\d*)");
        //regex headerMatch("^#");
        
        
        inputFile.open(filename);
        int mafPos = 0;
        for(int j = 0; getline(inputFile, line); j++)
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
        inputFile.close();
    }
}

void readInput::makePositionFile(string filename)
{
    string line;
    smatch match;
    
    if (!inputFile.is_open())
    {
        ofstream outFile("pos.txt", ofstream::out);
        inputFile.open(filename);
        for(int j = 0; getline(inputFile, line); j++)
        {
            if(!regex_search(line, match, headerMatch))
            {
                if(regex_search(line, match, posMatch))
                {
                    outFile << line << "\n";
                }
            }
        }
        
        outFile.close();
        inputFile.close();
    }
    //filename here should be the background file we use for unaffected.
    string subsetCommand = "vcftools -" + vcfType + " " + backgroundVcf + " --positions pos.txt --recode";
    system(subsetCommand.c_str());
    
}








