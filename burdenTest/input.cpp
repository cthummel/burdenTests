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
    caseCount = 0;
    vcfType = inputVcfType;
    
    //string line;
    string genoCommand = "vcftools -" + vcfType + " "  + vcfFile1 + " --extract-FORMAT-info GT";
    string mafCommand = "vcftools -" + vcfType + " "  + vcfFile1 + " --get-INFO AF";
    string subsetCommand = "vcftools -" + vcfType + " " + vcfFile1 + " --positions " + vcfFile2 + " --recode";
    //smatch match;
    
    //gMatch = regex("([0-9])(\\/|\\|)([0-9])");
    gMatch = regex("(.\\/.)|(\\d)(\\/|\\|)(\\d)");
    mafMatch = regex(";AF=(0\\.\\d*)");
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
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
        readVcfInitialInfo(vcfFile1);
        genotypeGslMatrix = gsl_matrix_calloc(variantCount, subjectCount);
        readGenotype(vcfFile1, genotypeGslMatrix);
        
        
        /*
        //Read in user data.
        readVcfInitialInfo(vcfFile1);
        gsl_matrix *affectedGenotype = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Subset background data to match user data and initialize background matrix.
        makePositionFile(vcfFile1);
        readVcfInitialInfo("out.recode.vcf");
        gsl_matrix *unaffectedGenotype = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Its probably going to be useful to only look through one chromosome at a time.
        //Im pretty sure this will cut down on the time to subset the background data.
        //For now I am assuming that all the background data is in a single file though.
        
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
         */
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
    regex caseMatch("HG\\d*");
    
    string statsFileName = filename.substr(0, filename.length() - vcfType.length());
    string summaryCommand = "bcftools stats " + filename + " > " + statsFileName + ".stats";
    system(summaryCommand.c_str());
    
    
    if(!inputFile.is_open())
    {
        inputFile.open(statsFileName + ".stats");
        for(int j = 0; getline(inputFile, line); j++)
        {
            if (subjectCount == 0)
            {
                if(regex_search(line, match, subjectCountMatch))
                {
                    cout << "subjectCount match 1: " << match[1] << endl;
                    subjectCount = stoi(match[1]);
                    //cout << "subjectCount: "  << subjectCount << endl;
                }
            }
            /*
            if(caseCount == 0)
            {
                while(regex_search(line, match, caseMatch))
                {
                    caseCount++;
                    line = match.suffix();
                    //cout << "caseCount: "  << caseCount << endl;
                }
            }
             */
            if(regex_search(line, match, variantCountMatch))
            {
                cout << "variantCount match 1: " << match[1] << endl;
                variantCount = stoi(match[1]);
                cout << "variantCount: "  << variantCount << endl;
            }
            if(variantCount != 0 && subjectCount !=0)
            {
                break;
            }
        }
        inputFile.close();
        
        inputFile.open(filename);
        for(int j = 0; getline(inputFile, line); j++)
        {
            if(caseCount == 0)
            {
                caseCount = 4;
                /*
                while(regex_search(line, match, caseMatch))
                {
                    caseCount++;
                    line = match.suffix();
                    //cout << "caseCount: "  << caseCount << endl;
                }
                 */
            }
            else
            {
                break;
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
                if(match[1] != "")
                {
                    gsl_matrix_set(inputMatrix, j-1, i, -1);
                }
                else
                {
                    gsl_matrix_set(inputMatrix, j-1, i, stoi(match[2]) + stoi(match[4]));
                }
                
                line = match.suffix();
            }
            //cout << "Line of GenoData done: " << j << endl;
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
    //string subsetCommand = "vcftools -" + vcfType + " " + backgroundVcf + " --positions pos.txt --recode";
    //system(subsetCommand.c_str());
    
}








