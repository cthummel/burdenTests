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
    gMatch = regex("(\\.\\/\\.)|(\\d)(\\/|\\|)(\\d)");
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

        double progress = 0.0;
        int barWidth = 70;
        for(int j = 0; getline(inputFile, line); j++)
        {
            for(int i = 0; regex_search(line, match, gMatch); i++)
            {
                //Since the first line of the file is just headers, we input into j-1
                gsl_matrix_set(genotypeGslMatrix, j-1, i, stoi(match[1]) + stoi(match[3]));
                genotypeMatrix[i][j-1] = stoi(match[1]) + stoi(match[3]);
                line = match.suffix();
            }
            
            
            std::cout << "[";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
            
            progress += j/variantCount;
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
        
        double progress = 0.0;
        int barWidth = 50;
        
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
            
            //Will remove this loading bar stuff later.
            //It was copy and pasted from the internet and is just there for testing sanity.
            std::cout << "Loading Data: [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
            double top = j * 1.0;
            progress = top / variantCount;
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








