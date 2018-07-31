//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"
#include <fstream>
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
    
    gMatch = regex("(\\.\\/\\.)|(\\d)(\\/|\\|)(\\d)");
    altAlleleCountMatch = regex("[ATGCN],[ATGCN]");
    //altAlleleCountMatch = regex("(?<=[ATGCN]),(?=[ACTGN])");
    mafMatch = regex(";AF=(0\\.*\\d*)");
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    headerMatch = regex("^#");
    posMatch = regex("^\\d{1,2}\\s\\d*");
    
    //Switching on test type to get proper input.
    if(testType == "wsbt")
    {
        readVcfInitialInfo(vcfFile1);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
        readGenotype(vcfFile1, genotypeGslMatrix);
    }
    else if(testType == "burden")
    {
        //Get basic info from VCF
        readVcfInitialInfo(vcfFile1);
        
        //Initilizing genotype data structures.
        genotypeMatrix = vector<vector<int> >(subjectCount, vector<int>(variantCount, 1));
        genotypeGslMatrix = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Initilizing maf data structure.
        maf = gsl_vector_alloc(variantCount);
        
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
    regex caseMatch("(\\t(\\d[^\\s]+))+");
    

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
                    subjectCount = stoi(match[1]);
                }
            }
            if(regex_search(line, match, variantCountMatch))
            {
                variantCount = stoi(match[1]);
            }
            if(variantCount != 0 && subjectCount !=0)
            {
                break;
            }
        }
        inputFile.close();
        
        cout << "SubjectCount is: " << subjectCount << endl;
        cout << "VariantCount is: " << variantCount << endl;
        
        inputFile.open(filename);
        for(int j = 0; getline(inputFile, line); j++)
        {
            if(caseCount == 0)
            {
                string templine;
                if(regex_search(line, match, caseMatch))
                {
                    templine = match[0];
                    for (int i = 0; regex_search(templine, match, regex("\\t(\\d[^\\s]+)")); i++)
                    {
                        caseCount++;
                        templine = match.suffix();
                    }
                }
            }
            else
            {
                break;
            }
        }
        inputFile.close();
        //Initilize the maf vector in case its used.
        maf = gsl_vector_alloc(variantCount);
        cout << "CaseCount is: " << caseCount << endl;
    }
}


void readInput::readGenotype(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;
        regex posMatch("POS\\t");
        regex refMatch("\\tHG\\d+");

        string genoCommand = "vcftools -" + vcfType + " "  + filename + " --extract-FORMAT-info GT";
        system(genoCommand.c_str());
        inputFile.open("out.GT.FORMAT");
        
        getline(inputFile, line);
        string templine = line;
        if (regex_search(line, match, posMatch))
        {
            line = match.suffix();
        }
        if (regex_search(line, match, posMatch))
        {
            templine = line.substr(0, match.suffix().length());
        }
        
        double progress = 0.0;
        int barWidth = 70;
        for(int j = 0; getline(inputFile, line); j++)
        {
            for(int i = 0; regex_search(line, match, gMatch); i++)
            {
                //Since the first line of the file is just headers, we input into j-1
                gsl_matrix_set(genotypeGslMatrix, j, i, stoi(match[2]) + stoi(match[4]));
                genotypeMatrix[i][j] = stoi(match[2]) + stoi(match[4]);
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
        regex_token_iterator<string::iterator> lineEnd;
        
        if(vcfType == "-gzvcf")
        {
            string genoCommand = "vcftools -" + vcfType + " "  + filename + " --extract-FORMAT-info GT";
            system(genoCommand.c_str());
            inputFile.open("out.GT.FORMAT");
        }
        else
        {
            inputFile.open(filename);
        }
        
        
        double progress = 0.0;
        int barWidth = 50;
        int mafPos = 0;
        int matrixInputLine = 0;
        
        for(int j = 0; getline(inputFile, line); j++)
        {
            //string templine = line;
            unsigned long int altAlleleCount = 1;
            if (regex_search(line, match, altAlleleCountMatch))
            {
                altAlleleCount = match.size() + 1;
            }
            if(!regex_search(line, match, headerMatch))
            {
                
                if(regex_search(line, match, mafMatch))
                {
                    gsl_vector_set(maf, mafPos, stod(match[1]));
                    mafPos++;
                }
                //Submatches:
                //1: ./.
                //2: left
                //4: right     of #|#.
                int submatches[] = {1,2,4};
                regex_token_iterator<string::iterator> genoParser(line.begin(), line.end(), gMatch, submatches);
                for(int i = 0; genoParser != lineEnd; i++)
                {
                    //We need to move the iterator (genoParser) forward 3 times to find a new match.
                    if(*genoParser == "./.")
                    {
                        gsl_matrix_set(inputMatrix, matrixInputLine, i, -1);
                        advance(genoParser, 3);
                    }
                    else if(*genoParser++ == "")
                    {
                        int left = 0;
                        int right = 0;
                        
                        left = stoi(*genoParser++);
                        if(left > 0)
                        {
                            left = 1;
                        }
                        right = stoi(*genoParser++);
                        if(right > 0)
                        {
                            right = 1;
                        }
                        gsl_matrix_set(inputMatrix, matrixInputLine, i, left + right);
                    }
                }
                matrixInputLine++;
            }
            
            //Will remove this loading bar stuff later.
            //It was copy and pasted from the internet and is just there for testing sanity.
            
            std::cout << "Loading Data: [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i)
            {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " % (" << matrixInputLine << "/" << variantCount << ") \r";
            std::cout.flush();
            progress = (matrixInputLine * 1.0) / variantCount;
             
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
        maf = gsl_vector_alloc(variantCount);
        
        inputFile.open(filename);
        int mafPos = 0;
        for(int j = 0; getline(inputFile, line); j++)
        {
            if(!regex_search(line, match, headerMatch))
            {
                if(regex_search(line, match, mafMatch))
                {
                    gsl_vector_set(maf, mafPos, stod(match[1]));
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








