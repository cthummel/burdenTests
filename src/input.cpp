//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"
#include <fstream>
#include <sstream>

using namespace std;

//Sets up initial test paramters from user information.
readInput::readInput(string dir, string tType, string geneRegionIndicator, string userVcf, string region, string phenoFile, string covFile)
{
    caseCount = 0;
    variantRegion = region;
    testType = tType;
    codingRegionType = geneRegionIndicator;
    testDir = dir;
    
    //Switching on test type to get proper input.
    if (testType == "wsbt")
    {
        if (variantRegion == "--g")
        {
            readVcfInitialInfo(userVcf);
            readCaseCount(userVcf);
            readSampleNames(userVcf);
            buildPosMap(userVcf);
            buildGeneInfo("orderedRefFlat.txt");
            //buildGeneInfo("nameOrderedRefFlat.txt");
            if(codingRegionType == "transcript")
            {
                cout << "--- Running initial gene check in transcript mode ---" << endl;
                matchGenesOnTranscript();
            }
            else if (codingRegionType == "exon")
            {
                cout << "--- Running initial gene check in exon mode ---" << endl;
                matchGenesOnExons();
            }
            
        }
    }
    else if (testType == "skat")
    {
        readVcfInitialInfo(userVcf);
        buildPosMap(userVcf);
        readCaseCount(userVcf);
        readSampleNames(userVcf);
        readPhenotype(phenoFile);
        readCovariates(covFile);
        if (variantRegion == "--g")
        {
            buildGeneInfo("orderedrefFlat.txt");
            matchGenesOnTranscript();
        }
    }
    else if (testType == "skato")
    {
        readVcfInitialInfo(userVcf);
        buildPosMap(userVcf);
        readCaseCount(userVcf);
        readSampleNames(userVcf);
        readPhenotype(phenoFile);
        readCovariates(covFile);
        if (variantRegion == "--g")
        {
            buildGeneInfo("orderedrefFlat.txt");
            matchGenesOnTranscript();
        }
    }
}


readInput::~readInput()
{
    //As a readInput object goes out of scope it should delete the gsl objects.
    if(covariates != nullptr)
    {
        gsl_matrix_free(covariates);
    }
    if(pheno != nullptr)
    {
        gsl_vector_free(pheno);
    }
    
}

void readInput::readVcfInitialInfo(string filename, string region, int thread_ID)
{
    subjectCount = -1;
    variantCount = -1;
    string line;
    string statsFile = "tmp/data" + to_string(thread_ID) + ".stats";
    ifstream in;
    smatch match;
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    
    //Build our stats file and GT file for later parsing.
    if(preMerged)
    {
        string command = externals_loc + "bcftools stats -r " + region + " " + filename + " > " + statsFile;
        system(command.c_str());
    }
    else
    {
        string command = externals_loc + "bcftools merge -r " + region + " " + userFile + " " + backFile +  " | " 
                       + externals_loc + "bcftools stats - > " + statsFile;
        system(command.c_str());
    }
    

    //Parsing background info
    in.open(statsFile);
    for (int j = 0; getline(in, line); j++)
    {
        if (subjectCount < 0)
        {
            if (regex_search(line, match, subjectCountMatch))
            {
                subjectCount = stoi(match[1]);
            }
        }
        if (regex_search(line, match, variantCountMatch))
        {
            variantCount = stoi(match[1]);
        }
        if (variantCount > -1 && subjectCount > -1)
        {
            break;
        }
    }
    in.close();

    if(variantCount < 1 || subjectCount < 1)
    {
        if(variantCount < 1 && subjectCount > 0)
        {
            //throw std::invalid_argument("Unable to find any variants in region " + region + ". ----Skipping----");
        }
        else if (subjectCount < 1 && variantCount > 0)
        {
            //throw std::invalid_argument("Unable to find any subjects in region " + region + ". ----Skipping----");
        }
        else
        {
            //throw std::invalid_argument("Unable to find any subjects or variants in region " + region + ". ----Skipping----");
        }
    }
    
}


void readInput::readVcfInitialInfo(string filename)
{
    string line;
    smatch match;
    ifstream in;
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    subjectCount = 0;
    variantCount = 0;

    string summaryCommand = externals_loc + "bcftools stats " + filename + " > tmp/user.stats";
    system(summaryCommand.c_str());
    
    if(!in.is_open())
    {
        in.open("tmp/user.stats");
        for(int j = 0; getline(in, line); j++)
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
        in.close();

	    cout << "Subject count: " << subjectCount << endl;
	    cout << "Variant count: " << variantCount << endl;
    }
}

//Needs fixing.
void readInput::readCaseCount(string filename)
{
    if(subjectCount > 2504)
    {
        caseCount = subjectCount - 2504;
    }
    else
    {
        caseCount = subjectCount;
    }
}

//If user has phenotype data we save it in a vector otherwise we assume binary phenotype.
void readInput::readPhenotype(string phenoFile)
{
    ifstream in;
    if(subjectCount < 2504)
    {   
        pheno = gsl_vector_calloc(2504 + caseCount);
    }
    else
    {
        pheno = gsl_vector_calloc(subjectCount);
    }
    
    if(phenoFile.compare("") == 0)
    {
        for(int i = 0; i < caseCount; i++)
        {
            gsl_vector_set(pheno, i, 1);
        }
    }
    else
    {
        //Parse phenotype and add to vector;
        
        in.open(phenoFile);
        string line;
        //Parse the header
        getline(in, line);
        vector<string> fid;
        vector<string> iid;
        vector<string> fatid;
        vector<string> matid;
        vector<int> sex;
        
        for(int i = 0; getline(in, line); i++)
        {
            size_t last = 0;
            size_t current = line.find('\t', last);
            for (int j = 0; current != string::npos; j++)
            {
                string token = line.substr(last, current - last);
                if (j == 0)
                {
                    fid.push_back(token);
                }
                if (j == 1)
                {
                    iid.push_back(token);
                }
                if (j == 2)
                {
                    fatid.push_back(token);
                }
                if (j == 3)
                {
                    matid.push_back(token);
                }
                if (j == 4)
                {
                    sex.push_back(stoi(token));
                }
                if (j == 5)
                {
                    gsl_vector_set(pheno, i, stoi(token) - 1);
                }
                
                last = current + 1;
                current = line.find('\t', last);
            }
        }
        in.close();
    }
}

//Need a good method for reading covariates
void readInput::readCovariates(string filename)
{
    ifstream in;
    if(filename.compare("") == 0)
    {
        covariates = gsl_matrix_calloc(subjectCount + 2504, 1);
        gsl_matrix_set_all(covariates, 1);
    }
}


//Builds an ordered (by chromosome then starting transcription region) list of genes.
void readInput::buildGeneInfo(string filename)
{
    string line;
    ifstream in;
    in.open(filename);
    for (int i = 0; getline(in, line); i++)
    {
        size_t last = 0;
        size_t current = line.find('\t', last);
        geneId currentGene;
        for (int j = 0; last != string::npos; j++)
        {
            
            string token = line.substr(last, current - last);
            //cout << i << "," << j << endl;
            if (j == 0)
            {
                currentGene.geneName = token;
            }
            if (j == 1)
            {
                currentGene.transcriptName = token;
            }
            if (j == 2)
            {
                //If coded as chr#
                if (token.length() > 3)
                {
                    currentGene.geneChrom = token.substr(3);
                }
                //If coded as just the number.
                else
                {
                    currentGene.geneChrom = token;
                }
            }
            if (j == 4)
            {
                currentGene.txStartPos = stoi(token);
            }
            if (j == 5)
            {
                currentGene.txEndPos = stoi(token);
            }
            if (j == 6)
            {
                currentGene.codingStartPos = stoi(token);
            }
            if (j == 7)
            {
                currentGene.codingEndPos = stoi(token);
            }
            //Exon Count
            if (j == 8 && codingRegionType == "exon")
            {
                currentGene.exonCount = stoi(token);
                currentGene.exonStarts = vector<int>(currentGene.exonCount);
                currentGene.exonEnds = vector<int>(currentGene.exonCount);
            }
            //Exon Start Position
            if(j == 9 && codingRegionType == "exon")
            {
                size_t exonStart = 0;
                for(int k = 0; k < currentGene.exonCount; k++)
                {
                    size_t exonEnd = token.find(',', exonStart);
                    string exon = token.substr(exonStart, exonEnd - exonStart);
                    currentGene.exonStarts[k] = stoi(exon);
                    exonStart = exonEnd + 1;
                }
            }
            //Exon end position
            if(j == 10 && codingRegionType == "exon")
            {
                size_t exonStart = 0;
                for(int k = 0; k < currentGene.exonCount; k++)
                {
                    size_t exonEnd = token.find(',', exonStart);
                    
                    currentGene.exonEnds[k] = stoi(token.substr(exonStart, exonEnd - exonStart));
                    exonStart = exonEnd + 1;
                }
            }
            if(current == string::npos)
            {
                //Breaks out of parsing this entry.
                last = string::npos;
            }
            else
            {
                last = current + 1;
                current = line.find('\t', last);
            }
        }
        //Setup the region for the current gene
        if (codingRegionType == "exon")
        {
            stringstream temp;
            for (int k = 0; k < currentGene.exonCount; k++)
            {
                temp << currentGene.geneChrom << ":" << currentGene.exonStarts[k] << "-" << currentGene.exonEnds[k];
                if (k + 1 < currentGene.exonCount)
                {
                    temp << ",";
                }
            }
            currentGene.region = temp.str();
        }
        else
        {
            currentGene.region = currentGene.geneChrom + ":" + to_string(currentGene.txStartPos) + "-" + to_string(currentGene.txEndPos);
        }
        info.push_back(currentGene);
    }
    in.close();
}

//Builds a map of variant positions to chromosome in user data.
void readInput::buildPosMap(string filename)
{
    ifstream in;
    string line;
    vector<int> vcfPos;
    string posFile = "tmp/user.posdata.txt";
    string command = externals_loc + "bcftools query -f '%CHROM,%POS\\n' " + filename + " > " + posFile;
    system(command.c_str());

    string currentChrom;
    in.open(posFile);
    for(int i = 0; getline(in, line); i++)
    {
        int comma = line.find_first_of(',');
        if (i == 0)
        {
            currentChrom = line.substr(0, comma);
        }
        if (currentChrom != line.substr(0, comma))
        {
            posMap.insert(pair<string, vector<int> >(currentChrom, vcfPos));
            cout << vcfPos.size() << " variants in chromosome " << currentChrom << endl;
            vcfPos.clear();
            currentChrom = line.substr(0, comma);
        }

        vcfPos.push_back(stoi(line.substr(comma + 1)));
    }
    posMap.insert(pair<string, vector<int> >(currentChrom, vcfPos));
    cout << vcfPos.size() << " variants in chromosome " << currentChrom << endl;

    in.close();
}


void readInput::matchGenesOnTranscript()
{
    cout << endl << "Looking through " << info.size() << " gene entries." << endl;
    for (int i = 0; i < info.size(); i++)
    {
        bool geneFound = false;
        //Check if the chromosome of current gene is even in the data set.
        if (posMap.find(info[i].geneChrom) == posMap.end())
        {
            continue;
        }

        //Iterate through all the variants in a chromosome to find gene matches.
        for (int j = 0; j < posMap[info[i].geneChrom].size(); j++)
        {
            //Consider not searching the early positions multiple times. Could initiate j to the last seen useful position.
            if (posMap[info[i].geneChrom][j] >= info[i].txStartPos && posMap[info[i].geneChrom][j] <= info[i].txEndPos)
            {
                //Try and enter the gene into the list of present genes.
                pair<map<string, string>::iterator, bool> duplicate;
                duplicate = regions.insert(pair<string, string>(info[i].geneName, info[i].region));
                //This gene was a duplicate of another already entered in. Take the largest transcript region.
                if(duplicate.second == false)
                {
                    string oldRegion = regions[info[i].geneName];
                    int oldStart = genePosMap[info[i].geneName].first;
                    int oldEnd = genePosMap[info[i].geneName].second;
                    int newWidth = info[i].txEndPos - info[i].txStartPos;
                    if(newWidth > (oldEnd - oldStart))
                    {
                        regions[info[i].geneName] = info[i].region;
                        genePosMap[info[i].geneName] = pair<int, int>(info[i].txStartPos, info[i].txEndPos);
                    }
                }
                else
                {
                    genePosMap.insert(pair<string, pair<int, int>>(info[i].geneName, pair<int,int>(info[i].txStartPos,info[i].txEndPos)));
                }
                geneFound = true;
            }
            if (posMap[info[i].geneChrom][j] > info[i].txEndPos || geneFound)
            {
                break;
            }
        }
    }
    cout << "Matched variant(s) to " << regions.size() << " unique genes." << endl;
    for(map<string,string>::iterator it = regions.begin(); it != regions.end(); it++)
    {
        cout << "Matched variant(s) to gene " << it->first << " in region " << it->second << endl;
    }
    
    //Dont need the posMap anymore.
    posMap.clear();
}

void readInput::matchGenesOnExons()
{
    cout << endl << "Looking through " << info.size() << " gene entries." << endl;
    for (int i = 0; i < info.size(); i++)
    {
        bool geneFound = false;
        //Check if the chromosome of current gene is even in the data set.
        if (posMap.find(info[i].geneChrom) == posMap.end())
        {
            continue;
        }

        //Iterate through all the variants in a chromosome to find gene matches.
        for (int j = 0; j < posMap[info[i].geneChrom].size(); j++)
        {
            for(int k = 0; k < info[i].exonCount; k++)
            {
                //Consider not searching the early positions multiple times. Could initiate j to the last seen useful position.
                if (posMap[info[i].geneChrom][j] >= info[i].exonStarts[k] && posMap[info[i].geneChrom][j] <= info[i].exonEnds[k])
                {
                    //Try and enter the gene into the list of present genes.
                    pair<map<string, string>::iterator, bool> duplicate;
                    duplicate = regions.insert(pair<string, string>(info[i].geneName, info[i].region));
                    //This gene was a duplicate of another already entered in. Take the largest transcript region.
                    if (duplicate.second == false)
                    {
                        string oldRegion = regions[info[i].geneName];
                        int oldStart = genePosMap[info[i].geneName].first;
                        int oldEnd = genePosMap[info[i].geneName].second;
                        int newWidth = info[i].txEndPos - info[i].txStartPos;
                        if (newWidth > (oldEnd - oldStart))
                        {
                            regions[info[i].geneName] = info[i].region;
                            genePosMap[info[i].geneName] = pair<int, int>(info[i].txStartPos, info[i].txEndPos);
                        }
                    }
                    else
                    {
                        genePosMap.insert(pair<string, pair<int, int>>(info[i].geneName, pair<int, int>(info[i].txStartPos, info[i].txEndPos)));
                    }
                    geneFound = true;
                    break;
                }
                
                if (posMap[info[i].geneChrom][j] > info[i].txEndPos)
                {
                    break;
                }
                
            }
            if(geneFound)
            {
                break;
            }
        }
    }
}

//Requires caseCount to be a known value.
void readInput::readSampleNames(string filename)
{
    string outputFile = "tmp/sampleNames.txt";

    string command = externals_loc + "bcftools query -l " + filename + " > " + outputFile;
    system(command.c_str());

    ifstream in(outputFile);
    string line;
    for(int i = 0; i < caseCount; i++)
    {
        getline(in,line);
        sampleNames.push_back(line);
    }
    in.close();
}





