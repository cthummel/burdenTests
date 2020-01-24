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

    //Currently we are ignoring X and Y chrom.
    skipXChrom = true;
    skipYChrom = true;
    
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
            ////buildGeneInfo("nameOrderedRefFlat.txt");
            //readAllGenes("orderedRefFlat.txt");
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

        if(skipXChrom)
        {
            //We are skipping X chromosomes
            if(currentGene.geneChrom == "X")
            {
                continue;
            }
        }
        if(skipYChrom)
        {
            //We are skipping Y chromosomes.
            if(currentGene.geneChrom == "Y")
            {
                continue;
            }
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

        //Honestly this take forever on 87 Million variants because we have to search 40Mill+ each time.
        //Lets try a better search method.
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
                    cout << "Matched variant(s) to gene " << info[i].geneName << " in region " << info[i].geneChrom << ":" << info[i].txStartPos << "-" << info[i].txEndPos << endl;
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
        //cout << "Matched variant(s) to gene " << it->first << " in region " << it->second << endl;
    }
    
    //Dont need the posMap anymore.
    posMap.clear();
}


//Not working yet. Need to finalize the skip forward step.
void readInput::matchGenesOnTranscriptv2()
{
    //cout << endl << "Looking through " << info.size() << " gene entries." << endl;
    
    for (std::map<std::string, vector<int>>::iterator it = posMap.begin(); it != posMap.end(); ++it)
    {
        vector<int> *currentPositionsInChromosome = &it->second;
        for(int i = 0; i < currentPositionsInChromosome->size(); ++i)
        {
            int positionsToMoveRight = 0;
            for(int j = 0; j < info.size(); ++j)
            {
                if(info[j].geneChrom == it->first)
                {
                    if(currentPositionsInChromosome->at(i) >= info[j].txStartPos && currentPositionsInChromosome->at(i) <= info[j].txEndPos)
                    {
                        //Try and enter the gene into the list of present genes.
                        pair<map<string, string>::iterator, bool> duplicate;
                        duplicate = regions.insert(pair<string, string>(info[j].geneName, info[j].region));
                        //This gene was a duplicate of another already entered in. Take the largest transcript region.
                        if (duplicate.second == false)
                        {
                            string oldRegion = regions[info[j].geneName];
                            int oldStart = genePosMap[info[j].geneName].first;
                            int oldEnd = genePosMap[info[j].geneName].second;
                            int newWidth = info[j].txEndPos - info[j].txStartPos;
                            if (newWidth > (oldEnd - oldStart))
                            {
                                regions[info[j].geneName] = info[j].region;
                                genePosMap[info[j].geneName] = pair<int, int>(info[j].txStartPos, info[j].txEndPos);
                            }
                            if(newWidth > positionsToMoveRight)
                            {
                                positionsToMoveRight = newWidth;
                            }
                        }
                        else
                        {
                            cout << "Matched variant(s) to gene " << info[j].geneName << " in region " << info[j].geneChrom << ":" << info[j].txStartPos << "-" << info[j].txEndPos << endl;
                            if (info[j].txEndPos - info[j].txStartPos > positionsToMoveRight)
                            {
                                positionsToMoveRight = info[j].txEndPos - info[j].txStartPos;
                            }
                            genePosMap.insert(pair<string, pair<int, int>>(info[j].geneName, pair<int, int>(info[j].txStartPos, info[j].txEndPos)));
                        }
                    }

                    //Increment the position count until we have a position not in the current gene.
                }
            }
        }
    }




    for (int i = 0; i < info.size(); i++)
    {
        bool geneFound = false;
        //Check if the chromosome of current gene is even in the data set.
        if (posMap.find(info[i].geneChrom) == posMap.end())
        {
            continue;
        }

        //Honestly this take forever on 87 Million variants because we have to search 40Mill+ each time.
        //Lets try a better search method.
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
                    cout << "Matched variant(s) to gene " << info[i].geneName << " in region " << info[i].geneChrom << ":" << info[i].txStartPos << "-" << info[i].txEndPos << endl;
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
        //cout << "Matched variant(s) to gene " << it->first << " in region " << it->second << endl;
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
                    //This gene was a duplicate of another already entered in. Merge exon regions to avoid variant double counts.
                    if (duplicate.second == false)
                    {
                        cout << "Duplicate record for " << info[i].geneName << ". Updating entry with combined range." << endl;
                        vector<int> oldStarts = geneExonPosMap[info[i].geneName].first;
                        vector<int> oldEnds = geneExonPosMap[info[i].geneName].second;
                        geneExonPosMap[info[i].geneName] = mergeExons(oldStarts, oldEnds, info[i].exonStarts, info[i].exonEnds);
                        stringstream newRegion;
                        // cout << "Matched Genes in " + info[i].geneName + " combined region ";
                        // for(int m = 0; m < geneExonPosMap[info[i].geneName].first.size(); m++)
                        // {
                        //     newRegion << geneExonPosMap[info[i].geneName].first[m] << "-" << geneExonPosMap[info[i].geneName].second[m] << ",";
                        // }
                        // cout << newRegion.str() << endl;
                        regions[info[i].geneName] = newRegion.str();
                    }
                    else
                    {
                        geneExonPosMap.insert(pair<string, pair<vector<int>, vector<int>>>(info[i].geneName, pair<vector<int>, vector<int>>(info[i].exonStarts, info[i].exonEnds)));
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

//Read all genes from base file and allow later checks to only run genes with a non-zero number of user variants.
void readInput::readAllGenes(string filename)
{
    cout << "Reading all genes for testing from " << filename << endl;

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

        if(skipXChrom)
        {
            //We are skipping X chromosomes
            if(currentGene.geneChrom == "X")
            {
                continue;
            }
        }
        if(skipYChrom)
        {
            //We are skipping Y chromosomes.
            if(currentGene.geneChrom == "Y")
            {
                continue;
            }
        }
        //Try and enter the gene into the list of present genes.
        pair<map<string, string>::iterator, bool> duplicate;
        duplicate = regions.insert(pair<string, string>(currentGene.geneName, currentGene.region));
        //This gene was a duplicate of another already entered in. Take the largest transcript region.
        if (duplicate.second == false)
        {
            string oldRegion = regions[currentGene.geneName];
            int oldStart = genePosMap[currentGene.geneName].first;
            int oldEnd = genePosMap[currentGene.geneName].second;
            int newWidth = currentGene.txEndPos - currentGene.txStartPos;
            if (newWidth > (oldEnd - oldStart))
            {
                regions[currentGene.geneName] = currentGene.region;
                genePosMap[currentGene.geneName] = pair<int, int>(currentGene.txStartPos, currentGene.txEndPos);
            }
        }
        else
        {
            //cout << "Entered potential matching gene " << currentGene.geneName << " in region " << currentGene.geneChrom << ":" << currentGene.txStartPos << "-" << currentGene.txEndPos << endl;
            ////genePosMap.insert(pair<string, pair<int, int>>(currentGene.geneName, pair<int, int>(currentGene.txStartPos, currentGene.txEndPos)));
        }
    }
    in.close();
}


pair<vector<int>, vector<int>> readInput::mergeExons(vector<int> oldStarts, vector<int>oldEnds, vector<int> newStarts, vector<int>newEnds)
{
    vector<int> mergedStarts;
    vector<int> mergedEnds;
    vector<int> combinedStarts = oldStarts;
    vector<int> combinedEnds = oldEnds;

    combinedStarts.insert(combinedStarts.end(), newStarts.begin(), newStarts.end());
    combinedEnds.insert(combinedEnds.end(), newEnds.begin(), newEnds.end());

    //Perm contains the sorted order by start position.
    size_t perm[combinedStarts.size()];
    gsl_sort_int_index(perm, combinedStarts.data(), 1, combinedStarts.size());
    
    //The first interval is pushed.
    int currentIndex = 0;
    mergedStarts.push_back(combinedStarts[perm[0]]);
    mergedEnds.push_back(combinedEnds[perm[0]]);

    for(int i = 1; i < combinedStarts.size(); i++)
    {
        //If the next interval lies beyond the current one.
        if(combinedStarts[perm[i]] > mergedEnds[currentIndex])
        {
            //We push it.
            mergedStarts.push_back(combinedStarts[perm[i]]);
            mergedEnds.push_back(combinedEnds[perm[i]]);
            currentIndex++;
        }
        //If the next interval starts before the end of the current one.
        else
        {
            //If it extends the interval we replace the old end.
            if(combinedEnds[perm[i]] > mergedEnds[currentIndex])
            {
                mergedEnds[currentIndex] = combinedEnds[perm[i]];
            }
        }
    }
    // cout << "Old Exons: ";
    // for (int i = 0; i < oldStarts.size(); i++)
    // {
    //     cout << oldStarts[i] << "-" << oldEnds[i] << ",";
    // }
    // cout << endl;
    // cout << "New Exons: ";
    // for (int i = 0; i < newStarts.size(); i++)
    // {
    //     cout << newStarts[i] << "-" << newEnds[i] << ",";
    // }
    // cout << endl;
    // cout << "Combined Exons: ";
    // for (int i = 0; i < mergedStarts.size(); i++)
    // {
    //     cout << mergedStarts[i] << "-" << mergedEnds[i] << ",";
    // }
    // cout << endl;

    return pair<vector<int>, vector<int>>(mergedStarts, mergedEnds);
}




