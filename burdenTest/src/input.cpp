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



readInput::readInput(string dir, string tType, string inputVcfType, string userVcf, string region, string phenoFile, string covFile)
{
    variantCount = 0;
    subjectCount = 0;
    caseCount = 0;
    variantRegion = region;
    vcfType = inputVcfType;
    testType = tType;
    testDir = dir;
    
    gMatch = regex("(\\.\\/\\.)|(\\d)(\\/|\\|)(\\d)");
    altAlleleCountMatch = regex("[ATGCN],[ATGCN]");
    mafMatch = regex(";AF=(0\\.*\\d*)");
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    headerMatch = regex("^#");
    posMatch = regex("^(\\d+)\\t(\\d+)");
    
    //Switching on test type to get proper input.
    if(testType == "wsbt")
    {
        readVcfInitialInfo(userVcf);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
        readMaf(userVcf);
        bcfInput(userVcf);
        readCaseCount(userVcf);
        //mergeData(userVcf);
    }
    else if(testType == "burden")
    {
        //Get basic info from VCF
        readVcfInitialInfo(userVcf);
        
        //Initilizing genotype data structures.
        genotypeGslMatrix = gsl_matrix_calloc(variantCount, subjectCount);
        
        //Initilizing maf data structure.
        maf = gsl_vector_alloc(variantCount);
        
        //Parse the maf
        readMaf(userVcf);

        //Input genotype data.
        bcfInput(userVcf);
        
    }
    else if (testType == "cast")
    {
        
    }
    else if (testType == "skat")
    {
        readVcfInitialInfo(userVcf);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
        pheno = gsl_vector_alloc(subjectCount);
        maf = gsl_vector_alloc(variantCount);

        readMaf(userVcf);
        bcfInput(userVcf);
        readCaseCount(userVcf);
        if(variantRegion == "-g")
        {
            readGenes("refFlat.txt");
        }
        

        if(phenoFile == "")
        {
            gsl_vector_set_zero(pheno);
            for(int i = 0; i < caseCount; i++)
            {
                gsl_vector_set(pheno, i, 1);
            }
        }
        else
        {
            readPhenotype(phenoFile);
        }
        
        if(covFile == "")
        {
            covariates = gsl_matrix_alloc(subjectCount, 1);
            gsl_matrix_set_all(covariates, 1);
        }
        else
        {

        }
    }
    else if (testType == "skato")
    {
        
    }
    
}

void readInput::readVcfInitialInfo(string filename)
{
    string line;
    string statsFileName = "";
    smatch match;
    subjectCount = 0;
    variantCount = 0;
    regex caseMatch("(\\t(\\d[^\\s]+))+");
    
    //So we can remove the directory that has the file and just keep the filename itself.
    int filePos;
    for(int i = filename.length(); i >= 0; i--)
    {
        if (filename[i] == '/')
        {
            filePos = i + 1;
            break;
        }
    }

    //Remove the filetype from filename.
    if(filename.substr(filename.length() - 7) == ".vcf.gz")
    {
        statsFileName = filename.substr(filePos, filename.length() - 7);
    }
    else if(filename.substr(filename.length() - 4) == ".vcf")
    {
        statsFileName = filename.substr(filePos, filename.length() - vcfType.length());
    }
    cout << "statsFileName: " << statsFileName << endl;
    string summaryCommand = bcftools_loc + " stats " + filename + " > " + testDir + "/tmp/" + statsFileName + ".stats";
    cout << "summaryCommand: " << summaryCommand << endl;
    system(summaryCommand.c_str());
    
    if(!inputFile.is_open())
    {
        inputFile.open(testDir + "/tmp/" + statsFileName + ".stats");
        for(int j = 0; getline(inputFile, line); j++)
        {
            if (subjectCount == 0)
            {
                if(regex_search(line, match, subjectCountMatch))
                {
		            //cout << "Match we are trying to stoi: " << match[1] << endl;
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

	    cout << "Subject count: " << subjectCount << endl;
	    cout << "Variant count: " << variantCount << endl;
        //Initilize the maf vector in case its used.
        maf = gsl_vector_alloc(variantCount);
    }
}

void readInput::readCaseCount(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;
        regex caseMatch("(\\t(\\d[^\\s]+))+");
        regex_token_iterator<string::iterator> lineEnd;
        
        if(vcfType == "-gzvcf")
        {
            string genoCommand = bgzip_loc + " -f -d " + filename;
            system(genoCommand.c_str());
            inputFile.open(filename.substr(0, filename.length() - 3));
        }
        else
        {
            inputFile.open(filename);
        }
        for(int j = 0; getline(inputFile, line); j++)
        {
            if(line[0] == '#')
            {
                //Since the headers dont have tabs, they will not match caseMatch prematurely.
                if(caseCount == 0)
                {
                    string templine;
                    if(regex_search(line, match, caseMatch))
                    {
                        templine = match[0];
                        regex tempMatch = regex("[^\\s]+");
                        regex_token_iterator<string::iterator> caseParser(templine.begin(), templine.end(), tempMatch);
                        while(caseParser != lineEnd)
                        {
                            caseCount++;
                            *caseParser++;
                        }
                        return;
                    }
                }
            }
        }
        inputFile.close();
    }
}

//If user has phenotype data we save it in a vector otherwise we assume binary phenotype.
void readInput::readPhenotype(string phenoFile)
{
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
        string line;
        inputFile.open(phenoFile);
        int phenoPos = 0;
        for(int i = 0; getline(inputFile, line, '\t'); i++)
        {
            if(((i + 1) % 6) == 0)
            {
                gsl_vector_set(pheno, phenoPos, stod(line));
                phenoPos++;
            }
        }
        ofstream outFile;
        outFile.open("tmp/inputpheno.txt");
        for (int i = 0; i < pheno->size; i++)
        {
            outFile << gsl_vector_get(pheno, i) << endl;
        }
        outFile.close();
        inputFile.close();
    }
}


void readInput::readMaf(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        string command = externals_loc + "bcftools query -f '%AF\\n' " + filename + " > mafdata.txt";
        system(command.c_str());
        maf = gsl_vector_alloc(variantCount);
        
        inputFile.open("mafdata.txt");
        for(int i = 0; getline(inputFile, line); i++)
        {
            gsl_vector_set(maf, i, stod(line));
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
    
}

//Merge the user vcf file with the right background set. The resulting vcf will be used in readGenotype().
//Merge Data is broken until we set it up with bcfInput().
void readInput::mergeData(string user_filename)
{
    gsl_vector* chr = gsl_vector_alloc(30);
    int currentChr = 0;
    string region = "gene";
    regex chrMatch("^\\d+");

    //Builds background with entire matching chromosome background data.
    if(strcmp(region.c_str(), "all") == 0)
    {
        readVcfInitialInfo(user_filename);
        int userSubjectCount = subjectCount;
        readVcfInitialInfo("background.vcf");
        subjectCount += userSubjectCount;
        cout << "Merge Data final subject count, variant count: " << subjectCount << ", " << variantCount << endl;
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);

        //Read in data.
        //readGenotype(user_filename, genotypeGslMatrix, 0);
        //readGenotype("background.vcf", genotypeGslMatrix, userSubjectCount);
    }
    //Builds background with matching gene region data.
    else if(strcmp(region.c_str(), "gene") == 0)
    {
        //Make the background data vcf


        if(vcfType == "-vcf")
        {
            string zip = bgzip_loc + " -f " + user_filename;
            system(zip.c_str());
            string tabix = externals_loc + "tabix -f -p vcf -s1 -b2 -e3 " + user_filename + ".gz";
            system(tabix.c_str());
            string merge = externals_loc + "bcftools merge -o mergetest.vcf " + user_filename + ".gz " + "fullbackground.vcf.gz";
            system(merge.c_str());
        }
        else
        {
            string tabix = externals_loc + "tabix -f -p vcf -s1 -b2 -e3 " + user_filename;
            system(tabix.c_str());
            string merge = externals_loc + "bcftools merge -o mergetest.vcf " + user_filename + " " + "fullbackground.vcf.gz";
            system(merge.c_str());
        }

        //Now read in from combined vcf
        readVcfInitialInfo("mergetest.vcf");
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
        //readGenotype("mergetest.vcf", genotypeGslMatrix, 0);
    }
    //Builds background with matching chromosome and position data.
    else if(strcmp(region.c_str(), "exact") == 0)
    {
        readVcfInitialInfo(user_filename);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount + 2504);
    }

}

//Reading in genotype data from a bcf should be much faster.
void readInput::bcfInput(string filename)
{
    string line;
    string geneFile = "tmp/" + filename + ".genodata.txt";
    string command = externals_loc + "bcftools query -f '%CHROM,%POS[ %GT]\\n' " + filename + " > " + geneFile;
    //cout << "Parse command: " << command << endl;
    system(command.c_str());

    string currentChrom;
    inputFile.open(geneFile);
    for(int i = 0; getline(inputFile, line); i++)
    {
        //Check for passing filter.
        /*
        if(line.find("PASS") == string::npos)
        {
            cout << "Variant did not pass VCF file's filter. Skipping." << endl;
            continue;
        }
        */
        //Pull out variant chomosome and base pair position.
        string::iterator it = line.begin();
        int comma = line.find_first_of(',');
        int pos = line.find_first_of(' ');
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

        vcfPos.push_back(stoi(line.substr(comma + 1, pos)));
        advance(it, pos);

        for(int j = 0; it != line.end(); j++)
        {
            //Skip space.
            it++;
            char leftAllele = *it++;
            char phase = *it++;
            char rightAllele = *it++;
            int left = 0;
            int right = 0;
            if(leftAllele == '.')
            {
                if(testType == "wsbt")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                }
                if (testType == "skat")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
                }
            }
            //Not set up to handle multiallelic sites.
            else if(leftAllele != '0')
            {
                left = 1;
            }
            if(rightAllele == '.')
            {
                //If we have already taken care of the missing data we dont need to again.
                if (leftAllele != '.')
                {
                    if (testType == "wsbt")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    }
                    if (testType == "skat")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
                    }
                }
                continue;
            }
            else if (rightAllele != '0')
            {
                right = 1;
            }  
            gsl_matrix_set(genotypeGslMatrix, i, j, left + right);
        }
    }
    posMap.insert(pair<string, vector<int> >(currentChrom, vcfPos));
    cout << vcfPos.size() << " variants in chromosome " << currentChrom << endl;
    inputFile.close();
}

//Only run after reading in combined data set.
//Generates a gsl_matrix_view for the variants in a given gene and saves it in a dictionary with the gene name as key.
void readInput::readGenes(string filename)
{
    string line;
    vector<string> geneName;        //Column 1
    vector<string> transcriptName;  //Column 2
    vector<string> geneChrom;       //Column 3
    vector<int> txStartPos;         //Column 5
    vector<int> txEndPos;           //Column 6
    vector<int> codingStartPos;     //Column 7
    vector<int> codingEndPos;       //Column 8

    //cout << "Reading gene data from " << filename << endl;

    //Honestly have no idea why the input stream is still in use.
    if (inputFile.is_open())
    {
        //cout << "inputFile already in use." << endl;
        inputFile.close();
    }
    inputFile.open(filename);
    for(int i = 0; getline(inputFile, line); i++)
    {
        size_t last = 0;
        size_t current = line.find('\t', last);
        for (int j = 0; current != string::npos; j++)
        {
            string token = line.substr(last, current - last);
            if (j == 0)
            {
                geneName.push_back(token);
            }
            if (j == 1)
            {
                transcriptName.push_back(token);
            }
            if (j == 2)
            {
                //Not chr#
                geneChrom.push_back(token.substr(3));
            }
            if (j == 4)
            {
                txStartPos.push_back(stoi(token));
            }
            if (j == 5)
            {
                txEndPos.push_back(stoi(token));
            }
            if (j == 6)
            {
                codingStartPos.push_back(stoi(token));
            }
            if (j == 7)
            {
                codingEndPos.push_back(stoi(token));
            }
            last = current + 1;
            current = line.find('\t', last);
        }
    }
    inputFile.close();
    
    string currentGene = geneName[0];
    vector<gsl_vector *> tempVariant;
    map<string, string> transcript;
    vector<double> tempMaf;
    cout << endl << "Looking through " << geneName.size() << " genes." << endl;
    //Subsets the user data so variants match with known genes.
    //Could improve speed by making sure we dont re-search over the same chromosome when looking for user variant matches.
    for(int i = 0; i < geneName.size(); i++)
    {
        //Once we have a new gene, gather all variants and then reset.
        if(currentGene != geneName[i] || i + 1 == geneName.size())
        {
            if (tempVariant.size() == 0)
            {
                //Gene chromosome matched some variants chromosome in user data set but positions didnt match.
            }
            else
            {
                gsl_matrix *genotypeSubset = gsl_matrix_alloc((int)tempVariant.size(), subjectCount);
                gsl_vector *mafSubset = gsl_vector_alloc((int)tempMaf.size());
                for (int j = 0; j < tempVariant.size(); j++)
                {
                    gsl_matrix_set_row(genotypeSubset, j, tempVariant[j]);
                    gsl_vector_set(mafSubset, j, tempMaf[j]);
                }
                pair<map<string,gsl_matrix *>::iterator,bool> duplicate;
                duplicate = genes.insert(pair<string, gsl_matrix *>(currentGene, genotypeSubset));
                //If we have already found the gene in a different transcript ID.
                //This is where we select the gene description with the most variants.
                if (duplicate.second == false)
                {
                    //If the new gene is bigger than the old one.
                    if (genes[currentGene]->size1 < genotypeSubset->size1)
                    {
                        genes[currentGene] = genotypeSubset;
                        geneMaf[currentGene] = mafSubset;
                        transcript[currentGene] = transcriptName[i-1];
                    }
                }
                else 
                {   
                    geneMaf.insert(pair<string, gsl_vector *>(currentGene, mafSubset));
                    transcript.insert(pair<string, string>(currentGene, transcriptName[i-1]));
                }
                
            }
            //Set up for next gene.
            currentGene = geneName[i];
            tempVariant.clear();
            tempMaf.clear();
        }
        //Check if the chromosome of current gene is even in the data set.
        if (posMap.find(geneChrom[i]) == posMap.end())
        {
            continue;
        }

        //Iterate through all the variants in a chromosome to find a match.
        for(int j = 0; j < posMap[geneChrom[i]].size(); j++)
        {
            //Consider not searching the early positions multiple times. Could initiate j to the last seen useful position.
            if (posMap[geneChrom[i]][j] >= txStartPos[i] && posMap[geneChrom[i]][j] <= txEndPos[i])
            {
                gsl_vector *tempGenoData = gsl_vector_alloc(subjectCount);
                gsl_matrix_get_row(tempGenoData, genotypeGslMatrix, j);
                tempVariant.push_back(tempGenoData);
                tempMaf.push_back(gsl_vector_get(maf, j));
            }
            if (posMap[geneChrom[i]][j] > txEndPos[i])
            {
                break;
            }
        }
    }
    /*
    //Handle data at the end of the geneList.
    if (tempVariant.size() == 0)
    {
        //Gene chromosome matched some variants chromosome in user data set but positions didnt match.
        //cout << geneName[i] << " was not found in data set. Skipping." << endl;
    }
    else
    {
        //cout << "Matched " << tempVariant.size() << " variant(s) to " << currentGene << endl;
        gsl_matrix *genotypeSubset = gsl_matrix_alloc((int)tempVariant.size(), subjectCount);
        gsl_vector *mafSubset = gsl_vector_alloc((int)tempMaf.size());
        for (int j = 0; j < tempVariant.size(); j++)
        {
            gsl_matrix_set_row(genotypeSubset, j, tempVariant[j]);
            gsl_vector_set(mafSubset, j, tempMaf[j]);
        }
        pair<map<string, gsl_matrix *>::iterator, bool> duplicate;
        duplicate = genes.insert(pair<string, gsl_matrix *>(currentGene, genotypeSubset));
        //If we have already found the gene in a different transcript ID.
        if (duplicate.second == false)
        {
            //If the new gene is bigger than the old one.
            if (genes[currentGene]->size1 < genotypeSubset->size1)
            {
                genes[currentGene] = genotypeSubset;
                geneMaf[currentGene] = mafSubset;
                transcript[currentGene] = transcriptName[transcriptName.size() - 1];
            }
        }
        else
        {
            geneMaf.insert(pair<string, gsl_vector *>(currentGene, mafSubset));
            transcript.insert(pair<string, string>(currentGene, transcriptName[transcriptName.size() - 1]));
        }
    }
    */

    map<string, gsl_matrix*>::iterator it;
    for (it = genes.begin(); it != genes.end(); it++)
    {
        cout << "Matched " << it->second->size1 << " variant(s) to " << it->first << " at transcript ID " << transcript[it->first] << endl;
    }

}


