//
//  input.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/6/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "input.hpp"
#include <fstream>
//#include <sstream>

using namespace std;

//Sets up initial test paramters from user information.
readInput::readInput(string dir, string tType, string inputVcfType, string userVcf, string region, string phenoFile, string covFile)
{
    variantCount = 0;
    subjectCount = 0;
    caseCount = 0;
    variantRegion = region;
    vcfType = inputVcfType;
    testType = tType;
    testDir = dir;
    
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    
    //Switching on test type to get proper input.
    if(testType == "wsbt")
    {
        bool newMethod = true;
        if(newMethod)
        {
            readVcfInitialInfo(userVcf);
            readCaseCount(userVcf);
            buildPosMap(userVcf);
            if (variantRegion == "-g")
            {
                buildGeneInfo("orderedRefFlat.txt");
                matchGenes();
            }
        }
        else
        {
            readVcfInitialInfo(userVcf);
            genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
            readMaf(userVcf);
            bcfInput(userVcf);
            readCaseCount(userVcf);
            if (variantRegion == "-g")
            {
                buildGeneInfo("orderedRefFlat.txt");
                readGenes("orderedRefFlat.txt");
            }
        }

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
            buildGeneInfo("orderedrefFlat.txt");
            readGenes("orderedrefFlat.txt");
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

//This is used to read in the data set for each gene when we need them.
//Only works on pre-merged vcf files where user has attached their data set to 1000Genomes.
readInput::readInput(bool userBackgroundIncluded, string user, string region, int count, int thread_ID)
{
    subjectCount = 0;
    variantCount = 0;
    caseCount = count;
    
    if(userBackgroundIncluded)
    {
        readVcfInitialInfo(user, region, thread_ID);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
        maf = gsl_vector_alloc(variantCount);
        //Parsing geno data.
        bcfInput(user, region, "data" + to_string(thread_ID) + ".txt");
        readMaf(user, region, "maf" + to_string(thread_ID) + ".txt");
    }
    else
    {
        string backgroundFilename = region.substr(0, region.find_first_of(':'));

        readVcfInitialInfo(backgroundFilename, region, thread_ID);
        genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount + count);
        maf = gsl_vector_alloc(variantCount);
        mergeData(user, backgroundFilename, region, "merge" + to_string(thread_ID) + ".txt");
    }

    
}

void readInput::readVcfInitialInfo(string filename, string region, int thread_ID)
{
    string line;
    ifstream in;
    smatch match;
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    
    //Build our stats file and GT file for later parsing.
    string command = externals_loc + "bcftools stats -r " + region + " " + filename + " > data" + to_string(thread_ID) + ".stats";
    system(command.c_str());

    //Parsing background info
    in.open("data" + to_string(thread_ID) + ".stats");
    for (int j = 0; getline(in, line); j++)
    {
        if (subjectCount == 0)
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
        if (variantCount != 0 && subjectCount != 0)
        {
            break;
        }
    }
    in.close();

    cout << "Subject count: " << subjectCount << endl;
    cout << "Variant count: " << variantCount << endl;
}


void readInput::readVcfInitialInfo(string filename)
{
    string line;
    string statsFileName = "";
    smatch match;
    ifstream in;
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

    string summaryCommand = bcftools_loc + " stats " + filename + " > " + testDir + "/tmp/" + statsFileName + ".stats";
    system(summaryCommand.c_str());
    
    if(!in.is_open())
    {
        in.open(testDir + "/tmp/" + statsFileName + ".stats");
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
    bool mergeMode = false;
    if(mergeMode)
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
                    gsl_vector_set(pheno, i, stoi(token));
                }
                
                last = current + 1;
                current = line.find('\t', last);
            }
        }
        in.close();
    }
}

void readInput::readMaf(string filename)
{
    ifstream in;
    if(!in.is_open())
    {
        string line;
        string command = externals_loc + "bcftools query -f '%AF\\n' " + filename + " > mafdata.txt";
        system(command.c_str());
        maf = gsl_vector_alloc(variantCount);
        in.open("mafdata.txt");
        for(int i = 0; getline(in, line); i++)
        {
            int alleleCount = 1;
            if (line.find(',') == string::npos)
            {
                gsl_vector_set(maf, i, stod(line));
            }
            else
            {
                size_t last = 0;
                size_t current = line.find(',', last);
                vector<double> mafs;
                while(current != string::npos)
                {
                    string token = line.substr(last, current - last);
                    mafs.push_back(stod(token));
                    last = current + 1;
                    current = line.find(',', last);
                }
                double max = 0;
                for(int i = 0; i < mafs.size(); i++)
                {
                    if (max < mafs[i])
                    {
                        max = mafs[i];
                    }
                }
                gsl_vector_set(maf, i, max);
            }
        }
        in.close();
    }
}

//Gene version of reading MAF data.
void readInput::readMaf(string filename, string region, string outfile)
{
    ifstream in;
    if(!in.is_open())
    {
        string line;
        string command = externals_loc + "bcftools query -r " + region + " -f '%AF\\n' " + filename + " > " + outfile;
        system(command.c_str());
        in.open(outfile);
        for(int i = 0; getline(in, line); i++)
        {
            int alleleCount = 1;
            if (line.find(',') == string::npos)
            {
                gsl_vector_set(maf, i, stod(line));
            }
            else
            {
                size_t last = 0;
                size_t current = line.find(',', last);
                vector<double> mafs;
                while(current != string::npos)
                {
                    string token = line.substr(last, current - last);
                    mafs.push_back(stod(token));
                    last = current + 1;
                    current = line.find(',', last);
                }
                double max = 0;
                for(int i = 0; i < mafs.size(); i++)
                {
                    if (max < mafs[i])
                    {
                        max = mafs[i];
                    }
                }
                gsl_vector_set(maf, i, max);
            }
        }
        in.close();
    }
}

//Basically the exact same as bcfInput but this time it has the merge step in the console command.
//Maybe dont need this function to be separate and can just bake into bcfInput.
void readInput::mergeData(string user, string back, string region, string outfile)
{
    string mergeCommand = externals_loc + "bcftools merge -r " + region + " " + user + " " + back + " | " +
                          externals_loc + "bcftools query -f '%CHROM,%POS[ %GT]\\n' > " + outfile;
    system(mergeCommand.c_str());

    ifstream in;
    string line;
    in.open(outfile);
    for(int i = 0; getline(in, line); i++)
    {
        string::iterator it = line.begin();
        bool missingData = false;
        int missingCount = 0;
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
                missingData = true;
                missingCount++;
                if(testType == "wsbt")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                }
                if (testType == "skat")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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
                    missingData = true;
                    missingCount++;
                    if (testType == "wsbt")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    }
                    if (testType == "skat")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                        //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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

        //Fix missing data points by imputing their value with the mean of the geno data for the variant
        if(missingData && testType != "wsbt")
        {
            int alleleCount = 0;
            int denominator = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                //Getting counts of non-missing alleles
                if (gsl_matrix_get(genotypeGslMatrix, i, j) >= 0)
                {
                    alleleCount += gsl_matrix_get(genotypeGslMatrix, i, j);
                    denominator++;
                }
            }
            double mean = (1.0 * alleleCount) / denominator;
            int fixed = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                if (gsl_matrix_get(genotypeGslMatrix, i, j) < 0)
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, mean);
                    fixed++;
                    if (fixed == missingCount)
                    {
                        break;
                    }
                }
            }
        }
        
    }
    in.close();
}

//Reading in genotype data from a bcf should be much faster.
void readInput::bcfInput(string filename)
{
    ifstream in;
    string line;
    vector<int> vcfPos;
    string geneFile = "tmp/" + filename + ".genodata.txt";
    string command = externals_loc + "bcftools query -f '%CHROM,%POS[ %GT]\\n' " + filename + " > " + geneFile;
    //cout << "Parse command: " << command << endl;
    system(command.c_str());

    string currentChrom;
    in.open(geneFile);
    for(int i = 0; getline(in, line); i++)
    {
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

        bool missingData = false;
        int missingCount = 0;
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
                missingData = true;
                missingCount++;
                if(testType == "wsbt")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                }
                if (testType == "skat")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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
                    missingData = true;
                    missingCount++;
                    if (testType == "wsbt")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    }
                    if (testType == "skat")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                        //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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

        //Fix missing data points by imputing their value with the mean of the geno data for the variant
        if(missingData && testType != "wsbt")
        {
            int alleleCount = 0;
            int denominator = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                //Getting counts of non-missing alleles
                if (gsl_matrix_get(genotypeGslMatrix, i, j) >= 0)
                {
                    alleleCount += gsl_matrix_get(genotypeGslMatrix, i, j);
                    denominator++;
                }
            }
            double mean = (1.0 * alleleCount) / denominator;
            //cout << mean << " ";
            int fixed = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                if (gsl_matrix_get(genotypeGslMatrix, i, j) < 0)
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, mean);
                    fixed++;
                    if (fixed == missingCount)
                    {
                        break;
                    }
                }
            }
        }
        
    }
    posMap.insert(pair<string, vector<int> >(currentChrom, vcfPos));
    cout << vcfPos.size() << " variants in chromosome " << currentChrom << endl;
    in.close();
}

//Gene version of input.
void readInput::bcfInput(string filename, string region, string outfile)
{
    ifstream in;
    string line;
    string command = externals_loc + "bcftools query -r " + region + " -f '[ %GT]\\n' " + filename + " > " + outfile;
    system(command.c_str());

    string currentChrom;
    in.open(outfile);
    for(int i = 0; getline(in, line); i++)
    {
        string::iterator it = line.begin();
        bool missingData = false;
        int missingCount = 0;
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
                missingData = true;
                missingCount++;
                if(testType == "wsbt")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                }
                if (testType == "skat")
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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
                    missingData = true;
                    missingCount++;
                    if (testType == "wsbt")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                    }
                    if (testType == "skat")
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, -1);
                        //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
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

        //Fix missing data points by imputing their value with the mean of the geno data for the variant
        if(missingData && testType != "wsbt")
        {
            int alleleCount = 0;
            int denominator = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                //Getting counts of non-missing alleles
                if (gsl_matrix_get(genotypeGslMatrix, i, j) >= 0)
                {
                    alleleCount += gsl_matrix_get(genotypeGslMatrix, i, j);
                    denominator++;
                }
            }
            double mean = (1.0 * alleleCount) / denominator;
            int fixed = 0;
            for(int j = 0; j < subjectCount; j++)
            {
                if (gsl_matrix_get(genotypeGslMatrix, i, j) < 0)
                {
                    gsl_matrix_set(genotypeGslMatrix, i, j, mean);
                    fixed++;
                    if (fixed == missingCount)
                    {
                        break;
                    }
                }
            }
        }
        
    }
    in.close();
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
        for (int j = 0; current != string::npos; j++)
        {
            string token = line.substr(last, current - last);
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
            last = current + 1;
            current = line.find('\t', last);
        }
        info.push_back(currentGene);
    }
    in.close();
}

//Only run after reading in combined data set.
//Generates a gsl_matrix_view for the variants in a given gene and saves it in a dictionary with the gene name as key.
void readInput::readGenes(string filename)
{
    string line;
    string currentGene = info[0].geneName;
    vector<gsl_vector *> tempVariant;
    map<string, string> transcript;
    vector<double> tempMaf;
    cout << endl << "Looking through " << info.size() << " gene entries." << endl;
    //Subsets the user data so variants match with known genes.
    //Could improve speed by making sure we dont re-search over the same chromosome when looking for user variant matches.
    for (int i = 0; i < info.size(); i++)
    {
        //Once we have a new gene, gather all variants and then reset.
        //if (currentGene != info[i].geneName || i + 1 == info.size())
        //{
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
                pair<map<string, gsl_matrix *>::iterator, bool> duplicate;
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
                        transcript[currentGene] = info[i-1].transcriptName;
                        string temp = info[i-1].geneChrom + ":" + to_string(info[i-1].txStartPos) + "-" + to_string(info[i-1].txEndPos);
                        region[currentGene] = temp;
                    }
                }
                else
                {
                    geneMaf.insert(pair<string, gsl_vector *>(currentGene, mafSubset));
                    transcript.insert(pair<string, string>(currentGene, info[i-1].transcriptName));
                    string temp = info[i-1].geneChrom + ":" + to_string(info[i-1].txStartPos) + "-" + to_string(info[i-1].txEndPos);
                    region.insert(pair<string, string>(currentGene, temp));
                }
            }
            //Set up for next gene.
            currentGene = info[i].geneName;
            tempVariant.clear();
            tempMaf.clear();
        //}
        //Check if the chromosome of current gene is even in the data set.
        if (posMap.find(info[i].geneChrom) == posMap.end())
        {
            continue;
        }

        //Iterate through all the variants in a chromosome to find a match.
        for (int j = 0; j < posMap[info[i].geneChrom].size(); j++)
        {
            //Consider not searching the early positions multiple times. Could initiate j to the last seen useful position.
            if (posMap[info[i].geneChrom][j] >= info[i].txStartPos && posMap[info[i].geneChrom][j] <= info[i].txEndPos)
            {
                gsl_vector *tempGenoData = gsl_vector_alloc(subjectCount);
                gsl_matrix_get_row(tempGenoData, genotypeGslMatrix, j);
                tempVariant.push_back(tempGenoData);
                tempMaf.push_back(gsl_vector_get(maf, j));
            }
            if (posMap[info[i].geneChrom][j] > info[i].txEndPos)
            {
                break;
            }
        }
    }
    map<string, gsl_matrix*>::iterator it;
    for (it = genes.begin(); it != genes.end(); it++)
    {
        cout << "Matched " << it->second->size1 << " variant(s) to " << it->first << " at transcript ID " << transcript[it->first] << " with region " << region[it->first] << endl;
    }

}

//Builds a map of variant positions to chromosome in user data.
void readInput::buildPosMap(string filename)
{
    ifstream in;
    string line;
    vector<int> vcfPos;
    string posFile = "tmp/" + filename + ".posdata.txt";
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


void readInput::matchGenes()
{
    string line;
    cout << endl << "Looking through " << info.size() << " gene entries." << endl;
    for (int i = 0; i < info.size(); i++)
    {
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
                string temp = info[i].geneChrom + ":" + to_string(info[i].txStartPos) + "-" + to_string(info[i].txEndPos);
                pair<map<string, string>::iterator, bool> duplicate;
                duplicate = region.insert(pair<string, string>(info[i].geneName, temp));
                //This gene was a duplicate of another already entered in. Take the largest transcript region.
                if(duplicate.second == false)
                {
                    string oldRegion = region[info[i].geneName];
                    int oldStart = genePosMap[info[i].geneName].first;
                    int oldEnd = genePosMap[info[i].geneName].second;
                    int newWidth = info[i].txEndPos - info[i].txStartPos;
                    if(newWidth > (oldEnd - oldStart))
                    {
                        region[info[i].geneName] = temp;
                        genePosMap[info[i].geneName] = pair<int, int>(info[i].txStartPos, info[i].txEndPos);
                    }
                }
                else
                {
                    genePosMap.insert(pair<string, pair<int, int>>(info[i].geneName, pair<int,int>(info[i].txStartPos,info[i].txEndPos)));
                }
            }
            if (posMap[info[i].geneChrom][j] > info[i].txEndPos)
            {
                break;
            }
        }
    }
    for(map<string,string>::iterator it = region.begin(); it != region.end(); it++)
    {
        cout << "Matched variant(s) to gene " << it->first << " in region " << it->second << endl;
    }
}




void readInput::geneWiseInput(string filename)
{
    ifstream in;
    string line;
    string geneFile = "tmp/" + filename + ".genodata.txt";
    string command = externals_loc + "bcftools query -f '%CHROM,%POS[ %GT]\\n' " + filename + " > " + geneFile;
    //cout << "Parse command: " << command << endl;
    system(command.c_str());

    geneId currentGene = info[0];
    int geneListPos = 0;
    vector<pair<int, gsl_vector*> > genodata;
    in.open(geneFile);
    for (int i = 0; getline(in, line); i++)
    {
        //Pull out variant chomosome and base pair position.
        string::iterator it = line.begin();
        int comma = line.find_first_of(',');
        int pos = line.find_first_of(' ');
        string currentChrom = line.substr(0, comma);
        int currentPos = stoi(line.substr(comma + 1, pos));
        bool update = false;

        if(info[geneListPos].geneChrom != currentChrom)
        {
            for(; geneListPos < info.size(); geneListPos++)
            {
                if(info[geneListPos].geneChrom == currentChrom)
                {
                    update = true;
                    break;
                }
            }   
        }
        if(info[geneListPos].txEndPos < currentPos)
        {
            for(; geneListPos < info.size(); geneListPos++)
            {
                if(info[geneListPos].txEndPos >= currentPos && info[geneListPos].txStartPos <= currentPos)
                {
                    update = true;
                    break;
                }
            }
        }
        //This shouldnt ever occur so it wastes time checking.
        /*
        if(geneListPos == info.size())
        {
            cout << "Variant in input file does not fall within any known genes." << endl;
            return;
        }
        */
        if(update)
        {
            //loadGene(currentGene, genodata);
            genodata.clear();
            currentGene = info[geneListPos];
            update = false;
        }
        
        bool missingData = false;
        int missingCount = 0;
        gsl_vector* varData = gsl_vector_alloc(subjectCount);
        for (int j = 0; it != line.end(); j++)
        {
            //Skip space.
            it++;
            char leftAllele = *it++;
            char phase = *it++;
            char rightAllele = *it++;
            int left = 0;
            int right = 0;
            if (leftAllele == '.')
            {
                missingData = true;
                missingCount++;
                if (testType == "wsbt")
                {
                    gsl_vector_set(varData, j, -1);
                }
                if (testType == "skat")
                {
                    gsl_vector_set(varData, j, -1);
                    //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
                }
            }
            //Not set up to handle multiallelic sites.
            else if (leftAllele != '0')
            {
                left = 1;
            }
            if (rightAllele == '.')
            {
                //If we have already taken care of the missing data we dont need to again.
                if (leftAllele != '.')
                {
                    missingData = true;
                    missingCount++;
                    if (testType == "wsbt")
                    {
                        gsl_vector_set(varData, j, -1);
                    }
                    if (testType == "skat")
                    {
                        gsl_vector_set(varData, j, -1);
                        //gsl_matrix_set(genotypeGslMatrix, i, j, 2 * gsl_vector_get(maf, i));
                    }
                }
                continue;
            }
            else if (rightAllele != '0')
            {
                right = 1;
            }
            gsl_vector_set(varData, j, left + right);
        }
        genodata.push_back(pair<int, gsl_vector *>(currentPos, varData));
    }
    in.close();
}


/*
void readInput::test(string filename)
{

    string line;
    
    vector<string> geneName;       //Column 1
    vector<string> transcriptName; //Column 2
    vector<string> geneChrom;      //Column 3
    vector<int> txStartPos;        //Column 5
    vector<int> txEndPos;          //Column 6
    vector<int> codingStartPos;    //Column 7
    vector<int> codingEndPos;      //Column 8
    

    //Honestly have no idea why the input stream is still in use. Be careful here if trying to multithread this code.
    if (in.is_open())
    {
        //cout << "in already in use." << endl;
        in.close();
    }
    //in.open("orderedRefFlat.txt");
    in.open("refFlat.txt");
    for (int i = 0; getline(in, line); i++)
    {
        size_t last = 0;
        size_t current = line.find('\t', last);
        geneId currentGene;
        for (int j = 0; current != string::npos; j++)
        {
            string token = line.substr(last, current - last);
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
            last = current + 1;
            current = line.find('\t', last);
        }
        stringstream region;
        stringstream in;
        region << currentGene.geneChrom << ":" << currentGene.txStartPos << "-" << currentGene.txEndPos;
        //string command = externals_loc + "bcftools query -r " + region.str() + " -f '%CHROM\\n' " + filename + " | head -c1 | wc -c";
        string command = externals_loc + "bcftools query -r " + region.str() + " -f '%CHROM,%POS[ %GT]\\n' " + filename;
        in = exec(command.c_str());
        if(in.str() != "")
        {
            cout << in.str() << endl;
        }
        
        //cout << currentGene.geneName << " " << in.str()[in.str().length() - 2] << endl;
        if(in.str()[in.str().length() - 2] == 1)
        {
            cout << "We found gene " << currentGene.geneName << endl;
            info.push_back(currentGene);
        }
        
    }
    in.close();
}

*/

/*

bool readInput::testReadFromStream(string filename, string region)
{
    string temp;
    stringstream in;
    //string command = externals_loc + "bcftools query -r " + region + " -f '%CHROM,%POS[ %GT]\\n' " + filename;
    string command = "bcftools query -r " + region + " -f '%CHROM\n' " + filename + " | head -c1 | wc -c";
    temp = exec(command.c_str());
    if (temp[0] == '1')
    {
        return true;
    }
    else
    {
        return false;
    }
}
*/


/*

//Copy and pasted from stackexchange. (modified to use a string stream rather than string for speed)
//https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
string readInput::exec(const char* cmd) {
    char buffer[64];
    string result;
    shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 64, pipe.get()) != nullptr)
            result += buffer;
    }
    return result;
}
*/
