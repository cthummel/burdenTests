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



readInput::readInput(string tType, string inputVcfType, string userVcf, string region, string phenoFile, string covFile)
{
    variantCount = 0;
    subjectCount = 0;
    caseCount = 0;
    vcfType = inputVcfType;
    testType = tType;
    
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
        //readGenotype(userVcf, genotypeGslMatrix, 0);
        readMaf(userVcf);
        bcfInput(userVcf);
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
        
        //Runs vcftools command to separate genotype data.
        //system(genoCommand.c_str());
        readGenotype(userVcf, genotypeGslMatrix, 0);
        
        //Parse the maf
        readMaf(userVcf);
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
        //readGenotype(userVcf, genotypeGslMatrix, 0);
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
    
    //Remove the filetype from filename.
    if(filename.substr(filename.length() - 7) == ".vcf.gz")
    {
        statsFileName = filename.substr(0, filename.length() - vcfType.length() - 1);
    }
    else if(filename.substr(filename.length() - 4) == ".vcf")
    {
        statsFileName = filename.substr(0, filename.length() - vcfType.length());
    }
    
    string summaryCommand = bcftools_loc + " stats " + filename + " > " + statsFileName + ".stats";
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

void readInput::readGenotype(string filename, gsl_matrix *inputMatrix, int subjectOffset)
{
    if(!inputFile.is_open())
    {
        string line;
        smatch match;
        regex caseMatch("(\\t(\\d[^\\s]+))+");
        regex_token_iterator<string::iterator> lineEnd;
        
        double progress = 0.0;
        int barWidth = 50;
        int mafPos = 0;
        int matrixInputLine = 0;
        
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
        bool binRead = false;
        if(binRead)
        {
            inputFile.close();
            FILE *input;
            input = fopen("SKAT_GenoData_Binary", "r");
            gsl_matrix_fread(input, inputMatrix);
            fclose(input);
            return;
        }
        
        for(int j = 0; getline(inputFile, line); j++)
        {
            //string templine = line;
            
            //This is useful if we want to catch the number of ALTs listed for a variant in the vcf.
            //unsigned long int altAlleleCount = 1;
            /*
            if (regex_search(line, match, altAlleleCountMatch))
            {
                altAlleleCount = match.size() + 1;
            }
            
            //Cutting out the background data from the testcase
            regex background("[A-Z]{2}\\d{5}");
            if (regex_search(line, match, background))
            {
                regex_token_iterator<string::iterator> backParser(line.begin(), line.end(), background);
                ofstream outfile;
                outfile.open("background.txt");
                while(backParser != lineEnd)
                {
                    outfile << *backParser++ << endl;
                }
                outfile.close();
            }
            */
            
            if(line[0] != '#')
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
                        if(testType == "skat")
                        {
                            //gsl_matrix_set(inputMatrix, matrixInputLine, i + subjectOffset, 2 * gsl_vector_get(maf, matrixInputLine));
                            gsl_matrix_set(inputMatrix, matrixInputLine, i + subjectOffset, 9);
                        }
                        if (testType == "wsbt")
                        {
                            gsl_matrix_set(inputMatrix, matrixInputLine, i + subjectOffset, -1);
                        }
                        
                        advance(genoParser, 3);
                    }
                    else if(*genoParser++ == "")
                    {
                        int left = 0;
                        int right = 0;
                        string temp;
			
                        //left = stoi(*genoParser++);
                        temp = *genoParser++;
                        if(temp.compare("0") != 0)
                        {
                            left = 1;
                        }
                        //right = stoi(*genoParser++);
                        //if(right > 0)
                        temp = *genoParser++;
                        if(temp.compare("0") != 0)
                        {
                            right = 1;
                        }
                        gsl_matrix_set(inputMatrix, matrixInputLine, i + subjectOffset, left + right);
                    }
                }
                matrixInputLine++;
            }
            else
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
                    }
                }
            }
            
            //Will remove this loading bar stuff later.
            //It was copy and pasted from the internet and is just there for testing sanity.
            /*
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
            */
             
             
        }
        inputFile.close();
        if(testType == "skat")
        {
            FILE *outfile;
            outfile = fopen("SKAT_GenoData_Binary", "w");
            int pass = gsl_matrix_fwrite(outfile, inputMatrix);
            if(pass != 0)
            {
                cout << "writing matrix failed" << endl;
            }
            fclose(outfile);
        }
        
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
        for(int i = 0; getline(inputFile, line); i++)
        {

        }
        inputFile.close();
    }
}


void readInput::readMaf(string filename)
{
    if(!inputFile.is_open())
    {
        string line;
        string command = externals_loc + "bcftools query -f '%AF\\n' " + filename + " > mafdata.txt";
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
    //filename here should be the background file we use for unaffected.
    //string subsetCommand = "vcftools -" + vcfType + " " + backgroundVcf + " --positions pos.txt --recode";
    //system(subsetCommand.c_str());
    
}

//Merge the user vcf file with the right background set. The resulting vcf will be used in readGenotype().
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
        readGenotype(user_filename, genotypeGslMatrix, 0);
        readGenotype("background.vcf", genotypeGslMatrix, userSubjectCount);
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
        readGenotype("mergetest.vcf", genotypeGslMatrix, 0);
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
    string command = externals_loc + "bcftools query -f '[ %GT]\\n' " + filename + " > bcfgenodata.txt";
    cout << "Parse command: " << command << endl;
    system(command.c_str());
    inputFile.open("bcfgenodata.txt");

    caseCount = 4;

    for(int i = 0; getline(inputFile, line); i++)
    {
        int j = 0;
        for(string::iterator it = line.begin(); it != line.end(); j++)
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

    inputFile.close();
}








