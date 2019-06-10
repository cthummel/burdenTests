#include "dataCollector.hpp"
#include <string>
#include <iostream>


using namespace std;


//This is used to read in the data set for each gene when we need them.
//Only works on pre-merged vcf files where user has attached their data set to 1000Genomes.
dataCollector::dataCollector(bool userBackgroundIncluded, string user, string back, string region, string test_type, int thread_ID)
{
    testType = test_type;
    preMerged = userBackgroundIncluded;
    userFile = user;
    backFile = back;

    if (userBackgroundIncluded)
    {
        readVcfInitialInfo(user, region, "tmp/data" + to_string(thread_ID) + ".stats");
        if(variantCount != 0 && subjectCount != 0)
        {
            genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
            bcfInput(user, "", region, "tmp/data" + to_string(thread_ID) + ".txt");
            readMaf(user, region, "tmp/maf" + to_string(thread_ID) + ".txt");
        }
        else
        {
            cout << "Data set did not contain variants in " << region << "." << endl;
        }
        
    }
    else
    {
        if(test_type == "wsbt")
        {
            readVcfInitialInfo(backFile, region, "tmp/data" + to_string(thread_ID) + ".stats");
            if(variantCount != 0 && subjectCount != 0)
            {
                genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
                bcfInput(user, backFile, region, "tmp/merge" + to_string(thread_ID) + ".txt");
            }
            else
            {
                cout << "Data set did not contain variants in " << region << "." << endl;
            }
            
        }
        else if(test_type == "skat")
        {
            readVcfInitialInfo(backFile, region, "tmp/data" + to_string(thread_ID) + ".stats");
            if(variantCount != 0 && subjectCount != 0)
            {
                genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
                maf = gsl_vector_alloc(variantCount);
                readMaf(user, region, "tmp/maf" + to_string(thread_ID) + ".txt");
                bcfInput(user, backFile, region, "tmp/merge" + to_string(thread_ID) + ".txt");
            }
            else
            {
                cout << "Data set did not contain variants in " << region << "." << endl;
            }
            
        }
        
    }
}

//Destructor
dataCollector::~dataCollector()
{
    //As a dataCollector object goes out of scope it should delete the gsl objects.
    if(genotypeGslMatrix != nullptr)
    {
        gsl_matrix_free(genotypeGslMatrix);
    }
    if(maf != nullptr)
    {
        gsl_vector_free(maf);
    }
}

void dataCollector::readVcfInitialInfo(string filename, string region, string outfile)
{
    subjectCount = -1;
    variantCount = -1;
    string line;
    ifstream in;
    smatch match;
    subjectCountMatch = regex("number of samples:\\s(\\d*)");
    variantCountMatch = regex("number of records:\\s(\\d*)");
    
    //Build our stats file and GT file for later parsing.
    if(preMerged)
    {
        string command = externals_loc + "bcftools stats -r " + region + " " + filename + " > " + outfile;
        system(command.c_str());
    }
    else
    {
        string command = externals_loc + "bcftools merge -r " + region + " " + userFile + " " + backFile +  " | " 
                       + externals_loc + "bcftools stats - > " + outfile;
        system(command.c_str());
    }
    

    //Parsing background info
    in.open(outfile);
    while(getline(in, line))
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
    
}


//Gene version of reading MAF data.
void dataCollector::readMaf(string filename, string region, string outfile)
{
    ifstream in;
    string line;
    string command;
    if(preMerged)
    {
        command = externals_loc + "bcftools query -r " + region + " -f '%AF\\n' " + userFile + " > " + outfile;
    }
    else
    {
        command = externals_loc + "bcftools merge -r " + region + " " + userFile + " " + backFile +  " | " 
                + externals_loc + "bcftools query -r " + region + " -f '%AF\\n' " + userFile + " > " + outfile;
    }
    system(command.c_str());
    in.open(outfile);
    for (int i = 0; i < maf->size; i++)
    {
        getline(in, line);
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
            while (current != string::npos)
            {
                string token = line.substr(last, current - last);
                mafs.push_back(stod(token));
                last = current + 1;
                current = line.find(',', last);
            }
            double max = 0;
            for (int i = 0; i < mafs.size(); i++)
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

void dataCollector::bcfInput(string filename, string back, string region, string outfile)
{
    ifstream in;
    string line;
    if(preMerged)
    {
        string command = externals_loc + "bcftools query -r " + region + " -f '[ %GT]\\n' " + filename + " > " + outfile;
        system(command.c_str());
    }
    else
    {
        string mergeCommand = externals_loc + "bcftools merge -r " + region + " " + filename + " " + back + " | " 
                            + externals_loc + "bcftools query -f '[ %GT]\\n' - > " + outfile;
        system(mergeCommand.c_str());
    }

    string currentChrom;
    in.open(outfile);
    for(int i = 0; i < genotypeGslMatrix->size1; i++)
    {
        getline(in, line);
        string::iterator it = line.begin();
        if(line.length() / 4 != genotypeGslMatrix->size2)
        {
            cerr << "The vcf genotype data for gene " << region << " in line " << i << " has length/4=" << line.length() / 4 << " while j is " << genotypeGslMatrix->size2 << endl;
        }
        bool missingData = false;
        int missingCount = 0;
        for(int j = 0; j < genotypeGslMatrix->size2; j++)
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
                gsl_matrix_set(genotypeGslMatrix, i, j, -1);
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
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
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
        if(missingData)
        {
            if(testType == "wsbt")
            {
                int fixed = 0;
                for (int j = 0; j < genotypeGslMatrix->size2; j++)
                {
                    if (gsl_matrix_get(genotypeGslMatrix, i, j) < 0)
                    {
                        gsl_matrix_set(genotypeGslMatrix, i, j, 0);
                        fixed++;
                        if (fixed == missingCount)
                        {
                            break;
                        }
                    }
                }
            }
            else
            {
                int alleleCount = 0;
                int denominator = 0;
                for (int j = 0; j < genotypeGslMatrix->size2; j++)
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
                for (int j = 0; j < genotypeGslMatrix->size2; j++)
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
    }
    in.close();
}
