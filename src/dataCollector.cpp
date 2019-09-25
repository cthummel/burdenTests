#include "dataCollector.hpp"
#include <string>
#include <sstream>
#include <iostream>


using namespace std;


//This is used to read in the data set for each gene when we need them.
//Only works on pre-merged vcf files where user has attached their data set to 1000Genomes.
dataCollector::dataCollector(bool userBackgroundIncluded, bool includeCADDWeights, string user, string back, string region, string test_type, int userCaseCount, int thread_ID)
{
    testType = test_type;
    preMerged = userBackgroundIncluded;
    userFile = user;
    backFile = back;
    caseCount = userCaseCount;

    if (userBackgroundIncluded)
    {
        readVcfInitialInfo(user, region, "tmp/data" + to_string(thread_ID) + ".stats");
        if(variantCount != 0 && subjectCount != 0)
        {
            shortGenotypeGslMatrix = gsl_matrix_short_alloc(variantCount, subjectCount);
            doubleBcfInput(user, "", region, "tmp/data" + to_string(thread_ID) + ".txt");
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
                shortGenotypeGslMatrix = gsl_matrix_short_alloc(variantCount, subjectCount);
                shortBcfInput(user, backFile, region, "tmp/merge" + to_string(thread_ID) + ".txt");
                if(includeCADDWeights)
                {
                    readCADD(user, back, region, "tmp/annoRaw" + to_string(thread_ID) + ".tsv", "tmp/anno" + to_string(thread_ID) + ".tsv","tmp/CADD" + to_string(thread_ID)+ ".txt");
                }
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
                doubleBcfInput(user, backFile, region, "tmp/merge" + to_string(thread_ID) + ".txt");
            }
            else
            {
                cout << "Data set did not contain variants in " << region << "." << endl;
            }
            
        }
        else if(test_type == "skato")
        {
            readVcfInitialInfo(backFile, region, "tmp/data" + to_string(thread_ID) + ".stats");
            if (variantCount != 0 && subjectCount != 0)
            {
                genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
                maf = gsl_vector_alloc(variantCount);
                readMaf(user, region, "tmp/maf" + to_string(thread_ID) + ".txt");
                doubleBcfInput(user, backFile, region, "tmp/merge" + to_string(thread_ID) + ".txt");
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
    if(shortGenotypeGslMatrix != nullptr)
    {
        gsl_matrix_short_free(shortGenotypeGslMatrix);
    }
    if(maf != nullptr)
    {
        gsl_vector_free(maf);
    }
    if(CADDWeights != nullptr)
    {
        gsl_vector_free(CADDWeights);
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

void dataCollector::shortBcfInput(string filename, string back, string region, string outfile)
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
    for(int i = 0; i < shortGenotypeGslMatrix->size1; i++)
    {
        getline(in, line);
        string::iterator it = line.begin();
        if(line.length() / 4 != shortGenotypeGslMatrix->size2)
        {
            cerr << "The vcf genotype data for gene " << region << " in line " << i << " has length=" << line.length() / 4 << " while subjectCount is " << shortGenotypeGslMatrix->size2 << endl;
        }
        //Used for SKAT test to impute missing genotypes.
        bool missingData = false;
        int missingCount = 0;

        //Will change to false if any non ./. records found.
        bool uniqueCaseVariant = true;
        bool uniqueBackVariant = true;

        for(int j = 0; j < shortGenotypeGslMatrix->size2; j++)
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
                gsl_matrix_short_set(shortGenotypeGslMatrix, i, j, -1);
            }
            else
            {
                //Not set up to handle multiallelic sites.
                if (leftAllele != '0')
                {
                    left = 1;
                }
                //Setting this as a non-unique variant.
                if(j < caseCount)
                {
                    uniqueBackVariant = false;
                }
                else
                {
                    uniqueCaseVariant = false;
                }
                
            }
            if(rightAllele == '.')
            {
                //If we have already taken care of the missing data we dont need to again.
                if (leftAllele != '.')
                {
                    missingData = true;
                    missingCount++;
                    gsl_matrix_short_set(shortGenotypeGslMatrix, i, j, -1);
                }
                continue;
            }
            else
            {
                if (rightAllele != '0')
                {
                    right = 1;
                }
                //Setting this as a non-unique variant.
                if(j < caseCount)
                {
                    uniqueBackVariant = false;
                }
                else
                {
                    uniqueCaseVariant = false;
                }
            }
            gsl_matrix_short_set(shortGenotypeGslMatrix, i, j, left + right);
        }

        if(uniqueCaseVariant)
        {
            caseUniqueVariantCount++;
        }
        if(uniqueBackVariant)
        {
            backgroundUniqueVariantCount++;
        }

        //Fix missing data points by imputing their value with the mean of the geno data for the variant
        if(missingData)
        {
            if(testType == "wsbt")
            {
                continue;
                int fixed = 0;
                for (int j = 0; j < shortGenotypeGslMatrix->size2; j++)
                {
                    if (gsl_matrix_short_get(shortGenotypeGslMatrix, i, j) < 0)
                    {
                        gsl_matrix_short_set(shortGenotypeGslMatrix, i, j, 0);
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
                for (int j = 0; j < shortGenotypeGslMatrix->size2; j++)
                {
                    //Getting counts of non-missing alleles
                    if (gsl_matrix_short_get(shortGenotypeGslMatrix, i, j) >= 0)
                    {
                        alleleCount += gsl_matrix_short_get(shortGenotypeGslMatrix, i, j);
                        denominator++;
                    }
                }
                double mean = (1.0 * alleleCount) / denominator;
                int fixed = 0;
                for (int j = 0; j < shortGenotypeGslMatrix->size2; j++)
                {
                    if (gsl_matrix_short_get(shortGenotypeGslMatrix, i, j) < 0)
                    {
                        gsl_matrix_short_set(shortGenotypeGslMatrix, i, j, mean);
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

void dataCollector::doubleBcfInput(string filename, string back, string region, string outfile)
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
            cerr << "The vcf genotype data for gene " << region << " in line " << i << " has length=" << line.length() / 4 << " while subjectCount is " << genotypeGslMatrix->size2 << endl;
        }
        //Used for SKAT test to impute missing genotypes.
        bool missingData = false;
        int missingCount = 0;

        //Will change to false if any non ./. records found.
        bool uniqueCaseVariant = true;
        bool uniqueBackVariant = true;

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
            else
            {
                //Not set up to handle multiallelic sites.
                if (leftAllele != '0')
                {
                    left = 1;
                }
                //Setting this as a non-unique variant.
                if(j < caseCount)
                {
                    uniqueBackVariant = false;
                }
                else
                {
                    uniqueCaseVariant = false;
                }
                
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
            else
            {
                if (rightAllele != '0')
                {
                    right = 1;
                }
                //Setting this as a non-unique variant.
                if(j < caseCount)
                {
                    uniqueBackVariant = false;
                }
                else
                {
                    uniqueCaseVariant = false;
                }
            }
            gsl_matrix_set(genotypeGslMatrix, i, j, left + right);
        }

        if(uniqueCaseVariant)
        {
            caseUniqueVariantCount++;
        }
        if(uniqueBackVariant)
        {
            backgroundUniqueVariantCount++;
        }

        //Fix missing data points by imputing their value with the mean of the geno data for the variant
        if(missingData)
        {
            if(testType == "wsbt")
            {
                continue;
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

void dataCollector::annotationParser(string filename, string back, string region, string outfile)
{
    string line;
    string command;
    if(preMerged)
    {
        command = externals_loc + "bcftools query -r " + region + " -f '%AF\\n' " + userFile + " > " + outfile;
    }
    else
    {
        command = externals_loc + "bcftools merge -r " + region + " " + userFile + " " + backFile +  " | " 
                //+ externals_loc + "bcftools annotate -a      -r " + region + " " + userFile + " | "
                + externals_loc + "bcftools query -r " + region + " -f '%INFO/CSQ\\n' " + userFile + " > " + outfile;
    }
    system(command.c_str());
    ifstream in(outfile);
    for (int i = 0; i < genotypeGslMatrix->size1; i++)
    {
        getline(in, line);
        size_t last = 0;
        size_t current = line.find('|', last);
        geneId currentGene;
        for (int j = 0; last != string::npos; j++)
        {
            
            string token = line.substr(last, current - last);


            if(current == string::npos)
            {
                //Breaks out of parsing this entry.
                last = string::npos;
            }
            else
            {
                last = current + 1;
                current = line.find('|', last);
            }
        }
    }
}


void dataCollector::weightImport(string region, string outfile)
{
    string command = externals_loc + "tabix ../../data/InDels_inclAnno.tsv.gz " + region + " > " + outfile;;
    system(command.c_str());
    ifstream in(outfile);
}

void dataCollector::readCADD(string filename, string back, string region, string rawFile, string annoFile, string outfile)
{
    if(preMerged)
    {
        //" | awk -v OFS='\\t' '{print $1,$2,$3,$4,$106,$107}' > "
        stringstream annotationFile, combineCommand, zipCommand, annotationCommand;
        annotationFile << externals_loc << "tabix ../../data/InDels_inclAnno.tsv.gz " << region << " | awk -v OFS='\\t' '{print $1,$2,$3,$4,$106,$107}' > " << rawFile;
        combineCommand << "cat tmp/CADDNames.txt " << rawFile << " > " << annoFile;
        zipCommand << externals_loc << "bgzip -f " << annoFile << " && tabix -f -p vcf " << annoFile << ".gz";
        annotationCommand << externals_loc << "bcftools annotate -a " << annoFile << ".gz -h ../../data/headerLines.txt -c CHROM,POS,REF,ALT,RawScore,PHRED -r " << region << " " << filename << " | " << 
                                   externals_loc << "bcftools query -f '%INFO/PHRED\\n' > " << outfile;

        system(annotationFile.str().c_str());
        system(combineCommand.str().c_str());
        system(zipCommand.str().c_str());
        system(annotationCommand.str().c_str());
    }
    else
    {
        //" | awk -v OFS='\\t' '{print $1,$2,$3,$4,$106,$107}' > "
        stringstream annotationFile, combineCommand, zipCommand, annotationCommand;
        annotationFile << externals_loc << "tabix ../../data/whole_genome_SNVs.tsv.gz " << region << " | awk -v OFS='\\t' '{print $1,$2,$3,$4,$5,$6}' > " << rawFile;
        combineCommand << "cat tmp/CADDNames.txt " << rawFile << " > " << annoFile;
        zipCommand << externals_loc << "bgzip -f " << annoFile << " && tabix -f -p vcf " << annoFile << ".gz";
        annotationCommand << externals_loc << "bcftools merge -r " << region << " " << userFile << " " << backFile <<  " | " <<
                                   externals_loc << "bcftools annotate -a " << annoFile << ".gz -h ../../data/headerLines.txt -c CHROM,POS,REF,ALT,RawScore,PHRED - | " << 
                                   externals_loc << "bcftools query -f '%INFO/PHRED\\n' > " << outfile;

        //cout << "annotationFile: " << annotationFile.str() << endl << endl;
        //cout << "combineFile: " << combineCommand.str() << endl << endl;
        //cout << "zipFile: " << zipCommand.str() << endl << endl;
        cout << "annotationCommand: " << annotationCommand.str() << endl << endl;

        system(annotationFile.str().c_str());
        system(combineCommand.str().c_str());
        system(zipCommand.str().c_str());
        system(annotationCommand.str().c_str());
    }
    string line;
    //ifstream in(outfile);
    ifstream in("anno.txt");
    CADDWeights = gsl_vector_calloc(variantCount);
    for(int i = 0; i < variantCount; i++)
    {
        getline(in, line);
        //If CADD has no weight for the variant.
        if(line == ".")
        {
            //Really should be dealt with properly. This is just to say that no score means no weight.
            gsl_vector_set(CADDWeights, i, 0.001);
        }
        else
        {
            gsl_vector_set(CADDWeights, i, stod(line));
        }
    }


}

