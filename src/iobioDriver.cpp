#include <string.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <gsl/gsl_matrix.h>
#include "dataCollector.cpp"
#include "wsbt.cpp"
#include "skat.cpp"


using namespace std;

gsl_matrix* readGeno(int subjectCount, int variantCount, string genoFile, string testType)
{
    gsl_matrix *genotypeGslMatrix = gsl_matrix_alloc(variantCount, subjectCount);
    string line = "";

    if(genoFile == "-")
    {
        for (int i = 0; getline(std::cin, line); i++)
        {
            string::iterator it = line.begin();
            bool missingData = false;
            int missingCount = 0;
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
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
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
            if (missingData && testType != "wsbt")
            {
                int alleleCount = 0;
                int denominator = 0;
                for (int j = 0; j < subjectCount; j++)
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
                for (int j = 0; j < subjectCount; j++)
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
    else
    {
        ifstream in(genoFile);
        for (int i = 0; getline(in, line); i++)
        {
            string::iterator it = line.begin();
            bool missingData = false;
            int missingCount = 0;
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
                    gsl_matrix_set(genotypeGslMatrix, i, j, -1);
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
            if (missingData && testType != "wsbt")
            {
                int alleleCount = 0;
                int denominator = 0;
                for (int j = 0; j < subjectCount; j++)
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
                for (int j = 0; j < subjectCount; j++)
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

    return genotypeGslMatrix;
}

gsl_vector *readMaf(int variantCount, string mafFile)
{
    gsl_vector *maf = gsl_vector_alloc(variantCount);
    string line = "";
    ifstream in(mafFile);
    for (int i = 0; getline(in, line); i++)
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
    return maf;
}

gsl_matrix* readCov(int variantCount, string covFile)
{
    gsl_matrix* cov = nullptr;
    ifstream in;
    if(strcmp(covFile.c_str(), "") == 0)
    {
        cov = gsl_matrix_alloc(variantCount, 1);
        gsl_matrix_set_all(cov, 1);
    }
    else
    {

    }

    return cov;
}

gsl_vector* readPheno(int subjectCount, int caseCount, string phenoFile)
{
    gsl_vector *pheno = gsl_vector_alloc(subjectCount);
    ifstream in;
    if(strcmp(phenoFile.c_str(), "") == 0)
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
    return pheno;
}

int main(int argc, char *argv[])
{
    int caseCount = 0;
    int variantCount = 0;
    int subjectCount = 0;
    bool testMode = false;
    string geneName = "";
    string testType = "";
    string dataFile = "";
    string mafFile = "";
    string phenoFile = "";
    string covFile = "";

    for (int i = 1; i < argc; i += 2)
    {
        if (strcmp(argv[i], "-c") == 0)
        {
            caseCount = stoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            variantCount = stoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            subjectCount = stoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-t") == 0)
        {
            testType = argv[i + 1];
        }
        else if (strcmp(argv[i], "-d") == 0)
        {
            dataFile = argv[i + 1];
        }
        else if (strcmp(argv[i], "-m") == 0)
        {
            mafFile = argv[i + 1];
        }
        else if (strcmp(argv[i], "-p") == 0)
        {
            phenoFile = argv[i + 1];
        }
        else if (strcmp(argv[i], "-c") == 0)
        {
            covFile = argv[i + 1];
        }
        else if (strcmp(argv[i], "-n") == 0)
        {
            geneName = argv[i + 1];
        }
        else if (strcmp(argv[i], "-testMode") == 0)
        {
            testMode = true;
        }
        else
        {
            cout << "Thanks for using the Iobio Burden Test Driver." << endl;
            cout << "./driver -c [caseCount] -v [variantCount] -s [subjectCount] -t [testType] -d [genotypeDataFileName] -n [geneName]" << endl;
            cout << "testType: wsbt, skat" << endl;
            return 0;
        }
    }

    if (testMode)
    {
        return 0;
    }
    
    if (strcmp(testType.c_str(), "wsbt") == 0)
    {
        if (caseCount == 0 || variantCount == 0 || subjectCount == 0)
        {
            cout << "Missing case count or variant count or subject count" << endl;
            return 2;
        }
        cout << "Subject Count: " << subjectCount << endl;
        cout << "Variant Count: " << variantCount << endl;
        cout << "Case Count: " << caseCount << endl;
        cout << "DataFile: " << dataFile << endl;
        cout << endl;
        gsl_matrix *data = readGeno(subjectCount, variantCount, dataFile, testType);
        auto runStartTime = chrono::high_resolution_clock::now();
        wsbt test = wsbt(data, caseCount, geneName, true);
        test.driverOutput();
        auto runEndTime = std::chrono::high_resolution_clock::now();
        cout << "Total runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(runEndTime - runStartTime).count() / 1000.0 << " seconds." << endl;

        gsl_matrix_free(data);
    }
    else if (strcmp(testType.c_str(), "skat") == 0)
    {
        cout << "Subject Count: " << subjectCount << endl;
        cout << "Variant Count: " << variantCount << endl;
        cout << "Case Count: " << caseCount << endl;
        cout << "DataFile: " << dataFile << endl;
        cout << "MafFile: " << mafFile << endl;
        cout << endl;
        gsl_matrix *data = readGeno(subjectCount, variantCount, dataFile, testType);
        gsl_vector *maf = readMaf(variantCount, mafFile);
        gsl_matrix *cov = readCov(variantCount, covFile);
        gsl_vector *pheno = readPheno(subjectCount, caseCount, phenoFile);
        auto runStartTime = chrono::high_resolution_clock::now();
        skat test = skat(data, maf, cov, pheno);
        auto runEndTime = std::chrono::high_resolution_clock::now();

        gsl_matrix_free(data);
        gsl_matrix_free(cov);
        gsl_vector_free(maf);
        gsl_vector_free(pheno);
    }
    else if (strcmp(testType.c_str(), "skato") == 0)
    {
    }
    else if (strcmp(testType.c_str(), "") == 0)
    {
        cout << "Type of test to be run was not specified." << endl;
        return 1;
    }
    else
    {
        cout << "Type of test to be run was specified incorrectly." << endl;
        return 2;
    }


    return 0;
}