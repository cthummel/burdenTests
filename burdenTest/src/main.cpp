//
//  main.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

//
#include <iostream>
#include <fstream>
#include <string.h>
#include <chrono>
#include "genericBurdenTest.cpp"
#include "wsbt.cpp"
#include "input.cpp"
#include "output.cpp"
//#include "cast.cpp"
#include "skat.cpp"

using namespace std;

int main(int argc, const char *argv[])
{
    string currentDir, vcffilename, vcfType, phenofilename, covfilename, filename, testType, region;
    bool geneBased = false;
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lasttime = startTime;

    currentDir = argv[0];
    if (currentDir[0] == '.')
    {
        currentDir = "";
    }
    else
    {
        currentDir = currentDir.substr(0, currentDir.length() - 5);
    }

    //Check for proper argument formatting and filetypes.
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-vcf") == 0)
        {
            vcfType = argv[i];
            if (argc > i)
            {
                vcffilename = argv[i + 1];
                if (vcffilename.substr(vcffilename.length() - 4) != ".vcf")
                {
                    //They input a filename that isnt a .vcf
                    cout << "burdenTest -vcf <filename.vcf>\n";
                    return 0;
                }
                i++;
            }
        }
        if (strcmp(argv[i], "-gzvcf") == 0)
        {
            vcfType = argv[i];
            if (argc > i)
            {
                vcffilename = argv[i + 1];
                if (vcffilename.substr(vcffilename.length() - 7) != ".vcf.gz")
                {
                    //They input a filename that isnt a .vcf.gz
                    cout << "burdenTest -gzvcf <filename.vcf.gz>\n";
                    return 0;
                }
                i++;
            }
        }
        if (strcmp(argv[i], "-pheno") == 0)
        {
            if (argc > i)
            {
                phenofilename = argv[i + 1];
                if (phenofilename.substr(phenofilename.length() - 4) != ".ped")
                {
                    //They input a filename that isnt a .pheno
                    cout << "burdenTest -pheno <filename.ped>\n";
                    return 0;
                }
                i++;
            }
        }
        if (strcmp(argv[i], "-cov") == 0)
        {
            if (argc > i)
            {
                covfilename = argv[i + 1];
                if (covfilename.substr(covfilename.length() - 4) != ".txt")
                {
                    //They input a filename that isnt a .txt
                    cout << "burdenTest [testtype] -cov <filename.txt>\n";
                    return 0;
                }
                i++;
            }
        }
        if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "-region") == 0)
        {
            if (strcmp(argv[i + 1], "all") == 0 || strcmp(argv[i + 1], "gene") == 0 || strcmp(argv[i + 1], "exact") == 0)
            {
                region = argv[i + 1];
                i++;
            }
            else
            {
                cout << "arguments error: incorrect region specified after indicator." << endl;
            }
        }
        if (strcmp(argv[i], "wsbt") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "burden") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "cast") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "skat") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "skato") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "-h") == 0)
        {
            testType = argv[i];
        }
        if (strcmp(argv[i], "-g") == 0)
        {
            geneBased = true;
            region = argv[i];
        }
    }

    //cout << "correct filename: " << vcffilename << "\n"; result;
    if (vcffilename != "" && vcfType != "")
    {
    }

    //Run input on the given test and save results in Input
    readInput result(currentDir, testType, vcfType, vcffilename, region, phenofilename, covfilename);

    currentTime = std::chrono::high_resolution_clock::now();
    lasttime = currentTime;
    cout << endl;
    cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
    cout << endl;

    if (testType == "wsbt")
    {
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "Before WSBT Test." << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() << endl;
        lasttime = currentTime;
        if (geneBased)
        {
            map<string, gsl_matrix *> subsets = result.getGeneSubsets();
            vector<double> pvalues;
            vector<double> permpvalues;
            vector<double> runTime;
            for (map<string, gsl_matrix *>::iterator iter = subsets.begin(); iter != subsets.end(); iter++)
            {
                cout << "Running WSBT test on gene " << iter->first << endl;
                wsbt test = wsbt(iter->second, result.getCaseCount(), result.getMaf(iter->first));
                pvalues.push_back(test.getPvalue());
                permpvalues.push_back(test.getPermPvalue());
                currentTime = std::chrono::high_resolution_clock::now();
                cout << "WSBT on " << iter->first << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 1000.0 << " minutes." << endl;
                runTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0);
                lasttime = currentTime;
            }
            int i = 0;
            for (map<string, gsl_matrix *>::iterator iter = subsets.begin(); iter != subsets.end(); i++, iter++)
            {
                cout << "Pvalue for gene " << iter->first << " is " << pvalues[i] << ", " << permpvalues[i];
                cout << " with a test run time of " << runTime[i] << " seconds." << endl;
            }
            double fisherStat = 0;
            int geneCount = 0;
            for (int i = 0; i < pvalues.size(); i++)
            {
                if(!gsl_isnan(pvalues[i]))
                {
                    fisherStat += -2 * log(pvalues[i]);
                    geneCount++;
                }
            }
            double fisherPvalue = gsl_cdf_chisq_P(fisherStat, 2 * geneCount);
            cout << "Fisher product test statistic is " << fisherStat << " with pvalue " << fisherPvalue << endl;
        }
        else
        {
            wsbt test = wsbt(result.getGslGenotype(), result.getCaseCount(), result.getMaf());
            currentTime = std::chrono::high_resolution_clock::now();
            cout << "WSBT Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0 << " minutes." << endl;
            lasttime = currentTime;
        }
        //File will be a vcf by this point whether it was zipped before or not.
        //writeOutput output(vcffilename, testType, test.getWeights());

        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 60000.0 << " minutes." << endl;
    }
    else if (testType == "burden")
    {
        gsl_vector *pheno = gsl_vector_alloc(result.getMaf()->size);
        genericBurdenTest test = genericBurdenTest(result.getGslGenotype(), result.getMaf(), pheno);
        cout << "Variant weights: ";
        for (int i = 0; i < test.getWeights()->size; i++)
        {
            cout << gsl_vector_get(test.getWeights(), i) << " ";
        }
        cout << endl;
    }
    else if (testType == "cast")
    {
        cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "skat")
    {
        currentTime = std::chrono::high_resolution_clock::now();
        //cout << "Before SKAT Test." << endl;
        lasttime = currentTime;
        if (geneBased)
        {
            map<string, gsl_matrix *> subsets = result.getGeneSubsets();
            vector<double> pvalues;
            vector<double> testStats;
            vector<double> rvTestStats;
            vector<double> runTime;
            for (map<string, gsl_matrix *>::iterator iter = subsets.begin(); iter != subsets.end(); iter++)
            {
                cout << "Running SKAT test on gene " << iter->first << endl;
                skat *test = new skat(iter->second, result.getMaf(iter->first), result.getCovariates(), result.getPheno());
                pvalues.push_back(test->getPvalue());
                testStats.push_back(test->getQ());
                rvTestStats.push_back(test->getRvQ());
                currentTime = std::chrono::high_resolution_clock::now();
                runTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0);
                cout << "SKAT Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0 << " minutes." << endl;
                cout << endl;
                lasttime = currentTime;
            }
            ofstream outFile;
            string output = vcffilename + ".SKAT.results";
            outFile.open(output);
            int i = 0;
            for (map<string, gsl_matrix *>::iterator iter = subsets.begin(); iter != subsets.end(); i++, iter++)
            {
                outFile << "Gene: " << iter->first << "\tpvalue:" << pvalues[i] << "\tQ:" << testStats[i] << "\truntime:" << runTime[i] << " minutes." << endl;
                //cout << "Gene: " << iter->first << "\tpvalue:" << pvalues[i] << "\tQ:" << testStats[i] << "\truntime:" << runTime[i] << " minutes." << endl;
                cout << "Gene: " << iter->first << "\tQ:" << testStats[i] << "\tpvalue:" << pvalues[i] << "\trvQ:" << rvTestStats[i] << endl;
            }
            
        }
        else
        {
            //Covariates and Phenotype not yet implemented. We will have to decide how thats done.
            skat test = skat(result.getGslGenotype(), result.getMaf(), result.getCovariates(), result.getPheno());
            currentTime = std::chrono::high_resolution_clock::now();
            cout << "SKAT Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0 << " minutes." << endl;
            lasttime = currentTime;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 60000.0 << " minutes." << endl;
    }
    else if (testType == "skato")
    {
        cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "-h")
    {
        cout << "Run tests with ./test [testname] [filetype filename] " << endl;

        cout << "That feel when even help isnt really implemented..." << endl;
    }
    else
    {
        cout << "Please specify desired test on command line. For help use -h" << endl;
    }

    //cout << "Hello, World!\n";
    return 0;
}
