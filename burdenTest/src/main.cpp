//
//  main.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//


#include <iostream>
#include <string.h>
#include <chrono>
#include "genericBurdenTest.cpp"
#include "wsbt.cpp"
#include "input.cpp"
#include "output.cpp"
#include "cast.cpp"
#include "skat.cpp"

using namespace std;

int main(int argc, const char * argv[])
{
    string vcffilename, vcfType, phenofilename, covfilename, filename, testType, region;
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lasttime = startTime;
    
    /*
    string line;
    cout << "Preloading background complete.";
    while(getline(cin, line))
    {

    }
    */

    //Check for proper argument formatting and filetypes.
    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-vcf") == 0)
        {
            vcfType = argv[i];
            if(argc > i)
            {
                vcffilename = argv[i+1];
                if (vcffilename.substr(vcffilename.length() - 4) != ".vcf")
                {
                    //They input a filename that isnt a .vcf
                    cout << "burdenTest -vcf <filename.vcf>\n";
                    return 0;
                }
                i++;
            }
        }
        if(strcmp(argv[i], "-gzvcf") == 0)
        {
            vcfType = argv[i];
            if(argc > i)
            {
                vcffilename = argv[i+1];
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
                phenofilename = argv[i+1];
                if (phenofilename.substr(phenofilename.length() - 6) != ".pheno")
                {
                    //They input a filename that isnt a .pheno
                    cout << "burdenTest -pheno <filename.pheno>\n";
                    return 0;
                }
                i++;
            }
        }
        if (strcmp(argv[i], "-cov") == 0)
        {
            if (argc > i)
            {
                covfilename = argv[i+1];
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
            if( strcmp(argv[i+1], "all") == 0 || strcmp(argv[i+1], "gene") == 0 || strcmp(argv[i+1], "exact") == 0)
            {
                region = argv[i+1];
                i++;
            }
            else
            {
                cout << "arguments error: incorrect region specified after indicator." << endl;
            }
        }
        if(strcmp(argv[i], "wsbt") == 0)
        {
            testType = argv[i];
        }
        if(strcmp(argv[i], "burden") == 0)
        {
            testType = argv[i];
        }
        if(strcmp(argv[i], "cast") == 0)
        {
            testType = argv[i];
        }
        if(strcmp(argv[i], "skat") == 0)
        {
            testType = argv[i];
        }
        if(strcmp(argv[i], "skato") == 0)
        {
            testType = argv[i];
        }
        if(strcmp(argv[i], "-h") == 0)
        {
            testType = argv[i];
        }
    }
    
    //cout << "correct filename: " << vcffilename << "\n"; result;
    if(vcffilename != "" && vcfType != "")
    {
        
    }
    
    //Run input on the given test and save results in Input
    readInput result(testType, vcfType, vcffilename, region, phenofilename, covfilename);
    
    
    currentTime = std::chrono::high_resolution_clock::now();
    lasttime = currentTime;
    cout << endl;
    cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-startTime).count() / 60000.0 << " minutes."<< endl;
    
    
    if(testType == "wsbt")
    {
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "Before WSBT Test." << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count() << endl;
        lasttime = currentTime;
        
        wsbt test = wsbt(result.getGslGenotype(), result.getCaseCount(), result.getMaf());
        
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "WSBT Took "  << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count()/60000.0 << " minutes."<< endl;
        lasttime = currentTime;
        
        //File will be a vcf by this point whether it was zipped before or not.
        writeOutput output(vcffilename, testType, test.getWeights());
        
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count()/60000.0 << " minutes." << endl;
    }
    else if (testType == "burden")
    {
        gsl_vector* pheno = gsl_vector_alloc(result.getMaf()->size);
        genericBurdenTest test = genericBurdenTest(result.getGslGenotype(), result.getMaf(), pheno);
        cout << "Variant weights: ";
        for(int i = 0; i < test.getWeights()->size; i++)
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
        cout << "Before SKAT Test." << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count() << endl;
        lasttime = currentTime;
        //Covariates and Phenotype not yet implemented. We will have to decide how thats done.
        skat test = skat(result.getGslGenotype(), result.getMaf(), result.getCovariates(), result.getPheno());
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "SKAT Took "  << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count()/60000.0 << " minutes."<< endl;
        lasttime = currentTime;
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count()/60000.0 << " minutes." << endl;
        cout << testType << " is not yet implemented." << endl;
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
