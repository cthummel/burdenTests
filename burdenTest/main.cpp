//
//  main.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include <iostream>
#include "genericBurdenTest.cpp"
#include "wsbt.cpp"
#include "input.cpp"
#include "cast.cpp"
#include "skat.cpp"

using namespace std;

int main(int argc, const char * argv[])
{
    std::string vcffilename, vcfType, phenofilename, filename, testType;
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lasttime = startTime;
    //int testType = 0;
    
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
                    //They input a filename that isnt a vcf
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
                    //They input a filename that isnt a vcf
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
    readInput result(testType, vcfType, vcffilename, vcffilename, phenofilename);
    currentTime = std::chrono::high_resolution_clock::now();
    lasttime = currentTime;
    cout << endl;
    cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-startTime).count() / 60000.0 << " minutes."<< endl;
    vector<double> pheno = vector<double>(result.getMaf()->size);
    
    if(testType == "wsbt")
    {
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "Before WSBT Test." << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count() << endl;
        lasttime = currentTime;
        wsbt test = wsbt(result.getGslGenotype(), result.getCaseCount(), result.getMaf());
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "Wsbt Took "  << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count()/60000.0 << " minutes."<< endl;
        lasttime = currentTime;
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count()/60000.0 << " minutes." << endl;
    }
    else if (testType == "burden")
    {
        genericBurdenTest test = genericBurdenTest(result.getGenotype(), result.getMaf(), pheno);
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
        
        //Covariates not yet implemented. We will have to decide how thats done.
        skat test = skat(result.getGslGenotype(), result.getMaf(), gsl_vector_alloc(10));
        currentTime = std::chrono::high_resolution_clock::now();
        cout << "SKAT Took "  << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime-lasttime).count()/60000.0 << " minutes."<< endl;
        lasttime = currentTime;
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count()/60000.0 << " minutes." << endl;
        //cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "skato")
    {
        cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "-h")
    {
        cout << "That feel when even help isnt implemented..." << endl;
    }
    else
    {
        cout << "Please specify desired test on command line." << endl;
    }
    
    //cout << "Hello, World!\n";
    return 0;
}
