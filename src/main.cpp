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
#include <vector>
#include <string.h>
#include <chrono>
#include <omp.h>
#include "genericBurdenTest.cpp"
#include "wsbt.cpp"
#include "input.cpp"
#include "output.cpp"
//#include "cast.cpp"
#include "skat.cpp"
#include "dataCollector.cpp"

using namespace std;

void handler(const char * reason, const char * file, int line, int gsl_errno)
{
    if (errno == GSL_EINVAL)
    {
        cout << "We out of bounds lads" << endl;
    }
}

int main(int argc, const char *argv[])
{
    string currentDir = "";
    string vcffilename = "";
    string vcfType = "";
    string backfilename = "";
    string phenofilename = "";
    string covfilename = "";
    string filename = "";
    string testType = ""; 
    string region = "";
    string regionFile = "";
    vector<string> geneList;
    bool geneBased = false;
    bool userBackgroundIncluded = true;
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lasttime = currentTime;

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
        //File parsing.
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
        if (strcmp(argv[i], "--d") == 0)
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
        //File option parsing.
        if (strcmp(argv[i], "--b") == 0 || strcmp(argv[i], "--back") == 0)
        {
            userBackgroundIncluded = false;
            if(argc > i + 1)
            {
                backfilename = argv[i+1];
                i++;
            }
            else
            {
                cout << "--back <backgroundFileName>" << endl;
            }
            //continue;
        }
        if (strcmp(argv[i], "--r") == 0 || strcmp(argv[i], "--region") == 0)
        {
            if(i + 1 < argc)
            {
                while(argv[i + 1][0] != '-')
                {
                    geneList.push_back(argv[i + 1]);
                    if(i + 2 < argc)
                    {
                        i++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            else
            {
                cout << "Did not provide a genetic region following the region comand." << endl;
                return 1;
            }
            continue;
        }
        if(strcmp(argv[i], "--R") == 0)
        {
            if(i + 1 < argc)
            {
                ifstream in(argv[i + 1]);
                string line;
                while(getline(in, line))
                {
                    if(line.find(':') != string::npos && line.find('-') != string::npos)
                    {
                        geneList.push_back(line);
                    }
                    //Manage genelist of names.
                    else
                    {
                        cout << "Regions in file " << argv[i + 1] << " should be coded as chromosome:start-stop" << endl;
                    }   
                }
                i++;
            }
            else
            {
                cout << "--R <RegionFileName.txt>" << endl;
            }
            continue;
        }
        //Test type parsing.
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
        //Test Parameter parsing.
        if (strcmp(argv[i], "--g") == 0)
        {
            geneBased = true;
            region = argv[i];
        }
        if (strcmp(argv[i], "--T") == 0)
        {
            if(argc > i)
            {
                omp_set_num_threads(stoi(argv[i+1]));
                i++;
            }
        }

    }

    //Run input on the given test and save results in Input
    if (testType == "wsbt")
    {
        if(geneList.size() > 0)
        {
            cout << "---Running WSBT in region mode---" << endl;

            size_t perm[geneList.size()];
            vector<string> genes(geneList.size());
            vector<double> pvalues(geneList.size());
            vector<double> permpvalues(geneList.size());
            vector<double> runTime(geneList.size());
            int index = 0;
            int skipped = 0;
            #pragma omp parallel for schedule(dynamic)
            for (int i = index; i < geneList.size(); i++)
            {
                //Thread safe method of iterating to the correct gene. Not elegant though.
                vector<string>::iterator iter = geneList.begin();
                advance(iter, i);
                auto runStartTime = chrono::high_resolution_clock::now();

                //Reads in genotype data from region.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, vcffilename, backfilename, *iter, testType, omp_get_thread_num());
                //Run the test.
                if(geneInput.getGslGenotype() == nullptr)
                {
                    cout << "Gene in region " << *iter << " has no variants in data sets. --Skipping--" << endl;
                    
                    genes[i] = *iter;
                    pvalues[i] = -1;
                    permpvalues[i] = -1;
                    skipped++; 
                }
                else
                {
                    wsbt test = wsbt(geneInput.getGslGenotype(), geneInput.getGslGenotype()->size2 - 2504, *iter);
                    genes[i] = *iter;
                    pvalues[i] = test.getPvalue();
                    permpvalues[i] = test.getPermPvalue();
                }
                
                //Check test speed.
                auto runCurrentTime = std::chrono::high_resolution_clock::now();
                runTime[i] = std::chrono::duration_cast<std::chrono::milliseconds>(runCurrentTime - runStartTime).count() / 60000.0;
            }
            gsl_sort_index(perm, pvalues.data(), 1, pvalues.size());
            ofstream out;
            out.open("wsbtResults.txt");
            for (int i = 0; i < pvalues.size(); i++)
            {
                //perm contains the sorted order.
                out << "Pvalue for gene " << genes[perm[i]] << " is " << pvalues[perm[i]] << ", " << permpvalues[perm[i]] << endl;
            }
            out.close();
            double fisherStat = 0;
            int geneCount = 0;
            for (int i = 0; i < pvalues.size(); i++)
            {
                if (pvalues[perm[i]] > 0)
                {
                    fisherStat += -2 * log(pvalues[perm[i]]);
                    geneCount++;
                }
                else
                {
                    cout << "Pvalue for gene " << genes[perm[i]] << " excluded from fisher product test because pvalue = " << pvalues[perm[i]] << endl;
                    if(pvalues[perm[i]] == 0)
                    {
                        skipped++;
                    }
                }
            }
            double fisherPvalue = gsl_cdf_chisq_Q(fisherStat, 2 * geneCount);
            cout << "Fisher product test statistic for " << geneList.size() - skipped << " gene(s) is " << fisherStat << " with pvalue " << fisherPvalue << endl;
        }
        else if (geneBased)
        {
            readInput result = readInput(currentDir, testType, vcfType, vcffilename, region, phenofilename, covfilename);
            currentTime = std::chrono::high_resolution_clock::now();
            cout << endl;
            cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
            cout << endl;
            lasttime = currentTime;

            map<string, string> regions = result.getRegions();
           // vector<geneId> *geneInfo = result.getGenes();

            size_t perm[regions.size()];
            vector<string> genes(regions.size());
            vector<double> pvalues(regions.size());
            vector<double> permpvalues(regions.size());

            gsl_set_error_handler_off();
            //gsl_set_error_handler(&handler);
            int index = 0;
            int skipped = 0;

            /*
            for(map<string, string>::iterator iter = regions.begin(); iter != regions.end(); iter++, index++)
            {
                if(iter->first[0] == 'P')
                {
                    break;
                }
            }
            */

            #pragma omp parallel for schedule(dynamic)
            for (int i = index; i < regions.size(); i++)
            {
                //Thread safe method of iterating to the correct gene. Not elegant though.
                map<string, string>::iterator iter = regions.begin();
                advance(iter, i);

                //Reads in genotype data from region.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, vcffilename, backfilename, iter->second, testType, omp_get_thread_num());
                //Run the test.
                wsbt test = wsbt(geneInput.getGslGenotype(), result.getCaseCount(), iter->first);
                genes[i] = iter->first;
                pvalues[i] = test.getPvalue();
                permpvalues[i] = test.getPermPvalue();
            }
            gsl_sort_index(perm, pvalues.data(), 1, pvalues.size());
            ofstream out("wsbtResults.csv");
            out << "Gene,Pvalue,Perm-Pvalue" << endl;
            for (int i = 0; i < pvalues.size(); i++)
            {
                //perm contains the sorted order.
                out << genes[perm[i]] << "," << pvalues[perm[i]] << "," << permpvalues[perm[i]] << endl;
            }
            out.close();
            if(out.is_open())
            {
                cout << "we aint done yet." << endl;
            }
            double fisherStat = 0;
            int geneCount = 0;
            for (int i = 0; i < pvalues.size(); i++)
            {
                if (pvalues[perm[i]] > 0)
                {
                    fisherStat += -2 * log(pvalues[perm[i]]);
                    geneCount++;
                }
                else
                {
                    cout << "Pvalue for gene " << genes[perm[i]] << " excluded from fisher product test because pvalue = " << pvalues[perm[i]] << endl;
                    if(pvalues[perm[i]] == 0)
                    {
                        skipped++;
                    }
                }
            }
            double fisherPvalue = gsl_cdf_chisq_Q(fisherStat, 2 * geneCount);
            cout << "Fisher product test statistic on " << pvalues.size() - skipped <<  " genes is " << fisherStat << " with pvalue " << fisherPvalue << endl;
        }
        else
        {
            /*
            readInput result(currentDir, testType, vcfType, vcffilename, region, phenofilename, covfilename);

            currentTime = std::chrono::high_resolution_clock::now();
            cout << endl;
            cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
            cout << endl;
            lasttime = currentTime;

            wsbt test = wsbt(result.getGslGenotype(), result.getCaseCount(), "All");
            currentTime = std::chrono::high_resolution_clock::now();
            cout << "WSBT Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0 << " minutes." << endl;
            lasttime = currentTime;
            */
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 60000.0 << " minutes." << endl;

        //Now manage output of results for user.
    }
    else if (testType == "burden")
    {
        /*
        readInput result(currentDir, testType, vcfType, vcffilename, region, phenofilename, covfilename);

        currentTime = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
        cout << endl;
        lasttime = currentTime;

        gsl_vector *pheno = gsl_vector_alloc(result.getMaf()->size);
        genericBurdenTest test = genericBurdenTest(result.getGslGenotype(), result.getMaf(), pheno);
        cout << "Variant weights: ";
        for (int i = 0; i < test.getWeights()->size; i++)
        {
            cout << gsl_vector_get(test.getWeights(), i) << " ";
        }
        cout << endl;
        */
    }
    else if (testType == "cast")
    {
        cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "skat")
    {
        readInput result(currentDir, testType, vcfType, vcffilename, region, phenofilename, covfilename);

        currentTime = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
        cout << endl;
        lasttime = currentTime;

        currentTime = std::chrono::high_resolution_clock::now();
        lasttime = currentTime;
        if (geneBased)
        {
            map<string, string> regions = result.getRegions();
            size_t perm[regions.size()];
            vector<string> genes(regions.size());
            vector<double> pvalues(regions.size());
            vector<double> testStats(regions.size());
            vector<double> rvTestStats(regions.size());
            vector<double> runTime(regions.size());
            
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < regions.size(); i++)
            {
                //Iterate to run the correct gene.
                map<string, string>::iterator iter = regions.begin();
                advance(iter, i);
                cout << "Running SKAT test on gene " << iter->first << endl;
                //Read in genotpe data from file for this gene.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, vcffilename, backfilename, iter->second, testType, omp_get_thread_num());
                //Run Test
                auto runStartTime = std::chrono::high_resolution_clock::now();
                skat test = skat(geneInput.getGslGenotype(), geneInput.getMaf(), result.getCovariates(), result.getPheno());
                auto runCurrentTime = std::chrono::high_resolution_clock::now();
                //Output and cleanup.
                genes[i] = iter->first;
                pvalues[i] = test.getPvalue();
                testStats[i] = test.getQ();
                rvTestStats[i] = test.getRvQ();
                runTime[i] = (std::chrono::duration_cast<std::chrono::milliseconds>(runCurrentTime - runStartTime).count() / 60000.0);
                cout << "SKAT Took " << runTime[i] << " minutes." << endl;
                cout << endl;
            }

            gsl_sort_index(perm, pvalues.data(), 1, pvalues.size());
            ofstream outFile;
            string output = vcffilename + ".SKAT.results";
            outFile.open(output);
            for (int i = 0; i < genes.size() ; i++)
            {
                outFile << "Gene: " << genes[perm[i]] << "\tpvalue:" << pvalues[perm[i]] << "\tQ:" << testStats[perm[i]] << "\truntime:" << runTime[perm[i]] << " minutes." << endl;
                cout << "Gene: " << genes[perm[i]] << "\tQ:" << testStats[perm[i]] << "\tpvalue:" << pvalues[perm[i]] << "\trvQ:" << rvTestStats[perm[i]] << endl;
            }
            
        }
        else
        {
            /*
            //Covariates and Phenotype not yet implemented. We will have to decide how thats done.
            skat test = skat(result.getGslGenotype(), result.getMaf(), result.getCovariates(), result.getPheno());
            currentTime = std::chrono::high_resolution_clock::now();
            cout << "SKAT Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lasttime).count() / 60000.0 << " minutes." << endl;
            lasttime = currentTime;
            */
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