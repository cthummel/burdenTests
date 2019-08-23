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
#include "skato.cpp"
#include "dataCollector.cpp"

using namespace std;

void handler(const char * reason, const char * file, int line, int gsl_errno)
{
    if (errno == GSL_EINVAL)
    {
        cout << "We out of bounds lads" << endl;
    }
}

//Returns true if all user data is missing. False otherwise.
bool checkMissingData(gsl_matrix_short *data, int caseCount)
{
    if(data->size1 == 0)
    {
        return true;
    }
    return false;
    for (int i = 0; i < data->size1; i++)
    {
        for (int j = 0; j < caseCount; j++)
        {
            if (gsl_matrix_short_get(data, i, j) > -1)
            {
                return false;
            }
        }
    }
    return true;
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
    string region = "--g";
    string regionFile = "";
    string outputFileName = "output";
    string variantRegion = "transcript";
    vector<string> geneList;
    bool agnosticGeneRun = true;
    bool userBackgroundIncluded = true;
    bool exactPvalueCalculation = true;
    bool useCADDWeights = false;
    auto startTime = chrono::high_resolution_clock::now();
    auto currentTime = startTime;
    auto lasttime = currentTime;

    omp_lock_t outputLock;
    omp_init_lock(&outputLock);

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
            if(argc > i)
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
            if(argc > i)
            {
                agnosticGeneRun = false;
                while(argv[i + 1][0] != '-')
                {
                    geneList.push_back(argv[i + 1]);
                    if(argc > i + 1)
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
            if(argc > i)
            {
                agnosticGeneRun = false;
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
                return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--o") == 0 || strcmp(argv[i], "--output") == 0) 
        {
            if (argc > i)
            {
                outputFileName = argv[i + 1];
                if(outputFileName.length() > 4)
                {
                    if (outputFileName.substr(outputFileName.length() - 4) == ".tsv")
                    {
                        outputFileName = outputFileName.substr(0, outputFileName.length() - 4);
                    }
                }
                
                i++;
            }
        }
        if (strcmp(argv[i], "--e") == 0 || strcmp(argv[i], "--exon") == 0) 
        {
            variantRegion = "exon";
            continue;
        }
        if (strcmp(argv[i], "--CADD") == 0)
        {
            useCADDWeights = true;
            string annotationFile = "echo '#Chr\\tPos\\tRef\\tAlt\\tRawScore\\tPHRED' > tmp/CADDNames.txt";
            system(annotationFile.c_str()); 
            continue;
        } 
        //Test type parsing.
        if (strcmp(argv[i], "wsbt") == 0)
        {
            testType = argv[i];
            continue;
        }
        if (strcmp(argv[i], "burden") == 0)
        {
            testType = argv[i];
            continue;
        }
        if (strcmp(argv[i], "cast") == 0)
        {
            testType = argv[i];
            continue;
        }
        if (strcmp(argv[i], "skat") == 0)
        {
            testType = argv[i];
            continue;
        }
        if (strcmp(argv[i], "skato") == 0)
        {
            testType = argv[i];
            continue;
        }
        if (strcmp(argv[i], "-h") == 0)
        {
            testType = argv[i];
            continue;
        }
        //Test Parameter parsing.
        if (strcmp(argv[i], "--g") == 0)
        {
            agnosticGeneRun = true;
            region = argv[i];
            continue;
        }
        if (strcmp(argv[i], "--p") == 0)
        {
            exactPvalueCalculation = false;
            continue;
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
            vector<int> effectSizes(geneList.size());
            vector<int> UserUniqueVariantCounts(geneList.size());
            vector<int> variantCounts(geneList.size());
            vector<double> pvalues(geneList.size());
            vector<double> scores(geneList.size());
            vector<double> testStats(geneList.size());
            vector<double> runTime(geneList.size());
            
            //int index = 0;
            int skipped = 0;

            ofstream out(outputFileName + ".tsv");
            out << "Region\tScore\tTestStat\tPvalue" << endl;

            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < geneList.size(); i++)
            {
                //Thread safe method of iterating to the correct gene. Not elegant though.
                vector<string>::iterator iter = geneList.begin();
                advance(iter, i);
                auto runStartTime = chrono::high_resolution_clock::now();

                //Reads in genotype data from region.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, useCADDWeights, vcffilename, backfilename, *iter, testType, 1, omp_get_thread_num());
                //Run the test.
                if(geneInput.getGslGenotype() == nullptr)
                {
                    cout << "Gene in region " << *iter << " has no variants in data sets. --Skipping--" << endl;
                    genes[i] = *iter;
                    pvalues[i] = -1;
                    scores[i] = -1;
                    testStats[i] = -9999;
                    skipped++; 
                }
                else
                {
                    gsl_vector* CADDWeights = nullptr;
                    if(useCADDWeights)
                    {
                        CADDWeights = geneInput.getCADDWeights();
                    }
                    wsbt test = wsbt(geneInput.getShortGslGenotype(), CADDWeights, geneInput.getShortGslGenotype()->size2 - 2504, *iter, exactPvalueCalculation);
                    test.driverOutput();
                    genes[i] = *iter;
                    pvalues[i] = test.getPvalue();
                    scores[i] = test.getScores()[0];
                    testStats[i] = test.getTestStat();
                }
                
                omp_set_lock(&outputLock);
                out << genes[i] << "\t" << scores[i] << "\t" << testStats[i] << "\t" << pvalues[i] << endl;
                omp_unset_lock(&outputLock);

                //Check test speed.
                auto runCurrentTime = std::chrono::high_resolution_clock::now();
                runTime[i] = std::chrono::duration_cast<std::chrono::milliseconds>(runCurrentTime - runStartTime).count() / 60000.0;
            }
            out.close();

            gsl_sort_index(perm, pvalues.data(), 1, pvalues.size());
            out.open("sorted_" + outputFileName + ".tsv");
            out << "Region\tScore\tTestStat\tPvalue\tEffectSize\tCaseUniqueVariants\tTotalVariants" << endl;
            for (int i = 0; i < pvalues.size(); i++)
            {
                //perm contains the sorted order.
                out << genes[perm[i]] << "\t" << scores[perm[i]] << "\t" << testStats[perm[i]] << "\t" << pvalues[perm[i]] 
                << "\t" << effectSizes[perm[i]] << "\t" << UserUniqueVariantCounts[perm[i]] << "\t" << variantCounts[perm[i]]
                << endl;
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
        else if (agnosticGeneRun)
        {
            cout << "Starting gene agnostic run of the Weighted Sums Burden Test." << endl;
            readInput result = readInput(currentDir, testType, variantRegion, vcffilename, region, phenofilename, covfilename);
            currentTime = std::chrono::high_resolution_clock::now();
            cout << endl;
            cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
            cout << endl;
            lasttime = currentTime;

            map<string, string> regions = result.getRegions();

            size_t perm[regions.size()];
            vector<string> genes(regions.size());
            vector<string> locations(regions.size());
            vector<int> effectSizes(regions.size());
            vector<int> UserUniqueVariantCounts(regions.size());
            vector<int> BackUniqueVariantCounts(regions.size());
            vector<int> variantCounts(regions.size());
            vector<double> pvalues(regions.size());
            vector<vector<double>> scores(regions.size(), vector<double>(result.getCaseCount()));
            vector<double> testStats(regions.size());

            //int index = 0;
            int skipped = 0;

            ofstream out(outputFileName + ".tsv");
            out << "Gene\tRegion\t";
            for(int i = 0; i < result.getCaseCount(); i++)
            {
                out << result.getSampleNames()->at(i) << "_Score\t";
            }
            out << "TestStat\tPvalue\tRegionSize\tCaseUniqueVariants\tBackgroundUniqueVariants\tTotalVariants" << endl;

            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < regions.size(); i++)
            {
                //Thread safe method of iterating to the correct gene. Not elegant though.
                map<string, string>::iterator iter = regions.begin();
                advance(iter, i);

                effectSizes[i] = stoi(iter->second.substr(iter->second.find("-")+1)) - stoi(iter->second.substr(iter->second.find(":")+1, iter->second.find("-")));

                //Reads in genotype data from region.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, useCADDWeights, vcffilename, backfilename, iter->second, testType, result.getCaseCount(), omp_get_thread_num());
                variantCounts[i] = geneInput.getShortGslGenotype()->size1;
                UserUniqueVariantCounts[i] = geneInput.getCaseUniqueVariantCount();
                BackUniqueVariantCounts[i] = geneInput.getBackgroundUniqueVariantCount();

                if(checkMissingData(geneInput.getShortGslGenotype(), result.getCaseCount()))
                {
                    cout << "User data for " << iter->first << " contains no genotype data. Skipping test.\n" << endl;
                    genes[i] = iter->first;
                    locations[i] = iter->second;
                    pvalues[i] = -1;
                    for(int j = 0; j < result.getCaseCount(); j++)
                    {
                        scores[i][j] = -9999;
                    }
                    testStats[i] = -9999;
                    if(outputFileName == "")
                    {
                        omp_set_lock(&outputLock);
                        cout << genes[i] << "\t" << locations[i] << "\t";
                        for (int j = 0; j < result.getCaseCount(); j++)
                        {
                            cout << scores[i][j] << "\t";
                        }
                        cout << testStats[i] << "\t" << pvalues[i] << "\t" << effectSizes[i] << "\t" << UserUniqueVariantCounts[i] << "\t"
                             << BackUniqueVariantCounts[i] << "\t" << variantCounts[i] << endl;
                        omp_unset_lock(&outputLock);
                    }
                    else
                    {
                        omp_set_lock(&outputLock);
                        out << genes[i] << "\t" << locations[i] << "\t";
                        for (int j = 0; j < result.getCaseCount(); j++)
                        {
                            out << scores[i][j] << "\t";
                        }
                        out << testStats[i] << "\t" << pvalues[i] << "\t" << effectSizes[i] << "\t" << UserUniqueVariantCounts[i] << "\t"
                            << BackUniqueVariantCounts[i] << "\t" << variantCounts[i] << endl;
                        omp_unset_lock(&outputLock);
                    }
                    
                }
                else
                {
                    //Run the test.
                    gsl_vector* CADDWeights = nullptr;
                    if(useCADDWeights)
                    {
                        CADDWeights = geneInput.getCADDWeights();
                    }
                    wsbt test = wsbt(geneInput.getShortGslGenotype(), CADDWeights, result.getCaseCount(), iter->first, exactPvalueCalculation);
                    if(exactPvalueCalculation)
                    {
                        omp_set_lock(&outputLock);
                        test.driverOutput();
                        omp_unset_lock(&outputLock);

                        genes[i] = iter->first;
                        locations[i] = iter->second;
                        scores[i] = test.getScores();
                        testStats[i] = test.getTestStat();
                        pvalues[i] = test.getPvalue();
                        

                        if (outputFileName == "")
                        {
                            omp_set_lock(&outputLock);
                            cout << genes[i] << "\t" << locations[i] << "\t";
                            for(int j = 0; j < result.getCaseCount(); j++)
                            {
                                cout << scores[i][j] << "\t";
                            } 
                            cout << testStats[i] << "\t" << pvalues[i] << "\t" << effectSizes[i] << "\t" << UserUniqueVariantCounts[i] << "\t" 
                            << BackUniqueVariantCounts[i] << "\t" << variantCounts[i] << endl;
                            omp_unset_lock(&outputLock);
                        }
                        else
                        {
                            omp_set_lock(&outputLock);
                            out << genes[i] << "\t" << locations[i] << "\t";
                            for(int j = 0; j < result.getCaseCount(); j++)
                            {
                                out << scores[i][j] << "\t";
                            } 
                            out << testStats[i] << "\t" << pvalues[i] << "\t" << effectSizes[i] << "\t" << UserUniqueVariantCounts[i] << "\t" 
                            << BackUniqueVariantCounts[i] << "\t" << variantCounts[i] << endl;
                            omp_unset_lock(&outputLock);
                        }
                    }
                    else
                    {
                        genes[i] = iter->first;
                        locations[i] = iter->second;
                        pvalues[i] = test.getPvalue();
                        testStats[i] = test.getTestStat();
                    }
                }
            }
            out.close();

            gsl_sort_index(perm, pvalues.data(), 1, pvalues.size());
            if (exactPvalueCalculation)
            {
                out.open("sorted_" + outputFileName + ".tsv");

                //Header
                out << "Gene\tRegion\t";
                for(int i = 0; i < result.getCaseCount(); i++)
                {
                    out << result.getSampleNames()->at(i) << "_Score\t";
                }
                out << "TestStat\tPvalue\tRegionSize\tCaseUniqueVariants\tBackgroundUniqueVariants\tTotalVariants" << endl;

                //Data
                for (int i = 0; i < pvalues.size(); i++)
                {
                    //perm contains the sorted order.
                    out << genes[perm[i]] << "\t" << locations[perm[i]] << "\t";
                    for(int j = 0; j < result.getCaseCount(); j++)
                    {
                        out << scores[perm[i]][j] << "\t"; 
                    } 
                    out << testStats[perm[i]] << "\t" << pvalues[perm[i]] << "\t" << effectSizes[perm[i]] << "\t" << UserUniqueVariantCounts[perm[i]] 
                    << "\t" << BackUniqueVariantCounts[perm[i]] << "\t" << variantCounts[perm[i]] << endl;
                }
                out.close();
            }
            else
            {
                out.open("sorted_" + outputFileName + ".tsv");
                out << "Gene\tRegion\tTestStat\tPvalue" << endl;
                for (int i = 0; i < pvalues.size(); i++)
                {
                    //perm contains the sorted order.
                    out << genes[perm[i]] << "\t" << locations[perm[i]] << "\t" << testStats[perm[i]] << "\t" << pvalues[perm[i]] << endl;
                }
                out.close();
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
                    skipped++;
                }
            }
            double fisherPvalue = gsl_cdf_chisq_Q(fisherStat, 2 * geneCount);
            cout << "Fisher product test statistic on " << pvalues.size() - skipped <<  " genes is " << fisherStat << " with pvalue " << fisherPvalue << endl;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 60000.0 << " minutes." << endl;

        //Now manage output of results for user.
    }
    else if (testType == "burden")
    {
       
    }
    else if (testType == "cast")
    {
        cout << testType << " is not yet implemented." << endl;
    }
    else if (testType == "skat")
    {
        readInput result(currentDir, testType, variantRegion, vcffilename, region, phenofilename, covfilename);

        currentTime = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
        cout << endl;
        lasttime = currentTime;

        currentTime = std::chrono::high_resolution_clock::now();
        lasttime = currentTime;
        if (agnosticGeneRun)
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
                dataCollector geneInput = dataCollector(userBackgroundIncluded, useCADDWeights, vcffilename, backfilename, iter->second, testType, result.getCaseCount(), omp_get_thread_num());
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
            outFile.open(outputFileName + ".SKAT.results");
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
        //cout << testType << " is not yet implemented." << endl;
        readInput result(currentDir, testType, variantRegion, vcffilename, region, phenofilename, covfilename);
        currentTime = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << "After Input. Took " << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.0 << " seconds." << endl;
        cout << endl;
        lasttime = currentTime;

        if(agnosticGeneRun)
        {
            map<string,string> regions = result.getRegions();
            size_t perm[regions.size()];
            vector<string> genes(regions.size());
            vector<string> locations(regions.size());
            vector<int> effectSizes(regions.size());
            vector<int> UserUniqueVariantCounts(regions.size());
            vector<int> BackUniqueVariantCounts(regions.size());
            vector<int> variantCounts(regions.size());
            vector<double> pvalues(regions.size());
            vector<vector<double>> scores(regions.size(), vector<double>(result.getCaseCount()));
            vector<double> testStats(regions.size());

            //int index = 0;
            int skipped = 0;

            ofstream out(outputFileName + ".tsv");
            out << "Gene\tRegion\t";
            for(int i = 0; i < result.getCaseCount(); i++)
            {
                out << result.getSampleNames()->at(i) << "_Score\t";
            }
            out << "TestStat\tPvalue\tRegionSize\tCaseUniqueVariants\tBackgroundUniqueVariants\tTotalVariants" << endl;

            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < regions.size(); i++)
            {
                //Iterate to run the correct gene.
                map<string, string>::iterator iter = regions.begin();
                advance(iter, i);
                cout << "Running SKATO test on gene " << iter->first << endl;
                //Read in genotpe data from file for this gene.
                dataCollector geneInput = dataCollector(userBackgroundIncluded, useCADDWeights, vcffilename, backfilename, iter->second, testType, result.getCaseCount(), omp_get_thread_num());
                //Run SKATO Test
                auto runStartTime = std::chrono::high_resolution_clock::now();
                skato skatoTest = skato(geneInput.getGslGenotype(), result.getCovariates(), geneInput.getMaf(), result.getPheno());
                auto runCurrentTime = std::chrono::high_resolution_clock::now();
                cout << "SKATO Took " << std::chrono::duration_cast<std::chrono::milliseconds>(runCurrentTime - runStartTime).count() / 60000.0 << " minutes." << endl;
                
                //After test Reporting
                effectSizes[i] = stoi(iter->second.substr(iter->second.find("-")+1)) - stoi(iter->second.substr(iter->second.find(":")+1, iter->second.find("-")));
                variantCounts[i] = geneInput.getGslGenotype()->size1;
                UserUniqueVariantCounts[i] = geneInput.getCaseUniqueVariantCount();
                BackUniqueVariantCounts[i] = geneInput.getBackgroundUniqueVariantCount();
                pvalues[i] = skatoTest.getPvalue();


                //Run WSBT
                //wsbt wsbtTest = wsbt(geneInput.getGslGenotype(), result.getCaseCount(), iter->first, exactPvalueCalculation);

            }
        }
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
