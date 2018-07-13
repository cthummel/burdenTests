//
//  main.cpp
//  burdenTest
//
//  Created by Corin Thummel on 6/29/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include <iostream>
#include "wsbt.cpp"
#include "input.cpp"
#include "cast.cpp"

using namespace std;

int main(int argc, const char * argv[])
{
    std::string vcffilename, vcfType, phenofilename, filename;
    
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
    }
    
    //cout << "correct filename: " << vcffilename << "\n";
    
    readInput result = readInput(vcffilename, vcfType, phenofilename);
    vector<double> pheno = vector<double>(result.getMaf().size());
    wsbt test = wsbt(result.getGenotype(), result.getMaf(), pheno);
    
    cout << "Variant weights: ";
    for(int i = 0; i < test.getWeights().size(); i++)
    {
        cout << test.getWeights()[1] << " ";
    }
    cout << endl;
    
    
    //cout << "Hello, World!\n";
    return 0;
}
