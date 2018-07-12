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

using namespace std;

int main(int argc, const char * argv[])
{
    std::string vcffilename, phenofilename, filename;
    
    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-vcf") == 0)
        {
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
    
    
    
    
    cout << "correct filename: " << vcffilename << "\n";
    
    readInput result = readInput(vcffilename, phenofilename);
    
    //std::vector<double> phenotype(result.getMaf.size());
    //wsbt test(result.getGenotype(), result.getMaf(), phenotype);
    
    
    
    /*
     if(argc < 1)
     {
        std::cout << "Please input a .vcf file";
        filename << std::cin;
     }
     else
     {
        filename = argv[1];
     }
     
     
     */
    
    
    
    
    
    cout << "Hello, World!\n";
    return 0;
}
