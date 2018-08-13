//
//  output.cpp
//  burdenTest
//
//  Created by Corin Thummel on 7/30/18.
//  Copyright © 2018 Corin Thummel. All rights reserved.
//

#include "output.hpp"


writeOutput::writeOutput(string filename, string test_type, gsl_vector* weights)
{
    string line, zip, testInfoTag, description;
    smatch match;
    regex posMatch("^(\\d{1,2})\\s(\\d*)");
    
    //Make sure the filename we are dealing with ends with ".vcf".
    if(filename.substr(filename.length() - 3) == ".gz")
    {
        filename = filename.substr(0, filename.length() - 3);
    }
    
    //If a test doesnt output weights then the vector will be empty.
    if(weights->size != 0)
    {
        outfile.open("anno.tab");
        infile.open(filename);
        for(int i = 0; getline(infile, line); )
        {
            if(regex_search(line, match, posMatch))
            {
                outfile << match[0] << "\t" << match[2] << "\t" << gsl_vector_get(weights, i) << endl;
                i++;
            }
        }
        outfile.close();
        infile.close();
    }
    
    if (test_type == "burden")
    {
        
    }
    else if(test_type == "cast")
    {
        
    }
    else if (test_type == "wsbt")
    {
        testInfoTag = "WSBT_WEIGHT";
        description = "\"Weight from Weighted Sums Burden Test. Low weight means high impact and high weight means low impact.\"";
        outfile.open("anno.hdr");
        outfile << "##INFO=<ID=" + testInfoTag + ",Number=1,Type=Float,Description=" + description + ">" << endl;
        outfile.close();
        
        
    }
    else if (test_type == "skat")
    {
        
    }
    else if (test_type == "skato")
    {
        
    }
    
    
    
    
    if(filename.substr(filename.length() - 4) == ".vcf")
    {
        zip = bgzip_loc + " " + filename;
        system(zip.c_str());
        filename = filename + ".gz";
    }
    
    
    
    zip = bgzip_loc + " anno.tab";
    system(zip.c_str());
    string tabix = tabix_loc + " -s1 -b2 -e3 anno.tab.gz";
    system(tabix.c_str());
    //zip = "bgzip test5.vcf";
    //system(zip.c_str());
    string annotateCommand = bcftools_loc + " annotate -a anno.tab.gz -h anno.hdr -c CHROM,FROM,TO," + testInfoTag + " " + filename + " > testOutput.vcf";
    system(annotateCommand.c_str());
    zip = bgzip_loc + " testOutput.vcf";
    system(zip.c_str());
}


