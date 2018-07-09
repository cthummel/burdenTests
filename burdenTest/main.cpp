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
    std::string filename;
    
    if(argc > 1)
    {
        filename = argv[1];
    }
    else
    {
        cout << "burdenTest <filename>\n";
        return 0;
    }
    if (filename.substr(filename.length() - 4) != ".vcf")
    {
        //They input a filename that isnt a vcf
        cout << "burdenTest <filename>\n";
        return 0;
    }
    
    cout << "correct filename: " << filename << "\n";
    
    readInput result = readInput(filename);
    
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
