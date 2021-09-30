// TESTING KARGA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <fstream>
#include <string>

#include <chrono>
#include <unordered_map>

//#include "Parsing.h"
#include "Tester.h"
using namespace std;
using namespace std::chrono;




int main()
{
    cout << "put the number of k" << endl;
    int k;
    cin >> k;
    
    Tester* parsing = new Tester(k);
    parsing->Aread("kargva_db_v5.fasta");
    auto start = high_resolution_clock::now();
    parsing->Qread("testdata_simul.fastq" , true);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "this is the second : " << duration.count() << endl;
}

//void readfastq(char* filename, int SRlength) {
//    _filelength = 0;
//    _SRlength = SRlength;
//
//    size_t bytes_read, bytes_expected;
//
//    FILE* fp;
//    fp = fopen(filename, "r");
//
//    fseek(fp, 0L, SEEK_END); //go to the end of file
//    bytes_expected = ftell(fp); //get filesize
//    fseek(fp, 0L, SEEK_SET); //go to the begining of the file
//
//    fclose(fp);
//
//    if ((_seqarray = (char*)malloc(bytes_expected / 2)) == NULL) //allocate space for file
//        err(EX_OSERR, "data malloc");
//
//
//    string name;
//    string seqtemp;
//    string garbage;
//    string phredtemp;
//
//    boost::iostreams::stream<boost::iostreams::file_source>file(filename);
//
//
//    while (std::getline(file, name)) {
//        std::getline(file, seqtemp);
//        std::getline(file, garbage);
//        std::getline(file, phredtemp);
//
//        if (seqtemp.size() != SRlength) {
//            if (seqtemp.size() != 0)
//                printf("Error on read in fastq: size is invalid\n");
//        }
//        else {
//            _names.push_back(name);
//
//            strncpy(&(_seqarray[SRlength * _filelength]), seqtemp.c_str(), seqtemp.length()); //do not handle special letters here, do on GPU
//
//            _filelength++;
//        }
//    }
//}
    


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
