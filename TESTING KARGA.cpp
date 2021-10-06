// TESTING KARGA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <fstream>
#include <string>

#include <chrono>
#include <unordered_map>

#include "Parsing.h"
//#include "Tester.h"
using namespace std;
using namespace std::chrono;




int main()
{
    string input[3];
    int finput;
    string fasta = "kargva_db_v5.fasta";
    string fastq = "testdata_simul.fastq";

    cout << "Please write the input elements: 1) number of K   2)FASTA file   3)FASTQ file." << endl;
    
    cin >> input[0] >> input[1] >>input[2];

    cout << "Do you need mapped genes informaiton in detail? Yes -> 1, No ->0" << endl;
    
    cin >> finput;
    cout << "the number of K: " << finput << endl;
    
    int Knum = stoi(input[0]);
    if (input[1].find("fasta") != string::npos) { fasta = input[1]; }
    else { cout << "This file is not FASTA file."; return (0); };
    if(input[1].find("fastq") != string::npos) fastq = input[2];
    else { cout << "This file is not FASTQ file."; return (0); };

    
   cout << " The number K is " << Knum << ", FASTA file name: " << fasta << ", FASTQ file name: " << fastq;

   // Replacement point 1 :start
    Parsing* parsing = new Parsing(Knum);
    parsing->Aread(fasta);
    auto start = high_resolution_clock::now();
    parsing->Qread(fastq , finput);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "this is the second : " << duration.count() << endl;
    //Replacement point 1 :end
}

// Tester includes zlib which can open/read/write zip file but it still has a few bugs (Details are written in the document.)
// Tester is expected to run faster than the current one. After fixing the bugs, the below codes can be replacec with block (Replacement point 1 :start - end).

/*Tester* parsing = new Tester(Knum);
parsing->Aread(fasta);
auto start = high_resolution_clock::now();
parsing->Qread(fastq, finput);
auto stop = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(stop - start);
cout << "this is the second : " << duration.count() << endl;
*/


// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
