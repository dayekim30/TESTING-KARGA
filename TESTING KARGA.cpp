// TESTING KARGA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <fstream>
#include <string>

#include <chrono>
#include <unordered_map>

#include "Parsing.h"

using namespace std;
using namespace std::chrono;




int main()
{
    cout << "put the number of k" << endl;
    int k;
    cin >> k;
    
    Parsing* parsing = new Parsing(k);
    parsing->Aread("kargva_db_v5.fasta");
    auto start = high_resolution_clock::now();
    parsing->Qread("testdata_simul.fastq" , true);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "this is the second : " << duration.count() << endl;
}
    //unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>> umap;
   
    //string filename("neclo.txt");
    //FILE* input_file = fopen(filename.c_str(), "r");
    //if (input_file == nullptr) {
    //    return EXIT_FAILURE;
    //}
    //
    //int c;
    //while ((c = fgetc(input_file)) != EOF) {
    //    int i = 1;
    //    array<int, FIXEDSIZE> tm = {c};
    //    while ((c = fgetc(input_file)) != EOF &&i<3) {
    //        tm[i] = c;
    //        i++;
    //    }
    //    umap[tm] = fgetc(input_file);
    //    fgetc(input_file);
    //}
    //std::cout << umap.size() << endl;
    //fclose( input_file);
    //
    //
    //{
    //    auto start = high_resolution_clock::now();

    //    string filename("testdata_simul.fastq");
    //    vector<char> bytes;

    //    FILE* input_file = fopen(filename.c_str(), "r");
    //    if (input_file == nullptr) {
    //        return EXIT_FAILURE;
    //    }
    //    int c;

    //    int i = 0;
    //    while ((c = fgetc(input_file)) != EOF) {
    //        int nlcount = 0;
    //       string id = "";
    //        vector<char> seq;

    //       
    //        id = id + (char)c;
    //        while ((c = fgetc(input_file)) != '\n') {
    //            id = id + (char)c;
    //        }
    //        array<int, FIXEDSIZE> tm;
    //        int a = 0;
    //        while ((c = fgetc(input_file)) != '\n') {
    //            tm[a++] = c;
    //            if (a == 3) {
    //                a = 0;
    //                if (umap.count(tm)) { seq.push_back(umap[tm]); }
    //                else seq.push_back('?');
    //            }
    //        }

    //        while ((c = fgetc(input_file)) != '\n') {
    //        }
    //        while ((c = fgetc(input_file)) != '\n' && (c = fgetc(input_file)) != EOF) {
    //        }
    //        i++;
    //       /* cout << "id : " << id << endl;
    //        string s(seq.begin(), seq.end());
    //        cout << "seq: " << s << endl;*/
    //            /*if ((c = fgetc(input_file)) == EOF) {
    //            break;
    //        }*/
    //    }


    //    //while ((c = fgetc(input_file)) != EOF) {
    //    //    //putchar(c);           
    //    //}
    //    cout << endl;
    //    fclose(input_file);
    //    auto stop = high_resolution_clock::now();
    //    auto duration = duration_cast<microseconds>(stop - start);

    //    // To get the value of duration use the count()
    //   // member function on the duration object
    //    cout << "this is the first : "<<duration.count() << endl;
    //    cout << "total i : " << i << endl;

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
