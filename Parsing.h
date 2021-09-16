#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <array>

#include "Sequence.h"

#include "RandomGenerator.h"
#define _CRT_SECURE_NO_DEPRECATE

typedef vector<string> heads;
//typedef unordered_set<string> sset;
typedef vector<vector<string>> mutation;

class Hashes {
	
	public:
	    std::size_t operator()(string const& str) const {
	        int p = 53;
	        int m = 1e9 + 9;
	        long long power_of_p = 1;
	        std::size_t hash_val = 0;
	
	        // Loop to calculate the hash value
	        // by iterating over the elements of string
	        for (int i = 0; i < str.length(); i++) {
	            hash_val
	               = (hash_val
	                    + (str[i] - 'a' + 1) * power_of_p)
	                % m;
	            power_of_p
	                = (power_of_p * p) % m;
	        }
	        return hash_val;
	    }
	
};


class Parsing {
#define FIXEDSIZE 3

	template<typename T, std::size_t N>
	class arrayHash {
	public:
		std::size_t operator()(std::array<T, N> const& arr) const {
			std::size_t sum(0);
			for (auto&& i : arr) sum += std::hash<T>()(i);
			return sum;
		}
	};
	struct AMRgene {
		unordered_map<string, int>* kmerFreq;
		unordered_map<string, float>* kmerMapped;
		AMRgene() {
			kmerFreq = new unordered_map<string, int>();
			kmerMapped = new unordered_map<string, float>();
		}
		~AMRgene() {}

	};
public:
	Parsing(const int& any) {
		kmerFromgene = new unordered_map<string, heads,Hashes>();
		//kmerFromgene = new unordered_map<string, heads>();
		mutFromkmer = new unordered_map<string, heads, Hashes>();
		AMRvar = new unordered_map<string, AMRgene>();
		umap = new unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>>();
		k = any;
	}
	~Parsing() {}

public:
	void Aread(const string& filename);
	void Qread(const string& filename, const bool& reportMultipleHits);
	string Frames(const string& input);
	void Nucio();

public:
	unordered_map<string, heads, Hashes>* kmerFromgene;
	unordered_map<string, heads, Hashes>* mutFromkmer;
	unordered_map<string, AMRgene>* AMRvar;
	unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>> *umap;
	//unordered_map<string, char>* nuciomap;
	int k;
};

