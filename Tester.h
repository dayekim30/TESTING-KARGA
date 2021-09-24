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
#include "kseq++.h"
#include "zlib.h"
#include "spp.h"

#include "Sequence.h"
#include "RandomGenerator.h"
using namespace klibpp;

#define _CRT_SECURE_NO_DEPRECATE
//#define _REMOVE_FPOS_SEEKPOS
#define _SILENCE_FPOS_SEEKPOS_DEPRECATION_WARNING
using spp::sparse_hash_map;
typedef vector<string> heads;
//typedef unordered_set<string> sset;
typedef vector<vector<string>> mutation;

class Hashes {

public:
	std::size_t operator()(string const& str) const {
		int p = 31;
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
struct AMRgene {
	unordered_map<string, int,Hashes>* kmerFreq;
	unordered_map<string, float, Hashes>* kmerMapped;
	AMRgene() {
		kmerFreq = new unordered_map<string, int, Hashes>();
		kmerMapped = new unordered_map<string, float, Hashes>();
	}
	~AMRgene() {}

};


class Tester {
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
	
public:
	Tester(const int& any) {
		//kmerFromgene = new unordered_map<string, heads, Hashes>();
		kmerFromgene = new unordered_map<string, heads, Hashes>();
		mutFromkmer = new unordered_map<string, heads, Hashes>();
		AMRvar = new unordered_map<string, AMRgene, Hashes>();
		//umap = new unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>>();
		umap = new unordered_map<string, char, Hashes>();
		k = any;
	}
	~Tester() {}

public:
	void Aread(const string& filename);
	void Qread(const string& filename, const bool& reportMultipleHits);
	string Frames(const string& input);
	void Nucio();

public:
	unordered_map<string, heads, Hashes>* kmerFromgene;
	unordered_map<string, heads, Hashes>* mutFromkmer;
	unordered_map<string, AMRgene, Hashes>* AMRvar;
	//unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>>* umap;
	unordered_map<string, char, Hashes>* umap;
	//unordered_map<string, char>* nuciomap;
	int k;
};

