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
		kmerFromgene = new unordered_map<string, heads>();
		//kmerFromgene = new unordered_map<string, heads>();
		mutFromkmer = new unordered_map<string, heads>();
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
	unordered_map<string, heads>* kmerFromgene;
	unordered_map<string, heads>* mutFromkmer;
	unordered_map<string, AMRgene>* AMRvar;
	unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>> *umap;
	//unordered_map<string, char>* nuciomap;
	int k;
};

