#include "Tester.h"
#include <string>

#define _CRT_SECURE_NO_DEPRECATE
#define FIXEDSIZE 3

void Tester::Aread(const string& filename)
{
	
	std::cout << "Reading Gene database, creating k - mer mapping (k = " << k << ")" << endl;

	KSeq record;
	const char *fname = filename.c_str();
	gzFile fp = gzopen(fname, "r");
	auto ks = make_kstream(fp, gzread, mode::in);
	string id;
	string seq;
	int i{ 0 };
	//testing
	

	while (ks >> record) {
		id= record.name;		
		seq =record.seq;
		if (seq.length() < k)continue;
		//Find mutation index number
		size_t found_f = id.find('|');
		size_t found_s = id.find('|', found_f + 1);
		string mutwin = id.substr(found_f + 1, found_s - found_f - 1);
		
		id = id.substr(0, found_f+1);
		

		vector <string> tokens;
		while (mutwin.find(';') != string::npos) {
			size_t com = mutwin.find(';');
			size_t ni = mutwin.find_first_of("0123456789");
			tokens.push_back(mutwin.substr(ni, com-ni-1));
			mutwin = (mutwin.substr(com + 1));
			//cout << "the editted tamp: " << tamp << endl;
		}
		if (mutwin.find(';') == string::npos) { 
			size_t ni = mutwin.find_first_of("0123456789");
			tokens.push_back(mutwin.substr(ni));
		}
		
		unique_ptr<vector<int>> idx(new vector<int>());
		for (string a : tokens) {
			idx->push_back(std::stoi(a)-1);
		}

		
		AMRgene amrgene;
		for (int j = 0; j < seq.size() - k + 1; j++) {
			
			string kmer = seq.substr(j, k);
			
			if (!(kmerFromgene->count(kmer))) {

				vector<string>* ai = new vector<string>();
				ai->push_back(id);

				kmerFromgene->insert(make_pair(kmer, *ai));
			}
			else {
				vector<string>& ai = kmerFromgene->at(kmer);
				ai.push_back(id);
				
				kmerFromgene->insert_or_assign(kmer, ai);
			}

			
			if (!amrgene.kmerFreq->count(kmer)) { amrgene.kmerFreq->insert(make_pair(kmer, 1)); }
			else { amrgene.kmerFreq->insert_or_assign(kmer, amrgene.kmerFreq->at(kmer) + 1); }
			
			//Mutation Collections
			for (int a : *idx) {

				if (j <= a && a - k + 1 <= j) {
					if (mutFromkmer->count(kmer)) {
						auto muth = mutFromkmer->at(kmer);
						muth.push_back(id);
						mutFromkmer->insert_or_assign(kmer, muth);
					}
					else { vector<string> mutHead{ id }; mutFromkmer->insert(make_pair(kmer, mutHead)); }
				}
			}
			


		}
		AMRvar->insert(make_pair(id, amrgene));
		i++;
		if (i % 1000 == 0) {
			std::cout << i << ".. ";
		}

	}
	std::cout << "\nthe number of genes: " << i << "\n";
	
	gzclose(fp);
	
	//HERE

	//fstream* fout = new fstream();
	//fout->open("amrgeneTEST1.csv", ios::out | ios::app);
	//*fout << "id, kmer/freq \n";
	//for (auto any = AMRvar->begin(); any != AMRvar->end(); ++any) {
	//	*fout << any->first;
	//	//amrgfreq
	//	AMRgene ata = any->second;
	//	for (auto snd = ata.kmerFreq->begin(); snd != ata.kmerFreq->end(); ++snd) {
	//		*fout << "," << snd->first<<"/"<<snd->second;
	//	}
	//	*fout << "\n";
	//}	//HERE
	//
	//fout->close();

		//fstream* fout = new fstream();
		//fout->open("mutFromKmerTEST1.csv", ios::out | ios::app);
		//*fout << "kmer, id \n";
		//for (auto any = mutFromkmer->begin(); any != mutFromkmer->end(); ++any) {
		//	*fout << any->first;
		//	//	//amrgfreq
		//	//	//AMRgene ata = any->second;
		//	for (auto snd = any->second.begin(); snd != any->second.end(); ++snd) {
		//		*fout << "," << snd->c_str();
		//	}
		//	*fout << "\n";
		//}
		//fout->close();
	
}

	
//#define _CRT_SECURE_NO_WARNINGS
void Tester::Nucio()
{
	ifstream file("neclo.txt");
	if (file.is_open()) {
		string line;
		while (getline(file, line)) {
			umap->insert(make_pair(line.substr(0, 3), line[4]));
			
		}
		file.close();
	}

}

#define _CRT_SECURE_NO_DEPRECATE
string Tester::Frames(const string& input)
{

	vector<char> sseq;
	for (int a = 0; a < input.length() - 2; a = a + 3) {
		char tm[3] = { input[a] ,input[a + 1] ,input[a + 2] };
		if (umap->count(tm)) { sseq.push_back(umap->at(tm)); }
		else { sseq.push_back('?'); }
	}
	string s(sseq.begin(), sseq.end());
	return s;

}

bool sortByVals(const pair<string, float>& a,
	const pair<string, float>& b)
{
	return (a.second < b.second);
}


void Tester::Qread(const string& filename, const bool& reportMultipleHits)
{
	
	const int numT = 12500;
	Nucio();

	double avg = 0;
	int tc = 0;

	KSeq record;
	const char* fname = filename.c_str();
	{
		gzFile qfp = gzopen(fname, "r");
		auto ks = make_kstream(qfp, gzread, mode::in);
		while (ks >> record) {
			avg = avg + (double)(record.seq.length());
			tc++;
		}
		gzclose(qfp); 
	}
	
	avg = (avg / (double)tc) + 0.5;
	std::cout << "average read length is " << (unsigned long int)avg << " bases" << endl;

	
	int md[numT];
	//int* matchDist = md;
	RanGen r;

	for (int m = 0; m < numT; m++) {
		Sequence* sq = new Sequence(r.stringGeratpr((int)avg));
		string fk = sq->forwardSeq();
		string rk = sq->backwardSeq();

		string frames[6] = { Frames(fk), Frames(fk.substr(1)),Frames(fk.substr(2)),Frames(rk), Frames(rk.substr(1)),Frames(rk.substr(2)) };
		md[m] = 0;

		for (string frame : frames) {
			int hits = 0;
			for (int g = 0; g < frame.length() - k + 1; g++) {
				string q = frame.substr(g, k);
				if (mutFromkmer->count(q)) hits++;
			}
			if (md[m] < hits) md[m] = hits;
		}
		if (m % (numT / 5) == 0) std::cout << m << "..";
		delete sq;
	}

	sort(md, md + numT);
	int pvalue = md[99 * numT / 100];
	int* pvalthres = &pvalue;
	std::cout << "99th percentile of random k-mers match distribution is " << *pvalthres << "(max is " << md[numT - 1] << " )" << endl;
	
	std::cout << "Reading file and mapping genes" << endl;

	

	fstream* fout = new fstream();
	fout->open("_KARGVA_mappedReadsTEST1.csv", ios::out | ios::app);
	*fout << "Idx, GeneAnnotation, GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal";
	if (reportMultipleHits) *fout << ",...";
	*fout << "\n";

	tc = 0;
	gzFile fp = gzopen(fname, "r");
	auto ks = make_kstream(fp, gzread, mode::in);
	string id;
	string seq;

	while (ks >> record) {
		id = record.name;
		seq = record.seq;
		tc++;
		if ((seq.length() / 3) > k) {
		
			std::unique_ptr<Sequence> sq(new Sequence(seq));
			string fk = sq->forwardSeq();
			string rk = sq->backwardSeq();
			string frames[6] = { Frames(fk), Frames(fk.substr(1)),Frames(fk.substr(2)),Frames(rk), Frames(rk.substr(1)),Frames(rk.substr(2)) };
		
			std::unique_ptr<vector<string>> kmerHitsMandatoryBest(new vector<string>);
			std::unique_ptr<vector<string>> kmerHitsBest(new vector<string>);
			std::unique_ptr<unordered_map<string, float>> geneHitsMandatoryWeightedBest(new unordered_map<string, float>());
			std::unique_ptr<unordered_map<string, float>> geneHitsWeightedBest(new unordered_map<string, float>());
			std::unique_ptr<unordered_map<string, int>> geneHitsMandatoryUnweightedBest(new unordered_map<string, int>());
			std::unique_ptr<unordered_map<string, int>> geneHitsUnweightedBest(new unordered_map<string, int>());
			std::unique_ptr<int> bestMandatory (new int(0));
			
			int ca = 0;
			for (string frame : frames) {
				
				vector<string> kmerHitsMandatory;
				vector<string> kmerHits;
				unordered_map<string, float> geneHitsMandatoryWeighted;
				unordered_map<string, float> geneHitsWeighted;
				unordered_map<string, int> geneHitsMandatoryUnweighted;
				unordered_map<string, int> geneHitsUnweighted;

				
				for (int l = 0; l < frame.length() - k + 1; l++) {
					string fk = frame.substr(l, k);

					if (mutFromkmer->count(fk)) {
						kmerHitsMandatory.push_back(fk);
						int size = mutFromkmer->at(fk).size();
						for (string head : mutFromkmer->at(fk)) {

							float frac = 1 / (float)size;
							
							if (!geneHitsMandatoryWeighted.count(head)) { geneHitsMandatoryWeighted[head] = frac; }
							else { geneHitsMandatoryWeighted[head] = geneHitsMandatoryWeighted[head] + frac; }

							if (!geneHitsMandatoryUnweighted.count(head)) { geneHitsMandatoryUnweighted[head] = 1; }
							else { geneHitsMandatoryUnweighted[head] = geneHitsMandatoryUnweighted[head] + 1; }

						}
					}

					if (kmerFromgene->count(fk)) {
						
						kmerHits.push_back(fk);
						int size = kmerFromgene->at(fk).size();
						for (string head : kmerFromgene->at(fk)) {

							float frac = 1 / (float)size;
							if (!geneHitsWeighted.count(head)) { geneHitsWeighted[head] = frac; }
							else { geneHitsWeighted[head] = geneHitsWeighted[head] + frac; }

							if (!geneHitsUnweighted.count(head)) { geneHitsUnweighted[head] = 1; }
							else { geneHitsUnweighted[head] = geneHitsUnweighted[head] + 1; }
						}

					}
				}
				if (kmerHitsMandatory.size() > kmerHitsMandatoryBest->size()) {
					*kmerHitsMandatoryBest = kmerHitsMandatory;
					*geneHitsMandatoryWeightedBest = geneHitsMandatoryWeighted;
					*geneHitsMandatoryUnweightedBest = geneHitsMandatoryUnweighted;
					*bestMandatory = ca;
					*kmerHitsBest = (kmerHits);
					*geneHitsWeightedBest = (geneHitsWeighted);
					*geneHitsUnweightedBest = (geneHitsUnweighted);
				}
				ca++;	
			}
			if (kmerHitsMandatoryBest->size() > *pvalthres) {
				unordered_map<string, float>* scores = new unordered_map<string, float>();

				for (auto sets : *geneHitsMandatoryWeightedBest) {
					string key = sets.first;

					float sc1 = (geneHitsWeightedBest->at(key)) / (float)(kmerHitsBest->size());
					float sc2 = (geneHitsUnweightedBest->at(key)) / (float)(kmerHitsMandatoryBest->size());
					float sc = sc1 + sc2 - (sc1 * sc2);
					scores->insert(make_pair(key, sc));
				}

				vector<pair<string, float>> sclist;
				for (auto scset : *scores) {
					sclist.push_back(scset);
				}

				sort(sclist.begin(), sclist.end(), sortByVals);
				*fout << id << ",";
				float ratio = sclist.at(0).second;
				int stp = 0;

				for (int b = 0; b < sclist.size(); b++) {
					*fout << sclist.at(b).first << ",";
					float fr = sclist.at(b).second;
					float fp = (float)(((int)(fr * 100 + 0.5)) / 100);
					*fout << fp << "/" << geneHitsMandatoryUnweightedBest->at(sclist.at(b).first) << "/" << geneHitsUnweightedBest->at(sclist.at(b).first) << "/" << kmerHitsBest->size() << "/" << frames[*bestMandatory].length() - k + 1;
					if (b > 19 || fr / ratio < 0.05 || !reportMultipleHits) { stp = b; break; }
					*fout << ",";
				}
				*fout << "\n";
				if (!reportMultipleHits) { stp = 1; }
				for (int b = 0; b < stp; b++) {
					AMRgene genehit = AMRvar->at(sclist.at(b).first);
					for (int c = 0; c < kmerHitsBest->size(); c++) {
						string kh = kmerHitsBest->at(c);
						if (genehit.kmerFreq->count(kh)) {
							if (!genehit.kmerMapped->count(kh)) { genehit.kmerMapped->insert(make_pair(kh, (int)1)); }
							else { genehit.kmerMapped->insert_or_assign(kh, (genehit.kmerMapped->at(kh)) + 1); }
						}

					}
				}
				delete scores;
			}
			else {

				*fout << id << "," << "?,?/?/?/?/?" << "\n";

			}
		}
		if (tc % 100000 == 0) {
			std::cout << tc << " reads processed." << endl;

		}
		
	}
	fout->close();
	gzclose(fp);


	fstream* f_ot = new fstream();
	f_ot->open("_KARGVA_mappedGenesTEST1.csv", ios::out | ios::app);
	*f_ot << "GeneIdx, PercentGeneCovered, AberageKMerDepth\n";

	vector<string> keyset;
	for (auto set : *AMRvar) {
		keyset.push_back(set.first);
	}
	sort(keyset.begin(), keyset.end());

	for (string key : keyset) {
		int totKmers = AMRvar->at(key).kmerFreq->size();
		double percCovered = 0;
		double kmerDepth = 0;
		for (auto kmer : *AMRvar->at(key).kmerFreq) {
			string fk = kmer.first;
			if (AMRvar->at(key).kmerMapped->count(fk)) {
				percCovered = percCovered + 1;
				double dd = 0;
				dd = dd + (double)AMRvar->at(key).kmerMapped->at(fk);
				dd = dd / (double)AMRvar->at(key).kmerFreq->at(fk);
				kmerDepth = kmerDepth + dd;

			}
		}
		percCovered = percCovered / (double)(totKmers);
		kmerDepth = kmerDepth / (double)(totKmers);
		if (percCovered > 0.001f) {
			*f_ot << key << "," << 100 * percCovered << "%," << kmerDepth << "\n";
		}
	}
	f_ot->close();
	delete f_ot;
	f_ot = NULL;

}


