#include "Parsing.h"
#include <chrono>


#define _CRT_SECURE_NO_DEPRECATE
#define FIXEDSIZE 3

using namespace std;
using namespace std::chrono;




void Parsing::Aread(const string& filename)
{

	fstream newfile;
	fstream* f_nf = &newfile;
	string id;
	string seq;

	f_nf->open(filename, ios::in);
	//newfile.open(filename, ios::in);

	std::cout << "Reading Gene database, creating k - mer mapping (k = " << k << ")" << endl;
	
	if (f_nf->is_open()) {
		string tp;
		int i = 0;
		//int l = 0;
		getline(*f_nf, tp);

		while (!f_nf->eof()) {

			if (tp[0] != '>') { std::cout << "wrong FASTA format" << endl; exit(0); }
			id = tp;
			seq = "";
			do {

				getline(*f_nf, tp);
				//l++;

				if (tp[0] == '>' || f_nf->eof())break;
				seq = seq + tp;

			} while (tp[0] != '>');

			i++;

			if (seq.size() < k)  continue;

			// triming -> id spliting |,; 
			size_t found = id.find('|');
			short first = found;
			string tamp = "";
			if (found != string::npos) {
				found = id.find('|', found + 1);
				if (found != string::npos) {
					tamp = id.substr(first + 1, found - first - 1);
				}
				id = id.substr(0, first + 1);
			}

			//cout << "tamp: " << tamp << endl;
			vector <string> tokens;

			while (tamp.find(';') != string::npos) {
				int com = tamp.find(';');
				tokens.push_back(tamp.substr(0, com));
				tamp = tamp.substr(com + 1);
				//cout << "the editted tamp: " << tamp << endl;
			}
			if (tamp.find(';') == string::npos) tokens.push_back(tamp);

			//getting mutation hash table from tokens


			vector<int>* positions = new vector<int>();
			for (string a : tokens) {
				int n;
				int st;
				//char c;
				for (int h = 0; h < a.length(); h++) {
					if (isdigit(a.at(h))) { st = h; break; }
				}
				size_t found = a.find("STOP");
				if (found != string::npos) {
					n = stoi(a.substr(st, found - 1));
					//c = '*';
				}
				else {
					n = stoi(a.substr(st, a.length() - 2));
					//c = a.at(a.length() - 1);
				}
				positions->push_back(n - 1);
				//cout << "target n:" << n << " tanget c: " << c << endl;
					// check there is targeted mutation in this sequence
				//if (seq.at(n - 1) == c) {
				//	//cout << "yes it is here" << endl;

				//	//make mutaion hashtable --- think about making this as hashset
				//	vector<string>* mt = new vector<string>();
				//	for (int x = 1; x < k + 1; x++) {
				//		int ft = (n - 1) - (k - x);
				//		if (ft < 0 || n + x - 1>seq.length()) continue;
				//		mt->push_back(seq.substr(ft, k));

				//	}
				//	muts->push_back(*mt);
				//	mutFromkmer->insert_or_assign(id, *muts);
				//}
			}

			////HERE
			//*fout << id;
			//for (int a : *positions) {
			//	//HERE
			//	*fout << "," << a;
			//}
			////HERE
			//*fout << "\n";
			
			AMRgene* am = new AMRgene();


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
					//kmerFromgene->insert(make_pair(kmer,ai));
					kmerFromgene->insert_or_assign(kmer, ai);
				}

				if (!am->kmerFreq->count(kmer)) { am->kmerFreq->insert(make_pair(kmer, 1)); }
				else { am->kmerFreq->insert_or_assign(kmer, am->kmerFreq->at(kmer) + 1); }

				//Collect Mutation 
				for (int a : *positions) {
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
			AMRvar->insert(make_pair(id, *am));
			delete positions;
			if (i % 1000 == 0) {
				std::cout << i << ".. ";
			}

		}
		std::cout << "\nthe number of genes: " << i << "\n";
		//cout << "the bucket size: " << kmerFromgene->size() << endl;
		//cout << " the bucket mut size: " << mutFromkmer->size() << endl;
		
		f_nf->close();
		f_nf = NULL;
		//fout->close();
	}

	
	
	



}
//#define _CRT_SECURE_NO_WARNINGS
void Parsing::Nucio()
{
	//unordered_map<std::array<int, FIXEDSIZE>, int, arrayHash<int, FIXEDSIZE>> umap;

	string filename("neclo.txt");
	FILE* input_file = fopen(filename.c_str(), "r");
	if (input_file == nullptr) {
		cout<< EXIT_FAILURE;
	}

	int c;
	while ((c = fgetc(input_file)) != EOF) {
		int i = 1;
		array<int, FIXEDSIZE> tm = { c };
		while ((c = fgetc(input_file)) != EOF && i < 3) {
			tm[i] = c;
			i++;
		}
		umap->insert(make_pair(tm,fgetc(input_file)));
		fgetc(input_file);
	}
	
	fclose(input_file);
}
#define _CRT_SECURE_NO_DEPRECATE
string Parsing::Frames(const string& input)
{

	vector<char> sseq;
	for (int a = 0; a < input.length() - 2; a = a + 3) {
		array<int, FIXEDSIZE> tm;
		tm[0] = input[a];
		tm[1] = input[a + 1];
		tm[2] = input[a + 2];
		if (umap->count(tm)) { sseq.push_back(umap->at(tm)); }
		else { sseq.push_back('?'); }
	}
	string s(sseq.begin(), sseq.end());
	return s;

}

bool sortByVal(const pair<string, float>& a,
	const pair<string, float>& b)
{
	return (a.second < b.second);
}


void Parsing::Qread(const string& filename, const bool& reportMultipleHits)
{   
	/*unsigned long long _filelength = 0;
	int _SRlength = 8192;

	size_t bytes_read, bytes_expected;*/

	const int numT = 12500;
	Nucio();
	
	//const char* filenames = filename.c_str();
	
	//fstream infile;

	fstream* i_nf = new fstream();
	i_nf->open(filename, ios::in);
	

	//cout << "Translation starts (k = " << k << ")" << endl;
	//FILE* i_nf;
	//i_nf = fopen(filenames, "r");
	//fseek(i_nf, 0L, SEEK_END); //go to the end of file
 //   bytes_expected = ftell(i_nf); //get filesize
 //   fseek(i_nf, 0L, SEEK_SET); //go to the begining of the file
	//fclose(i_nf);
	//if (((char*)malloc(bytes_expected / 2)) == NULL) //allocate space for file
	////        err(EX_OSERR, "data malloc");
	//	string name;
	//    string seqtemp;
	////    string garbage;
	////    string phredtemp;

	//boost::iostreams::stream<boost::iostreams::file_source>file(filename);



	double avg = 0;
	int t = 0;
	int& tc = t;
	if (i_nf->is_open()) {
		string tp;

		while (getline(*i_nf, tp)) {
			if (tp[0] != '@') {
				std::cout << "this is wrong FASTQ file format." << endl;
				exit(0);
			}
			getline(*i_nf, tp);
			avg = avg + (double)(tp.length());
			getline(*i_nf, tp);
			getline(*i_nf, tp);
			tc++;
		}
		i_nf->close();
		delete i_nf;
		i_nf = NULL;
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
		
		string frames[6] = {Frames(fk), Frames(fk.substr(1)),Frames(fk.substr(2)),Frames(rk), Frames(rk.substr(1)),Frames(rk.substr(2)) };
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
	//delete matchDist;
	//matchDist = NULL;
	std::cout << "Reading file and mapping genes" << endl;

	tc = 0;

	fstream* fout = new fstream();
	fout->open("_KARGVA_mappedReads.csv", ios::out | ios::app);
	*fout << "Idx, GeneAnnotation, GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal";
	if (reportMultipleHits) *fout << ",...";
	*fout << "\n";

	fstream* f_in = new fstream();
	f_in->open(filename, ios::in);
	string id;
	string seq;
	if (f_in->is_open()) {
		string tp;

		while (getline(*f_in, tp)) {
			if (tp[0] != '@') {
				std::cout << "this is wrong FASTQ file format." << endl;
				exit(0);
			}
			id = tp;
			getline(*f_in, tp);
			seq = tp;
			getline(*f_in, tp);
			getline(*f_in, tp);
			tc++;
			if ((seq.length() / 3) > k) {

				Sequence* sq = new Sequence(seq);
				string fk = sq->forwardSeq();
				string rk = sq->backwardSeq();
				string frames[6] = { Frames(fk), Frames(fk.substr(1)),Frames(fk.substr(2)),Frames(rk), Frames(rk.substr(1)),Frames(rk.substr(2)) };
				/*cout << "id: " << id << endl;
				cout << "seq: " << frames[0] << endl;*/
				/*vector<string>* kmerHitsMandatoryBest = new vector<string>();
				vector<string>* kmerHitsBest = new vector<string>();
				unordered_map<string, float>* geneHitsMandatoryWeightedBest = new unordered_map<string, float>();
				unordered_map<string, float>* geneHitsWeightedBest = new unordered_map<string, float>();
				unordered_map<string, int>* geneHitsMandatoryUnweightedBest = new unordered_map<string, int>();
				unordered_map<string, int>* geneHitsUnweightedBest = new unordered_map<string, int>();*/
				std::unique_ptr<vector<string>> kmerHitsMandatoryBest(new vector<string>);
				std::unique_ptr<vector<string>> kmerHitsBest(new vector<string>);// = nullptr;
				std::unique_ptr<unordered_map<string, float>> geneHitsMandatoryWeightedBest(new unordered_map<string, float>());// = nullptr;
				std::unique_ptr<unordered_map<string, float>> geneHitsWeightedBest(new unordered_map<string, float>());// = nullptr;
				std::unique_ptr<unordered_map<string, int>> geneHitsMandatoryUnweightedBest(new unordered_map<string, int>());// = nullptr;
				std::unique_ptr<unordered_map<string, int>> geneHitsUnweightedBest(new unordered_map<string, int>());// = nullptr;

				int* bestMandatory = new int(0);
				int ca = 0;
				for (string frame : frames) {
					//vector<string>* kmerHitsMandatory = new vector<string>();
					//vector<string>* kmerHits = new vector<string>();
					//unordered_map<string, float>* geneHitsMandatoryWeighted = new unordered_map<string, float>();
					//unordered_map<string, float>* geneHitsWeighted = new unordered_map<string, float>();
					//unordered_map<string, int>* geneHitsMandatoryUnweighted = new unordered_map<string, int>();
					//unordered_map<string, int>* geneHitsUnweighted = new unordered_map<string, int>();
					
					vector<string> kmerHitsMandatory;
					vector<string> kmerHits;
					unordered_map<string, float> geneHitsMandatoryWeighted;
					unordered_map<string, float> geneHitsWeighted;
					unordered_map<string, int> geneHitsMandatoryUnweighted;
					unordered_map<string, int> geneHitsUnweighted;			
					
					ca++;
					for (int l = 0; l < frame.length() - k + 1; l++) {
						string fk = frame.substr(l, k);

						if (mutFromkmer->count(fk)) {
							kmerHitsMandatory.push_back(fk);
							int size = mutFromkmer->at(fk).size();
							for (string head : mutFromkmer->at(fk)) {

								float frac = 1 / (float)size;
								/*if (!geneHitsMandatoryWeighted->count(head)) { geneHitsMandatoryWeighted->insert(make_pair(head, frac)); }
								else { geneHitsMandatoryWeighted->insert_or_assign(head, geneHitsMandatoryWeighted->at(head) + frac); }

								if (!geneHitsMandatoryUnweighted->count(head)) { geneHitsMandatoryUnweighted->insert(make_pair(head, 1)); }
								else { geneHitsMandatoryUnweighted->insert_or_assign(head, geneHitsMandatoryUnweighted->at(head) + 1); }*/
								
								if (!geneHitsMandatoryWeighted.count(head)) { geneHitsMandatoryWeighted[head] = frac; }
								else { geneHitsMandatoryWeighted[head]= geneHitsMandatoryWeighted[head] + frac; }

								if (!geneHitsMandatoryUnweighted.count(head)) { geneHitsMandatoryUnweighted[head] =  1; }
								else { geneHitsMandatoryUnweighted[head]= geneHitsMandatoryUnweighted[head] + 1; }

							}
						}

						if (kmerFromgene->count(fk)) {
							/*kmerHits.push_back(fk);
							int size = kmerFromgene->at(fk).size();
							for (string head : kmerFromgene->at(fk)) {

								float frac = 1 / (float)size;
								if (!geneHitsWeighted->count(head)) { geneHitsWeighted->insert(make_pair(head, frac)); }
								else { geneHitsWeighted->insert_or_assign(head, geneHitsWeighted->at(head) + frac); }

								if (!geneHitsUnweighted->count(head)) { geneHitsUnweighted->insert(make_pair(head, 1)); }
								else { geneHitsUnweighted->insert_or_assign(head, geneHitsUnweighted->at(head) + 1); }
							}*/

							kmerHits.push_back(fk);
							int size = kmerFromgene->at(fk).size();
							for (string head : kmerFromgene->at(fk)) {

								float frac = 1 / (float)size;
								if (!geneHitsWeighted.count(head)) { geneHitsWeighted[head]=frac; }
								else { geneHitsWeighted[head]= geneHitsWeighted[head] + frac; }

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
					
					//delete kmerHitsMandatory, kmerHits, geneHitsMandatoryWeighted, geneHitsWeighted, geneHitsMandatoryUnweighted, geneHitsUnweighted;

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

					sort(sclist.begin(), sclist.end(), sortByVal);
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
								if (!genehit.kmerMapped->count(kh)) { genehit.kmerMapped->insert(make_pair(kh, 1)); }
								else { genehit.kmerMapped->insert_or_assign(kh, genehit.kmerMapped->at(kh) + 1); }
							}

						}
					}
					delete scores;
				}
				else {

					*fout << id << "," << "?,?/?/?/?/?" << "\n";

				}

				delete sq, kmerHitsMandatoryBest, kmerHitsBest, geneHitsMandatoryWeightedBest, geneHitsWeightedBest, geneHitsMandatoryUnweightedBest, geneHitsUnweightedBest, bestMandatory;
			}
			if (tc % 100000 == 0) {
				std::cout << tc << " reads processed." << endl;

			}



		}
		f_in->close();
		fout->close();
	}
	delete fout, f_in;
	fout = NULL;
	f_in = NULL;



	fstream* f_ot = new fstream();
	f_ot->open("_KARGVA_mappedGenes.csv", ios::out | ios::app);
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


