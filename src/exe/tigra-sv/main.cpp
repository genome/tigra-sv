/***************************************************************************
 *   Copyright (C) 2010 by Ken Chen & Xian Fan   *
 *   kchen@genome.wustl.edu, xfan@genome.wustl.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include <version.h>
#include <iostream>
#include <fstream>
//#include <strstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <map>
#include <assert.h>
#include <sstream>
//#include "tigra.h"
#include <stdio.h>
#include <stdlib.h>
#include <locale>
#include <algorithm>
#include <inttypes.h>
#include  "sam.h"
#include "bam.h"
#include "faidx.h"
#include "ksort.h"
#include "khash.h"
#include "samtools.h"
#include <sys/stat.h> 

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif


namespace std
{
	using namespace __gnu_cxx;
}


using namespace std;

string options ("");

typedef struct eqstr
{
	bool operator()(string s1, string s2) const
	{
		return s1.compare(s2) == 0;
	}
};

typedef struct test_HH
{
  	int n;
  	int AI;
	int CI;
	int GI;
	int TI;
	int AO;
	int CO;
	int GO;
	int TO;
	int tag;
	int tag2;
};

typedef struct BD_data{
	string chr1;
	int pos1;
	string ori1;
	string chr2;
	int pos2;
	string ori2;
	string type;
	int size;
	int score;
	int nreads;
	string nreads_lib;
	vector<string> extra;
	string line;
	string cnstr;
	string gene;
	string database;
	string bam_related;
};

typedef struct test_data
{
	string key2;
	int value;
};

typedef struct read_data{
	int32_t NM;
	int32_t MF;
	string name;
	uint32_t flag;
	uint32_t mqual;
	string seq;
	float PercMapped;
};


string itos(int i){
    stringstream i_str_stream;
    i_str_stream << i;
    return i_str_stream.str();
}

string ftos(float i){
    stringstream i_str_stream;
    i_str_stream << i;
    return i_str_stream.str();
}

// tools to be used in the following classes
class tools {
public:
/*	// split rReads by the sep, sep could be extended continuously, results to be splitted
	tools(){}
	void split(string rReads, string sep, vector<string> &splitted);
	string chomp(string str);
	string reverse(string str);
	string tr(string str, string str_old, string str_new);
	void sort_index(map<int, int> &a, vector<int> &index, int begin, int end);
	void sort_index_random_key(map<int, int> &a, vector<int> &index)
	void sort_value(map<int, int> &a, vector<int> &value, int begin, int end);
	void sort_value_reverse_key(map<int, int> &h, vector<int> &sorted_key)
};*/

	void split(string rReads, string sep, vector<string> &splitted){
		size_t pos_start = 0;
		size_t pos;
		if(sep.length() > 0){
			do {
				pos = rReads.find(sep, pos_start);
				if(pos != string::npos){
					if(pos - pos_start > 0){
						string str_tmp = rReads.substr(pos_start, pos - pos_start);
						splitted.push_back(str_tmp);
					}
					rReads = rReads.substr(pos + 1);
				}
				else if(rReads.length() != 0)
					splitted.push_back(rReads);
			} while (pos != string::npos);
		}
		else{
			for(int i = 0; i < rReads.length(); i++){
				splitted.push_back(rReads.substr(i,1));
			}
		}
		
		return;
	}
	
	// move the last endl sign, and return back the string
	string chomp(string str){
                string ret;
                string::size_type pos = str.find_last_not_of("\n");
                if(pos != string::npos)
                    ret = str.substr(0, pos+1);
                else
                    ret = str;
                return ret;
	}
	
	// take the string and reverse everything
	string reverse(string str){
		string str_ret = "";
		for(int i = str.length()-1; i >= 0; i--){
			str_ret += (str.substr(i, 1)).c_str()[0];			
		}		
		return str_ret;
	}
	
	// take the string and replace ATGC with TACG
	string tr(string str, string str_old, string str_new){
		char str_[str.length()];
		strcpy(str_, str.c_str());
		for(int i = 0; i < str.length(); i++){
			for(int j = 0; j < str_old.length(); j++){
				//if(str.substr(i, 1).compare(str_old.substr(j, 1)) == 0){
				if(str[i] == str_old[j]){
					//str.replace(i,1,str_new.substr(j, 1));
					str_[i] = str_new[j];
					break;
				}
			}
		}
		string ret(str_);
		//return str;
		return ret;
	}		
	
	void sort_index(map<int, float> &a, vector<int> &index, int begin, int end){
		vector<float> b;
		for(int i = begin; i <= end; i++){
			b.push_back(a[i]);
			index.push_back(i);
		}
		for(int i = 0; i < b.size() - 1; i++){
			for(int j = i + 1; j < b.size(); j++){
				if(b[i] > b[j]){
					float tmp = b[i];
					b[i] = b[j];
					b[j] = tmp;
					int tmp_index = index[i];
					index[i] = index[j];
					index[j] = tmp_index;
				}
			}
		}
	}
	
	void sort_index_random_key(map<int, int> &a, vector<int> &index){
		vector<int> b;
		for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it++){
			b.push_back((*a_it).second);
			index.push_back((*a_it).first);
		}
		for(int i = 0; i < index.size() - 1; i++){
			for(int j = i+1; j < index.size(); j++){
				if(b[i] < b[j]){
					int tmp = b[i];
					b[i] = b[j];
					b[j] = tmp;
					int tmp_index = index[i];
					index[i] = index[j];
					index[j] = tmp_index;
				}
			}
		}
	}
	
	void sort_value(map<int, float> &a, vector<float> &value, int begin, int end){
		for(int i = begin; i <= end; i++){
			value.push_back(a[i]);
		}
		for(int i = 0; i < value.size() - 1; i++){
			for(int j = i + 1; j < value.size(); j++){
				if(value[i] > value[j]){
					float tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
				}
			}
		}
	}
	
	/*void sort_value_reverse_key_int(map<int, int> &h, vector<int> &sorted_key){
		vector<int> value;
		for(map<int, int>::iterator h_it = h.begin(); h_it != h.end(); h_it ++){
			value.push_back((*h_it).second);
			sorted_key.push_back((*h_it).first);
		}
		for(int i = 0; i < h.size() - 1; i ++){
			for(int j = 0; j < h.size(); j++){
				if(value[i] < value[j]){
					int tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
					int tmp_index = sorted_key[i];
					sorted_key[i] = sorted_key[j];
					sorted_key[j] = tmp_index;
				}
			}
		}
	}*/

	void sort_value_reverse_key_int(map<int, int> &h, vector<int> &sorted_key, int begin, int end){
		vector<int> value;
		if(begin < 0 || end < 0){
			for(map<int, int>::iterator h_it = h.begin(); h_it != h.end(); h_it ++){
				value.push_back((*h_it).second);
				sorted_key.push_back((*h_it).first);
			}
		}
		else{
			for(int i = begin; i <= end; i++){
				value.push_back(h[i]);
				sorted_key.push_back(i);
			}
		}
		for(int i = 0; i < value.size() - 1; i ++){
			for(int j = i+1; j < value.size(); j++){
				if(value[i] < value[j]){
					int tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
					int tmp_index = sorted_key[i];
					sorted_key[i] = sorted_key[j];
					sorted_key[j] = tmp_index;
				}
			}
		}
	}

	void sort_value_reverse_key_float(map<int, float> &h, vector<int> &sorted_key, int begin, int end){
		vector<float> value;
                if(h.size() == 1){
                    sorted_key.push_back((*(h.begin())).first);
                    return;
                }
                else if(h.size() == 0)
                    return;
                //cout << "begin: " << begin << "\tend: " << end << endl;
                //cout << "h size: " << h.size() << endl;
		if(begin < 0 || end <= 0){
			for(map<int, float>::iterator h_it = h.begin(); h_it != h.end(); h_it ++){
				value.push_back((*h_it).second);
				sorted_key.push_back((*h_it).first);
			}
		}
		else{
			for(int i = begin; i <= end; i++){
				value.push_back(h[i]);
				sorted_key.push_back(i);
			}
		}
//                cout << "value size: " << value.size()-1 << endl;
		for(int i =0; i < value.size()-1; i ++){
			for(int j = i+1; j <value.size(); j++){
				if(value[i] < value[j] || value[i] == value[j] && sorted_key[i] > sorted_key[j]){
					float tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
					int tmp_index = sorted_key[i];
					sorted_key[i] = sorted_key[j];
					sorted_key[j] = tmp_index;
				}
			}
		}
		/*for(int i = 0; i < value.size() - 1; i++){
			for(int j = i+1; j < value.size(); j++){
				if(value[i] == value[j] && sorted_key[i] > sorted_key[j]){
					float tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
					int tmp_index = sorted_key[i];
					sorted_key[i] = sorted_key[j];
					sorted_key[j] = tmp_index;
				}
			}
		}*/
		
	}
	
	void sort_value_reverse_key_str(map<string, int> &h, vector<string> &sorted_key){
		vector<int> value;
		for(map<string, int>::iterator h_it = h.begin(); h_it != h.end(); h_it ++){
			value.push_back((*h_it).second);
			sorted_key.push_back((*h_it).first);
		}
		for(int i = 0; i < h.size() - 1; i ++){
			for(int j = 0; j < h.size(); j++){
				if(value[i] < value[j]){
					int tmp = value[i];
					value[i] = value[j];
					value[j] = tmp;
					string tmp_index = sorted_key[i];
					sorted_key[i] = sorted_key[j];
					sorted_key[j] = tmp_index;
				}
			}
		}
	}
	
	void sort_by_longest_length(vector<string> &longest_uniq_paths){
	        if(longest_uniq_paths.size()<=0) return;
		vector<vector<string> > paths;
		// separate them to numbers by numbers first
		string sep = ".";
		for(int i = 0; i < longest_uniq_paths.size(); i++){
			vector<string> splitted;
			split(longest_uniq_paths[i], sep, splitted);
			paths.push_back(splitted);
		}
		// see which string contains more numbers
		for(int i = 0; i < longest_uniq_paths.size() - 1; i++){
			for(int j = i+1; j < longest_uniq_paths.size(); j++){
				if(paths[i].size() < paths[j].size()){
					string tmp = longest_uniq_paths[j];
					longest_uniq_paths[j] = longest_uniq_paths[i];
					longest_uniq_paths[i] = tmp;
					vector<string> tmp_vec = paths[i];
					paths[i] = paths[j];
					paths[j] = tmp_vec;
				}
			}
		}
	}
	
	string StringToLower(string strToConvert)
	{//change each element of the string to lower case
		for(unsigned int i=0;i<strToConvert.length();i++)
		{
			char tmp = strToConvert[i];
			strToConvert[i] = tolower(tmp);
		}
		return strToConvert;//return the converted string
	}

	void analyze_BAM_list(map<string, string> &bams, string file_list){
		ifstream BAM;
		BAM.open(file_list.c_str());
		char line_[1024];
		if(BAM.is_open()){
			while(BAM.good()){
				BAM.getline(line_, 1024);
				string line(line_);
				if(line.length() == 0)
					break;
				vector<string> tmp;
				split(line, ":", tmp);
				bams[tmp[0]] = tmp[1];
			}
			BAM.close();
		}
	}
	
};

// generate kmer by hh[kmer]=number_of_kmer, in which number_of_kmer < C && > c
class kmergen {
public:
	int k;// = 25;
	int c;// = 2;
	int C;// = 2e9;
	string filter;
	
	kmergen(){//int k, int c, int C){
		k = 25;//k;
		c = 2;//c;
		C = 2e9;//C;
		filter = "";
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	void set_filter(string filter_){
		/*tools tl;
		ifstream FASTAS;
		FASTAS.open(filter_.c_str());
		string header;
		char line_[1024*1024];
		if(FASTAS.is_open()){
			while(FASTAS.good()){
				FASTAS.getline(line_, 1024*1024);
				if(line_[0] != '>'){
					string line(line_);
					line = tl.chomp(line);
					filter += line;
				}
			}
			FASTAS.close();
			cout << "\n" << filter << endl;
		}*/
		filter = filter_;		
	}
	
	int doit(vector<string> &rReads, map<string, int> &hh){
		tools tl;
		int totalcount = 0;
		for(vector<string>::iterator rReads_it = rReads.begin(); rReads_it != rReads.end(); rReads_it++){			
//?? my @segemnts = split/N+/ -> split the string wherever there is "N": seems like that
			string sep = "N";
			vector<string> segments;
			tl.split(*rReads_it, sep, segments);	// need to work on split function, done
			for(vector<string>::iterator seg = segments.begin(); seg != segments.end(); seg ++){
				int kmers = (*seg).length() - k;
				if(kmers < 1)
					continue;
				for(int i = 0; i <= kmers; i++){
					string w = (*seg).substr(i, k);
					
					if(filter.length() != 0 && filter.find(w) != string::npos){ // need to see how to implement i: not specifying the case
						int test1 = 0;
						continue;
					}
					if(hh.find(w) != hh.end())
						hh[w] += 1;
					else {
						//string vw = tl.reverse(w);	// need to work on reverse function, done
						string vw = w;
						reverse(vw.begin(),vw.end());
						vw = tl.tr(vw, "ATGC", "TACG"); // need to work on tr function, done
						if(filter.length() != 0 && filter.find(vw) != string::npos){ // need to see how to implement i: not specifying the case
							int test1 = 0;
							continue;
						}
						if(hh.find(vw) != hh.end())
							hh[vw] += 1;
						else
							hh[w] = 1;
					}
					totalcount ++;
				}
			}
		}
	
		// the following code seems useless in the class 
		//int aa[1001];
		//for(int i = 1; i <= 1000; i++)
		//	aa[i] = 0;
	
		int *aa=(int*)calloc(1001,sizeof(int));
		map<string, int> hh_tmp;
		for(map<string, int>::iterator hh_it = hh.begin(); hh_it != hh.end(); hh_it ++){
			(*hh_it).second > 1000 ? aa[1000]+=1 : aa[(*hh_it).second]+=1;			
			if((*hh_it).second < c || (*hh_it).second > C)
				totalcount -= (*hh_it).second;			
			else
				hh_tmp[(*hh_it).first] = (*hh_it).second;
		}
		free(aa);
		hh = hh_tmp;

		return totalcount;
	}
	
	void printMer(string fout, map<string, int> &hh){
		char fout_[fout.length()];
		strcpy(fout_, fout.c_str());
		ofstream fh;
		fh.open(fout_);
		if(fh.is_open()){
			for(map<string, int>::iterator hh_it = hh.begin(); hh_it != hh.end(); hh_it ++){
				fh << (*hh_it).first << "\t" << (*hh_it).second << endl;
			}
		}
		else {
			cerr << "unable to open " << fout << endl;
		}
		//delete []fout_; // need to see if necessary to release this memory
		fh.close();
	}	
};

// with HH generated in graphgen, walk through with the input(I) and output(O) alleles, and connect them to form a contig, and store it to Contigs[node]=ContigsSeq
class walknodes {
public:
	int k;// = 25;
	int Thin;// = 2;
	int Walkcutoff;// = 4;
	
	map<string, test_HH > HH; // Xian: map<kmer, map<info, count> >; HH[kmer][tag] = contigID
	map<int, string> Contigs; // Xian: map<contigID, contig>
	map<int, int> Contiglens; // Xian: map<contigID, contigLength>
	map<int, float> Contigcovs; // Xian: map<contigID, contigCoverages>
	map<int, string> Contigtypes; // Xian: map<contigID, contigTypes>
	int Contignum;
	map<int, int> Contigtags; 
	
	walknodes(){//int k, int Thin, int Walkcutoff){
		k = 25;//k;
		Thin = 2;// Thin;
		Walkcutoff = 4;//Walkcutoff;
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	~walknodes(){
		HH = map<string,test_HH >();
		Contigs = map<int, string>();
		Contiglens = map<int, int>();
		Contigcovs = map<int, float>();
		Contigtypes = map<int, string>();
		Contigtags = map<int, int>();
	}
	
	// generate proto-contigs
	// Xian: this class is to maximally connect the root branch from one key (node) -> Contig
	void strictwalk(map<string, test_HH > &rHH){
//ofstream fh;
//fh.open("leftright");
		tools tl;
		Contignum = 0;
		HH = rHH;
int tmp = HH.size();
		for(map<string, test_HH >::iterator HH_it = HH.begin(); HH_it != HH.end(); HH_it++){
			HH[(*HH_it).first].tag = 0;
		}	/////// this loop seems unnecessary seeing the following loop
		//cerr << "Strictwalk nodes .. \n";
		for(map<string, test_HH >::iterator HH_it = HH.begin(); HH_it != HH.end(); HH_it++){
			// moved here since we are counting from 1, there are both signs
			/*map<string, int> HH_tmp_second = (*HH_it).second;
			if(HH_tmp_second["tag"] != 0)
				continue;*/
			if(HH[(*HH_it).first].tag != 0)
				continue;
			if((*HH_it).first.length() < k)
				continue;
			Contignum += 1;
			string key = (*HH_it).first; // Xian: key is kmer
//if(key.compare("AAAGTAGCAAGATAGGCTGCCCCTC") == 0){
//int a = 0;
//}
			HH[key].tag = Contignum; // Xian: contig id
			string right, left;
			int rightsum, leftsum, righttype, lefttype;
			walk(&right, &rightsum, &righttype, key, "O", Contignum);
			walk(&left, &leftsum, &lefttype, key, "I", Contignum);
			string lefttype_str = itos(lefttype);
			string righttype_str = itos(righttype);
			//left = tl.reverse(left);
			reverse(left.begin(), left.end());
			Contigs[Contignum] = left + key + right;
//fh << "left: " << left << "\nkey:  " << key << "\nright:" << right << "\n\n";
			Contiglens[Contignum] = Contigs[Contignum].length();
			int tmp_ = Contiglens[Contignum] - k + 1;
			if(tmp_ == 0)
				tmp_ = 1;
			Contigcovs[Contignum] = (float)((int)100*(rightsum+leftsum+HH[key].n)/tmp_)/100; /////// seems unnecessary to multiply 100 and divide by 100
			Contigtypes[Contignum] = lefttype_str + righttype_str; // Xian: left branches number and right branches number
		}
		//cerr << " Num of Contigs: " << Contignum << endl;
//fh.close();
	}
	
	// do not check existance of nodes, always use solid nodes
	// break when out-degree != 1 or the in-degree of the next node > 1
	void walk(string *side, int *sum, int *sidetype, string node, string dir, int contignum){
		tools tl;
		string antidir = "I";
		if(dir.compare("I") == 0)
			antidir = "O";
		string nbase = "";
		int count = 0;
		if(dir.compare("I") == 0){
			if(HH[node].AI > 0){
				nbase = "A";
				count ++;
			}
			if(HH[node].CI > 0){
				nbase = "C";
				count ++;
			}
			if(HH[node].GI > 0){
				nbase = "G";
				count ++;
			}
			if(HH[node].TI > 0){
				nbase = "T";
				count ++;
			}
		}
		else{
			if(HH[node].AO > 0){
				nbase = "A";
				count ++;
			}
			if(HH[node].CO > 0){
				nbase = "C";
				count ++;
			}
			if(HH[node].GO > 0){
				nbase = "G";
				count ++;
			}
			if(HH[node].TO > 0){
				nbase = "T";
				count ++;
			}
		}

		//char ACGT[4] = {'A','C','G','T'};
		/*for(int i = 0; i < 4; i++){
			if(HH[node][ACGT[i] + dir] > 0){
				nbase = ACGT[i];
				count ++;
			}
		}*/
		if(count == 0 || count > 1){
			*side = "";
			*sum = 0;
			*sidetype = count;
			return;
		}
		string nnode, vnnode;
		if(dir.compare("O") == 0)
			nnode = node.substr(1, k - 1) + nbase;
		else 
			nnode = nbase + node.substr(0, k - 1);
		vnnode = nnode;
		reverse(vnnode.begin(), vnnode.end());
		//vnnode = tl.reverse(nnode);
		vnnode = tl.tr(vnnode, "ATGC", "TACG");
		if(HH.find(nnode) != HH.end()){
			if(HH[nnode].tag != 0){
				*side = "";
				*sum = 0;
				*sidetype = count;
				return;
			}
			int vcount = 0;
			if(antidir.compare("I") == 0){
				if(HH[nnode].AI > 0){
					//nbase = "A";
					vcount ++;
				}
				if(HH[nnode].CI > 0){
					//nbase = "C";
					vcount ++;
				}
				if(HH[nnode].GI > 0){
					//nbase = "G";
					vcount ++;
				}
				if(HH[nnode].TI > 0){
					//nbase = "T";
					vcount ++;
				}
			}
			else{
				if(HH[nnode].AO > 0){
					//nbase = "A";
					vcount ++;
				}
				if(HH[nnode].CO > 0){
					//nbase = "C";
					vcount ++;
				}
				if(HH[nnode].GO > 0){
					//nbase = "G";
					vcount ++;
				}
				if(HH[nnode].TO > 0){
					//nbase = "T";
					vcount ++;
				}
			}

			/*for(int i = 0; i < 4; i++){
				if(HH[nnode][ACGT[i] + antidir] > 0)
					vcount ++; 
			}*/
			if(vcount != 1){			// or > 1?? // don't know why
				*side = "";
				*sum = 0;
				*sidetype = count;
				return;
			}
			HH[nnode].tag = contignum;
			string rnext;
			int sumnext, endtype;
			walk(&rnext, &sumnext, &endtype, nnode, dir, contignum);
			*side = nbase + rnext;
			*sum = HH[nnode].n + sumnext; // Xian: how many reads in all (matters to only one root branch)
			*sidetype = endtype; // Xian: how many branches for closest to the root
			return;
		}
		else if(HH.find(vnnode) != HH.end()){
			if(HH[vnnode].tag != 0){
				*side = "";
				*sum = 0;
				*sidetype = count;
				return;
			}
			int vcount = 0;
			if(dir.compare("I") == 0){
				if(HH[vnnode].AI > 0){
					//nbase = "A";
					vcount ++;
				}
				if(HH[vnnode].CI > 0){
					//nbase = "C";
					vcount ++;
				}
				if(HH[vnnode].GI > 0){
					//nbase = "G";
					vcount ++;
				}
				if(HH[vnnode].TI > 0){
					//nbase = "T";
					vcount ++;
				}
			}
			else{
				if(HH[vnnode].AO > 0){
					//nbase = "A";
					vcount ++;
				}
				if(HH[vnnode].CO > 0){
					//nbase = "C";
					vcount ++;
				}
				if(HH[vnnode].GO > 0){
					//nbase = "G";
					vcount ++;
				}
				if(HH[vnnode].TO > 0){
					//nbase = "T";
					vcount ++;
				}
			}

			/*for(int i = 0; i < 4; i++){
				if(HH[vnnode][ACGT[i] + dir] > 0)
					vcount ++;
			}*/
			if(vcount != 1){			// or >1 ??
				*side = "";
				*sum = 0;
				*sidetype = count;
				return;
			}
			HH[vnnode].tag = - contignum;
			string rnext;
			int sumnext, endtype;
			walk(&rnext, &sumnext, &endtype, vnnode, antidir, -contignum);
			rnext = tl.tr(rnext, "ATGC", "TACG");
			*side = nbase + rnext;
			*sum = HH[vnnode].n + sumnext;
			*sidetype = endtype;
			return;
		}
	}

	// with both the info in contigs and the second key of the contig
	void dump_protocontigs(map<int, map<string, string> > &protocontigs) { // map<contigID, map<info, info_details> > : summarizing contigs
				
		for(int ii = 1; ii <= Contigs.size(); ii++){	// should counting from 1 since there is negative, compare with 0 /////// contigs.size() or contigs. - 1?
			string node = Contigs[ii].substr(0, k); /////// why previous with -1 but not this one?
			int dir;
			string truenode = true_func(node, &dir);
			int contig = HH[truenode].tag * dir * -1;
			map<int, int> a;
			nextcontigs_warc(contig, a);
			string i("I");
			string o("O");
			for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it ++){
				int con = (*a_it).first;
				i += itos(con) + ":" + itos(a[con]) + ",";
			}
//int tmp = Contiglens[ii];
			if(Contiglens[ii]-k < 0)
				continue;
			node = Contigs[ii].substr(Contiglens[ii] - k, k);
			truenode = true_func(node, &dir);
			contig = HH[truenode].tag*dir;
			a = map<int, int>();
			nextcontigs_warc(contig, a);
			for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it ++){
				int con = (*a_it).first;
				o += itos(con) + ":" + itos(a[con]) + ",";
			}
			if(Contigtags.find(ii) == Contigtags.end())
				Contigtags[ii] = 0;
			map<string, string> protocontig;
			// be careful it may rewrite the previous information for the stringstream
			protocontig["id"] = itos(ii);
			
			protocontig["lens"] = itos(Contiglens[ii]);
			
			protocontig["covs"] = ftos(Contigcovs[ii]);
			
			protocontig["types"] = Contigtypes[ii];
			
			protocontig["tags"] = itos(Contigtags[ii]);
			protocontig["I"] = i;
			protocontig["O"] = o;
			protocontig["seq"] = Contigs[ii];
			protocontigs[ii] = protocontig;
		}
	}
				
	// nextcontigs_warc is to get the second key next to the key: a = map<direciton*contigID, # of next key>
	void nextcontigs_warc(int contignum, map<int, int> &a){
		// return a hash containing the connections between the query proto-contig and its connected proto-contigs with the thickness of the edges
		string node;
		if(contignum > 0)
			node = Contigs[contignum].substr(Contiglens[contignum] - k, k);
		else {
			node = Contigs[-contignum].substr(0, k);
			node = revcom(node);
		}
		//char ACGT[4] = {'A','C','G','T'};
		int dir;
		string O = "O";
		string I = "I";
		if(HH.find(node) != HH.end()){
			string node_tmp = node.substr(1, k-1);
			string truenode;
			string vnode;
			if(HH[node].AO != 0){
				vnode = node_tmp + "A";
				truenode = true_func(vnode, &dir);
				a[dir*HH[truenode].tag] = HH[node].AO;
			} 
			if(HH[node].CO != 0){
				vnode = node_tmp + "C";
				truenode = true_func(vnode, &dir);
				a[dir*HH[truenode].tag] = HH[node].CO;
			} 
			if(HH[node].GO != 0){
				vnode = node_tmp + "G";
				truenode = true_func(vnode, &dir);
				a[dir*HH[truenode].tag] = HH[node].GO;
			} 
			if(HH[node].TO != 0){
				vnode = node_tmp + "T";
				truenode = true_func(vnode, &dir);
				a[dir*HH[truenode].tag] = HH[node].TO;
			} 

			/*for(int i = 0; i < 4; i++){
				if(HH[node][ACGT[i] + O] == 0 || HH[node].find(ACGT[i] + O) == HH[node].end())
					continue;
				string vnode = node.substr(1, k - 1) + ACGT[i];
				string truenode = true_func(vnode, &dir);
				a[dir*HH[truenode]["tag"]] = HH[node][ACGT[i] + O];
			}*/
		}
		else {
			string rnode = revcom(node);
			string rnode_tmp = rnode.substr(0, k-1);
			string truenode;
			string vnode;
			if(HH[rnode].AI != 0){
				vnode = "A" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode].tag] = HH[rnode].AI;
			}
			if(HH[rnode].CI != 0){
				vnode = "C" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode].tag] = HH[rnode].CI;
			}
			if(HH[rnode].GI != 0){
				vnode = "G" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode].tag] = HH[rnode].GI;
			}
			if(HH[rnode].TI != 0){
				vnode = "T" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode].tag] = HH[rnode].TI;
			}

			/*for(int i = 0; i < 4; i++){
				if(HH[rnode][ACGT[i] + I] == 0 || HH[rnode].find(ACGT[i] + I) == HH[rnode].end())
					continue;
				string vnode = ACGT[i] + rnode.substr(0, k - 1);
				string truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode]["tag"]] = HH[rnode][ACGT[i] + I];
			}*/
		}
	}

	string revcom(string seq){
		tools tl;
		//seq = tl.reverse(seq);
		reverse(seq.begin(), seq.end());
		seq = tl.tr(seq, "ATGC", "TACG");
		return seq;
	}
	
	string true_func(string node, int *dir){
		if(HH.find(node) != HH.end()){
			*dir = 1;
			return node;
		}
		else {
			*dir = -1;
			return revcom(node);
		}
	}
};

// track from the end of a contig which does not have any branches, and see how long it could go: Contigtips[Contignum]=TipLength
class processtips {

public:
	
	int K;// = 25;
	int Tipcount;// = 0;
	int Thin;// = 4;
	int DefaultTip;// = 1000;
	map<string, test_HH > HH; // map<kmer, map<info, count> >
	map<int, map<string, string> > Contigs; // map<contigID, map<info, details> >
	map<int, int> Contiglens;
	map<int, float> Contigcovs;
	map<int, string> Contigtypes;
	int Contignum;// = 0;
	map<int, int> Contigtips; // map<contigID, contig
	map<int, map<string, string> > Contigs2;
	
	processtips(){//int K, int Thin, int Tipcount, int DefaultTip){
		K = 25;//K;
		Thin = 4;//Thin;
		Tipcount = 0;//Tipcount;
		DefaultTip = 1000;//DefaultTip;
		Contignum = 0;
	}
	
	void set_kmer_size(int k_){
		K = k_;
	}
	
	~processtips(){
		HH = map<string, test_HH >();
		Contigs = map<int, map<string, string> >();
		Contiglens = map<int, int>();
		Contigcovs = map<int, float>();
		Contigtypes = map<int, string>();
		Contigtips = map<int, int>();
		Contigs2 = map<int, map<string, string> >();
	}
	
	void tipwrap(int cutoff, map<int, map<string, string> > &rContigs, map<string, test_HH> &rHH){
		HH = rHH;
		Contigs = rContigs;
		//cerr << "labeling tips up to " << cutoff << " .. \n" ;
		for (int i = 1; i <= Contigs.size(); i++){
			map<string, string> contig = Contigs[i];
			assert(contig["types"].length()-3);// be careful that it might be three digits
			assert(contig["types"].length()-4);// be careful that it might be four digits
			int left = atoi(contig["types"].substr(0,1).c_str());
			int right = atoi(contig["types"].substr(1,1).c_str());
			if(left == 0)
				labeltip(-i, K-1, cutoff);
			if(right == 0)
				labeltip(i, K-1, cutoff);
		}
		//return Contigtips;
	}
	
	// to see the longest tip + contig's length
	void labeltip(int contignum, int dist, int tipcutoff){
		dist += atoi(Contigs[abs(contignum)]["lens"].c_str()) - K + 1;
		if(Contigtips.find(contignum) == Contigtips.end())
			Contigtips[contignum] = DefaultTip;
		if(dist >= tipcutoff || Contigtips[contignum] <= dist)
			return;		
		Contigtips[contignum] = dist;
		Tipcount ++;
		vector<int> nextcontigs;
		nextcontigs_func(-contignum, nextcontigs);
		for(int k = 0; k < nextcontigs.size(); k++){
			vector<int> ncon;
			nextcontigs_func(-nextcontigs[k], ncon);
			int indist = dist;
			for (int j = 0; j < ncon.size(); j++) {
				int i = ncon[j];
				if(Contigtips.find(i) == Contigtips.end())
					Contigtips[i] = DefaultTip;
				if(Contigtips[i] > indist)
					indist = Contigtips[i]; // get the longest contig => indist
			}
			if(indist < tipcutoff - 1) // to see when to stop label tips, depend on the tip length
				labeltip(-nextcontigs[k], indist, tipcutoff);
		}
	}	
	
	void breaktip(int cutoff, map<string, test_HH> &rHH, map<int, int> &rcontigtips, map<int, map<string, string> > &rcontigs, map<int, map<string, string> > &rcontigs2){
		HH = rHH;
		Contigtips = rcontigtips;
		Contigs = rcontigs;
		Contigs2 = rcontigs2;
		
		//cerr << "Breaking tips upto " << cutoff << " .. \n";
		Tipcount = 0;
		for (int i = 1; i <= Contigs2.size(); i++) {
			string end = contigend2(i);
			int dir;
			string trueend = true_func(end, &dir);
			int si = dir*HH[trueend].tag;
			vector<int> n;
			nextcontigs_func(si, n);
			int dist = K - 1;
			int tag = 0;
			if(n.size() > 0){
				for (int k = 0; k < n.size(); k++) {
					int j = n[k];
					if(Contigtips.find(j) == Contigtips.end())
						Contigtips[j] = DefaultTip;
					if (contigend(-j).compare(contigend2(-(j/abs(j))*atoi(Contigs[abs(j)]["tags"].c_str()))) == 0 ) {
						if(dist < Contigtips[j])
						   dist = Contigtips[j];
					}
					else 
						tag = 1;
				}
				if(tag == 0)
					continue;
				if(dist < cutoff - 1)
					labeltip(si, dist, cutoff);
			}
			end = contigend2(-i);
			trueend = true_func(end, &dir);
			si = dir*HH[trueend].tag;
			n = vector<int>();
			nextcontigs_func(si, n);
			dist = K - 1;
			tag = 0;
			if(n.size() > 0){
				for (int k = 0; k < n.size(); k++) {
					int j = n[k];
					if(Contigtips.find(j) == Contigtips.end())
						Contigtips[j] = DefaultTip;
					if (contigend(-j).compare(contigend2(-(j/abs(j))*atoi(Contigs[abs(j)]["tags"].c_str()))) == 0 ) {
						if(dist < Contigtips[j])
							dist = Contigtips[j];
					}
					else 
						tag = 1;
				}
				if(tag == 0)
					continue;
				if(dist < cutoff - 1)
					labeltip(si, dist, cutoff);
			}
		}
		//cerr << "Tip count: " << Tipcount << endl;
		//return &(Contigtips);
	}
	
	void thicken(float cutoff, map<string, test_HH > &rHH, map<int, int> &rcontigtips, map<int, map<string, string> > &rcontigs, map<int, map<string, string> > &rcontigs2){
		HH = rHH;
		Contigtips = rcontigtips;
		Contigs = rcontigs;
		Contigs2 = rcontigs2;
		//cerr << "Thicken with cutoff " << cutoff << " ..\n";
		int Smalltip = 100;
		string I = "I";
		string O = "O";
		int size = Contigs2.size();
		for(int i = -size; i <= size; i++){
			if(i == 0)
				continue;
			if(atof(Contigs2[abs(i)]["covs"].c_str()) < cutoff || Contigs2[abs(i)]["types"].compare("11") != 0)
				continue;
			string end = contigend2(i);
			int dir;
			string trueend = true_func(end, &dir);
			int si = dir*HH[trueend].tag;
			if(Contigtips.find(-si) == Contigtips.end())
				Contigtips[-si] = DefaultTip;
			if(Contigtips[-si] < Smalltip)
				continue;
			char ACGT[4] = {'A','C','G','T'};
			if(dir == 1){
				int big = 0;
				string base("");
				/*for(int j = 0; j < 4; j++){
					if(HH[trueend][ACGT[j] + O] > big){
						big = HH[trueend][ACGT[j] + O];
						base = ACGT[j];
					}
				}*/
				if(HH[trueend].AO > big){
					big = HH[trueend].AO;
					base = "A";
				}
				if(HH[trueend].CO > big){
					big = HH[trueend].CO;
					base = "C";
				}
				if(HH[trueend].GO > big){
					big = HH[trueend].GO;
					base = "G";
				}
				if(HH[trueend].TO > big){
					big = HH[trueend].TO;
					base = "T";
				}
				for(int j = big; j <= Thin; j++)
					addstring(trueend + base);
			}
			else {
				int big = 0;
				string base("");
				/*for(int j = 0; j < 4; j++){
					if(HH[trueend][ACGT[j] + I] > big){
						big = HH[trueend][ACGT[j] + I];
						base = ACGT[j];
					}
				}*/

				if(HH[trueend].AI > big){
					big = HH[trueend].AI;
					base = "A";
				}
				if(HH[trueend].CI > big){
					big = HH[trueend].CI;
					base = "C";
				}
				if(HH[trueend].GI > big){
					big = HH[trueend].GI;
					base = "G";
				}
				if(HH[trueend].TI > big){
					big = HH[trueend].TI;
					base = "T";
				}
				for(int j = big; j <= Thin; j++)
					addstring(base + trueend);
			}
		}
		// return map<string, map<string, int> > &HH;
	}

	void addstring(string rd){
		int l = rd.length();
		tools tl;
		string I = "I";
		string O = "O";
		for(int i = K; i <= l -1; i++){
			string w = rd.substr(i - K, K);
			string rw = w;
			reverse(rw.begin(), rw.end());
			//string rw = tl.reverse(w);
			rw = tl.tr(rw, "ATGC", "TACG");
			if(HH.find(w) != HH.end()){
				string u = rd.substr(i - K + 1, K);
				string ru = u;
				reverse(ru.begin(), ru.end());
				//string ru = tl.reverse(u);
				ru = tl.tr(ru, "ATGC", "TACG");
				if(HH.find(u) != HH.end()){
					//string next = u.substr(K-1, 1);
					char next = u[K-1];
					if(next == 'A')
						HH[w].AO ++;
					else if(next == 'C')
						HH[w].CO ++;
					else if(next == 'G')
						HH[w].GO ++;
					else if(next == 'T')
						HH[w].TO ++;
					//HH[w][next + O] ++;
					next = w[0];
					if(next == 'A')
						HH[u].AI ++;
					else if(next == 'C')
						HH[u].CI ++;
					else if(next == 'G')
						HH[u].GI ++;
					else if(next == 'T')
						HH[u].TI ++;
					//next = w.substr(0, 1);
					//HH[u][next + I] ++;
				}
				else if(HH.find(ru) != HH.end()){
					//string next = u.substr(K-1, 1);
					char next = u[K-1];
					if(next == 'A')
						HH[w].AO ++;
					else if(next == 'C')
						HH[w].CO ++;
					else if(next == 'G')
						HH[w].GO ++;
					else if(next == 'T')
						HH[w].TO ++;
					//HH[w][next + O] ++;
					next = rw[K-1];
					if(next == 'A')
						HH[ru].AO ++;
					else if(next == 'C')
						HH[ru].CO ++;
					else if(next == 'G')
						HH[ru].GO ++;
					else if(next == 'T')
						HH[ru].TO ++;
					//next = rw.substr(K - 1, 1);
					//HH[ru][next + O] ++;
				}
				else {
					i++;
				}
			}
			else if(HH.find(rw) != HH.end()){
				string u = rd.substr(i-K+1, K);
				string ru = u;
				reverse(ru.begin(), ru.end());
				//string ru = tl.reverse(u);
				ru = tl.tr(ru, "ATGC", "TACG");
				if(HH.find(u) != HH.end()){
					//string next = ru.substr(0,1);
					char next = ru[0];
					if(next == 'A')
						HH[rw].AI ++;
					else if(next == 'C')
						HH[rw].CI ++;
					else if(next == 'G')
						HH[rw].GI ++;
					else if(next == 'T')
						HH[rw].TI ++;
					//HH[rw][next + I] ++;
					next = w[0];
					if(next == 'A')
						HH[u].AI ++;
					else if(next == 'C')
						HH[u].CI ++;
					else if(next == 'G')
						HH[u].GI ++;
					else if(next == 'T')
						HH[u].TI ++;
					//next = w.substr(0,1);
					//HH[u][next + I] ++;
				}
				else if(HH.find(ru) != HH.end()){
					//string next = ru.substr(0,1);
					char next = ru[0];
					if(next == 'A')
						HH[rw].AI ++;
					else if(next == 'C')
						HH[rw].CI ++;
					else if(next == 'G')
						HH[rw].GI ++;
					else if(next == 'T')
						HH[rw].TI ++;
					//HH[rw][next + I] ++;
					next = rw[K-1];
					if(next == 'A')
						HH[ru].AO ++;
					else if(next == 'C')
						HH[ru].CO ++;
					else if(next == 'G')
						HH[ru].GO ++;
					else if(next == 'T')
						HH[ru].TO ++;
					//next = rw.substr(K-1,1);
					//HH[ru][next + O] ++;
				}
				else {
					i++;
				}
			}
		}
	}
	
	// Get the head or tail of the contig, get the next node of it, push the number of the next node into the vector
	void nextcontigs_func(int contignum, vector<int> &a){
		string node;
		if(contignum > 0)
			node = Contigs[contignum]["seq"].substr(atoi(Contigs[contignum]["lens"].c_str()) - K, K);
		else {
			node = Contigs[-contignum]["seq"].substr(0, K);
			node = revcom(node);
		}
		//char ACGT[4] = {'A','C','G','T'};
		string I = "I";
		string O = "O";
		if(HH.find(node) != HH.end()){
			string node_tmp = node.substr(1, K-1);
			int dir;
			string truenode;
			if(HH[node].AO != 0){
				string vnode = node_tmp + "A";
				truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			if(HH[node].CO != 0){
				string vnode = node_tmp + "C";
				truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			if(HH[node].GO != 0){
				string vnode = node_tmp + "G";
				truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			if(HH[node].TO != 0){
				string vnode = node_tmp + "T";
				truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			/*for(int i = 0; i < 4; i++){
				if (HH[node][ACGT[i] + O] == 0) // may mean it's empty??
					continue;
				string vnode = node.substr(1, K - 1) + ACGT[i];
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}*/
		}
		else {
			string rnode = revcom(node);
			string rnode_tmp = rnode.substr(0, K-1);
			int dir;
			string truenode;
			if(HH[rnode].AI != 0){
				string vnode = "A" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}
			if(HH[rnode].CI != 0){
				string vnode = "C" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}
			if(HH[rnode].GI != 0){
				string vnode = "G" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}
			if(HH[rnode].TI != 0){
				string vnode = "T" + rnode_tmp;
				truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}

			/*for (int i = 0; i < 4; i++) {
				if(HH[rnode][ACGT[i] + I] == 0) // may mean it's empty?
					continue;
				string vnode = ACGT[i] + rnode.substr(0, K-1);
				int dir;
				string truenode = true_func(vnode, &dir);
//int tmp = HH[truenode]["tag"];
				a.push_back(-dir*HH[truenode].tag);
			}*/
		}
	}


	string revcom(string seq){
		tools tl;
		reverse(seq.begin(), seq.end());
		//seq = tl.reverse(seq);
		seq = tl.tr(seq, "ATGC", "TACG");
		return seq;
	}
	
	string true_func(string node, int *i){
		string rtn;
		if(HH.find(node) != HH.end()){
			*i = 1;
			rtn = node;
		}
		else{
			*i = -1;
			rtn = revcom(node);
		}
		return rtn;
	}
	
	string contigend(int contig){
		string rtn;
		if(contig > 0)
			rtn = Contigs[contig]["seq"].substr(atoi(Contigs[contig]["lens"].c_str()) - K, K);
		else 
			rtn = revcom(Contigs[-contig]["seq"].substr(0, K));
		return rtn;
	}
	
	string contigend2(int contig){
		string rtn;
		if(contig > 0){
			if(Contigs2.find(contig) != Contigs2.end())
				rtn = Contigs2[contig]["seq"].substr(atoi(Contigs2[contig]["lens"].c_str()) - K, K);
			else
				rtn = "";
		}
		else {
			if(Contigs2.find(-contig) != Contigs2.end())
				rtn = revcom(Contigs2[-contig]["seq"].substr(0, K));
			else
				rtn = "";
		}
		return rtn;
	}
};


// recover low frequency kmers in high quality reads that bridge separated non-tip proto-contig graphs
class addbridgekmer {
	
public:
	
	//map<string, map<string, int> > HH;
	map<string, test_HH> HH;
	int Rdnum;// = 0;
	map<int, map<string, string> > Contigs;
	int Contignum;// = 0;
	map<int, int> Contigtips;
	vector<string> newKmers;
	int k;// = 25;
	int Bridgeanchor;// = 1;
	int DefaultTip;// = 1000;
	
	addbridgekmer(){//, int Bridgeanchor, int DefaultTip, map<string, map<string, int> > &graph){
		k = 25;
		Bridgeanchor = 1;//Bridgeanchor;
		DefaultTip = 1000;//DefaultTip;
		Rdnum = 0;
		Contignum = 0;
		//if(graph.size() > 0)
		//HH = graph;
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	void set_HH(map<string, test_HH > &graph){
		if(graph.size() > 0)
			HH = graph;
	}
	
	~addbridgekmer(){
		//HH = map<string, map<string, int> >();
		HH = map<string, test_HH>();
		Contigs = map<int, map<string, string> >();
		Contigtips = map<int, int>();
		newKmers = vector<string>();
	}
	
	//void doit(map<string, map<string, int> > &rHH, map<int, map<string, string> > &rContigs, map<int, int> &rContigtips, vector<string> &rPR, vector<string> &rReads){
	void doit(map<string, test_HH > &rHH, map<int, map<string, string> > &rContigs, map<int, int> &rContigtips, vector<string> &rPR, vector<string> &rReads){
		HH = rHH;
		Contigs = rContigs;
		Contigtips = rContigtips;
		vector<string> PR = rPR;
		for(int i = 0; i < rReads.size(); i++){
			if(i >= PR.size())
                            continue;
                        string nums = PR[i];
			addbridge(rReads[i], nums);
		}
		// return map<string, map<string, int> > HH; vector<string> newKmers;
	}
	
	void addbridge(string rd, string nums){
		tools tl;
		vector<string> aqual;
		/*int begin = 0;
		int stop = 0;
		while (stop == 0) {
			size_t end = nums.find(" ", begin);
			if(end == string::npos){
				end = nums.length()-1;
				stop = 1;
			}
			int nums_ = atoi(nums.substr(begin, end - begin));
			aqual.push_back(nums_);
			begin = end + 1;
		}*/
		string sep = " ";
		tl.split(nums, sep, aqual);
		int lefts = -1;
		int lefte = 0;
		int leftsum = 0;
		int leftbig = 0;
		int leftbigpos = 0;
		for(int i = 0; i < aqual.size(); i++){
			int aqual_ = atoi(aqual[i].c_str());
			if(aqual_ > 1){
				if(lefts == -1)
					lefts = i;
				leftsum += aqual_;
				if(leftbig < aqual_){
					leftbig = aqual_;
					leftbigpos = i;
				}
			}
			else if(lefts > -1){ // enter into non breaking bridge before
				if(leftsum > Bridgeanchor){ // the previous bridge strength is strong
					lefte = i; // this was broken, now it may be amended
					break;
				}
				else {
					lefts = -1;
					lefte = 0;
					leftsum = 0;
					leftbig = 0;
					leftbigpos = 0;
				}
			}
		}
		if(lefte == 0)
			return;
		
		int rights = -1;
		int righte = 0;
		int rightsum = 0;
		int rightbig = 0;
		int rightbigpos = 0;
		for(int i = aqual.size() - 1; i >= lefte; i--){
			int aqual_ = atoi(aqual[i].c_str());
			if(aqual_ > 1){
				if(rights == -1)
					rights = i;
				rightsum += aqual_;
				if(rightbig < aqual_){
					rightbig = aqual_;
					rightbigpos = i;
				}
			}
			else if(rights > -1){
				if(rightsum > Bridgeanchor){
					righte = i;
					break;
				}
				else{
					rights = -1;
					righte = 0;
					rightsum = 0;
					rightbig = 0;
					rightbigpos = 0;
				}
			}
		}
		if(righte == 0)
			return;
		
		if(righte - lefte > k -2)	// do not allow real single cov bridge
			return; // try to add the kmer back if the left and right are strong, and they can be connected with this kmer
			
		if(leftbigpos >= rd.length() || leftbigpos < 0) // to protect substr
			return;
		string node = rd.substr(leftbigpos, k);
		int dir;
		string truenode = true_func(node, &dir);
		int contigl = dir * HH[truenode].tag;
		if(rightbigpos >= rd.length() || rightbigpos < 0) // to protect substr
			return;
		node = rd.substr(rightbigpos, k);
		node = revcom(node);
		truenode = true_func(node, &dir);
		int contigr = dir * HH[truenode].tag;
		if(Contigtips.find(contigl) == Contigtips.end())
			Contigtips[contigl] = DefaultTip;
		if(Contigtips.find(contigr) == Contigtips.end())
			Contigtips[contigr] = DefaultTip;
		if(Contigtips[contigl] - atoi(Contigs[abs(contigl)]["lens"].c_str()) < 100 || Contigtips[contigr] - atoi(Contigs[abs(contigr)]["lens"].c_str()) < 100) ///// need to write a atoi override function from string
			addback(rd, lefte, righte);
	}

	void addback(string rd, int le, int re){
		if(le-1 < 0 || le-1 >= rd.length()) // to protect substr
			return;
		string u = rd.substr(le-1, k);
		string w;
		///// is le always < re? Right!
		for(int i = le; i <= re + 1; i++){
			w = u;
			if(i < 0 || i >= rd.length()) // to protect substr
				continue;
			u = rd.substr(i, k);
			if(u.find("N") != string::npos || w.find("N") != string::npos)
				continue;
			if(! (HH.find(u) != HH.end() || HH.find(revcom(u)) != HH.end())){
				HH[u].n = 1;
				HH[u].AI = 0;
				HH[u].CI = 0;
				HH[u].GI = 0;
				HH[u].TI = 0;
				HH[u].AO = 0;
				HH[u].CO = 0;
				HH[u].GO = 0;
				HH[u].TO = 0;
				HH[u].tag = 0;
				HH[u].tag2 = 0;
				newKmers.push_back(u); // include all those kmers spanning [le, re]

			}
			addedge(w,u);
		}
	}

	void addedge(string a, string b){
		int dira, dirb;
		string truea = true_func(a, &dira);
		string trueb = true_func(b, &dirb);
		if(dira == 1){
			if(dirb == 1){
				if(b.length() < k)
					return;
				char next = b[k-1];
				if(next == 'A')
					HH[a].AO ++;
				else if(next == 'C')
					HH[a].CO ++;
				else if(next == 'G')
					HH[a].GO ++;
				else if(next == 'T')
					HH[a].TO ++;
				//HH[a][next + "O"] ++;
				next = a[0];				
				if(next == 'A')
					HH[b].AI ++;
				else if(next == 'C')
					HH[b].CI ++;
				else if(next == 'G')
					HH[b].GI ++;
				else if(next == 'T')
					HH[b].TI ++;
				//HH[b][next + "I"] ++;
			}
			else{
				char next = b[k-1];
				if(b.length() < k)
					return;
				if(next == 'A')
					HH[a].AO ++;
				else if(next == 'C')
					HH[a].CO ++;
				else if(next == 'G')
					HH[a].GO ++;
				else if(next == 'T')
					HH[a].TO ++;
				//HH[a][next + "O"] ++;
				string tmp = revcom(a);
				next = tmp[k-1];
				if(next == 'A')
					HH[trueb].AO ++;
				else if(next == 'C')
					HH[trueb].CO ++;
				else if(next == 'G')
					HH[trueb].GO ++;
				else if(next == 'T')
					HH[trueb].TO ++;
				//HH[trueb][next + "O"] ++;
			}
		}
		else{
			if(dirb == 1){
				string tmp = revcom(b);
				char next = tmp[0];
				if(next == 'A')
					HH[truea].AI ++;
				else if(next == 'C')
					HH[truea].CI ++;
				else if(next == 'G')
					HH[truea].GI ++;
				else if(next == 'T')
					HH[truea].TI ++;
				//HH[truea][next + "I"] ++;
				next = a[0];				
				if(next == 'A')
					HH[b].AI ++;
				else if(next == 'C')
					HH[b].CI ++;
				else if(next == 'G')
					HH[b].GI ++;
				else if(next == 'T')
					HH[b].TI ++;
				//HH[b][next + "I"] ++;
			}
			else{
				string tmp = revcom(b);
				char next = tmp[0];
				if(next == 'A')
					HH[truea].AI ++;
				else if(next == 'C')
					HH[truea].CI ++;
				else if(next == 'G')
					HH[truea].GI ++;
				else if(next == 'T')
					HH[truea].TI ++;
				//HH[truea][next + "I"] ++;
				tmp = revcom(a);
				if(tmp.length() < k)
					return;
				next = tmp[k-1];			
				if(next == 'A')
					HH[trueb].AO ++;
				else if(next == 'C')
					HH[trueb].CO ++;
				else if(next == 'G')
					HH[trueb].GO ++;
				else if(next == 'T')
					HH[trueb].TO ++;
				//HH[trueb][next + "O"] ++;
			}
		}
	}
	
	string revcom(string seq){
		tools tl;
		reverse(seq.begin(), seq.end());
		//seq = tl.reverse(seq);
		seq = tl.tr(seq, "ATGC", "TACG");
		return seq;
	}
	
	string true_func(string node, int *i){
		*i = 0;
		string rtn("0");
		if(HH.find(node) != HH.end()){
			*i = 1;
			rtn = node;
		}
		else{
			string rnode = revcom(node);
			if(HH.find(rnode) != HH.end()){
				*i = -1;
				rtn = rnode;
			}
			else {
				*i = 0;
				rtn = "0";
			}
		}
		return rtn;	// not guarantee this is what is wanted for *i
	}

	void nextbase(string node, int *A, int *T, int *G, int *C){
		int dir;
		string truenode = true_func(node, &dir);
		if(dir == 1){
			*A = HH[truenode].AO;
			*T = HH[truenode].TO;
			*G = HH[truenode].GO;
			*C = HH[truenode].CO;
		}
		else {
			*A = HH[truenode].TO;
			*T = HH[truenode].AO;
			*G = HH[truenode].CO;
			*C = HH[truenode].GO;
		}
	}
	
	int min_vector(vector<int> &list){
		int min = list[list.size() - 1];
		list.pop_back();
		for(int i = 0; i < list.size(); i++){
			if(list[i] < min)
				min = list[i];
		}
		return min;
	}
};

// extend proto-contigs to contigs by removing tips and collapse bubbles with heuristic cut-offs
class walkcontig {

public:

	//map<string, map<string, int> > HH;
	map<string, test_HH> HH;
	map<int, int> Contigtips;
	map<int, string> Contigs;
	map<int, int> Contiglens;
	map<int, float> Contigcovs;
	map<int, string> Contigtypes;
	uint32_t Contignum;
	map<int, int> Contigtags;

	map<int, string> Contigs2;
	map<int, int> Contiglens2;
	map<int, float> Contigcovs2;
	map<int, string> Contigtypes2;
	int Contignum2;
	float Ratiocutoff;// = 0.3;
	int Convergestepwall;// = 12;
	int ConvergeLengthwall;// = 150;

	int k;// = 25;
	int Thin;// = 2;
	int Smalltip;// = 100;
	int Walkcutoff;// = 3;
	int Tipcount;// = 0;
	int DefaultTip;// = 1000;


	walkcontig(){//, int Thin, int Smalltip, int Walkcutoff, int Tipcount, int DefaultTip){
		k = 25;
		Thin = 2;
		Smalltip = 100;
		Walkcutoff = 3;
		Tipcount = 0;
		DefaultTip = 1000;
		Ratiocutoff = 0.3;
		Convergestepwall = 12;
		ConvergeLengthwall = 150;

	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	~walkcontig(){
		//HH = map<string, map<string, int> >();
		HH = map<string, test_HH>();
		Contigtips = map<int, int>();
		Contigs = map<int, string>();
		Contiglens = map<int, int>();
		Contigcovs = map<int, float>();
		Contigtypes = map<int, string>();
		Contigtags = map<int, int>();
		Contigs2 = map<int, string>();
		Contiglens2 = map<int, int>();
		Contigcovs2 = map<int, float>();
		Contigtypes2 = map<int, string>();
		// Contignum, Contignum2, Ratiocutoff, Convergestepwall, Convergelenthwall
	}
	
	void walkcontigwrap(int cutoff_, float ratiocutoff_, map<string, test_HH> &HH_, map<int, int> &tips, map<int, map<string, string> > &in_contigs){
//            cout << "Before doing walkcontig:\n";
		tools tl;
		int cutoff;
		float ratiocutoff;
		if(cutoff_ != -1)
			cutoff = cutoff_;
		else 
			cutoff = Walkcutoff;

		if(ratiocutoff_ != -1)
			ratiocutoff = ratiocutoff_;
		else
			ratiocutoff = Ratiocutoff; // need to deal with this better by separating the assignment
		
		//cout << "Walking contigs with cutoff " << cutoff << ", ratio_cutoff " << ratiocutoff << "..\n";
		HH = HH_;
		Contigtips = tips;
		
		int i;
		for(i = 1; i <= in_contigs.size(); i++){
			Contiglens[i] = atoi(in_contigs[i]["lens"].c_str());
			Contigcovs[i] = atof(in_contigs[i]["covs"].c_str());
			Contigtypes[i] = in_contigs[i]["types"];
			Contigtags[i] = atoi(in_contigs[i]["tags"].c_str());
			Contigs[i] = in_contigs[i]["seq"];
		}
		
		Contignum2 = 0;
		vector<int> index;
//                cout << "Problem in reverse key:\n";
		tl.sort_value_reverse_key_float(Contigcovs, index, 1, Contigs.size()); // need to write this in tools
//                cout << "Before walkcontig loop" << index.size() << endl;
		for(int ii = 0; ii < index.size(); ii++){
//                    cout << ii << " out of " << index.size() << endl;
			int i = index[ii];
			if(Contigtags[i] != 0)
				continue;

			Contignum2 ++;
			Contigtags[i] = Contignum2;
			string left, right;
			int leftsum, rightsum, lefttype, righttype;
			walkcontig_doit(i, Contignum2, cutoff, ratiocutoff, &right, &rightsum, &righttype);

			walkcontig_doit(-i, -Contignum2, cutoff, ratiocutoff, &left, &leftsum, &lefttype);
			left = revcom(left);
			Contigs2[Contignum2] = left + Contigs[i] + right;
			Contiglens2[Contignum2] = Contigs2[Contignum2].length();
			Contigcovs2[Contignum2] =(float)(int(100*(rightsum + leftsum + (Contiglens[i] - k + 1)*Contigcovs[i]))/int(Contiglens2[Contignum2] - k +1))/(float)100;
			Contigtypes2[Contignum2] = itos(lefttype) + itos(righttype);
			int j;
			for (j = k; j <= Contiglens2[Contignum2]; j++) {
				
				string node = Contigs2[Contignum2].substr(j-k, k);
				int dir;
				string truenode = true_func(node, &dir);
				HH[truenode].tag2 = dir*Contignum2; 
			}
		}

		for (i = 1; i <= Contigs.size(); i++)
			in_contigs[i]["tags"] = itos(Contigtags[i]);
		
		//cerr << "Done walking contigs, number of contigs2: " << Contignum2 << endl;
	}

	void walkcontig_doit(int contig, int contiglab, int cutoff, float ratiocutoff, string *side, int *sidesum, int *sidetype){
		vector<int> ncontigs;
		step(contig, cutoff, ratiocutoff, ncontigs);
		int nx = 0;
		if(ncontigs.size() != 0){
			nx = ncontigs.back();
			ncontigs.pop_back();
		}
		if(nx == 0 || Contigtags[abs(ncontigs[0])] != 0){
			*side = "";
			*sidesum = 0;
			*sidetype = nx;
			return;
		}
		int ncontig = 0;
		if(ncontigs.size() > 1){
			for(int i = 0; i < ncontigs.size(); i++){
				int ncontig = ncontigs[i];
				vector<int> vcontigs;
				step(-ncontig, cutoff, ratiocutoff, vcontigs);
				vcontigs.pop_back();
				if(vcontigs.size() != 1 || vcontigs[0]!= - contig){
					*side = "";
					*sidesum = 0;
					*sidetype = nx;
					return;
				}
			}
			int y;
			converge(ncontigs, &ncontig, &y);
		}
		else {
			int pseudo;
			int convergeon;
			pseudowalkcontig(-ncontigs[0], cutoff, ratiocutoff, &pseudo, &convergeon);
			if(pseudo == -contig || convergeon != 0 && abs(convergeon)/convergeon*Contigtags[abs(convergeon)] == -contiglab)
				ncontig = ncontigs[0];
		}
		if(ncontig == 0){
			*side = "";
			*sidesum = 0;
			*sidetype = nx;
			return;
		}
		string seq;
		if(ncontig > 0){
			Contigtags[ncontig] = contiglab;
			if(k-1 < Contigs[ncontig].length() && k-1 >= 0) // to protect substr
				seq = Contigs[ncontig].substr(k - 1, Contiglens[ncontig] - k + 1);
		}
		else {
			Contigtags[-ncontig] = - contiglab;
			seq = Contigs[-ncontig].substr(0, Contiglens[-ncontig] - k + 1);
			seq = revcom(seq);										   
		}
		string rnext;
		int sumnext;
		int endtype;
		walkcontig_doit(ncontig, contiglab, cutoff, ratiocutoff, &rnext, &sumnext, &endtype);
		*side = seq + rnext;
		*sidesum = (Contiglens[abs(ncontig)] - k + 1)*Contigcovs[abs(ncontig)]+sumnext;
		*sidetype = endtype;
		return;
	}
	
	void pseudowalkcontig(int contig, int cutoff, float ratiocutoff, int *ncontig, int *convergeon){
		vector<int> ncontigs;
		step(contig, cutoff, ratiocutoff, ncontigs);
		int nx = 0;
		if(ncontigs.size() != 0){
		    nx = ncontigs.back();
			ncontigs.pop_back();
		}
		if(nx == 0){
			*ncontig = 0;
			*convergeon = 0;
			return;
		}
		for (vector<int>::iterator ncontigs_it = ncontigs.begin(); ncontigs_it != ncontigs.end(); ncontigs_it ++) {
			*ncontig = *ncontigs_it;
			vector<int> vcontigs;
			step(-(*ncontig), cutoff, ratiocutoff, vcontigs);
			vcontigs.pop_back();
			if(vcontigs.size() != 1 || vcontigs[0] != -contig){
				*ncontig = 0;
				*convergeon = 0;
				return;
			}
		}
		converge(ncontigs, ncontig, convergeon);
		return;
	}
		
	void converge(vector<int> &contigs, int *ncontig, int *convergeon){
		int stepcutoff = Walkcutoff;
		float stepratiocutoff = Ratiocutoff;
		if(contigs.size() == 0 || contigs[0] == 0){
			*ncontig = 0;
			*convergeon = 0;
			return;
		}
		if(contigs.size() == 1){
			*ncontig = contigs[0];
			*convergeon = contigs[0];
			return;
		}
		if(contigs.size() > 2){
			vector<int> a;
			a.push_back(contigs[0]);
			a.push_back(contigs[1]);
			contigs.erase(contigs.begin(), contigs.begin()+2);
			int x;
			int y;
			converge(a, &x, &y);
			contigs.insert(contigs.begin(), x);
			converge(contigs, ncontig, convergeon);
			return;
		}
		else {
			vector<int> a;
			vector<int> b;
			a.push_back(contigs[0]);
			b.push_back(contigs[1]);
			int lengtha = Contiglens[abs(a[0])] - k + 1;
			int lengthb = Contiglens[abs(b[0])] - k + 1;
			int length = lengtha < lengthb ? lengtha : lengthb;
			int step_ = 0;
			
			while(length < ConvergeLengthwall && step_ <= Convergestepwall){
				step_ ++;
				if(a[a.size() - 1] != 0){
					vector<int> temp;
//if(temp.size() == 0)
//	continue;
					step(a[a.size() - 1], stepcutoff, stepratiocutoff, temp);
					int x = temp[0];
					if(x != 0){
						temp = vector<int>();
						step(-x, stepcutoff, stepratiocutoff, temp);
						if(temp[0] != - a[a.size() - 1])
							x = 0;
						
						else {
							int b_contains_x = 0;
							for(int i = 0; i < b.size(); i++){
								if(b[i] == x){
									b_contains_x = 1;
									break;
								}
							}
							if(b_contains_x == 1){
								lengthb = 0;
								for(int i = 0; i < b.size(); i++){
									int b_ = b[i];
									if(b_ == x)
										break;
									lengthb += Contiglens[abs(b_)] - k + 1;
								}
								int xx = abs(lengtha - lengthb);
								if(xx > 3 || lengthb == 0){
									*ncontig = 0;
									*convergeon = 0;
									return;	//make sure the other branch not 0
								}
								*ncontig = a[0];
								*convergeon = x;
								return;
							}
							lengtha += Contiglens[abs(x)] - k + 1;
							a.push_back(x);
						}
					}
				}
				if(b[b.size() - 1] != 0){
					vector<int> temp;
//if(temp.size() == 0)
//	continue;
					step(b[b.size() - 1], stepcutoff, stepratiocutoff, temp);
					int x = temp[0];
					if(x != 0){
						int a_contains_x = 0;
						for(int i = 0; i < a.size(); i++){
							if(a[i] == x){
								a_contains_x = 1;
								break;
							}
						}
						if(a_contains_x == 1){
							lengtha = 0;
							for(int i = 0; i < a.size(); i++){
								int a_ = a[i];
								if(a_ == x)
									break;
								lengtha += Contiglens[abs(a_)] - k + 1;
							}
							int xx = abs(lengtha - lengthb);
							if(xx > 3 || lengtha == 0){
								*ncontig = 0;
								*convergeon = 0;
								return;
							}
							*ncontig = a[0];
							*convergeon = x;
							return;
						}
						lengthb += Contiglens[abs(x)] - k + 1;
					}
					b.push_back(x);
				}
				length = lengtha < lengthb ? lengtha:lengthb;
			}
			*ncontig = 0;
			*convergeon = 0;
			return;			  
		}
		return;
	}
	
	void step(int contig, int cutoff, float ratiocutoff,vector<int> &aa){
		tools tl;
		map<int, int> n;
		nextcontigs_warc(contig, n);
		int nx = n.size();
		if(nx == 0){
			aa.push_back(0);
			return;
		}
		map<int , int> long_;
		int bigpro = 0;
		int bigprocontig = 0;
		int biglong = 0;
		int biglongcontig = 0;
		map<int, int> ncontigh;
		for (map<int, int>::iterator n_it = n.begin(); n_it != n.end(); n_it ++) {
			int n_key = (*n_it).first; // id in HH
			int n_value = (*n_it).second; // count in HH
			if(Contigtips.find(n_key) == Contigtips.end())
				Contigtips[n_key] = DefaultTip;
			if(n_value * Contigtips[n_key] > bigpro){
				bigpro = n_value * Contigtips[n_key];
				bigprocontig = n_key;
			}
			if(Contigtips[n_key] > Smalltip && n_value > Thin){ // long path must also be thick
				long_[n_key] = n_value;
				if(long_[n_key] > biglong){
					biglong = long_[n_key];
					biglongcontig = n_key;
				}
			}
		}
		if(long_.size() == 0)
			ncontigh[bigprocontig] = n[bigprocontig];
		else{
			for (map<int, int>::iterator long_it = long_.begin(); long_it != long_.end(); long_it ++) {
				int i = (*long_it).first;
				if(float(long_[i]) > ratiocutoff*(float)biglong || long_[i] > cutoff)
					ncontigh[i] = long_[i];
			}
		}
		tl.sort_value_reverse_key_int(ncontigh, aa, -1, -1);
		aa.push_back(nx);
	}
	
	void dump_contigs2(map<int, map<string, string> > &dContigs2){
		for(int ii = 1; ii <= Contigs2.size(); ii++){
			string node = Contigs2[ii].substr(0, k);
			int dir;
			string truenode = true_func(node, &dir);
			int contig = (HH[truenode].tag)*dir*(-1);
			map<int, int> a;
			nextcontigs_warc(contig, a);
			string i = "";
			string o = "";
			for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it++){
				int con = (*a_it).first;
				i = i + itos(Contigtags[abs(con)]*con/abs(con));
				i = i + ":" + itos(a[con]) + ","; // ID:count in HH
			}
			if(Contiglens2[ii]-k < 0 || Contiglens2[ii]-k > Contigs2[ii].length()) // protect substr
				continue;
			node = Contigs2[ii].substr(Contiglens2[ii]-k, k);
			truenode = true_func(node, &dir);
			contig = (HH[truenode].tag)*dir;
			a = map<int, int>();
			nextcontigs_warc(contig, a);
			for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it++){
				int con = (*a_it).first;
				o = o + itos(Contigtags[abs(con)]*con/abs(con));
				o = o + ":" + itos(a[con]) + ",";
			}
			map<string, string> contig2;
			contig2["id"] = itos(ii);
			contig2["lens"] = itos(Contiglens2[ii]);
			contig2["covs"] = ftos(Contigcovs2[ii]);
			contig2["types"] = Contigtypes2[ii];
			contig2["I"] = i;
			contig2["O"] = o;
			contig2["seq"] = Contigs2[ii];
			dContigs2[ii] = contig2;
		}
		return;
	}
			
			
	// find the next kmer of the contig, and count how many they are: map<id, count> n
	void nextcontigs_warc(int contignum, map<int, int> &n){
		string node;
		if(contignum > 0){
			if(Contiglens[contignum]-k < 0) // protect substr
				return;
			node = Contigs[contignum].substr(Contiglens[contignum] - k, k);
		}

		else {
			node = Contigs[-contignum].substr(0, k);
			node = revcom(node);
		}
		string I = "I";
		string O = "O";
		//char ACGT[4] = {'A','C','G','T'};
		if(HH.find(node) != HH.end()){
			if(node.length() <= 1) // protect substr
				return;
			if(HH[node].AO != 0){
				string vnode = (node.substr(1, k - 1)) + "A";
				int dir;
				string truenode = true_func(vnode, &dir);
				n[dir*(HH[truenode].tag)] = HH[node].AO;
			}
			
			if(HH[node].CO != 0){
				string vnode = (node.substr(1, k - 1)) + "C";
				int dir;
				string truenode = true_func(vnode, &dir);
				n[dir*(HH[truenode].tag)] = HH[node].CO;
			}
			
			if(HH[node].GO != 0){
				string vnode = (node.substr(1, k - 1)) + "G";
				int dir;
				string truenode = true_func(vnode, &dir);
				n[dir*(HH[truenode].tag)] = HH[node].GO;
			}

			if(HH[node].TO != 0){
				string vnode = (node.substr(1, k - 1)) + "T";
				int dir;
				string truenode = true_func(vnode, &dir);
				n[dir*(HH[truenode].tag)] = HH[node].TO;
			}
		}
		else {
			string rnode = revcom(node);
			if(HH[rnode].AI != 0){
				string vnode = "A" + (rnode.substr(0, k - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				n[-dir*(HH[truenode].tag)] = HH[rnode].AI;
			}
			if(HH[rnode].CI != 0){
				string vnode = "C" + (rnode.substr(0, k - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				n[-dir*(HH[truenode].tag)] = HH[rnode].CI;
			}

			if(HH[rnode].GI != 0){
				string vnode = "G" + (rnode.substr(0, k - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				n[-dir*(HH[truenode].tag)] = HH[rnode].GI;
			}

			if(HH[rnode].TI != 0){
				string vnode = "T" + (rnode.substr(0, k - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				n[-dir*(HH[truenode].tag)] = HH[rnode].TI;
			}
		}
	}
	
	string true_func(string node, int *dir){
		if(HH.find(node) != HH.end()){
			*dir = 1;
			return node;
		}
		else {
			*dir = -1;
			return revcom(node);
		}
	}
	
	string revcom(string seq){
		tools tl;
		reverse(seq.begin(), seq.end());
		//seq = tl.reverse(seq);
		seq = tl.tr(seq, "ATGC", "TACG");
		return seq;
	}
};

class maprdtocontig{
	
public:
	//map<string, map<string, int> > HH;
	map<string, test_HH> HH;
	int Rdnum;// = 0;
	int K;
	map<int, string > Contigs;
	map<int, int> Contiglens;
	map<int, float> Contigcovs;
	map<int, string> Contigtypes;
	int Contignum;// = 0;
	map<int, int> Contigtips;
	int DefaultTip;// = 1000;
	
	map<int, int> Contigtags;
	map<int, string > Contigs2;
	map<int, int> Contiglens2;
	map<int, float> Contigcovs2;
	map<int, string> Contigtypes2;
	int Contignum2;// = 0;
	int Smalltip;// = 100;
	int Thin;// = 2; // 4 for trep; tricky here
	float Ratiocutoff;// = 0.3;
	int Walkcutoff; // 10 fro trep
	int Convergestepwall;// = 12;
	int Convergelengthwall;// = 150;
	
	map<int, vector<int> > Interest2_rd;
	map<int, vector<int> > Interest2;
	map<int, vector<int> > Interest;
	int Insertwall;// = 35;
	int Interestlength;// = 100;
	float Scaffold_factor;// = 1.2;
	map<int, int> Contigtags2;
	vector<string> Reads;
	
	maprdtocontig(){
		Contignum = 0;
		DefaultTip = 1000;
		Contignum2 = 0;
		Smalltip = 100;
		Thin = 2;
		Ratiocutoff = 0.3;
		Walkcutoff = 10;
		Convergestepwall = 12;
		Convergelengthwall = 150;
		Insertwall = 35;
		Interestlength = 100;
		Scaffold_factor = 1.2;
		Walkcutoff = Thin * 3 + 1; // 10 for trep
	}
	
	~maprdtocontig(){
		//HH = map<string, map<string, int> >();
		HH = map<string, test_HH>();
		Contigtips = map<int, int>();
		Contigs = map<int, string >();
		Contiglens = map<int, int>();
		Contigcovs = map<int, float>();
		Contigtypes = map<int, string>();
		Contigtags = map<int, int>();
		Contigs2 = map<int, string >();
		Contiglens2 = map<int, int>();
		Contigcovs2 = map<int, float>();
		Contigtypes = map<int, string>();
		Interest2_rd = map<int, vector<int> >();
		Interest2 = map<int, vector<int> >();
	}
	
	void set_kmer_size(int k){
		K = k;
	}
	
	void doit(vector<string> &reads, map<string, test_HH > &rHH, map<int, int> &tips, map<int, map<string, string> > &contig, map<int, map<string, string> > &contig2){
		Reads = reads;
		HH = rHH;
		Contigtips = tips;
		map<int, map<string, string> > in_contigs = contig;
		map<int, map<string, string> > in_contigs2 = contig2;
		
		for(int i = 1; i <= in_contigs.size(); i++){
			Contiglens[i] = atoi(in_contigs[i]["lens"].c_str());
			Contigcovs[i] = atof(in_contigs[i]["covs"].c_str());
			Contigtypes[i] = in_contigs[i]["types"];
			Contigtags[i] = atoi(in_contigs[i]["tags"].c_str());
			Contigs[i] = in_contigs[i]["seq"];
		}	
		
		for(int i = 1; i <= in_contigs2.size(); i++){
			Contiglens2[i] = atoi(in_contigs2[i]["lens"].c_str());
			Contigcovs2[i] = atof(in_contigs2[i]["covs"].c_str());
			Contigtypes2[i] = in_contigs2[i]["types"];
			Contigs2[i] = in_contigs2[i]["seq"];
		}
		
		for(int i = 1; i <= Contigs2.size(); i++)
			Contigtags2[i] = i;
		
		interestgen();
		mapreads();
		scaffold1();
	}
	
	void scaffold1(){
		tools tl;
		//cerr << "Scaffolding .. \n";
		map<int, int> contig2check;
		for(map<int, vector<int> >::iterator Interest2_rd_it = Interest2_rd.begin(); Interest2_rd_it != Interest2_rd.end(); Interest2_rd_it++)
			contig2check[(*Interest2_rd_it).first] = 0;
		for(map<int, vector<int> >::iterator Interest2_rd_it = Interest2_rd.begin(); Interest2_rd_it != Interest2_rd.end(); Interest2_rd_it++){
			int key_Interest2_rd = (*Interest2_rd_it).first;
			if(contig2check[key_Interest2_rd] == 1)
			   continue;
			int xcontig2;
			string path;
			checkIcontig(key_Interest2_rd, &xcontig2, &path);
			contig2check[key_Interest2_rd] = 1;
			if(xcontig2 == 0)
				continue;
			int xxcontig2;
			string xpath;
			checkIcontig(-xcontig2, &xxcontig2, &xpath);
			if(key_Interest2_rd == -xxcontig2){
				contig2check[-xcontig2] = 1; // Xian: check to see if x should be xx
				string sep = " ";
				vector<string> p;
				tl.split(path, sep, p); // Xian: need to write this function
				vector<float> sort;
				tl.sort_value(Contigcovs, sort, abs(atoi(p[0].c_str())), abs(atoi(p[p.size() - 1].c_str()))); // Xian: done. need to write this function, may be float
				//cerr << "bingo " << key_Interest2_rd << "\t" << xcontig2 << "\t" << path << "\tsort " << sort[0] << " " << sort[1] << " "<< sort[1] - sort[0] << " "<< Walkcutoff*Scaffold_factor<< " " << Ratiocutoff*Scaffold_factor*sort[1] << "\n";
				if((sort[1] - sort[0]) < float(Walkcutoff)*Scaffold_factor && (sort[1] - sort[0]) < Ratiocutoff*Scaffold_factor*float(sort[1])){
					//cerr << "Bingo " << key_Interest2_rd << "\t" << xcontig2 << "\t" << path << "\n";
					
					
					if(atoi(p[0].c_str()) == Interest2[key_Interest2_rd][0] && atoi(p[p.size() - 1].c_str()) == - Interest2[-xcontig2][0]){ // merge go here
						//cerr << "BINGO "<< key_Interest2_rd << "\t"<< xcontig2 << "\t"<< path<< "\n";
						string pathseq = atoi(p[0].c_str()) > 0 ? Contigs[abs(atoi(p[0].c_str()))] : revcom(Contigs[abs(atoi(p[0].c_str()))]);
						if(pathseq.length()-K + 1 < 0) // protect subseq
							continue;
						pathseq = pathseq.substr(pathseq.length()-K + 1, K - 1);
						for(int i = 1; i <= p.size() - 2; i++){
							string tmp = atoi(p[i].c_str()) > 0 ? Contigs[abs(atoi(p[i].c_str()))]:revcom(Contigs[abs(atoi(p[i].c_str()))]);
							pathseq = tmp.substr(K - 1, tmp.length() - K + 1);
						}
						merge(key_Interest2_rd, xcontig2, pathseq);
					}
				}
			}
		}
	}
	
	void merge(int contigx, int contigy, string path){
		tools tl;
		int contigxreal = realcontig2(contigx);
		int contigyreal = realcontig2(contigy);
		if(abs(contigxreal) == abs(contigyreal))
			return;
		if(contigxreal > 0){
			string seqy = Contigs2[abs(contigyreal)];
			if(abs(contigyreal)< 0)
				seqy = revcom(seqy);
			Contigs[contigxreal] = Contigs2[contigxreal].substr(0, Contiglens2[contigxreal] - K + 1) + path + seqy.substr(K - 1, Contiglens2[abs(contigyreal)] - K + 1);
			Contigcovs2[contigxreal] = (Contigcovs2[contigxreal]*(float)Contiglens2[contigxreal] + (float)Contigcovs2[abs(contigyreal)]*Contiglens2[abs(contigyreal)])/(float)(Contiglens2[contigxreal] + Contiglens2[abs(contigyreal)]);
			Contiglens2[contigxreal] = (Contigs2[contigxreal]).length();
			string x = Contigtypes2[abs(contigyreal)];
			if(contigyreal < 0)
				reverse(x.begin(), x.end());
				//x = tl.reverse(x);
			string sep = "";
			vector<string> tmp;
			tl.split(Contigtypes2[contigxreal]+x, sep, tmp);
			Contigtypes2[contigxreal] = tmp[0] + tmp[tmp.size() - 1];
			Contigs2[abs(contigyreal)] = "";
		}
		else if(contigxreal < 0){
			string seqy = Contigs2[abs(contigyreal)];
			if(contigyreal > 0)
				seqy = revcom(seqy);
			Contigs2[-contigxreal] = (seqy.substr(0, Contiglens2[abs(contigyreal)] - K + 1)) + revcom(path)+ (Contigs2[-contigxreal].substr(K - 1, Contiglens2[-contigxreal] - K + 1));
			Contigcovs2[-contigxreal] = (Contigcovs2[-contigxreal] * float(Contiglens2[-contigxreal]) + Contigcovs2[abs(contigyreal)] * float(Contiglens2[abs(contigyreal)]))/(float(Contiglens[-contigxreal]) + float(Contiglens2[abs(contigyreal)]));
			Contiglens2[-contigxreal] = (Contigs2[-contigxreal]).length();
			string x = Contigtypes2[abs(contigyreal)];
			if(contigyreal > 0)
				reverse(x.begin(), x.end());
				//x = tl.reverse(x);
			string sep = "";
			vector<string> tmp;
			tl.split(x + Contigtypes2[-contigxreal], sep, tmp);
			Contigtypes2[-contigxreal] = tmp[0] + tmp[tmp.size() - 1];
			Contigs2[abs(contigyreal)] = "";
		}
		Contigtags2[abs(contigyreal)] = contigxreal *contigyreal / abs(contigyreal);
	}
	
	int realcontig2(int contig2){
		if(contig2 == 0)
			return 0;
		if(Contigtags2[abs(contig2)] == abs(contig2))
			return contig2;
		else 
			return realcontig2(Contigtags2[abs(contig2)]*contig2/abs(contig2));
	}

	void checkIcontig(int contig2, int *xcontig2, string *path){
		map<int, int> h;
		map<int, map<string, int> > ha;
		int n = 0;
		for(int ii = 0; ii != Interest2_rd[contig2].size(); ii++){
			int i = Interest2_rd[contig2][ii];
			string rd = Reads[abs(i) - 1]; // be careful here, -1
			if(i < 0)
				rd = revcom(rd);
			int xcontig2;
			string path;
			walkrd(rd, contig2, &xcontig2, &path);
			if(xcontig2 == 0)
				continue;
			h[xcontig2]++;
			if(ha.find(xcontig2) == ha.end() || ha[xcontig2].find(path) == ha[xcontig2].end())
				ha[xcontig2][path] = 0; // changed here for cpp
			else 
				ha[xcontig2][path] ++;
			n++;
		}
		tools tl;
		vector<int> sorted;
		tl.sort_value_reverse_key_int(h, sorted, -1, -1);//Xian: need to write this function, reverse the value, get the key
		if(h[sorted[1]] < Walkcutoff && float(h[sorted[0]])*Ratiocutoff > h[sorted[1]] && h[sorted[0]] >= Thin){
			vector<string> sortedpath;
			tl.sort_value_reverse_key_str(ha[sorted[0]], sortedpath);// Xian: same thing
			if(ha[sorted[0]][sortedpath[0]] > 0.5*h[sorted[0]]){
				*xcontig2 = sorted[0];
				*path = sortedpath[0];
				return;
			}
		}
		*xcontig2 = 0;
		*path = "0";
	}
	
	void walkrd(string rd, int contig2, int *xcontig2, string *path){
		vector<int> s;
		vector<int> s2;
		vector<int> p;
		float tag = 0;
		string end = "";
		int endtag = 0;
		for(int i = 0; i <= rd.length() - K; i++){
			string w = rd.substr(i, K);
			if(HH.find(w)!=HH.end() || HH.find(revcom(w)) != HH.end()){
				int dir;
				string truew = true_func(w, &dir);
				if(dir*(HH[truew].tag) != s[s.size() - 1] || endtag == 1){
					endtag = 0;
					int x = Contigtags[abs(HH[truew].tag)]*HH[truew].tag/abs(HH[truew].tag)*dir;
					s.push_back(dir*HH[truew].tag);
					if(Interest2_rd.find(x) != Interest2_rd.end() && Interest2_rd[x].size() != 0){
						if(x!= s2[s2.size() - 1]){
							s2.push_back(x);
							if(tag == 1.5)
								tag = 2;
						}
						if(x==contig2)
							tag = 1;
						else if(tag == 1){
							p = vector<int>();
							p.push_back(s[s.size() - 2]);
							p.push_back(s[s.size() - 1]);
							tag = 1.5;
						}
						else if(tag == 1.5)
							p.push_back(s[s.size() - 1]);
						else if(tag == 2){
							p.push_back(s[s.size() - 1]);
							tag = 3;
							break;
						}					
					}
					else{
						if(tag == 1){
							p = vector<int>();
							p.push_back(s[s.size() - 2]);
							p.push_back(s[s.size() - 1]);
							tag = 2;
						}
						else if(tag == 2 || tag == 1.5){
							p.push_back(s[s.size() - 1]);
							if(tag == 1.5)
								tag = 2;
						}
						end = Contigs[abs(s[s.size() - 1])];
						if(s[s.size() - 1] < 0)
							end = revcom(end);
					}
				}
				if(w.compare(end.substr(end.length() - K, K)) == 0)
					endtag = 1;
			}
			else{
				tag = 0;
			}
		}
		vector<int> tmp = Interest2[contig2];
		/*cerr << "debugx contig2 " << contig2 << " tmp ";
		for(int i = 0; i < tmp.size(); i++)
			cerr << tmp[i] << " ";
		cerr << " s ";
		for(int i = 0; i < s.size(); i++)
			cerr << s[i] << " ";
		cerr << " s2 ";
		for(int i = 0; i < s2.size(); i++)
			cerr << s2[i] << " ";
		cerr << " p ";
		for(int i = 0; i < p.size(); i++)
			cerr << p[i] << " ";
		cerr << "\n";*/
		if(tag == 3){
			*xcontig2 = s2[s2.size() - 1];
			*path = "@p";
			return;
		}
		*xcontig2 = 0;
		*path = "0";
	}
	
	void test(){
		
	}
		
	/*void mapreads(){
		tools tl;
		//cerr << "Mapping reads to contigs .. \n";
		int index = 0; // reads start from 1, negative mean revcom
		for(int j = 0; j < Reads.size(); j++){
			string rd = Reads[j];
			rd = tl.chomp(rd);
			index ++;
			//if(index % 100000 == 0)
			//	cerr << "Mapped " << index << " reads\n";
			map<int, int> a;
			for(int i = 0; i <= rd.length() - K; i+=2){
				string w = rd.substr(i, K);
				if(!(HH.find(w)!= HH.end() && HH[w].size() != 0 || HH.find(revcom(w))!=HH.end() && HH[revcom(w)].size() != 0))
					continue;
				int dir;
				string truew = true_func(w, &dir);
				if(HH.find(truew) != HH.end() && HH[truew].find("n") != HH[truew].end() && HH[truew]["n"] <= 1)
					continue;
				int contig = dir*HH[truew]["tag"];
				if(! (Interest.find(abs(contig))!=Interest.end()))
					continue;
				if(a.find(contig) == a.end() || a[contig] == 0){
					Interest[abs(contig)].push_back(contig/abs(contig)*index);
					a[contig] = 1;
				}
			}
		}
	}*/

	void mapreads(){
		tools tl;
		//cerr << "Mapping reads to contigs .. \n";
		int index = 0; // reads start from 1, negative mean revcom
		for(int j = 0; j < Reads.size(); j++){
			string rd = Reads[j];
			if(rd.length() < K)
				continue;
			rd = tl.chomp(rd);
			index ++;
			//if(index % 100000 == 0)
			//	cerr << "Mapped " << index << " reads\n";
			map<int, int> a;
			for(int i = 0; i <= rd.length() - K; i+=2){
				string w = rd.substr(i, K);
				if(!(HH.find(w)!= HH.end() || HH.find(revcom(w))!=HH.end()))
					continue;
				int dir;
				string truew = true_func(w, &dir);
				if(HH.find(truew) != HH.end() && HH[truew].n <= 1)
					continue;
				int contig = dir*HH[truew].tag;
				if(! (Interest.find(abs(contig))!=Interest.end()))
					continue;
				if(a.find(contig) == a.end() || a[contig] == 0){
					Interest[abs(contig)].push_back(contig/abs(contig)*index);
					a[contig] = 1;
				}
			}
		}
	}
	
	void printmapping(string myfile){
		ofstream MYF;
		char myfile_[myfile.length()];
		strcpy(myfile_, myfile.c_str());
		MYF.open(myfile_);
		//cerr << "Print reads mapping to file " << myfile << " .. \n";
		for(map<int, vector<int> >::iterator i_it = Interest2.begin(); i_it != Interest2.end(); i_it++){
			int k = (*i_it).first;
			map<int, int> a;
			vector<int> temp = Interest2[k];
			MYF << "BC " << k;
			for(int i = 0; i < temp.size(); i++)
				MYF << " " << temp[i]; // Big Contig
			MYF << "\n";
			Interest2_rd[k] = vector<int>();
			for(int ii = 0; ii < temp.size(); ii++){
				int i = temp[ii];
				int dir = i/abs(i);
				vector<int> tmp = Interest[abs(i)];
				for(int jj = 0; jj < tmp.size(); jj++){
					int j = tmp[jj];
					int rd = j*dir;
					if(! (a.find(rd) != a.end() && a[rd] != 0)){
						Interest2_rd[k].push_back(rd);
						a[rd] = 1;
					}
				}
				MYF << "SC " << abs(i);
				for(int jj = 0; jj < tmp.size(); jj++)
					MYF << " " << tmp[jj]; // Small (proto) Contig
				MYF << "\n";
			}
			
			vector<int> t = Interest2_rd[k];
			MYF << "AC " << k;
			for(int ii = 0; ii < t.size(); ii++)
				MYF << " " << t[ii];
			MYF << "\n";			
		}
		MYF.close();
	}
	
	void interestgen(){// no int cutoff as the input used
		//cerr << "Generating contigs list ..\n";
		for(int i = 1; i <= Contigs2.size(); i++){
			if(Contiglens2[i] < Interestlength)
				continue;
			Interest2[i] = vector<int>();
			Interest2[-i] = vector<int>();
			int wall = Insertwall;
			if(Contiglens2[i]/2 < wall)
				wall = Contiglens2[i]/2;
			int length = K - 1;
			while(1){
				string node = Contigs2[i].substr(length - K + 1, K);
				int dir;
				string truenode = true_func(node, &dir);
				int contig = HH[truenode].tag;
				Interest2[-i].push_back(contig*-dir);
				if(!(Interest.find(abs(contig)) != Interest.end() && Interest[abs(contig)].size() != 0))
				   Interest[abs(contig)] = vector<int>();
				length += Contiglens[abs(contig)] - K + 1;
				if(length >= wall)
					break;
			}
			
			length = K - 1;
			while(1){
				string node = Contigs2[i].substr(Contiglens2[i] - length - 1, K);
				int dir;
				string truenode = true_func(node, &dir);
				int contig = HH[truenode].tag;
				Interest2[i].push_back(contig*dir);
				if(!(Interest.find(abs(contig)) != Interest.end() && Interest[abs(contig)].size() != 0))
					Interest[abs(contig)] = vector<int>();
				length += Contiglens[abs(contig)] - K + 1;
				if(length >= wall)
					break;
			}
		}
		int temp = Interest.size();
		//cerr << temp << " interesting contigs found.\n";
	}
	
	string revcom(string seq){
		tools tl;
		reverse(seq.begin(), seq.end());
		//seq = tl.reverse(seq);
		seq = tl.tr(seq, "ATGC", "TACG");
		return seq;
	}	
	
	string true_func(string node, int *dir){
		if(HH.find(node) != HH.end()){
			*dir = 1;
			return node;
		}
		else {
			*dir = -1;
			return revcom(node);
		}
	} // slightly different, need to check if HH[rnode] exist. If neither exist?
	
	void nextcontigs(int contignum, vector<int> &a){ // this function has never been used in this class
		string node;
		if(contignum > 0)
			node = Contigs[contignum].substr(Contiglens[contignum] - K, K);
		else{
			node = Contigs[-contignum].substr(0, K);
			node = revcom(node);
		}
		string I = "I";
		string O = "O";
		//char ACGT[4] = {'A','C','G','T'};
		if(HH.find(node) != HH.end()){
			if(HH[node].AO != 0){
				string vnode = (node.substr(1, K - 1)) + "A";
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			
			if(HH[node].CO != 0){
				string vnode = (node.substr(1, K - 1)) + "C";
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			
			if(HH[node].GO != 0){
				string vnode = (node.substr(1, K - 1)) + "G";
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}

			if(HH[node].TO != 0){
				string vnode = (node.substr(1, K - 1)) + "T";
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(dir*HH[truenode].tag);
			}
			
		}
		else {// be careful if rnode is not the key
			string rnode = revcom(node);
			if(HH[rnode].AI != 0){
				string vnode = "A" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}
			if(HH[rnode].CI != 0){
				string vnode = "C" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}

			if(HH[rnode].GI != 0){
				string vnode = "G" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}

			if(HH[rnode].TI != 0){
				string vnode = "T" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a.push_back(-dir*HH[truenode].tag);
			}

		}
		return;
	}
		
	void nextbase(string node, int *A, int *T, int *G, int *C){// this function has never been used in this class. copied from addbridgekmer
		int dir;
		string truenode = true_func(node, &dir);
		if(dir == 1){
			*A = HH[truenode].AO;
			*T = HH[truenode].TO;
			*G = HH[truenode].GO;
			*C = HH[truenode].CO;
		}
		else {
			*A = HH[truenode].TI;
			*T = HH[truenode].AI;
			*G = HH[truenode].CI;
			*C = HH[truenode].GI;
		}
	}
		
	void dump_contigs2_m(map<int, map<string, string> > &dContigs2){
		for(int ii = 1; ii <= Contigs2.size(); ii++){
			map<string, string> contig2;
			if(Contigtags2[ii] == ii){
				string node = Contigs2[ii].substr(0, K);
				int dir;
				string truenode = true_func(node, &dir);
				int contig = (HH[truenode].tag)*dir*(-1);
				map<int, int> a;
				nextcontigs_warc(contig, a);
				string i = "";
				string o = "";
				for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it ++){
					int con = (*a_it).first;
					int tmp = realcontig2(Contigtags[abs(con)]*con/abs(con));
					string conhead = con > 0 ? Contigs[con].substr(0, K):revcom(Contigs[-con].substr(Contiglens[-con]-K, K));
					string tmphead = tmp > 0 ? Contigs2[tmp].substr(0, K):revcom(Contigs2[-tmp].substr(Contiglens2[-tmp]-K, K));
					string tmp_str;
					if(tmphead.compare(conhead) != 0)
						tmp_str = itos(tmp) + "M";
					else
						tmp_str = itos(tmp);
					i = i + tmp_str;
					i = i + ":" + itos(a[con]) + ",";
				}
				node = Contigs2[ii].substr(Contiglens2[ii]-K, K);
				truenode = true_func(node, &dir);
				contig = (HH[truenode].tag)*dir;
				a = map<int, int>();
				nextcontigs_warc(contig, a);
				for(map<int, int>::iterator a_it = a.begin(); a_it != a.end(); a_it ++){
					int con = (*a_it).first;
					int tmp = realcontig2(Contigtags[abs(con)]*con/abs(con));
					string conhead = con > 0 ? Contigs[con].substr(0, K):revcom(Contigs[-con].substr(Contiglens[-con]-K, K));
					string tmphead = tmp > 0 ? Contigs2[tmp].substr(0, K):revcom(Contigs2[-tmp].substr(Contiglens2[-tmp]-K, K));
					string tmp_str;
					if(tmphead.compare(conhead) != 0)						
						tmp_str = itos(tmp) + "M";
					else
						tmp_str = itos(tmp);
					o = o + tmp_str;
					o = o + ":" + itos(a[con]) + ",";
				}				
				contig2["id"] = itos(ii);
				contig2["lens"] = itos(Contiglens2[ii]);
				contig2["covs"] = ftos(Contigcovs2[ii]);
				contig2["types"] = Contigtypes2[ii];
				contig2["I"] = i;
				contig2["O"] = o;
				contig2["tags"] = itos(Contigtags2[ii]);
			}
			else{
				contig2["id"] = itos(ii);
				contig2["lens"] = itos(Contiglens2[ii]);
				contig2["covs"] = ftos(Contigcovs2[ii]);
				contig2["types"] = Contigtypes2[ii];
				contig2["I"] = "";
				contig2["O"] = "";
				contig2["tags"] = itos(Contigtags2[ii]) + " *";
			}
			contig2["seq"] = Contigs2[ii];
			dContigs2[ii] = contig2;
		}
		return;
	}
				
	void nextcontigs_warc(int contignum, map<int, int> &a){ // this function has never been used in this class
		string node;
		if(contignum > 0)
			node = Contigs[contignum].substr(Contiglens[contignum] - K, K);
		else{
			node = Contigs[-contignum].substr(0, K);
			node = revcom(node);
		}
		string O = "O";
		string I = "I";
		//char ACGT[4] = {'A','C','G','T'};
		if(HH.find(node) != HH.end()){
			if(HH[node].AO != 0){
				string vnode = (node.substr(1, K - 1)) + "A";
				int dir;
				string truenode = true_func(vnode, &dir);
				a[dir*(HH[truenode].tag)] = HH[node].AO;
			}
			
			if(HH[node].CO != 0){
				string vnode = (node.substr(1, K - 1)) + "C";
				int dir;
				string truenode = true_func(vnode, &dir);
				a[dir*(HH[truenode].tag)] = HH[node].CO;
			}
			
			if(HH[node].GO != 0){
				string vnode = (node.substr(1, K - 1)) + "G";
				int dir;
				string truenode = true_func(vnode, &dir);
				a[dir*(HH[truenode].tag)] = HH[node].GO;
			}

			if(HH[node].TO != 0){
				string vnode = (node.substr(1, K - 1)) + "T";
				int dir;
				string truenode = true_func(vnode, &dir);
				a[dir*(HH[truenode].tag)] = HH[node].TO;
			}


			/*for(int i = 0; i < 4; i++){
				if(HH[node].find(ACGT[i]+O) == HH[node].end() || HH[node].find(ACGT[i]+O) != HH[node].end() && HH[node][ACGT[i]+O] == 0)
					continue;
				string vnode = (node.substr(1, K - 1)) + ACGT[i];
				int dir;
				string truenode = true_func(vnode, &dir);
				a[dir*HH[truenode]["tag"]] = HH[node][ACGT[i]+O];
			}*/
		}
		else {
			string rnode = revcom(node);
			if(HH[rnode].AI != 0){
				string vnode = "A" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a[-dir*(HH[truenode].tag)] = HH[rnode].AI;
			}
			if(HH[rnode].CI != 0){
				string vnode = "C" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a[-dir*(HH[truenode].tag)] = HH[rnode].CI;
			}

			if(HH[rnode].GI != 0){
				string vnode = "G" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a[-dir*(HH[truenode].tag)] = HH[rnode].GI;
			}

			if(HH[rnode].TI != 0){
				string vnode = "T" + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a[-dir*(HH[truenode].tag)] = HH[rnode].TI;
			}
			/*for(int i = 0; i < 4; i++){
				if(HH[rnode].find(ACGT[i]+I) == HH[rnode].end() || HH[rnode].find(ACGT[i]+I) != HH[rnode].end() && HH[rnode][ACGT[i]+I] == 0)
					continue;
				string vnode = ACGT[i] + (rnode.substr(0, K - 1));
				int dir;
				string truenode = true_func(vnode, &dir);
				a[-dir*HH[truenode]["tag"]] = HH[rnode][ACGT[i]+I];
			}*/
		}
		return;
	}	
};
		
// generate a graph of the kmer created in kmergen, HH, so that we know of the numbers of each input allele (on the left) and the output allele (on the right) of each kmer in hh. Will not connect the two at different status: one on a branch, the other as a unique one.
class graphgen{
public:
	
	int k;// = 25;
	map<string, test_HH > HH; // map<contig, map<info, count> >
	
	graphgen(){
		k = 25;
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}

	void set_HH(map<string, int> &kmers){
		for(map<string, int>::iterator kmers_it = kmers.begin(); kmers_it != kmers.end(); kmers_it ++){
			string key = (*kmers_it).first;
			int value = (*kmers_it).second;
			test_HH t;
			t.n = value; t.AI = 0; t.CI = 0; t.GI = 0; t.TI = 0; t.AO = 0; t.CO = 0; t.GO = 0; t.TO = 0; t.tag = 0; t.tag2 = 0;
			HH[key] = t;
		}
	}

	map<string, test_HH > get_HH(){
		/*map<string, map<string, int> > rHH;
		for(map<string, test_HH>::iterator it = HH.begin(); it != HH.end(); it++){
			map<string, int> rHH_second;
			rHH_second["n"] = (*it).second.n;
			rHH_second["AI"] = (*it).second.AI;
			rHH_second["CI"] = (*it).second.CI;
			rHH_second["TI"] = (*it).second.TI;
			rHH_second["GI"] = (*it).second.GI;
			rHH_second["AO"] = (*it).second.AO;
			rHH_second["CO"] = (*it).second.CO;
			rHH_second["GO"] = (*it).second.GO;
			rHH_second["TO"] = (*it).second.TO;
			rHH_second["tag"] = (*it).second.tag;
			rHH_second["tag2"] = (*it).second.tag2;
			rHH[(*it).first] = rHH_second;
		}
		return rHH;*/
		return HH;
	}
	
	~graphgen(){
		HH = map<string, test_HH >();
	}
	
	// generate HH, which take the infos in hh, and also the next kmer (O) or the previous kmer (I) number for each kmer 
	int doit(vector<string> &rReads, vector<string> &PR){
		tools tl;
		int RDnum = 0;
		for(int ii = 0; ii < rReads.size(); ii++){

			string rd = rReads[ii];
			string PRstr("");
			int l = rd.length();
if(l < k)
	continue;

// efficient code:
string rrd = rd;
reverse(rrd.begin(), rrd.end());
rrd = tl.tr(rrd, "ATGC", "TACG");
//
			for(int i = k; i <= l-1; i++){ // ??? Need to ask Xian: need to leave space for the next kmer
				string w = rd.substr(i-k, k);
				string rw = rrd.substr(l-i, k);
				//string rw = w;			
				//reverse(rw.begin(), rw.end());
				//string rw = tl.reverse(w);
				//rw = tl.tr(rw, "ATGC", "TACG");
				int occur = 0;
				int tag = 0; // skip bases when kmer do not exit
				if(HH.find(w) != HH.end()){
					occur = HH[w].n;
					string u = rd.substr(i-k+1, k);
					string ru = rrd.substr(l-i-1, k);
					//string ru = u;
					//reverse(ru.begin(), ru.end());
					//string ru = tl.reverse(u);
					//string old_str = "ATGC";
					//string new_str = "TACG";
					//ru = tl.tr(ru, old_str, new_str);
					if(HH.find(u) != HH.end()){ // $w $u  // Xian: notice; the next contig should also exit in HH ###

						char next = u[k-1];
						if(next == 'A')
							HH[w].AO ++;
						else if(next == 'C')
							HH[w].CO ++;
						else if(next == 'G')
							HH[w].GO ++;
						else if(next == 'T')
							HH[w].TO ++;
						//HH[w][next+"O"] ++;  // go out count up
						next = w[0];
						if(next == 'A')
							HH[u].AI ++;
						else if(next == 'C')
							HH[u].CI ++;
						else if(next == 'G')
							HH[u].GI ++;
						else if(next == 'T')
							HH[u].TI ++;
						
						//HH[u][next+"I"] ++; // go in for u count up
					}
					else if(HH.find(ru) != HH.end()){ // $w $ru
						char next = u[k-1];
						if(next == 'A')
							HH[w].AO ++;
						else if(next == 'C')
							HH[w].CO ++;
						else if(next == 'G')
							HH[w].GO ++;
						else if(next == 'T')
							HH[w].TO ++;
						//HH[w][next+"O"] ++; // Xian: cannot understand this
						next = rw[k-1];
						if(next == 'A')
							HH[ru].AO ++;
						else if(next == 'C')
							HH[ru].CO ++;
						else if(next == 'G')
							HH[ru].GO ++;
						else if(next == 'T')
							HH[ru].TO ++;
						//HH[ru][next+"O"] ++;
					}
					else{
						i++;
						tag = 1;
					}
				}
				else if(HH.find(rw) != HH.end()){
					occur = HH[rw].n;
					string u = rd.substr(i-k+1, k);
					string ru = rrd.substr(l-i-1, k);
					//string ru = u;
					//reverse(ru.begin(), ru.end());
					//string ru = tl.reverse(u);
					//string old_str = "ATGC";
					//string new_str = "TACG";
					//ru = tl.tr(ru, old_str, new_str);
					if(HH.find(u) != HH.end()){ // rw, u
						char next = ru[0];
						if(next == 'A')
							HH[rw].AI ++;
						else if(next == 'C')
							HH[rw].CI ++;
						else if(next == 'G')
							HH[rw].GI ++;
						else if(next == 'T')
							HH[rw].TI ++;
						//HH[rw][next+"I"] ++;
						next = w[0];
						if(next == 'A')
							HH[u].AI ++;
						else if(next == 'C')
							HH[u].CI ++;
						else if(next == 'G')
							HH[u].GI ++;
						else if(next == 'T')
							HH[u].TI ++;
						//HH[u][next+"I"] ++;
					}
					else if(HH.find(ru) != HH.end()){ // rw ru
						char next = ru[0];
						if(next == 'A')
							HH[rw].AI ++;
						else if(next == 'C')
							HH[rw].CI ++;
						else if(next == 'G')
							HH[rw].GI ++;
						else if(next == 'T')
							HH[rw].TI ++;
						//HH[rw][next+"I"] ++;
						next = rw[k-1];
						if(next == 'A')
							HH[ru].AO ++;
						else if(next == 'C')
							HH[ru].CO ++;
						else if(next == 'G')
							HH[ru].GO ++;
						else if(next == 'T')
							HH[ru].TO ++;
						//HH[ru][next+"O"] ++;
					}
					else{
						i++;
						tag = 1;
					}
				}
				else
					occur = 1; // no such kmer
				
				PRstr = PRstr + itos(occur) + " "; // Xian: numbers of the kmers
				if(tag == 1 && i < l)
					PRstr = PRstr + "1 ";
			}
			
			string w = rd.substr(l-k, k); // start starting from 0 or 1?
			string rw = w;
			reverse(rw.begin(), rw.end());
			//string rw = tl.reverse(w);
			rw = tl.tr(rw, "ATGC", "TACG");

			if(HH.find(w) != HH.end())
				PRstr = PRstr + itos(HH[w].n);
			else if(HH.find(rw) != HH.end())
				PRstr = PRstr + itos(HH[rw].n);
			else 
				PRstr = PRstr + "1";
			PR.push_back(PRstr);
			
			RDnum ++;
			if(RDnum % 100000 == 0)
				cerr << "Added " << RDnum << " reads\n";
		}
		return RDnum;
	}

	void printPR(string file, vector<string> &rPR){
		ofstream PRF;
		char file_char[file.length()];
		strcpy(file_char, file.c_str());
		PRF.open(file_char);
		PRF.close();
		PRF.open(file_char, ofstream::app);
		for(int i = 0; i < rPR.size(); i++)
			PRF << rPR[i] << "\n";
		PRF.close();
	}
};
		
/*class graphgen{
public:
	
	int k;// = 25;
	map<string, map<string, int> > HH; // map<contig, map<info, count> >
	
	graphgen(){
		k = 25;
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	void set_HH(map<string, int> &kmers){
		for(map<string, int>::iterator kmers_it = kmers.begin(); kmers_it != kmers.end(); kmers_it ++){
			string key = (*kmers_it).first;
			int value = (*kmers_it).second;
			HH[key]["n"] = value;
			HH[key]["AI"] = 0;
			HH[key]["CI"] = 0;
			HH[key]["GI"] = 0;
			HH[key]["TI"] = 0;		
			HH[key]["AO"] = 0;		
			HH[key]["CO"] = 0;			
			HH[key]["GO"] = 0;					
			HH[key]["TO"] = 0;			
			HH[key]["tag"] = 0;			
			HH[key]["tag2"] = 0;							 
		}
	}

	map<string, map<string, int> > &get_HH(){
		return HH;
	}
	
	~graphgen(){
		HH = map<string, map<string, int> >();
	}
	
	// generate HH, which take the infos in hh, and also the next kmer (O) or the previous kmer (I) number for each kmer 
	int doit(vector<string> &rReads, vector<string> &PR){
		tools tl;
		int RDnum = 0;
		for(int ii = 0; ii < rReads.size(); ii++){

			string rd = rReads[ii];
			string PRstr("");
			int l = rd.length();

// efficient code:
string rrd = rd;
reverse(rrd.begin(), rrd.end());
rrd = tl.tr(rrd, "ATGC", "TACG");
//
			for(int i = k; i <= l-1; i++){ // ??? Need to ask Xian: need to leave space for the next kmer
				string w = rd.substr(i-k, k);
				string rw = rrd.substr(l-i, k);
				//string rw = w;			
				//reverse(rw.begin(), rw.end());
				//string rw = tl.reverse(w);
				//rw = tl.tr(rw, "ATGC", "TACG");
				int occur = 0;
				int tag = 0; // skip bases when kmer do not exit
				if(HH.find(w) != HH.end()){
					occur = HH[w]["n"];
					string u = rd.substr(i-k+1, k);
					string ru = rrd.substr(l-i-1, k);
					//string ru = u;
					//reverse(ru.begin(), ru.end());
					//string ru = tl.reverse(u);
					//string old_str = "ATGC";
					//string new_str = "TACG";
					//ru = tl.tr(ru, old_str, new_str);
					if(HH.find(u) != HH.end()){ // $w $u  // Xian: notice; the next contig should also exit in HH ###

						string next = u.substr(k-1, 1);
						HH[w][next+"O"] ++;  // go out count up
						next = w.substr(0,1);
						HH[u][next+"I"] ++; // go in for u count up
					}
					else if(HH.find(ru) != HH.end()){ // $w $ru
						string next = u.substr(k-1, 1);
						HH[w][next+"O"] ++; // Xian: cannot understand this
						next = rw.substr(k-1,1);
						HH[ru][next+"O"] ++;
					}
					else{
						i++;
						tag = 1;
					}
				}
				else if(HH.find(rw) != HH.end()){
					occur = HH[rw]["n"];
					string u = rd.substr(i-k+1, k);
					string ru = rrd.substr(l-i-1, k);
					//string ru = u;
					//reverse(ru.begin(), ru.end());
					//string ru = tl.reverse(u);
					//string old_str = "ATGC";
					//string new_str = "TACG";
					//ru = tl.tr(ru, old_str, new_str);
					if(HH.find(u) != HH.end()){ // rw, u
						string next = ru.substr(0,1);
						HH[rw][next+"I"] ++;
						next = w.substr(0,1);
						HH[u][next+"I"] ++;
					}
					else if(HH.find(ru) != HH.end()){ // rw ru
						string next = ru.substr(0,1);
						HH[rw][next+"I"] ++;
						next = rw.substr(k-1,1);
						HH[ru][next+"O"] ++;
					}
					else{
						i++;
						tag = 1;
					}
				}
				else
					occur = 1; // no such kmer
				
				PRstr = PRstr + itos(occur) + " "; // Xian: numbers of the kmers
				if(tag == 1 && i < l)
					PRstr = PRstr + "1 ";
			}
			
			string w = rd.substr(l-k, k); // start starting from 0 or 1?
			string rw = w;
			reverse(rw.begin(), rw.end());
			//string rw = tl.reverse(w);
			rw = tl.tr(rw, "ATGC", "TACG");

			if(HH.find(w) != HH.end())
				PRstr = PRstr + itos(HH[w]["n"]);
			else if(HH.find(rw) != HH.end())
				PRstr = PRstr + itos(HH[rw]["n"]);
			else 
				PRstr = PRstr + "1";
			PR.push_back(PRstr);
			
			RDnum ++;
			if(RDnum % 100000 == 0)
				cerr << "Added " << RDnum << " reads\n";
		}
		return RDnum;
	}

	void printPR(string file, vector<string> &rPR){
		ofstream PRF;
		char file_char[file.length()];
		strcpy(file_char, file.c_str());
		PRF.open(file_char);
		PRF.close();
		PRF.open(file_char, ofstream::app);
		for(int i = 0; i < rPR.size(); i++)
			PRF << rPR[i] << "\n";
		PRF.close();
	}
};*/
	
// connect the contigs and make a graph
class allpaths{
	
public:
	
	int k;
	int n;
	float cov;
	map<int, map<string, string> > rcontigs;
	string f_contig;
	int min_degree;
	map<int, map<string, string> > nodes;
	vector<int> bnodes;
	int graph;// 0 or 1 for printing the graph
	
	allpaths(){
		min_degree = 2;
		k = 25;
		n = 100;
		cov = 3;
		graph = 0;
		f_contig = "";
	}
	
	~allpaths(){
		rcontigs = map<int, map<string, string> > ();
		nodes = map<int, map<string, string> >();
		bnodes = vector<int>();
	}
	
	void set_kmer_size(int k_){
		k = k_;
	}
	
	void set_n(int n_){
		n = n_;
	}
	
	void set_cov(float cov_){
		cov = cov_;
	}
	
	void set_rcontig(map<int, map<string, string> > &contig){
		rcontigs = contig;
	}
	
	void set_file_name(string filename_){
		f_contig = filename_;
	}
	
	void set_degree(int degree_){
		min_degree = degree_;
	}
	
	void set_graph(){
		graph = 1;
	}
	
	void doit(map<int, map<string, string> > &Contigs, map<int, map<string, string> > &retContigs){
		tools tl;
		rcontigs = Contigs;
		for(map<int, map<string, string> >::iterator rcontigs_it = rcontigs.begin(); rcontigs_it != rcontigs.end(); rcontigs_it ++){
			map<string, string> contig = rcontigs[(*rcontigs_it).first];
			int id = atoi((contig["id"]).c_str());
			nodes[id] = contig;
			if(contig["types"].length() != 0 && atof(contig["covs"].c_str()) > cov){
				string sep = "";
				vector<string> degree;
				tl.split(contig["types"], sep, degree);
				if(atoi(degree[0].c_str()) >= min_degree)
					bnodes.push_back(-id);
				if(atoi(degree[1].c_str()) >= min_degree)
					bnodes.push_back(id);
			}
		}
		
		if(bnodes.size() -1 >= n){
			cerr << "Skip Graph size too large! \n";
			return;
		}
		
		vector<string> longest_uniq_paths;
		vector<string> allpaths;
		string pathstr = "";
		for(int i = 0; i < bnodes.size(); i++){
			int id = bnodes[i];
			map<int, map<int, string> > edges;
			CreateGraph(id, f_contig, graph, edges);
			map<int, map<int, int> > visited;
			vector<string> paths;
			getPath(id, edges, visited, paths);
			for(int i = 0; i < paths.size(); i++){
				string p = paths[i];
				allpaths.push_back(p);
			}
		}
		
		string sep = ".";
		string minus = "-";
		tl.sort_by_longest_length(allpaths);
		for(int i = 0; i < allpaths.size(); i++){
			// Mirror path
			string p = allpaths[i];
			vector<string> rpns;
			tl.split(p, sep, rpns);
			for(int i = 0; i <= rpns.size() -1; i++)
				rpns[i] = itos(-atoi(rpns[i].c_str()));
				//rpns[i] = rpns[i] > 0 ? -rpns[i]:abs(rpns[i]);
			string ap = rpns[rpns.size() -1];
			for(int i = rpns.size() - 2; i >= 0; i--)
				ap += "." + rpns[i];
			size_t rindex1 = pathstr.rfind(p);
			size_t rindex2 = pathstr.rfind(ap);
			if(rindex1 != string::npos){ // not sure if this is the case
				string char_ = pathstr.substr(rindex1+p.length(), 1);
				if(atoi(char_.c_str()) == 0)
					continue;
			}
			if(rindex2 != string::npos){
				string char_ = pathstr.substr(rindex2+ap.length(), 1);
				if(atoi(char_.c_str()) == 0)
					continue;
			}
			pathstr += "|" + p;
			longest_uniq_paths.push_back(p);
		}
		
		int ii = 1;
		tl.sort_by_longest_length(longest_uniq_paths);
		for(int i = 0; i < longest_uniq_paths.size(); i++){
			// edit fasta
			string p = longest_uniq_paths[i];
			vector<string> pns;
			tl.split(p, sep, pns);
			string fasta = "";
			float kmercovsum = 0;
			for(int i = 0; i <= pns.size() -1; i++){
				string nid = pns[i];
				string fa = nodes[abs(atoi(nid.c_str()))]["seq"];
				if(atoi(nid.c_str()) < 0){
					fa = tl.tr(fa, "ACGT", "TGCA");
					reverse(fa.begin(), fa.end());
					//fa = tl.reverse(fa);
				}
				if(fasta.length() == 0)
					fasta = fa;
				else {
					fa = fa.substr(k-1);
					fasta += fa;
				}
				kmercovsum += atof(nodes[abs(atoi(nid.c_str()))]["covs"].c_str()) * (atof(nodes[abs(atoi(nid.c_str()))]["lens"].c_str()) - k + 1);
			}
			map<string, string> contig;
			int lens = fasta.length();
			float avgkmercov = lens > 0 ? float(int(kmercovsum*100/lens))/float(100):0;
			contig["id"] = p;
			contig["seq"] = fasta;
			contig["lens"] = itos(lens);
			contig["covs"] = ftos(avgkmercov);
			retContigs[ii] = contig;
			ii ++;
		}
	}

	
	void getPath(int id, map<int, map<int, string> > &edges, map<int, map<int, int> > &visited, vector<string> &path){
		// recursion
		string sep = ".";
		if(edges[id].size() > 0){
			for(map<int, string>::iterator edges_id_it = edges[id].begin(); edges_id_it != edges[id].end(); edges_id_it ++){
				int nd = (*edges_id_it).first;
				if(visited.find(id) != visited.end() && visited[id].find(nd) != visited[id].end())
					continue; // avoid repeat
				visited[id][nd] = 1;
				vector<string> spath;
				getPath(nd, edges, visited, spath);
				for(int i = 0; i < spath.size(); i++)
					path.push_back(itos(id) + sep + spath[i]);
			}
		}
		else{
			path.push_back(itos(id));
		}
	}	
	
	void CreateGraph(int seed, string f_contig, int f_graph, map<int, map<int, string> > &dedges){
		vector<int> tails;
		tails.push_back(seed);
		//map<> cn;
		map<int, map<int, string> > edges;
		map<int, int> snodes;
		while(tails.size() > 0){
			map<int, int> newtails;
			for(int i = 0; i < tails.size(); i++){ // breadth first search
				int t = tails[i];
				map<int, string> neighbor;
				InNodes(t, neighbor);
				OutNodes(t, neighbor);
				snodes[t] = 1;
				for(map<int, string>::iterator neighbor_it = neighbor.begin(); neighbor_it != neighbor.end(); neighbor_it ++){
					int n = (*neighbor_it).first;
					if(n == t)
						continue;	// skip looping back to itself
					if(edges[t].find(n) != edges[t].end())
						continue;
					edges[t][n] = neighbor[n];
					
					if(snodes.find(n) == snodes.end()){ // a new node
						newtails[n] = 1;
					}
				}
			}
			tails = vector<int>();
			for(map<int, int>::iterator newtails_it = newtails.begin(); newtails_it != newtails.end(); newtails_it++)
				tails.push_back((*newtails_it).first);
		}
		
		if(f_graph == 1){
			// visualization; skip temporarily
		}

		map<int, int> visited;
		tails = vector<int>();
		tails.push_back(seed);
		int nedges = 0;
		while(tails.size() > 0){
			map<int, int> newtails;
			for(int i = 0; i < tails.size(); i++){ // breadth first search
				int t = tails[i];
				visited[t] = 1;
				map<int, string> outnodes;
				OutNodes(t, outnodes);
				for(map<int, string>::iterator edges_t_it = edges[t].begin(); edges_t_it != edges[t].end(); edges_t_it ++){
					int n = (*edges_t_it).first;
					if(edges[n].find(t) != edges[n].end() && outnodes.find(n) != outnodes.end()){ // only care end to end connections
						// g->add_edge
						dedges[t][n] = edges[t][n];
						if(visited.find(n) == visited.end())
							newtails[n] = 1;
						nedges ++;
					}
				}
			}
			tails = vector<int>();
			for(map<int, int>::iterator newtails_it = newtails.begin(); newtails_it != newtails.end(); newtails_it++)
				tails.push_back((*newtails_it).first);
		}

		if(f_graph == 1 && nedges > 0){
			// visualization: save image; skip temporarily
		}
	}
		
		
	void InNodes(int id, map<int, string> &neighbor){
		tools tl;
		map<string, string> node = nodes[abs(id)];
		string sep = ",";
		string sep1 = ":";
		string minus = "-";
		if(id > 0){
			if(node.find("I") != node.end()){
				vector<string> nis;
				tl.split(node["I"], sep, nis);
				for(int i = 0; i < nis.size(); i++){
					string ni = nis[i];
					vector<string> no_split;
					tl.split(ni, sep1, no_split);
					string id_split = no_split[0];
					string nreads_split = no_split[1];
					if(id_split.find("M") != string::npos)
						continue;
					int id = atoi(id_split.c_str());
					id = (id<0)?abs(id):-abs(id);
					neighbor[id] = nreads_split;
				}
			}
		}
		else{
			if(node.find("O") != node.end()){
				vector<string> nis;
string tmp = node["O"];
				tl.split(node["O"], sep, nis);
				for(int i = 0; i < nis.size(); i++){
					string ni = nis[i];
					vector<string> no_split;
					tl.split(ni, sep1, no_split);
					string id_split = no_split[0];
					string nreads_split = no_split[1];
					if(id_split.find("M") != string::npos)
						continue;
					int id = atoi(id_split.c_str());
					int end = (id>0)?1:0;  //inbound reverse in-node orientation
					// 1: the start end of a contig, 0: the end of a contig
					id = (id<0)?abs(id):-abs(id);
					neighbor[id] = nreads_split;
				}
			}
		}
	}
			
			
				
		
	void OutNodes(int id, map<int, string> &neighbor){
		tools tl;
		map<string, string> node = nodes[abs(id)];
		string sep = ",";
		string sep1 = ":";

		if(id > 0){
			vector<string> nos;
			tl.split(node["O"], sep, nos);
			for(int i = 0; i < nos.size(); i++){
				string no = nos[i];
				vector<string> no_split;
				tl.split(no, sep1, no_split);
				string id_split = no_split[0];
				string nreads_split = no_split[1];
				if(id_split.find("M") != string::npos)
					continue;
				int id = atoi(id_split.c_str());
				neighbor[id] = nreads_split;
			}
		}
		else {
			vector<string> nos;
			tl.split(node["I"], sep, nos);
			for(int i = 0; i < nos.size(); i++){
				string no = nos[i];
				vector<string> no_split;
				tl.split(no, sep1, no_split);
				string id_split = no_split[0];
				string nreads_split = no_split[1];
				if(id_split.find("M") != string::npos)
					continue;
				int id = atoi(id_split.c_str());
				neighbor[id] = nreads_split;
			}
			
		}
	}

		
	/*int byLongestLength(string a, string b){
		tools tl;
		string sep = ".";
		vector<string> na = tl.split(a, sep);
		vector<string> nb = tl.split(b, sep);
		if(nb.size() -1 > na.size() - 1)
			return 1;
		else if(nb.size() - 1 < na.size() - 1)
			return -1;
		else
			return 0;
		return 0;
	}*/
};

class main_functions{
	
public:
	
	int debug;// = 0;
	int low_kmer;
	int high_kmer;
	int kmer_size;// = 0;
	int tip_number;
	int total;
	string estimate_STR;
	
	map<int, map<string, string> > contigs;
	map<int, int> contigtips;
    	vector<string> fastas;
	vector<string> Reads;
	map<int, map<string, string> > to_print_contigs;
	map<string, int> hh;
	
	main_functions(int debug_, int low_kmer_, int high_kmer_, map<int, map<string, string> > &contigs_, map<int, int> &contigtips_, vector<string> &fastas_, int tip_number_, string estimate_STR_){
		debug = debug_;
		low_kmer = low_kmer_;
		high_kmer = high_kmer_;
		contigs = contigs_;
		contigtips = contigtips_;
		fastas = fastas_;
		tip_number = tip_number_;
		estimate_STR = estimate_STR_;
		total = 0;
	}
	
	void set_kmer_size(int kmer_size){
		kmer_size = kmer_size;
		return;
	}
	
	/*void getReads(int fastas_, int contigs_, vector<string> &Reads, string fasta_str, vector<string> &pReads){
		tools tl;
		vector<string> fastas_local;
		//vector<string> pReads;
		if(fastas_ == 1)
			fastas_local = fastas;
		if(fasta_str.length() != 0)
			fastas_local.push_back(fasta_str);			
			
		if(fastas_ == 1){
			for(int i = 0; i < fastas_local.size(); i++){
				string fasta = fastas_local[i];
				ifstream FASTAS;
				FASTAS.open(fasta.c_str());
				string header;
				char line_[500];
				if(FASTAS.is_open()){
					while(FASTAS.good()){
						FASTAS.getline(line_, 500);
						string line(line_);
						line = tl.chomp(line);
						if(line.substr(0,1).compare(">")==0)
							header = line;
						FASTAS.getline(line_, 500);
						line = line_;
						line = tl.chomp(line);
						if(line.length() > 0)
							Reads.push_back(line); // Xian: Reads from fasta
						if(estimate_STR.length() > 0 && header.find(estimate_STR) != string::npos && line.length() > 0)
							pReads.push_back(line);						
					}
					FASTAS.close();
				}
			}
		}
		
		if(contigs_ == 1){
			for(int i = 1; i <= contigs.size(); i++){
				map<string, string> contig = contigs[i];
				if(! (atoi(contig["lens"].c_str()) > kmer_size + 5))
					continue;
				Reads.push_back(contig["seq"]);
			}
		}
		return;
	}*/

	void getReads(int fastas_, int contigs_, vector<string> &Reads, string fasta_str, vector<string> &pReads){
		if(fastas_ == 1){
			string header;
			for(int i = 0; i < fastas.size(); i+=2){
				if(fastas[i][0] == '>'){
					header = fastas[i];
					if(i+1 < fastas.size() && fastas[i+1].length() > 0){
						Reads.push_back(fastas[i+1]);
						if(estimate_STR.length() > 0 && header.find(estimate_STR) != string::npos)
							pReads.push_back(fastas[i+1]);
					}
				}
			}
		}
		fastas = vector<string>();
		if(contigs_ == 1){
			for(int i = 1; i <= contigs.size(); i++){
				map<string, string> contig = contigs[i];
				if(! (atoi(contig["lens"].c_str()) > kmer_size + 5))
					continue;
				Reads.push_back(contig["seq"]);
			}
		}
		return;
	}
	
	void initial_iteration(int kmersize){
		/////////////////////////////////// asm_1_kmergen2.pl	
		kmergen kg;
		kg.set_kmer_size(kmersize);
		kg.c = low_kmer;
		kg.C = high_kmer;
		map<string, int> hh; // Xian: map<kmer, count>
		int total = kg.doit(Reads, hh);
		if(debug == 1)
			kg.printMer("Mer", hh);
		
//                cout << "Before graphgen:\n";
		graphgen gg;
		gg.set_kmer_size(kmersize);
		gg.set_HH(hh);
		map<string, test_HH > HH;
		vector<string> PR;
		int RDnum = gg.doit(Reads, PR);
		HH = gg.get_HH();
		if(debug == 1){
			printnodes("Mynodes", HH);
			gg.printPR("MyREADnum", PR);
		}
		
  //              cout << "Before walknode:\n";
		// from de bruigin graph to proto-contig graph
		walknodes wn;
		wn.set_kmer_size(kmersize);
		wn.strictwalk(HH);
		HH = wn.HH;
		map<int, map<string, string> > protocontigs;
		wn.dump_protocontigs(protocontigs);
		if(debug == 1){
			printnodes("Mynodes.strwlk1", HH);
			int switch_ = 1;
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.strwlk1");
			outputcontigs(ofs, switch_, 0);
			//outputcontigs("Mycontigs.strwlk1", switch_);
			ofs.close();
		}
	
    //    cout << "Before processtips:\n";        
		// compute tip value: max{nucleotide distances to leaves} for all proto-contigs and the anti-proto-contigs, start from the leaves
		processtips pt;
		pt.set_kmer_size(kmersize);
		pt.tipwrap(tip_number, protocontigs, HH);
		contigtips = pt.Contigtips;
		if(debug == 1){
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.tiplabel1");
			outputcontigs(ofs, 0, 0);
			//outputcontigs("Mycontigs.tiplabel1", 0);
			ofs.close();
		}

//	cout << "Before addbridgekmer:\n";	
		//recover low frequency kmers in high quality reads that bridge separated non-tip proto-contig graphs
		addbridgekmer ab;
		ab.set_kmer_size(kmersize);
		vector<string> newkmers;
		ab.doit(HH, protocontigs, contigtips, PR, Reads);
		HH = ab.HH;
		newkmers = ab.newKmers;
		if(debug == 1){
			printnodes("Mynodes.add", HH);
			ofstream FOUT;
			FOUT.open("Newmer", ofstream::app);
			if(FOUT.is_open()){
				for(int i = 0; i < newkmers.size(); i++)
					FOUT << newkmers[i] << "\n";
			}
			FOUT.close();
		}
		
//                cout << "Before walknode:\n";
		contigtips = map<int, int>();
		wn.strictwalk(HH);
		HH = wn.HH;
		wn.dump_protocontigs(protocontigs);
		if(debug==1){
			printnodes("Mynodes.strwlk2", HH);
			int switch_ = 1;
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.strwlk2");
			outputcontigs(ofs, switch_, 0);
			//outputcontigs("Mycontigs.strwlk2", switch_);
			ofs.close();
		}
		
  //              cout << "Before second processtips:\n";
		//updated the set of proto-contigs with the expanded hash
		processtips pt1;
		pt1.set_kmer_size(kmersize);
		pt1.tipwrap(tip_number, protocontigs, HH);
		contigtips = pt1.Contigtips;
		if(debug == 1){
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.tiplabel2");				
			outputcontigs(ofs, 0, 0);
			//outputcontigs("Mycontigs.tiplabel2", 0);
			ofs.close();
		}
			
//		cout << "Before walkcontig:\n";
		//extend proto-contigs to contigs by removing tips and collapse bubbles with heuristic cut-offs
		walkcontig wc;
		wc.set_kmer_size(kmersize);
		wc.walkcontigwrap(3, 0.3, HH, contigtips, protocontigs);
		HH = wc.HH;
//		contigtips = wc.Contigtips;
		contigs = protocontigs;
		map<int , map<string, string> > contigs2;
		wc.dump_contigs2(contigs2);
		if(debug==1){
			printnodes("Mynodes.wlkcon", HH);
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.wlkcon");
			outputcontigs(ofs, 3, 0);
			ofs.close();
			//outputcontigs("Mycontigs.wlkcon", 3);
			to_print_contigs = contigs2;
			ofstream ofs1("Mycontigs2.wlkcon");
			outputcontigs(ofs1, 2, 0);
			ofs1.close();
			//outputcontigs("Mycontigs2.wlkcon", 2);
		}
		
//                cout << "Before breaktip:\n";
		//relabel tips on proto-contigs connected to the middle of a contig
		pt1.breaktip(tip_number, HH, contigtips, contigs, contigs2);
		HH = pt1.HH;
		contigs = pt1.Contigs;
		contigs2 = pt1.Contigs2;
		contigtips = pt1.Contigtips;
		if(debug == 1){
			printnodes("Mynodes.btip", HH);
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.btip");
			outputcontigs(ofs, 3, 0);
			//outputcontigs("Mycontigs.btip", 3);
			ofs.close();
		}
		
//                cout << "Before second time walkcontig:\n";
		//similar to 6, different param, more sensitive to weak branches
		walkcontig wc1;
		wc1.set_kmer_size(kmersize);
		wc1.walkcontigwrap(3, 0.2, HH, contigtips, protocontigs);
		contigs = protocontigs;
		HH = wc1.HH;
		if(debug == 1){
			printnodes("Mynodes.wlkcon2", HH);
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.wlkcon2");
			outputcontigs(ofs, 1, 0);
			//outputcontigs("Mycontigs.wlkcon2", 1);
			ofs.close();
			to_print_contigs = contigs2;
			ofstream ofs1("Mycontigs2.wlkcon2");
			outputcontigs(ofs1, 2, 0);
			//outputcontigs("Mycontigs2.wlkcon2", 2);
			ofs1.close();
		}
		
//                cout << "Before maprdtocontig:\n";
		//use entire read length to resolve small repeats
		maprdtocontig mr;
		mr.set_kmer_size(kmersize);
		mr.doit(Reads, HH, contigtips, contigs, contigs2);
		contigs = map<int, map<string, string> >();
		mr.dump_contigs2_m(contigs);
		if(debug == 1){
			mr.printmapping("Mapfile");
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs2.scaf");
			outputcontigs(ofs, 4, 0);
			//outputcontigs("Mycontigs2.scaf", 4);
			ofs.close();
		}
		to_print_contigs = contigs;
		return;
	}
	
	void iteration(int kmersize, map<int, map<string, string> > &fakereads, vector<string> &r_Reads){
		Reads = r_Reads;
		vector<string> fakeReads;
		contigs = fakereads;
		int contigs_ = 1;
		int fastas_ = 0;
		kmer_size = kmersize;
		vector<string> pReads;
		string fasta_ = "";
		getReads(fastas_, contigs_, fakeReads, fasta_, pReads);
		
		kmergen kg;
		kg.set_kmer_size(kmersize);
		kg.c = low_kmer;
		kg.C = high_kmer;
		vector<string> totalReads = Reads;
		for(int ii = 0; ii < fakeReads.size(); ii++)
			totalReads.push_back(fakeReads[ii]);
		for(int ii = 0; ii < fakeReads.size(); ii++)
			totalReads.push_back(fakeReads[ii]);
		map<string, int> hh;
		int total = kg.doit(totalReads, hh);
		if(debug == 1)
			kg.printMer("Mer", hh);
		
		totalReads = Reads;
		for(int ii = 0; ii < fakeReads.size(); ii++)
			totalReads.push_back(fakeReads[ii]);
		graphgen gg;
		gg.set_kmer_size(kmersize);
		gg.set_HH(hh);
		map<string, test_HH > HH;
		vector<string> PR;
		int RDnum = gg.doit(totalReads, PR);
		HH = gg.get_HH();
		if(debug == 1){
			printnodes("Mynodes", HH);
			gg.printPR("MyREADnum", PR);
		}
		
		//from de bruigin graph to proto-contig graph
		walknodes wn;
		wn.set_kmer_size(kmersize);
		wn.strictwalk(HH);
		HH = wn.HH;
		map<int, map<string, string> > protocontigs;
		wn.dump_protocontigs(protocontigs);
		if(debug == 1){
			printnodes("Mynodes.strwlk1", HH);
			int switch_ = 1;
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.strwlk1");
			outputcontigs(ofs, switch_, 0);
			//outputcontigs("Mycontigs.strwlk1", switch_);
			ofs.close();
		}
		
		// compute tip value: max{nucleotide distances to leaves} for all proto-contigs and the anti-proto-contigs, start from the leaves
		processtips pt;
		pt.set_kmer_size(kmersize);
		pt.tipwrap(tip_number, protocontigs, HH);
		contigtips = pt.Contigtips;
		if(debug == 1){
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.tiplabel1");
			outputcontigs(ofs, 0, 0);
			//outputcontigs("Mycontigs.tiplabel1", 0);
			ofs.close();
		}
		
		
		//recover low frequency kmers in high quality reads that bridge separated non-tip proto-contig graphs
		addbridgekmer ab;
		ab.set_kmer_size(kmersize);
		vector<string> newkmers;
		//vector<string> PR;
		ab.doit(HH, protocontigs, contigtips, PR, Reads);
		HH = ab.HH;
		newkmers = ab.newKmers;
		if(debug == 1){
			printnodes("Mynodes.add", HH);
			ofstream FOUT;
			FOUT.open("Newmer");
			if(FOUT.is_open()){
				for(int i = 0; i < newkmers.size(); i++)
					FOUT << newkmers[i] <<"\n";
			}
			FOUT.close();
		}
		
		contigtips = map<int, int>();
		wn.strictwalk(HH);
		HH = wn.HH;
		wn.dump_protocontigs(protocontigs);
		if(debug==1){
			printnodes("Mynodes.strwlk2", HH);
			int switch_ = 1;
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.strwlk2");
			outputcontigs(ofs, switch_, 0);
			//outputcontigs("Mycontigs.strwlk2", switch_);
			ofs.close();
		}
		
		//updated the set of proto-contigs with the expanded hash
		processtips pt1;
		pt1.set_kmer_size(kmersize);
		pt1.Thin = 2;
		pt1.tipwrap(tip_number, protocontigs, HH);
		contigtips = pt1.Contigtips;
		if(debug == 1){
			to_print_contigs = protocontigs;
			ofstream ofs("Mycontigs.tiplabel2");				
			outputcontigs(ofs, 0, 0);
			//outputcontigs("Mycontigs.tiplabel2", 0);
			ofs.close();
		}
		
		//extend proto-contigs to contigs by removing tips and collapse bubbles with heuristic cut-offs
		walkcontig wc;
		wc.set_kmer_size(kmersize);
		wc.walkcontigwrap(3, 0.3, HH, contigtips, protocontigs);
		HH = wc.HH;
		map<int, map<string, string> > contigs2;
		wc.dump_contigs2(contigs2);
		map<int, map<string, string> > contigs = protocontigs;
		if(debug==1){
			printnodes("Mynodes.wlkcon", HH);
			int switch_ = 3;
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.wlkcon");
			outputcontigs(ofs, 3, 0);
			ofs.close();
			//outputcontigs("Mycontigs.wlkcon", 3);
			to_print_contigs = contigs2;
			ofstream ofs1("Mycontigs2.wlkcon");
			outputcontigs(ofs1, 2, 0);
			ofs1.close();
			//outputcontigs("Mycontigs2.wlkcon", 2);
		}
		
		//relabel tips on proto-contigs connected to the middle of a contig
		pt1.thicken(2.5, HH, contigtips, contigs, contigs2);
		HH = pt1.HH;
		contigs = pt1.Contigs;
		contigs2 = pt1.Contigs2;
		contigtips = pt1.Contigtips;
		if(debug == 1){
			printnodes("Mynodes.btip", HH);
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.btip");
			outputcontigs(ofs, 3, 0);
			//outputcontigs("Mycontigs.btip", 3);
			ofs.close();
		}
		
		//similar to 6, different param, more sensitive to weak branches
		walkcontig wc1;
		wc1.set_kmer_size(kmersize);
		wc1.walkcontigwrap(3, 0.2, HH, contigtips, protocontigs);
		if(debug == 1){
			printnodes("Mynodes.wlkcon2", HH);
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs.wlkcon2");
			outputcontigs(ofs, 1, 0);
			//outputcontigs("Mycontigs.wlkcon2", 1);
			ofs.close();
			to_print_contigs = contigs2;
			ofstream ofs1("Mycontigs2.wlkcon2");
			outputcontigs(ofs1, 2, 0);
			//outputcontigs("Mycontigs2.wlkcon2", 2);
			ofs1.close();
		}
		
		//use entire read length to resolve small repeats
		maprdtocontig mr;
		mr.set_kmer_size(kmer_size);
		mr.doit(Reads, HH, contigtips, contigs, contigs2);
		contigs = map<int, map<string, string> >();
		mr.dump_contigs2_m(contigs);
		if(debug == 1){
			mr.printmapping("Mapfile");
			to_print_contigs = contigs;
			ofstream ofs("Mycontigs2.scaf");
			outputcontigs(ofs, 4, 0);
			//outputcontigs("Mycontigs2.scaf", 4);
			ofs.close();			
		}
		to_print_contigs = contigs;
		return;
		
	}
	
	void printnodes(string filename, map<string, test_HH > &HH){
		ofstream fh;
		fh.open(filename.c_str());
		for(map<string, test_HH >::iterator HH_it = HH.begin(); HH_it != HH.end(); HH_it++){
			string key = (*HH_it).first;
			fh << key << " " << HH[key].n << " " << HH[key].AI << " " << HH[key].TI<< " " << HH[key].GI<< " " << HH[key].CI<< " " << HH[key].AO<< " " << HH[key].TO<< " " << HH[key].GO<< " " << HH[key].CO<< " " << HH[key].tag<< " " << HH[key].tag2<<"\n";
		}
		fh.close();
	}
	
	void outputcontigs(std::ostream &fh, int switch_, int hh_turn){
		//ofstream fh;
		//fh.open(filename.c_str());
		tools tl;
		for(map<int, map<string, string> >::iterator contig_ = to_print_contigs.begin(); contig_ != to_print_contigs.end(); contig_++){
			map<string, string> contig = (*contig_).second;
			if(contig.size() == 0)
				continue;
			int contig_kmerUtil = 0;
			string case_seq = "";
			
			if(hh.size() != 0 && hh_turn){
				string seq = contig["seq"];
				case_seq = KmerUtility(kmer_size, seq, hh, total, &contig_kmerUtil);
			}
			else {
				case_seq = tl.StringToLower(contig["seq"]);
			}
			
			if(switch_ == 1){
				fh << ">Contig" << contig["id"] << " " << contig["lens"] << " " << contig["covs"] << " " << contig["types"] << " ";
				if(contigtips.find(atoi(contig["id"].c_str())) == contigtips.end())
					fh << tip_number << " ";
				else
					fh <<contigtips[atoi(contig["id"].c_str())]<<" ";
				if(contigtips.find(-atoi(contig["id"].c_str())) == contigtips.end())
					fh << tip_number << " ";
				else
					fh <<contigtips[-atoi(contig["id"].c_str())] << " ";
 				fh << contig["tags"] << " I" << contig["I"] << " O" << contig["O"] << " " << contig_kmerUtil;
			}
			else if(switch_ == 2)
				fh << ">Contig" << contig["id"] << " " << contig["lens"] << " " << contig["covs"] << " " << contig["types"] << " I" << contig["I"]<< " O" << contig["O"] << " " << contig_kmerUtil;
			else if(switch_ == 3){
				fh << ">Contig"<<contig["id"]<<" "<<contig["lens"]<<" "<<contig["covs"]<<" "<<contig["types"]<<" ";
				if(contigtips.find(atoi(contig["id"].c_str())) == contigtips.end())
					fh<< tip_number << " ";
				else
					fh << contigtips[atoi(contig["id"].c_str())] << " ";
				if(contigtips.find(-atoi(contig["id"].c_str())) == contigtips.end())
					fh<< tip_number << " ";
				else
					fh << contigtips[-atoi(contig["id"].c_str())] << " ";
				fh << contig["tags"] << " " << contig_kmerUtil;
			}
			else if(switch_ == 4)
				fh << ">Contig"<<contig["id"]<<" "<<contig["lens"]<<" "<<contig["covs"]<<" "<<contig["types"]<<" I"<<contig["I"]<<" O"<<contig["O"]<<" "<<contig["tags"] << " " << contig_kmerUtil;
			else if(switch_ == 5)
				fh << ">Contig"<<contig["id"]<<" "<<contig["lens"]<<" "<<contig["covs"]<<" "<<contig_kmerUtil; 
			else{
				fh << ">Contig"<<contig["id"]<<" "<<contig["lens"]<<" "<<contig["covs"]<<" "<<contig["types"]<<" ";
				if(contigtips.find(atoi(contig["id"].c_str())) == contigtips.end())
					fh << tip_number << " ";
				else
					fh <<contigtips[atoi(contig["id"].c_str())]<<" ";
				if(contigtips.find(-atoi(contig["id"].c_str())) == contigtips.end())
					fh << tip_number;
				else
					fh <<contigtips[-atoi(contig["id"].c_str())];
				fh << " "<<contig_kmerUtil;
			}
			fh << endl;// << contig["seq"] << endl;	
			if(case_seq.length() > 0)
				fh << case_seq << endl;
			else 
				fh << contig["seq"] << endl;

		}
		//fh.close();
	}
	
	string KmerUtility(int kmersize, string seq, map<string, int> &rhh, int total, int *utility){
		tools tl;
		int occur = 0;
		vector<int> uniqpos;
		for(int i = kmersize; i <= seq.length(); i++){
			string w = seq.substr(i-kmersize, kmersize);
			string rw = w;
			reverse(rw.begin(), rw.end());
			//string rw = tl.reverse(w);
			rw = tl.tr(rw, "ATGC", "TACG");
			
			if(rhh.find(w) != rhh.end()){
				occur += rhh[w];
				uniqpos.push_back(i);
			}
			else if(rhh.find(rw) != rhh.end()){
				occur += rhh[rw];
				uniqpos.push_back(i);
			}
		}
		
		seq = tl.StringToLower(seq);
		if(uniqpos.size() > 0){
			//string sep = "";
			//vector<string> bases;
			//tl.split(seq, sep, bases);
			//string seq_ = "";
			for(int j = 0; j < uniqpos.size(); j++){
				int pos = uniqpos[j];
				for(int i = pos - kmersize; i < pos; i++){
					//bases[i] = (toupper(bases[i]);
					//char tmp = toupper(seq[i]);
					//seq_ += tmp;
					seq.replace(i, 1, 1, toupper(seq[i]));
				}
			}
			//seq = seq_;
			// seq = "";
			//for(int i = 0; i < bases.size(); i++)
			//	seq = seq + bases[i];
			// Xian: Implemnt this in my own way
		}
		*utility = total > 0 ? int(occur*100/total) : 0;
		
		return seq;
	}
};						
			
class tigra{
public:
	
	string kmer_size; // -k
	float min_kmer_cov; // -c
	int low_kmer; // -m
	int high_kmer; // -M
	string assembly_file; // -o
	int tip_number; // -t
	string alternative_haplotype; // -h
	string estimate_STR; // -p 
	string reference_for_screen; // -r // now the string is the sequence
	int max_node; // -n
	int min_degree; // -N
	string graph_file; // -g
	int debug; // -d
	
	tigra(){
		kmer_size = "25";
		min_kmer_cov = 3;
		low_kmer = 2;
		high_kmer = 2e9;
		tip_number = 1000;
		max_node = 100;
		min_degree = 2;
		debug = 0;
	}
	
	int run_tigra(vector<string> &fastas){
		if(fastas.size() == 0){
			fprintf(stderr, "\n");
			fprintf(stderr, "tigra <fasta.fa>\n\n");
			fprintf(stderr, "Options: \n");
			fprintf(stderr, "       -k STRING		Specify Kmer sizes, Use comma to delimit multiple kmers\n");
			fprintf(stderr, "	-c FLOAT		Minimal average kmer coverage [%f] for a contig to be considered in the alternative haplotypes\n", min_kmer_cov);
			fprintf(stderr, "       -m INT          Lowest Kmer frequency in initial hashing [%d]\n", low_kmer);		 
			fprintf(stderr, "       -M INT          Highest Kmer frequency in initial hashing [%d]\n", high_kmer);		
			fprintf(stderr, "       -o STRING       File to save primary assembly contigs\n");		 
			fprintf(stderr, "	-p STRING		Estimate the utility of the subset of reads containing STR in header\n");
			fprintf(stderr, "	-r STRING		Provide a reference to screen for non-reference kmers\n");
			fprintf(stderr, "       -t INT          Tip number [%d]\n", tip_number);	
			fprintf(stderr, "       -h STRING       File to save alternative haplotypes\n");	
			fprintf(stderr, "	-n INT			Maxmimal number of nodes allowed for constructing alternative haplotypes [%d]\n", max_node); 
			fprintf(stderr, "	-N INT			Minimal number of degrees for a node to consider as start of alternative haplotype [%d]\n", min_degree);
			fprintf(stderr,	"	-g STRING		Save assembly graph(png) in FILE\n");
			fprintf(stderr, "       -d INT          Turn on debug mode, generate intermediate files\n");
			fprintf(stderr, "Version: %s (commit %s)\n", __g_prog_version, __g_commit_hash);

			fprintf(stderr, "\n");
			return 1;
		}
		
		
		// print options
		char options_[500];	
		sprintf(options_,"-c %f -m %d -M %d -t %d -n %d -N %d -d %d", min_kmer_cov, low_kmer, high_kmer, tip_number, max_node, min_degree, debug);
		options = options_;
		options += "-k" + kmer_size;
		options += "-o" + assembly_file;
		options += "-h" + alternative_haplotype;
		options += "-p" + estimate_STR;
		//options += "-r" + reference_for_screen;
		options += "-g" + graph_file;
		
		map<int, int> contigtips;
		map<int, map<string, string> > contigs;	
		//vector<string> fastas;
		
		/*int optindex = optind;
		 while(optindex < argc){
		 fastas.push_back(argv[optindex]);
		 optindex ++;
		 }*/
	
//        cout << "Before Main Function:\n";        
		main_functions mf(debug, low_kmer, high_kmer, contigs, contigtips, fastas, tip_number, estimate_STR);
		int fastas_ = 1;
		int contigs_ = 0;
		vector<string> Reads;
		string fasta_ = "";
		vector<string> pReads;
		mf.getReads(fastas_, contigs_, Reads, fasta_, pReads);
		mf.Reads = Reads;
		if(debug == 1){
			string fileName = "getReads";
			ofstream fh;
			fh.open(fileName.c_str());
			for(int i = 0; i < Reads.size(); i++)
				fh << Reads[i] << "\n";
			fh.close();
		}
		
		vector<string> kmers;
		string sep = ",";
		tools tl;
		tl.split(kmer_size, sep, kmers);
		string kmer_debug = kmers[0];
//                cout << "Before tigra:\n";
		mf.initial_iteration(atoi(kmers[0].c_str()));
  //              cout << "After first tigra:\n";
		mf.contigtips = map<int, int>();
		for(int i = 1; i < kmers.size(); i++){
			int kmer = atoi(kmers[i].c_str());
			mf.iteration(kmer, mf.contigs, Reads);
		}
//		cout << "After second tigra:\n";
		/*string refseq = "";
		 if(reference_for_screen.length() != 0){
		 refseq = reference_for_screen;// read that in kmergen rather than here
		 }*/
		
		int total_;
		map<string, int> phh;
		if(int(pReads.size()) - 1 >= 0){
			kmergen pkg;
			pkg.set_kmer_size(atoi(kmers[kmers.size() - 1].c_str()));
			pkg.c = low_kmer;
			pkg.C = high_kmer;
			//map<string, int> hh;
			pkg.set_filter(reference_for_screen); // need to go through filter in kmergen again
			total_ = pkg.doit(pReads, phh);
			if(debug == 1)
				pkg.printMer("MerP", phh);
		}
		
		// Output assembled contigs
		int switch_ = 4;
		mf.total = total_;
		mf.hh = phh;
		mf.kmer_size = atoi(kmers[kmers.size()-1].c_str());
		if(assembly_file.length() != 0){
			ofstream ofs(assembly_file.c_str());
			mf.outputcontigs(ofs, switch_, 1);
			ofs.close();
		}
		else {
			mf.outputcontigs(std::cout, switch_, 1);
		}
		
		
		// do the following when allpaths is ready
		if (alternative_haplotype.length() > 0 || graph_file.length() > 0) {
			allpaths ap;
			ap.set_kmer_size(atoi(kmers[kmers.size()-1].c_str()));
			ap.set_cov(min_kmer_cov);
			ap.set_n(max_node);
			ap.set_rcontig(mf.to_print_contigs);
			
			ap.set_degree(min_degree);
			ap.set_file_name(fastas[0]);
			ap.set_graph();
			map<int, map<string, string> > ret_contigs;
			ap.doit(mf.to_print_contigs, ret_contigs);
			mf.to_print_contigs = ret_contigs;
			int switch_ = 5;
			if(alternative_haplotype.length() != 0){
				ofstream ofs(alternative_haplotype.c_str());
				mf.outputcontigs(ofs, switch_, 5);
				ofs.close();
			}
			else {
				mf.outputcontigs(std::cout, switch_, 5);
			}
			
		}
	}
};		

class assemble {
public:
	string version;
	vector<string> bams;
	map<string, string> bamsmap;
	string BreakDancer_file;
	string PCR_file;
	int estimate_max_ins; // A
	int flanking_size; // l
	int pad_local_ref; // w
	int assemble_read_qual; // q
	int num_mismatch_poor_map; // N
	int high_depth_skip; // p
	int qual_threshold; // Q
	string library_to_skip_; // L
	string datadir; // I
	int write_to_ref; // r
	int write_to_read; // d
	string reference_file; // R
	string chromosome; // c
	int skip_call; // z
        int min_size_threshold; // M
        int max_node; // h
        string kmers;
	
        assemble(){		
		qual_threshold = 0;
	}
	
	
	string get_string(uint8_t *pt, int32_t length){
		char seq_[length + 1];
		for(int i = 0; i < length ; i++)
			seq_[i] = bam_nt16_rev_table[bam1_seqi(pt, i)];
		//seq_[i] = *(pt+i);
		seq_[length] = '\0';
		string seq(seq_);// the good way to turn char * to string
		return seq;
	}
	
	// prepare: read bam file by a chromosome
	pair64_t * ReadBamChr_prep(string chr_str, string bam_name, int *tid, int *beg, int *end, samfile_t *in, int *n_off){	
		char *bam_name_;
		bam_name_ = new char[bam_name.length()+1];
		strcpy(bam_name_, bam_name.c_str());

		pair64_t *off=NULL;
		bam_index_t *idx = 0;
		idx = bam_index_load(bam_name_);// index   		
		delete [] bam_name_;

		//if(idx == 0){
		//	off = (pair64_t*)calloc(1, 16);
		//	off[0].u = -1;
		//	off[0].v = -1;
		//	cout << "Error: should do sort and index first if specifying chromosome!" << endl;
		//	return off;
		//}
		if(idx){
		  char *chr_str_;
		  chr_str_ = new char[chr_str.length()+1];
		  strcpy(chr_str_, chr_str.c_str());
		  bam_parse_region(in->header, chr_str_, tid, beg, end);// parse
		  delete []chr_str_;


		  // return the file handle for handle
		  //*fp = in->x.bam;
		  //bamFile fp = in->x.bam;
		  off = get_chunk_coordinates(idx, *tid, *beg, *end, n_off);
		}
		bam_index_destroy(idx);
		return off;
	}
	
	// read bam file by a chromosome by one line; fp will track where we are
	int ReadBamChr(bam1_t *b, bamFile fp, int tid, int beg, int end, uint64_t *curr_off, int *i, int *n_seeks, pair64_t *off, int n_off){
		
		if (off == 0) return 0;
		
		if (*curr_off == 0 || (*i>=0 && *curr_off >= off[*i].v)) { // then jump to the next chunk
			if (*i == n_off - 1) return 0; // no more chunks
			if (*i >= 0 && *curr_off != off[*i].v) return 0;//assert(*curr_off == off[*i].v); // otherwise bug
			if (*i < 0 || (*i>=0 && off[*i].v != off[*i+1].u)) { // not adjacent chunks; then seek
				if(*i + 1 >= 0){
					bam_seek(fp, off[*i+1].u, SEEK_SET);
					*curr_off = bam_tell(fp);
					++(*n_seeks);
				}
			}
			++(*i);
		}
		if (bam_read1(fp, b) > 0) {
			*curr_off = bam_tell(fp);
			if (b->core.tid != tid || b->core.pos >= end) return 0; // no need to proceed, will jump out
			else if (is_overlap(beg, end, b)) return 1;
			else return -1; // not yet arrived
		} 
		else 
			return 0; // will jump out
		return 1;
	}
	
	// read PCR coor; need bam_related in the data type
	void ReadPCRCoor(string PCR_file, vector<BD_data> &coor){
		tools tl;
		ifstream PCR;
		PCR.open(PCR_file.c_str());
		char line_[1024];
		if(PCR.is_open()){
			while(PCR.good()){
				PCR.getline(line_, 1024);
				if(strcmp(line_, "CHR") == 0)
					continue;
				string line(line_);
				if(line.length() == 0)
					continue;
				line = tl.chomp(line);
				BD_data cr;
				vector<string> fields;
				tl.split(line, "\t", fields);
				cr.chr1 = fields[0];
				cr.chr2 = fields[0];
				cr.pos1 = atoi(fields[2].c_str());
				cr.pos2 = atoi(fields[3].c_str());
				cr.type = fields[5];
				cr.size = abs(atoi(fields[6].c_str()));
				cr.bam_related = fields[9];
				coor.push_back(cr);
			}
		}
		PCR.close();
	}
	
	void ReadBDCoor(string BreakDancer_file, vector<BD_data> &coor){
		tools tl;
		//vector<BD_data> coor;
		vector<int> col_cn;
		int col_gene = -1, col_database = -1;
		vector<string> headfields;
		vector<string> nlibs;
		if(library_to_skip_.length() != 0)
			tl.split(library_to_skip_, ",", nlibs);
		ifstream BD;
		BD.open(BreakDancer_file.c_str());
		char line_[30000];
                int call_num = 0;
		if(BD.is_open()){
			while (BD.good()) {
				BD.getline(line_, 30000);
                                call_num ++;
				if(line_[0] == '#'){
					// next unless(/Chr1\s/); skip
					// Parse the header line
					vector<string> headfields;
					tl.split(line_, "\t", headfields);
					int col_normalcn = -1;
					for(int i = 0; i < headfields.size(); i++){
						if(headfields[i].find(".bam") != string::npos){
							if(headfields[i].find("normal") != string::npos){
								col_normalcn = i;
							}
							else {
								col_cn.push_back(i);
							}
						}
						else if(headfields[i].find("Gene") != string::npos){
							col_gene = i;
						}
						else if(headfields[i].find("DataBases") != string::npos){
							col_database = i;
						}
					}
					if(col_normalcn != -1)
						col_cn.push_back(col_normalcn);
				}
				else {
                                    if(call_num < skip_call)
                                        continue;
					string line(line_);
					if(line.length() == 0)
						continue;
					line = tl.chomp(line);
					BD_data cr;
					vector<string> fields;
					tl.split(line, "\t", fields);
					/*
					if(fields.size()<8){
					  cerr<<"Must have at last 8 columns:Chr1,Pos1,Ori1,Chr2,Pos2,Ori2,Type,Size"<<endl;
					  exit EXIT_FAILURE;
					  }*/
					cr.chr1 = fields[0];
					cr.pos1 = atoi(fields[1].c_str());
					cr.ori1 = fields[2];
					cr.chr2 = fields[3];
					cr.pos2 = atoi(fields[4].c_str());
					cr.ori2 = fields[5];
					cr.type = fields[6];
					cr.size = abs(atoi(fields[7].c_str()));
                                        
					if(fields.size()>8) cr.score = atoi(fields[8].c_str()) ;
					if(fields.size()>9) cr.nreads = atoi(fields[9].c_str()) ;
					if(fields.size()>10) cr.nreads_lib = fields[10].c_str() ;
					if(fields.size() > 11){
                                            for(int j = 11; j < fields.size(); j++)
    					        cr.extra.push_back(fields[j]);
                                        }
					cr.line = line;
					
					if(cr.chr1.find("NT")!=string::npos || cr.chr1.find("RIB") != string::npos || cr.chr2.find("NT")!=string::npos || cr.chr2.find("RIB") != string::npos)
						continue;
					
					// skip two next
					// insert -Q
					if(cr.score < qual_threshold || cr.size < min_size_threshold)
						continue;
					
					// Ignore events detected in a library
					int ignore = 0;
					for(int i = 0; i < nlibs.size(); i++){
						string nlib = nlibs[i];
						if(cr.nreads_lib.find(nlib) != string::npos) {
							ignore = 1;
							break;
						}
					}
					
					// Include Copy Number Altered events if available
					/*if(col_cn.size() != 0){
					 for(int i = 0; i < col_cn.size() - 1; i++){
					 if(fields[col_cn[i]].find("NA") != string::npos || fields[col_cn[col_cn.size()-1]].find("NA") != string::npos)
					 continue;
					 float diff_cn = abs(atof(fields[col_cn[i]].c_str()) - atof(fields[col_cn[col_cn.size()-1]].c_str()));
					 //	$ignore=0 if($diff_cn>$opts{C}); skip
					 }
					 //vector<string> cns;
					 string tmp;
					 tmp = headfields[col_cn[0]] + ":" + fields[col_cn[0]];
					 for (int i = 1; i < col_cn.size(); i++) {
					 tmp += "," + headfields[col_cn[i]] + ":" + fields[col_cn[i]];
					 }
					 cr.cnstr = tmp;
					 }*/
					if(col_gene != -1 && fields.size() > col_gene)
						cr.gene = fields[col_gene];
					if(col_database != -1 && fields.size() > col_database)
						cr.database = fields[col_database];
					
					if(cr.line.find("cancer") != string::npos || cr.line.find("coding") != string::npos)
						ignore = 0;
					if(cr.line.find("gene") != string::npos || cr.type.find("ctx") != string::npos)
						ignore = 0;
					if(ignore > 0)
						continue;
					coor.push_back(cr);
					
				}
			}
		}
		BD.close();
	}
	
	int AssembleBestSV(string const& chr1, int start, string const& chr2, int end, string const& type, int size, string const& ori, BD_data SV, int a, int b){
		string a_str(itos(a));
		string b_str(itos(b));
		BD_data maxSV;
		tools tl;
		
		// skip
		int seqlen = 0;
		int nreads = 0;
		int nSVreads = 0;
		//map<> readhash;
		int start1, end1, start2, end2, regionsize, refsize;
		vector<string> refs;
		string posstr;

        stringstream prefixss;
        prefixss << chr1 << "."
            << start << "."
            << chr2 << "."
            << end << "."
            << type << "."
            << size << "."
            << ori;
        string prefix(prefixss.str());

		cerr << prefix << "\ta: " << a << "\tb: " << b << endl;
		int makeup_size = 0;
		int concatenated_pos = 0;
		vector<string> samtools;
		string pos_1, pos_2;
		for (int i = 0; i < bams.size(); i++) {
			string fbam = bams[i];
			if(type.compare("ITX") == 0){
				start1 = start - b;
				end1 = start + estimate_max_ins;
				start2 = end - estimate_max_ins;
				end2 = end + a;
			}
			else if(type.compare("INV") == 0){
				start1 = start - flanking_size;
				end1 = start + flanking_size;
				start2 = end - flanking_size;
				end2 = end + flanking_size;
			}
			else {
				start1 = start - flanking_size;
				end1 = start + a;
				start2 = end - b;
				end2 = end + flanking_size;
			}
			
			if(ori.compare("+-")  == 0 && chr1.compare(chr2) == 0 && start2 < end1){
				if(refs.size() == 0)
                                {
                                        int tmp1 = start - flanking_size - pad_local_ref;
                                        tmp1 = tmp1 > 0 ? tmp1 : 0;
                                        int tmp2 = end + flanking_size + pad_local_ref;
                                        tmp2 = tmp2 > 0 ? tmp2 : 0;
					refs.push_back(chr1 + ":" + itos(tmp1) + ":" + itos(tmp2)); //itos(start - flanking_size - pad_local_ref) + ":" + itos(end + flanking_size + pad_local_ref));
					refsize = end - start + 1 + 2*flanking_size + 2*pad_local_ref;
				}
				posstr = chr1 + "_" + itos(start1 - pad_local_ref) + "_" + chr1 + "_" + itos(start1 - pad_local_ref) + "_" + type + "_" + itos(size) + "_" + ori; // this may contain a bug
				pos_1 = itos(max(start1 - pad_local_ref,0));
				pos_2 = itos(max(start1 - pad_local_ref,0));
				// samtools view bam chr1:start1-end2;
				//string tmp = "samtools view " + fbam + " " + chr1 + ":" + itos(start1) + "-" + itos(end2);
				//samtools.push_back(tmp);
				
				samtools.push_back(fbam);
				samtools.push_back(chr1 + ":" + itos(start1) + "-" + itos(end2));
				regionsize = end2 - start1 + 1;
			}
			else {
				if(type.compare("CTX") == 0){
					if(refs.size() == 0){
                                                int tmp1 = start - flanking_size - pad_local_ref;
                                                tmp1 = tmp1 > 0 ? tmp1 : 0;
                                                int tmp2 = start + flanking_size + pad_local_ref;
                                                tmp2 = tmp2 > 0 ? tmp2 : 0;
						refs.push_back(chr1 + ":" + itos(tmp1) + ":" + itos(tmp2)); //itos(start - flanking_size - pad_local_ref) + ":" + itos(start + flanking_size + pad_local_ref));
                                                tmp1 = end - flanking_size - pad_local_ref;
                                                tmp1 = tmp1 > 0 ? tmp1 : 0;
                                                tmp2 = end + flanking_size + pad_local_ref;
                                                tmp2 = tmp2 > 0 ? tmp2 : 0;
						refs.push_back(chr2 + ":" + itos(tmp1) + ":" + itos(tmp2));//itos(end - flanking_size - pad_local_ref) + ":" + itos(end + flanking_size + pad_local_ref));
						refsize = 2*(flanking_size + pad_local_ref) + 1;
					}					
					posstr = chr1 + "_" + itos(start - flanking_size - pad_local_ref) + "_" + chr2 + "_" + itos(end - flanking_size - pad_local_ref) + "_" + type + "_" + itos(size) + "_" + ori;
					pos_1 = itos(max(start - flanking_size - pad_local_ref,0));
					pos_2 = itos(max(end - flanking_size - pad_local_ref,0));
				}
				else if(size > 99999){
					if(refs.size() == 0){
                                                
                                                int tmp1 = start - flanking_size - pad_local_ref;
                                                tmp1 = tmp1 > 0 ? tmp1 : 0;
                                                int tmp2 = start + flanking_size;
                                                tmp2 = tmp2 > 0 ? tmp2 : 0;
						refs.push_back(chr1 + ":" + itos(tmp1) + ":" + itos(tmp2)); // itos(start - flanking_size - pad_local_ref) + ":" + itos(start + flanking_size));
						
                                                tmp1 = end - flanking_size;
                                                tmp1 = tmp1 > 0 ? tmp1 : 0;
                                                tmp2 = end + flanking_size + pad_local_ref;
                                                tmp2 = tmp2 > 0 ? tmp2 : 0;
                                                refs.push_back(chr2 + ":" + itos(tmp1) + ":" + itos(tmp2)); //itos(end - flanking_size) + ":" + itos(end + flanking_size + pad_local_ref));
						refsize = 2*flanking_size + pad_local_ref + 1;
					}
					
					makeup_size = (end - flanking_size) - (start + flanking_size) - 1;
					concatenated_pos = 2*flanking_size + pad_local_ref;
					posstr = chr1 + "_" + itos(start - flanking_size - pad_local_ref) + "_" + chr2 + "_" + itos(end - 3*flanking_size - pad_local_ref) + "_" + type + "_" + itos(size) + "_" + ori;
					pos_1 = itos(max(start - flanking_size - pad_local_ref,0));
					pos_2 = itos(max(end - 3*flanking_size - pad_local_ref,0));
				}
				else {
					if(refs.size() == 0){

                                                int tmp1 = start - flanking_size - pad_local_ref;
                                                tmp1 = tmp1 > 0 ? tmp1 : 0;
                                                int tmp2 = end + flanking_size + pad_local_ref;
                                                tmp2 = tmp2 > 0 ? tmp2 : 0;
						refs.push_back(chr1 + ":" + itos(tmp1) + ":" + itos(tmp2));// itos(start - flanking_size - pad_local_ref) + ":" + itos(end + flanking_size + pad_local_ref));
						refsize = end - start + 1 + 2*(flanking_size + pad_local_ref);	
					}
					posstr = chr1 + "_" + itos(start1 - pad_local_ref) + "_" + chr1 + "_" + itos(start1 - pad_local_ref) + "_" + type + "_" + itos(size) + "_" + ori;
					pos_1 = itos(max((start1 - pad_local_ref),0));
					pos_2 = itos(max(start1 - pad_local_ref,0));
				}
				
				string reg1, reg2;
				if(ori.compare("+-") == 0){
					reg1 = chr1 + ":" + itos(start1) + "-" + itos(end1);
					reg2 = chr2 + ":" + itos(start2) + "-" + itos(end2);
				}
				else if(ori.compare("-+") == 0){
					reg1 = chr2 + ":" + itos(end - flanking_size) + "-" + itos(end + a);
					reg2 = chr1 + ":" + itos(start - b) + "-" + itos(start + flanking_size);
				}
				else if(ori.compare("++") == 0){
					reg1 = chr1 + ":" + itos(start - flanking_size) + "-" + itos(start + a);
					reg2 = chr2 + ":" + itos(end - flanking_size) + "-" + itos(end + b);
				}
				else if(ori.compare("--") == 0){
					reg1 = chr1 + ":" + itos(start - a) + "-" + itos(start + flanking_size);
					reg2 = chr2 + ":" + itos(end - b) + "-" + itos(end + flanking_size);
				}
				regionsize = a + b + 2*flanking_size + 1;
				
				samtools.push_back(fbam);
				samtools.push_back(reg1);
				//string tmp = "samtools view " + fbam + " " + reg1;
				//samtools.push_back(tmp);
				//tmp = "samtools view " + fbam + " " + reg2;
				//samtools.push_back(tmp);
				samtools.push_back(fbam);
				samtools.push_back(reg2);
			}
			
		}
		
		if(refsize > 1e6){
			cerr << "reference size " << refsize << " too large, skip ... " << endl;
			return 0;
		}
		
		string cmd;
		// skip;
		// create reference
		cerr << "Use samtools faidx to fetch fai\n";
		string ref_string(""); // variable to be passed to tigra
		if(reference_file.length() != 0){
			string ref_file = reference_file;
			ofstream REF;
			if(write_to_ref){
				// move here to handle open and append
				string write_to_ref_file = datadir + "/" + prefix + ".ref.fa";
				REF.open(write_to_ref_file.c_str());
			}
                        int ii = 0; // for the track of the time to write to ref
			for (int i = 0; i < refs.size(); i++) {
				string ref = refs[i];
				//cout << "ref region: " << ref << endl;
                                vector<string> tmp;
				tl.split(ref, ":", tmp);
				string chr_ref = tmp[0];
				string start_ref = tmp[1];
				string end_ref = tmp[2];
				
				//string ref_file = "/gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome." + chr_ref + ".fa";
				string region = chr_ref + ":" + start_ref + "-" + end_ref;
				// skip
				/*if(i == 0)
				 cmd = "expiece " + start_ref + " " + end_ref + " /gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome." + chr_ref + ".fa > " + datadir + "/" + prefix + ".ref.fa";
				 else
				 cmd = "expiece " + start_ref + " " + end_ref + " /gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome." + chr_ref + ".fa >> " + datadir + "/" + prefix + ".ref.fa";*/
				// skip
				//system(cmd.c_str());
				//cerr << cmd << endl;
				// use samtools faidx rather than expiece
				char *s=NULL;
				int l;
				faidx_t *fai;
				fai = fai_load(ref_file.c_str());
				if(fai == 0){
					cerr << "Error in reading reference file! Will skip this step. \n";
					continue;
				}
				cerr << "Reference: " << fai << endl;
				s = fai_fetch(fai, region.c_str(), &l);
				if(s==NULL || strlen(s)<=0){
				  cerr << "Failed to load the reference.  Please check chromosome naming between the reference fasta header and the SV calls." << endl;
				}
				else{
				  ref_string += s;
				  if(write_to_ref){
				    int j, k;
				    ii++;
				    if(REF.is_open()){
                                        if(makeup_size <= 0 || ii == 1)
                                        REF << ">" << region << "\n";// write the first header only when makeup_size > 0
				      for(j = 0; j < l; j+= 60){
					for(k = 0; k < 60 && k < l-j; ++k)
					  REF << s[j+k];
    				        REF << endl;
                                      }
				    }					
				  }
				}
				free(s);
				fai_destroy(fai);
			}
			if(REF.is_open())
				REF.close();
		}
		else
			cerr << "Skip Reading reference. \n";
		
/*		if((makeup_size > 0) && write_to_ref){ // piece together 2 refs as one
			string tmp = "head -n 1 " + datadir + "/" + prefix + ".ref.fa > " + datadir + "/" + prefix + ".ref.fa.tmp";
			system(tmp.c_str());
			tmp = "grep -v : " + datadir + "/" + prefix + ".ref.fa >> " + datadir + "/" + prefix + ".ref.fa.tmp";
			system(tmp.c_str());
			tmp = "mv " + datadir + "/" + prefix + ".ref.fa.tmp " + datadir + "/" + prefix + ".ref.fa";
			system(tmp.c_str());
		}*/
		
		string freads = datadir + "/" + prefix + /*".a" + a_str + ".b" + b_str +*/ ".fa";
		float avgdepth = 0;
		samfile_t *in = 0; // get ready
		char in_mode[5] = "rb";
		char *fn_list = 0;
		char *fn_ref = 0;
		char *fn_rg = 0;
		
		vector<string> buffer;
		ifstream ifile(freads.c_str());
		if(1){//! ifile || ifile){ // no I skip
			//list<string> buffer;
			map<string, int> uniqreadnames;
			for (int i = 0; i < samtools.size(); i+=2) {
				string scmd = samtools[i];
				vector<read_data> reads;
				map<string, uint32_t> max_mapqual;
				map<string, int32_t> min_NM;
				cerr << scmd << endl;
				//system(scmd);
				char *bam_name_;
				bam_name_ = new char[samtools[i].length()+1];
				strcpy(bam_name_, samtools[i].c_str());
				if ((in = samopen(bam_name_, in_mode, fn_list)) == 0) {
					fprintf(stderr, "[main_samview] fail to open file for reading.\n");
					return 0;
				}
				delete [] bam_name_;
				if (in->header == 0) {
					fprintf(stderr, "[main_samview] fail to read the header.\n");
					return 0;
				}
				bamFile fp = in->x.bam;
				bam1_t *b = bam_init1();
				// chromosome defined
				bam_index_t *idx;		
				int tid, beg, end, n_off;		
				string chr_str = samtools[i+1];
				pair64_t *off;
				//bamFile fp;
				cerr << samtools[i] << chr_str << endl;
				off = ReadBamChr_prep(chr_str, samtools[i], &tid, &beg, &end, in, &n_off);
				//if(off[0].u == -1 && off[0].v == -1){
				if(off==NULL){
				  cerr << "Error: can't load indexed bam for " << chr_str << endl;
				  return 0;
				}
				//bamFile fp = in->x.bam;
				uint64_t curr_off;
				int ii, ret, n_seeks;
				n_seeks = 0; ii = -1; curr_off = 0;
				while(int while_index = ReadBamChr(b, fp, tid, beg, end, &curr_off, &ii, &n_seeks, off, n_off)){
					if(while_index == -1)
						continue;
					//char *mtid;
					//mtid = in->header->target_name[b->core.tid];
					if(b->core.tid < 0 || b->core.tid != tid || b->core.pos >= end)
						continue;
					// define read
					read_data read;
					
					if(uint8_t *tmp = bam_aux_get(b, "NM"))
						read.NM = bam_aux2i(tmp);// int32_t
					if(uint8_t *tmp = bam_aux_get(b, "MF"))
						read.MF = bam_aux2i(tmp);
					int NBaseMapped = 0;
					uint32_t *cigar = bam1_cigar(b);
					for(int k = 0; k < b->core.n_cigar; ++k){
						int op = cigar[k] & BAM_CIGAR_MASK;
						if(op == BAM_CMATCH)
							NBaseMapped += cigar[k] >> BAM_CIGAR_SHIFT;
					}
					read.name = bam1_qname(b);
					read.flag = b->core.flag;
					read.mqual = b->core.qual;
					read.seq = get_string(bam1_seq(b), b->core.l_qseq);
					int readlen = read.seq.length();
					read.PercMapped = readlen > 0 ? NBaseMapped/readlen : 0;
					
					seqlen += readlen;
                                        if(regionsize != 0 && seqlen/regionsize > high_depth_skip*2){
                                            cerr << "High depth, average seq coverage: " << seqlen/regionsize << ", skip..\n";
                                            samclose(in);
                                            return 0;
                                        }
					if(read.name.find("/1") != string::npos)
						read.name.replace(read.name.find("/1"),2,"");
					if(read.name.find("/2") != string::npos)
						read.name.replace(read.name.find("/2"),2,"");
					if(max_mapqual.find(read.name) == max_mapqual.end() || max_mapqual[read.name] < read.mqual)
						max_mapqual[read.name] = read.mqual;
					if(min_NM.find(read.name) == min_NM.end() || min_NM[read.name]>read.NM)
						min_NM[read.name] = read.NM;
					reads.push_back(read);
				}
				samclose(in);
				for(int i = 0; i < reads.size(); i++){
					read_data read = reads[i];
					if(read.flag & 0x0200 || read.flag & 0x0400)
						continue;
					if(max_mapqual.find(read.name) == max_mapqual.end() || max_mapqual[read.name] < assemble_read_qual)
						continue;
					if(max_mapqual[read.name] >= 40 
					   && read.flag & 0x0001 
					   && read.flag & 0x0004 // one end unmapped reads
					   || read.NM > num_mismatch_poor_map && min_NM[read.name] < 3 // one-end poorly mapped
					   || read.MF == 130 // Maq Smith-Waterman aligned
					   || read.PercMapped < 0.9 // more than 10% of bases are not mapped
					   ){
						read.name += ",SV";
						nSVreads ++;
					}
					if(uniqreadnames.find(read.name + itos(read.flag)) != uniqreadnames.end())
						continue;
					string readname_tmp = ">" + read.name;
					buffer.push_back(readname_tmp);
					buffer.push_back(read.seq);
					uniqreadnames[read.name + itos(read.flag)] ++;
					nreads++;
				}
		
				free(off);
				bam_destroy1(b);
			}
			if(nreads <= 0){
				cerr << "no qualified reads from bam, skip ... \n";
				return 0;
			}
			float avgseqlen = float(seqlen)/float(nreads);

			if(write_to_read){
				ofstream OUT;
				string tmp = datadir + "/" + prefix + /*".a" + a_str + ".b" + b_str +*/ ".fa";
				OUT.open(tmp.c_str());
				/*while (buffer.size() != 0) {
					OUT << ">" << buffer.front() << endl;
					buffer.pop_front();
					string sequence = buffer.front();
					buffer.pop_front();
					OUT << sequence << endl;
					// skip
				}*/

				// change buffer back to vector
				for(int bi = 0; bi < buffer.size(); bi += 2){
					OUT << buffer[bi] << endl;
					if(bi + 1 < buffer.size())
						OUT << buffer[bi+1] << endl;
				}
			
				// skip -h
				// skip -D
				OUT.close();
			}
			if(regionsize <= 0){
				cerr << "region size <= 0, skip ... \n";
				return 0;
			}
			regionsize += int(2*avgseqlen);
			avgdepth = regionsize > 0 ? (float)seqlen/(float)regionsize: 0;
			if(avgdepth <= 0 || avgdepth > high_depth_skip){ // skip high depth region
				if(avgdepth <= 0)
					cerr << "no coverage, skip ...\n";
				if(avgdepth > high_depth_skip)
				  cerr << "coverage: " << avgdepth << " > " << high_depth_skip << ", skip ... \n";
				return 0;
			}
		}
		cerr << "#Reads:"<< nreads << "\t#SVReads:" << nSVreads << "\tRegionSize:" << regionsize << "\tAvgCoverage:" << avgdepth << endl;
		// Assemble
		tigra TG;
		TG.alternative_haplotype = datadir + "/" + prefix + /*".a" + a_str + ".b" + b_str +*/  "." + itos(regionsize) + "." + pos_1 + "." + pos_2 + ".fa.contigs.het.fa";
		TG.assembly_file = datadir + "/" + prefix + /*".a" + a_str + ".b" + b_str +*/ "." + itos(regionsize)  + "." + pos_1 + "." + pos_2 + ".fa.contigs.fa";
		//TG.reference_for_screen = datadir + "/" + prefix + ".ref.fa";
		TG.reference_for_screen = ref_string;
		TG.estimate_STR = "SV";
        TG.kmer_size = kmers;
                TG.max_node = max_node;
                float tmp1 = avgdepth/100.0;
                tmp1 = (tmp1 < 2) ? 2:tmp1;
		TG.low_kmer = tmp1;
		//vector<string> fastas;
		//fastas.push_back(datadir + "/" + prefix + /*".a" + a_str + ".b" + b_str +*/ ".fa");
		// now to read from buffer directly;

		if(type.compare("ITX") == 0 || type.compare("INS") == 0){
		  TG.min_degree=1;
		}
		//TG.run_tigra(fastas);
		TG.run_tigra(buffer);
		// done with tigra
	}	
	
	int write_to_bams(BD_data SV){
		int ret = 0;
		if(SV.bam_related.length() != 0){
			tools tl;
			vector<string> splitted;
			tl.split(SV.bam_related, ",", splitted);
			for(int i = 0; i < splitted.size(); i++){
				string bam_first_name = splitted[i];
				if(bamsmap.find(bam_first_name) != bamsmap.end()){
					string bam = bamsmap[bam_first_name];
					bams.push_back(bam);
				}
			}
			if(bams.size() != 0)
				ret = 1;
			else 
				cerr << "Didn't find the bam in the bam list file.\n";
			
		}
		else 
			cerr << "Cannot find bam related column or string!\n";
		
		return ret;
	}
	
	void assemble_func(){//string BreakDancer_file, vector<string> bams){
		struct stat stFileInfo; 
		int intStat = stat(reference_file.c_str(),&stFileInfo); 
		if(intStat != 0){
			cerr << "Cannot find the reference file: " << reference_file << ", please check if it exists. " << endl;
			reference_file = "";
		}
		vector<BD_data> SVs;
		if(BreakDancer_file.length() != 0)
			ReadBDCoor(BreakDancer_file, SVs); // skip library_to_skip_ at the moment
		else // support pcr format
			ReadPCRCoor(PCR_file, SVs);
		cout << "#" << SVs.size() << " SVs to be assembled from" << endl;
		cout << "#Bams: ";
		if(bams.size() != 0){
			cout << bams[0];
			for(int i = 1; i < bams.size(); i++)
				cout << ", " << bams[i];
		}
		else{
			for(map<string, string>::iterator bamsmap_it = bamsmap.begin(); bamsmap_it != bamsmap.end(); bamsmap_it++)
				cout << (*bamsmap_it).first << ":" << (*bamsmap_it).second << "\n";
		}
		cout << endl;
		
		cout << "#Version-" << version;
		cout << "\tParameters: \n";
		// foreach my $opt(keys %opts){printf "\t%s",join(':',$opt,$opts{$opt});} print "\n"; skip
		//if($opts{k}){ skip
		//cout << "#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tAllele_frequency\tVersion\tRun_Param\tAsmChr1\tAsmStart1\tAsmChr2\tAsmStart2\tAsmOri\tAsmSize\tAsmHet\tAsmScore\tAlnScore\twAsmScore\n";
                cout << "#CHR1\tPOS1\tCHR2\tPOS2\tORI\tSIZE\tTYPE\tHET\twASMSCORE\tTRIMMED_CONTIG_SIZE\tALIGNED\%\tNUM_SEG\tNUM_FSUB\tNUM_FINDEL\tBP_FINDEL\tMicroHomology\tMicroInsertion\tPREFIX\tASMPARM\tCopyNumber\tGene\tKnown\n";
		
		//else{
		//#printf "%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%d\t%s\ta%d.b%d\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$maxSV->{scarsize},$prefix,$maxSV->{a},$maxSV->{b};
		//			print "\#CHR1\tPOS1\tCHR2\tPOS2\tORI\tSIZE\tTYPE\tHET\twASMSCORE\tTRIMMED_CONTIG_SIZE\tALIGNED\%\tNUM_SEG\tNUM_FSUB\tNUM_FINDEL\tBP_FINDEL\tMicroHomology\tMicroInsertion\tPREFIX\tASMPARM\tCopyNumber\tGene\tKnown\n";
		//		} skip
		
		/*if($opts{v}){open(FN,">$opts{v}") || die "unable to open $opts{v}\n";}
		 if($opts{r}){open(ALNOUT,">$opts{r}") || die "unable to open $opts{r}\n";}
		 
		 my @as; my @bs;
		 @as=(defined $opts{a})?($opts{a}):(50,100,150);
		 @bs=(defined $opts{b})?($opts{b}):(50,100,150);
		 if($opts{H}){@as=(50); @bs=(100);$opts{p}=10000;}  #parameteres for capture data
		 if($opts{e}){my $stepsize=int($opts{l}/2); @as=($stepsize); @bs=($stepsize);} skip */
		int as = 50;
		int bs = 100;
		
		//srand(time ^ $$); skip
		string prefix;
		for(int i = 0; i < SVs.size(); i++){
			BD_data SV = SVs[i];
			string chr1 = SV.chr1;
			//if(chr1.find("chr") != string::npos)
			//	chr1.replace(chr1.find("chr"), 3, "");
			int start = SV.pos1;
			string chr2 = SV.chr2;
			//if(chr2.find("chr") != string::npos)
			//	chr2.replace(chr2.find("chr"), 3, "");
                        if(chromosome.length() > 0 && chr2.compare(chromosome) != 0)
                            continue;
			int end = SV.pos2;
			string type = SV.type;
                        //if(type.compare("INV") == 0){
			//	cerr << "Skip " << chr1 << "." << start << "." << chr2 << "." << end << "." << type << ", do not support INV assembly now. " << endl;
			//	continue;
			//}
			int size = SV.size;
			//$chr1=~s/chr//; $chr2=~s/chr//;skip
			//next unless ($start=~/^\d/ && $end=~/^\d/ && $size=~/^\d/);skip
			//string datadir; // skip data directory
			if(chr1.compare(chr2) == 0 && start > end){
				int tmp = start;
				start = end;
				end = tmp;
			}
			string start_end;
			/*if($opts{e}){
			 for(my $s=$start-$opts{e};$s<=$start+$opts{e};$s+=$opts{l}){
			 for(my $e=$end-$opts{e}; $e<=$end+$opts{e};$e+=$opts{l}){
			 push @start_ends,join(',',$s,$e);
			 }
			 }
			 }
			 else{
			 @start_ends=join(',',$start,$end);
			 }
			 skip*/
			
			/*foreach (@start_ends){
			 ($start,$end)=split /\,/; skip */
			vector<string> orientations;
            vector<string> oris;
			oris.push_back("+-");
			
			if(chr1.compare(chr2) != 0){
				oris.push_back("++");
				oris.push_back("--");
				oris.push_back("-+");
			}
			
			if(bamsmap.size() != 0){
				bams.clear();
				int ret = write_to_bams(SV); // write the related bams to the vector bams
				if(ret == 1){
					for(int j = 0; j < oris.size(); j++){
						AssembleBestSV(chr1, start, chr2, end, type, size, oris[j], SV, as, bs);
					}
				}
			}
			else {
                for(int j = 0; j < oris.size(); j++){
                    AssembleBestSV(chr1, start, chr2, end, type, size, oris[j], SV, as, bs);
				}
			}
		}
		cerr << "AllDone\n";
	}	
};

namespace {
	const int DEFAULT_FLANKING_SIZE = 500; // l
	const int DEFAULT_ASSEMBLE_READ_QUAL = 1; // q
	const int DEFAULT_NUM_MISMATCH_POOR_MAP = 5; // N
	const int DEFAULT_HIGH_DEPTH_SKIP = 1000; // p
	const int DEFAULT_PAD_LOCAL_REF = 200; // w
    const int DEFAULT_MIN_SIZE_THRESHOLD = 3; // M skip those with input size smaller than 3
    const int DEFAULT_MAX_NODE = 100; // h
    const string DEFAULT_KMERS = "15,25";

void usage() {
        fprintf(stderr, "\n./tigra_sv <SV file> <a.bam> <b.bam> ...\n\n");
        fprintf(stderr, "\n Or: ./tigra_sv <SV file> <bam_list_file>\n\nOptions: \n");
        //
        //fprintf(stderr, "    -A INT     Esimated maximal insert size [%d]\n", estimate_max_ins);
        //
        fprintf(stderr, "    -l INT     Flanking size for assembly [%d] bp\n", DEFAULT_FLANKING_SIZE);
        fprintf(stderr, "    -c STR     Only assemble calls on specified chromosome\n");
        fprintf(stderr, "    -R STR     Reference file location with the full path\n");
        fprintf(stderr, "    -q INT     Only assemble reads with mapping quality > [%d]\n", DEFAULT_ASSEMBLE_READ_QUAL);
        fprintf(stderr, "    -N INT     Number of mismatches required to be tagged as poorly mapped [%d]\n", DEFAULT_NUM_MISMATCH_POOR_MAP);
        fprintf(stderr, "    -p INT     Ignore cases that have average read depth greater than [%d]\n", DEFAULT_HIGH_DEPTH_SKIP);
        fprintf(stderr, "    -r         Write local reference to a file with .ref.fa as the suffix\n");
        fprintf(stderr, "    -d         Dump reads to fasta files\n");
        fprintf(stderr, "    -I STR     Save output files into an existing directory\n");
        fprintf(stderr, "    -w INT     Pad local reference by additional [%d] bp on both ends\n", DEFAULT_PAD_LOCAL_REF);
        fprintf(stderr, "    -b         Check when the input format is breakdancer\n");
        fprintf(stderr, "    -M INT     Skip those calls with input size smaller than [%d]\n", DEFAULT_MIN_SIZE_THRESHOLD);
        fprintf(stderr, "    -h INT     Maximum node to assemble, by default [%d]\n", DEFAULT_MAX_NODE);
        fprintf(stderr, "    -k STR     List of kmer sizes to use as a comma delimited string [%s]\n", DEFAULT_KMERS.c_str());
        //
        //fprintf(stderr, "    -Q INT     Minimal BreakDancer score required for analysis [%d]\n", qual_threshold);
        //fprintf(stderr, "    -L STRING  Ignore calls supported by libraries that contains (comma separated) STRING\n");
        //
        fprintf(stderr, "Version: %s (commit %s)\n", __g_prog_version, __g_commit_hash);
}
}

int main(int argc, char *argv[])
{
	int c;
	vector<string> bams;
	string BreakDancer_file;
	int estimate_max_ins = 500; // A
	int flanking_size = DEFAULT_FLANKING_SIZE; // l
	int assemble_read_qual = DEFAULT_ASSEMBLE_READ_QUAL; // q
	int num_mismatch_poor_map = DEFAULT_NUM_MISMATCH_POOR_MAP; // N
	int high_depth_skip = DEFAULT_HIGH_DEPTH_SKIP; // p
	int pad_local_ref = DEFAULT_PAD_LOCAL_REF; // w
    int min_size_threshold = DEFAULT_MIN_SIZE_THRESHOLD; // M skip those with input size smaller than 3
    int max_node = DEFAULT_MAX_NODE; // h
    string kmers = DEFAULT_KMERS;
	int qual_threshold = 0; // Q
	int breakdancer_format = 0; // b by default off
	string library_to_skip(""); // L
	string datadir = "."; // I mandatory
	int write_to_ref = 0; // r by default off
	int write_to_read = 0; // d by default off
    string chr = ""; // specify chromosome to parallelizing
    int skip_call = 0; // z skip running calls before this number
	string reference_file;
	
	while((c = getopt(argc, argv, "A:l:w:q:N:p:I:Q:L:brdR:c:M:h:k:")) >= 0){
		switch(c) {
			case 'A': estimate_max_ins = atoi(optarg); break;
			case 'l': flanking_size = atoi(optarg); break;
			case 'w': pad_local_ref = atoi(optarg); break;
			case 'q': assemble_read_qual = atoi(optarg); break;
			case 'N': num_mismatch_poor_map = atoi(optarg); break;
			case 'p': high_depth_skip = atoi(optarg); break;
			case 'I': datadir = optarg; break;
			case 'Q': qual_threshold = atoi(optarg); break;
			case 'L': library_to_skip = optarg; break;
			case 'b': breakdancer_format = 1; break;
			case 'r': write_to_ref = 1; break;
			case 'd': write_to_read = 1; break;
			case 'R': reference_file = optarg; break;
                        case 'c': chr = optarg; break;
                        case 'z': skip_call = atoi(optarg); break;
                        case 'M': min_size_threshold = atoi(optarg); break;          
                        case 'h': max_node = atoi(optarg); break;
                        case 'k': kmers = string(optarg); break;
			default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
				return 1;
		}
	}

    if (reference_file.empty()) {
        cerr << "No value specified for required argument -R! Abort\n\n";
        usage();
        return 1;
    }
	
	if(optind == argc){
        usage();
		return 1;
	}
	
	assemble AS;
	AS.estimate_max_ins = estimate_max_ins;
	AS.flanking_size = flanking_size;
	AS.pad_local_ref = pad_local_ref;
	AS.assemble_read_qual = assemble_read_qual;
	AS.num_mismatch_poor_map = num_mismatch_poor_map;
	AS.high_depth_skip = high_depth_skip;
	AS.datadir = datadir;
	AS.qual_threshold = qual_threshold;
	AS.library_to_skip_ = library_to_skip;
	AS.write_to_ref = write_to_ref;
	AS.write_to_read = write_to_read;
	AS.reference_file = reference_file;
        AS.chromosome = chr;
        AS.skip_call = skip_call;
        AS.min_size_threshold = min_size_threshold;
        AS.max_node = max_node;
        AS.kmers = kmers;
	
	if(breakdancer_format == 1)
		AS.BreakDancer_file = argv[optind];
	else
		AS.PCR_file = argv[optind];
	int optindex = optind + 1;
	string arg1(argv[optindex]);
	if(arg1.find(".bam") == string::npos){
		map<string, string> bams;
		tools tl;
		tl.analyze_BAM_list(bams, arg1);
		AS.bamsmap = bams;
	}
	else{
		while(optindex < argc){
			AS.bams.push_back(argv[optindex]);
			optindex ++;
		}
	}
	AS.assemble_func();
	return EXIT_SUCCESS;
}
	
/*int main(int argc, char *argv[])
{
	int c;
	string kmer_size("15");
	float min_kmer_cov = 3;
	int low_kmer = 2;
	int high_kmer = 2e9;
	string assembly_file;
	int tip_number = 1000;
	string alternative_haplotype;
	string estimate_STR;
	string reference_for_screen;
	int max_node = 100;
	int min_degree = 2;
	string graph_file;
	int debug = 0;
	cout << "PID of this process: " << getpid() << endl;
	while((c = getopt(argc, argv, "k:c:m:M:o:t:h:p:r:n:N:g:d")) >= 0){
		switch(c) {
			case 'k': kmer_size = strdup(optarg); break;
			case 'c': min_kmer_cov = atof(optarg); break;
			case 'm': low_kmer = atoi(optarg); break;
			case 'M': high_kmer = atoi(optarg); break;
			case 'o': assembly_file = strdup(optarg); break;
			case 't': tip_number = atoi(optarg); break;
			case 'h': alternative_haplotype = strdup(optarg); break;
			case 'p': estimate_STR = strdup(optarg); break;
			case 'r': reference_for_screen = strdup(optarg); break;
			case 'n': max_node = atoi(optarg); break;
			case 'N': min_degree = atoi(optarg); break;
			case 'g': graph_file = strdup(optarg); break;
			case 'd': debug = 1; break;
				
			default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
				return 1;
		}
	}
	if(optind == argc){
		fprintf(stderr, "\n");
		fprintf(stderr, "tigra <fasta.fa>\n\n");
		fprintf(stderr, "Options: \n");
		if(optind == argc){
			fprintf(stderr, "\n");
			//fprintf(stderr, "./tigra <fasta.fa>\n\n");
			//fprintf(stderr, "Options: \n");
			fprintf(stderr, "       -k STRING		Specify Kmer sizes, Use comma to delimit multiple kmers\n");
			fprintf(stderr, "	-c FLOAT		Minimal average kmer coverage [%f] for a contig to be considered in the alternative haplotypes\n", min_kmer_cov);
			fprintf(stderr, "       -m INT          Lowest Kmer frequency in initial hashing [%d]\n", low_kmer);		 
			fprintf(stderr, "       -M INT          Highest Kmer frequency in initial hashing [%d]\n", high_kmer);		
			fprintf(stderr, "       -o STRING       File to save primary assembly contigs\n");		 
			fprintf(stderr, "	-p STRING		Estimate the utility of the subset of reads containing STR in header\n");
			fprintf(stderr, "	-r STRING		Provide a reference to screen for non-reference kmers\n");
			fprintf(stderr, "       -t INT          Tip number [%d]\n", tip_number);	
			fprintf(stderr, "       -h STRING       File to save alternative haplotypes\n");	
			fprintf(stderr, "	-n INT			Maxmimal number of nodes allowed for constructing alternative haplotypes [%d]\n", max_node); 
			fprintf(stderr, "	-N INT			Minimal number of degrees for a node to consider as start of alternative haplotype [%d]\n", min_degree);
			fprintf(stderr,	"	-g STRING		Save assembly graph(png) in FILE\n");
			fprintf(stderr, "       -d INT          Turn on debug mode, generate intermediate files\n");
			fprintf(stderr, "Version: %s\n", version);
			fprintf(stderr, "\n");
			return 1;
		}
		fprintf(stderr, "\n");
		return 1;
	}
	
	// print options
	char options_[500];	
	sprintf(options_,"-c %f -m %d -M %d -t %d -n %d -N %d -d %d", min_kmer_cov, low_kmer, high_kmer, tip_number, max_node, min_degree, debug);
	options = options_;
	options += "-k" + kmer_size;
	options += "-o" + assembly_file;
	options += "-h" + alternative_haplotype;
	options += "-p" + estimate_STR;
	options += "-r" + reference_for_screen;
	options += "-g" + graph_file;
	
	map<int, int> contigtips;
	map<int, map<string, string> > contigs;	
	vector<string> fastas;
	
	int optindex = optind;
	while(optindex < argc){
		fastas.push_back(argv[optindex]);
		optindex ++;
	}
	
	main_functions mf(debug, low_kmer, high_kmer, contigs, contigtips, fastas, tip_number, estimate_STR);
	int fastas_ = 1;
	int contigs_ = 0;
	vector<string> Reads;
	string fasta_ = "";
	vector<string> pReads;
	mf.getReads(fastas_, contigs_, Reads, fasta_, pReads);
	mf.Reads = Reads;
	if(debug == 1){
		string fileName = "getReads";
		ofstream fh;
		fh.open(fileName.c_str());
		for(int i = 0; i < Reads.size(); i++)
			fh << Reads[i] << "\n";
		fh.close();
	}
	
	vector<string> kmers;
	string sep = ",";
	tools tl;
	tl.split(kmer_size, sep, kmers);
	string kmer_debug = kmers[0];
	mf.initial_iteration(atoi(kmers[0].c_str()));
	for(int i = 1; i < kmers.size(); i++){
		mf.contigtips = map<int, int>();
		int kmer = atoi(kmers[i].c_str());
		mf.iteration(kmer, mf.contigs, Reads);
	}
	
	/*string refseq = "";
	if(reference_for_screen.length() != 0){
		refseq = reference_for_screen;// read that in kmergen rather than here
	}*/
	
	/*int total_;
	map<string, int> phh;
	if(int(pReads.size()) - 1 >= 0){
		kmergen pkg;
		pkg.set_kmer_size(atoi(kmers[kmers.size() - 1].c_str()));
		pkg.c = low_kmer;
		pkg.C = high_kmer;
		//map<string, int> hh;
		pkg.set_filter(reference_for_screen); // need to go through filter in kmergen again
		total_ = pkg.doit(pReads, phh);
		if(debug == 1)
			pkg.printMer("MerP", phh);
	}
	
	// Output assembled contigs
	int switch_ = 4;
	mf.total = total_;
	mf.hh = phh;
	mf.kmer_size = atoi(kmers[kmers.size()-1].c_str());
	if(assembly_file.length() != 0){
		ofstream ofs(assembly_file.c_str());
		mf.outputcontigs(ofs, switch_, 1);
		ofs.close();
	}
	else {
		mf.outputcontigs(std::cout, switch_, 1);
	}
	
	
	// do the following when allpaths is ready
	if (alternative_haplotype.length() > 0 || graph_file.length() > 0) {
		allpaths ap;
		ap.set_kmer_size(atoi(kmers[kmers.size()-1].c_str()));
		ap.set_cov(min_kmer_cov);
		ap.set_n(max_node);
		ap.set_rcontig(mf.to_print_contigs);
		
		ap.set_degree(min_degree);
		ap.set_file_name(argv[optind]);
		ap.set_graph();
		map<int, map<string, string> > ret_contigs;
		ap.doit(mf.to_print_contigs, ret_contigs);
		mf.to_print_contigs = ret_contigs;
		int switch_ = 5;
		if(alternative_haplotype.length() != 0){
			ofstream ofs(alternative_haplotype.c_str());
			mf.outputcontigs(ofs, switch_, 5);
			ofs.close();
		}
		else {
			mf.outputcontigs(std::cout, switch_, 5);
		}

	 }
	
}		*/
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
