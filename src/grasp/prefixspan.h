#ifndef PREFIXSPAN_H_
#define PREFIXSPAN_H_

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "tree.h"
#include "logreg.h"

#include <cmath>

using std::set;
using std::map;
using std::pair;
using std::vector;
using std::string;
using std::multimap;
using std::ostringstream;
using std::unordered_map;

typedef map< int, vector<int> > occurence;

struct OutObject {
public:
  string substr;
  string subgraph;
  set<int> instances;
  Vector pos, neg;
  string notes;
  int sup;
  int hit_g, hit_s;
  int num_g, num_s;
};

class Prefixspan {
 private:
  map<char,occurence> scan_prefixspan(occurence&);
  unsigned int support(occurence&);
 public:
  bool verbose;
  bool classify;
  string current_pref;
  multimap<double,OutObject> ranking;
  unsigned maxout;
  unsigned num_out;
  int num_links;
  int neg_num_links;
  unordered_map<int, unsigned> degree;
  unordered_map<int, unsigned> neg_degree;
  unordered_map<string,int> mapper;
  unsigned int minsup;
  unsigned int delta;
  vector<string> sdata;
  template<typename T> void read(T& is);
  void traverse(string&, occurence&, unordered_map<int,vector<Destpair> >&, tree<vertex>&, int, unsigned);
  void prefixspan(string& pref, occurence& oc, unordered_map<int,vector<Destpair> >&, Gspan&);
  void co_mine(unordered_map<int,vector<Destpair> >&, Gspan&);
};

template<typename T> void Prefixspan::read(T& is){
  string line, name;
  string seq = "";
  bool noread = true;
  int file = 0;
  while (getline(is, line)) {
    string str = "";
    std::istringstream stream(line);
    if(line[0] == '>'){
      if (noread){
	noread = false;
      }else{
	//std::cout << seq << std::endl;
	sdata.push_back(seq);
	seq = "";
      }
      stream >> name;
      name.erase(0,1);
      //std::cout << name << "->" << file << std::endl;
      mapper[name] = file;
      file += 1;
      continue;
    }
    stream >> str;
    seq += str;
  }
  //std::cout << seq << std::endl;
  sdata.push_back(seq);
}

map<char,occurence> Prefixspan::scan_prefixspan(occurence& oc){
  map<char,occurence> m;
  for(occurence::iterator p=oc.begin(); p != oc.end(); ++p){
    int id = p->first;
    for(vector<int>::iterator q=p->second.begin(); q != p->second.end(); ++q){
      //hash_set<int> checked;
      unsigned int ulimit = sdata[id].length();
      if(*q+delta+1<ulimit) ulimit = *q+delta+1;
      for(unsigned int pos=*q+1; pos<ulimit; ++pos){
	char key = sdata[id][pos];
	//if(checked.find(key)==checked.end()){
	  m[key][id].push_back(pos);
	  //}
	  //checked.insert(key);
      }
    }
  }
  return m;
}

unsigned int Prefixspan::support(occurence& oc){
  int sup = 0;
  for(occurence::iterator x=oc.begin(); x != oc.end(); ++x){
    sup += degree[x->first];
  }
  return sup;
}

inline void Prefixspan::traverse(string& pref, occurence& oc, 
				 unordered_map<int,vector<Destpair> >& pdata,
				 tree<vertex>& stree, int n_id, unsigned gdata_size){
  vertex& gnode = stree.node[n_id];
  unsigned int psup=0, nsup=0;
  unsigned int psup_neg=0, nsup_neg=0;

  OutObject encode;
  encode.substr = pref;
  encode.subgraph = gnode.name;

  if(classify){
    for(occurence::iterator x=oc.begin(); x != oc.end(); ++x){
      for(vector<Destpair>::iterator y=pdata[x->first].begin(); y!=pdata[x->first].end(); ++y){
	if(gnode.exist.find(y->graph) != gnode.exist.end()){
	  if(y->is_positive){
	    psup++;
	  }else{
	    psup_neg++;
	  }
	  encode.instances.insert(y->id);
	}else{
	  if(y->is_positive){
	    nsup++;
	  }else{
	    nsup_neg++;
	  }
	}
      }
    }
  }else{
    for(occurence::iterator x=oc.begin(); x != oc.end(); ++x){
      for(vector<Destpair>::iterator y=pdata[x->first].begin(); y!=pdata[x->first].end(); ++y){
	if(gnode.exist.find(y->graph) != gnode.exist.end()){
	  psup++;
	  encode.instances.insert(y->id);
	}else{
	  nsup++;
	}
      }
    }
  }

  if(psup < minsup) return;

  num_out++;

  int gYes_sNo, gNo_sYes, gYes_sYes, gNo_sNo;

  encode.pos.resize(4);
  encode.neg.resize(4);

  //Vector pos(4), neg(4);

  if (classify){
    // g:No  s:Yes
    encode.pos(2) = nsup;
    encode.neg(2) = nsup_neg; 

    // g:Yes s:Yes
    encode.pos(3) = psup;
    encode.neg(3) = psup_neg; 
    
    // g:Yes s:No
    encode.pos(1) = gnode.support - encode.pos(3); 
    encode.neg(1) = gnode.neg_support - encode.neg(3); 

    // g:No s:No
    encode.pos(0) = num_links - encode.pos(1) - encode.pos(2) - encode.pos(3);
    encode.neg(0) = neg_num_links - encode.neg(1) - encode.neg(2) - encode.neg(3);

  }else{
    int pos_g, pos_s, neg_g, neg_s;

    pos_g = gnode.num;
    pos_s = oc.size();
    neg_g = gdata_size - pos_g;
    neg_s = sdata.size() - pos_s;

    gYes_sNo  = gnode.support - psup;
    gNo_sYes  = nsup;
    gYes_sYes = psup;
    gNo_sNo   = num_links - gYes_sNo - gNo_sYes - gYes_sYes;

    encode.hit_g = pos_g;
    encode.hit_s = pos_s;
    encode.num_g = gdata_size;
    encode.num_s = sdata.size();

    encode.pos(0) = gNo_sNo;   encode.neg(0) = neg_g*neg_s-gNo_sNo;
    encode.pos(1) = gYes_sNo;  encode.neg(1) = pos_g*neg_s-gYes_sNo;
    encode.pos(2) = gNo_sYes;  encode.neg(2) = neg_g*pos_s-gNo_sYes;
    encode.pos(3) = gYes_sYes; encode.neg(3) = pos_g*pos_s-gYes_sYes;
  }

  double score = log10_p(encode.pos,encode.neg);
  
  if(score < 0){
    ostringstream memo;
    memo << "(g)" << gnode.num << "/" << gdata_size << " (s)" << oc.size()
	 << "/" << sdata.size();
    encode.notes = memo.str();
    encode.sup = psup;
    ranking.insert(std::pair<double,OutObject>(score,encode));
    if(maxout != 0 && ranking.size()>maxout){
      multimap<double,OutObject>::iterator rit = ranking.end();
      ranking.erase(--rit);
    }

    if (verbose){
      std::cout << "log10(p)=" << score << " str=" << pref << "[" << oc.size() <<
	"] dfs=" << gnode.name << "[" << gnode.num << "]" << std::endl ;
      std::cout << "    No  No => pos:" << encode.pos(0) << " neg:" << encode.neg(0) << std::endl;
      std::cout << "   Yes  No => pos:" << encode.pos(1) << " neg:" << encode.neg(1) << std::endl;
      std::cout << "    No Yes => pos:" << encode.pos(2) << " neg:" << encode.neg(2) << std::endl;
      std::cout << "   Yes Yes => pos:" << encode.pos(3) << " neg:" << encode.neg(3) << std::endl;
    }
  }

  for(list<int>::iterator it=stree.children[n_id].begin();
      it != stree.children[n_id].end(); ++it){
    traverse(pref,oc,pdata,stree,*it,gdata_size);
  }
}

void Prefixspan::prefixspan(string& pref, occurence& oc, 
			    unordered_map<int,vector<Destpair> >& pdata, Gspan& gspan){
  if (support(oc) < minsup) return;
  if (current_pref != pref){
    std::cerr << "\r" << pref << "  ";
    current_pref = pref;
  }

  tree<vertex>& stree = gspan.search_tree;
  for(list<int>::iterator it=stree.children[stree.root].begin();
      it != stree.children[stree.root].end(); ++it){
    traverse(pref,oc,pdata,stree,*it,gspan.gdata.size());
  }
	
  map<char,occurence> extend = scan_prefixspan(oc);
	
  for(map<char,occurence>::iterator it=extend.begin(); it != extend.end(); ++it){
    string s = pref+it->first;
    prefixspan(s,it->second,pdata,gspan);
  }
}


void Prefixspan::co_mine(unordered_map<int,vector<Destpair> >& pdata, Gspan& gspan){
  map<char,occurence> m;
  current_pref = "";

  for(unsigned int id=0; id<sdata.size(); ++id){
    for(unsigned int pos=0; pos<sdata[id].length(); ++pos){
      char key = sdata[id][pos];
      m[key][id].push_back(pos);
    }
  }
  for(map<char,occurence>::iterator it=m.begin(); it != m.end(); ++it){
    string prefix(1,it->first);
    prefixspan(prefix,it->second,pdata,gspan);
  }
}

#endif /*PREFIXSPAN_H_*/
