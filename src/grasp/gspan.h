#ifndef GSPAN_H_
#define GSPAN_H_


#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include "tree.h"

using std::set;
using std::map;
using std::string;
using std::vector;
using std::greater;
using std::unordered_map;

struct Triplet {
public:
  int x, y, z;
  Triplet reverse();
  explicit Triplet(){}
  explicit Triplet(int _x,int _y,int _z): x(_x),y(_y),z(_z){}
};

struct Pair { public: int a, b; void set(int, int); };
inline void Pair::set(int _a, int _b){ a = _a; b = _b; }

struct VertexPair: public Pair { public: int id; };
struct Edge   { public: int to, id; Triplet labels; };
struct DFSCode { public: Triplet labels; Pair time; };

typedef vector< vector<Edge> > AdjacentList;

class Graph: public AdjacentList {
 public:
  vector<int> label;
  int num_of_edges;
  void resize(unsigned s){ AdjacentList::resize(s); label.resize(s); }
};

class EdgeTracer {
 public:
  VertexPair  vpair;
  EdgeTracer* predec;
  explicit EdgeTracer(){};
  void set(int,int,int,EdgeTracer*);
};

inline void EdgeTracer::set(int a,int b,int id, EdgeTracer* pr){
  vpair.a = a;
  vpair.b = b;
  vpair.id = id;
  predec = pr;
}

typedef vector<EdgeTracer> Tracers;
typedef map<int,Tracers> GraphToTracers;
typedef map<Pair,GraphToTracers> PairSorter;

inline Triplet Triplet::reverse(){ return Triplet(z,y,x); }

inline bool operator< (const Triplet& left, const Triplet& right){
  if (left.x!=-1 && right.x!=-1 && left.x != right.x) return (left.x < right.x);
  if (left.y!=-1 && right.y!=-1 && left.y != right.y) return (left.y < right.y);
  return (left.z < right.z);
}

inline bool operator<= (const Triplet& left, const Triplet& right){
  return !(right < left);
}

inline bool operator< (const Pair& left, const Pair& right){
  if (left.a != right.a) return (left.a < right.a);
  return (left.b < right.b);
}

inline bool operator== (const DFSCode& left, const DFSCode& right){
  if(left.time.a != right.time.a) return false;
  if(left.time.b != right.time.b) return false;	
  if(left.labels.x != right.labels.x) return false;
  if(left.labels.y != right.labels.y) return false;
  return (left.labels.z == right.labels.z);	
}

inline bool operator!= (const DFSCode& x, const DFSCode& y){
  return !(x==y);
}

inline std::ostream& operator<< (std::ostream& os, const vector<DFSCode> pattern){
  if(pattern.empty()) return os;
  os << "(" << pattern[0].labels.x << ") " << pattern[0].labels.y << " (0f" << pattern[0].labels.z << ")";
  for(unsigned i=1; i<pattern.size(); ++i){
    if(pattern[i].time.a < pattern[i].time.b){
      os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
    }else{
      os << " " << pattern[i].labels.y << " (b" << pattern[i].time.b << ")";
    }
  }
  return os;
}

const vector<Graph> readGraphs(std::istream&);

class Gspan {
 private:
  bool is_min();
  Graph toGraph(vector<DFSCode>&);
  unsigned support(GraphToTracers&);
  void scan_rm(vector<DFSCode>&, vector<int>&);
  bool min_checker(vector<DFSCode>&, Graph&, Tracers&);
  int  scan_gspan(GraphToTracers&, PairSorter&, map<int,PairSorter,greater<int> >&);
 public:
  bool classify;
  tree<vertex> search_tree;
  unordered_map<int, unsigned> degree;
  unordered_map<int, unsigned> pos_degree;
  unordered_map<int, unsigned> neg_degree;
  unordered_map<string,int> mapper;
  unsigned minsup;
  unsigned minpat;
  vector<Graph>   gdata;
  vector<DFSCode> pattern;
  template<typename T>
    void read(T& ins, unordered_map<string,int>& node_dict);
  void edge_grow(GraphToTracers&, int);
  int report(GraphToTracers&, int);
  void build_tree();
};

struct Destpair { public: int graph; bool is_positive; int id; };

struct Spair { public: string key; unsigned value; };
inline bool operator< (const Spair& left, const Spair& right){
  if (left.value != right.value) return (left.value > right.value);
  return (left.key < right.key);
}

template<typename T>
void make_dict(T& ins, unsigned minsup, unordered_map<string,int>& node_dict,string filename){
  std::string line;
  int i = 0;
  int num_v = 0;
  int num_e = 0;
  map<string,int> nodes;
  int file = 0;

  while ( getline(ins,line) ){
    if (3 < i && i <= 3 + num_v) {
      //std::cout << i << ": " << line << std::endl;
      nodes[line.substr(31,3)] += 1;
    }else if (i == 3){
      if(line.substr(33,6) != " V2000"){
	std::cerr << "SDF file error: Not V2000 Format." << std::endl;
	std::cerr << "line" << i << ":" << line << std::endl;
	return;
      }
      num_v = atoi(line.substr(0,3).c_str());
      num_e = atoi(line.substr(3,3).c_str());
    }
    //    else if (3+num_v < i && i <= 3 + num_v + num_e){
    //      edges[line.substr(6,3)] += 1;
//       std::cout << "Line " << i << " : " << line.substr(0,3)
// 		<< "|" << line.substr(3,3) 
// 		<< "|" << line.substr(6,3) << std::endl;
//    }
    //std::cout << "Line " << i << " : " << line << std::endl;
    i += 1;
    //if (i > 10) break;
    //if (line.find("$$$$") < line.length()) break;
    if (line.find("$$$$") < line.length()){
      file += 1;
      i = 0;
    }
  }

  std::cout << ">" << file << " SDF files read." << std::endl;

  std::cout << ">relabel the nodes in descending frequency." << std::endl;

  Spair p;
  set<Spair> out;
  for(map<string,int>::iterator it = nodes.begin();
      it != nodes.end(); ++it){
    p.key   = it->first;
    p.value = it->second;
    out.insert(p);
  }

  std::cout << ">remove infrequent edges (and hydrogen edges)." << std::endl;

  std::ofstream ofs(filename.c_str());
  int id = 0;
  for(set<Spair>::iterator it = out.begin();
      it != out.end(); ++it){
    if (it->value > minsup && it->key != "H  ") {
      //if (it->value > minsup) {
      ofs << "v " << id << " " << it->key << std::endl;
      //ofs << " : " << it->value << std::endl;
      node_dict[it->key] = id;
      id += 1;
    }else{
      node_dict[it->key] = -1;
    }
  }

  ofs << "e 0 1" << std::endl;
  ofs << "e 1 2" << std::endl;
  ofs << "e 2 3" << std::endl;
  ofs << "e 3 4" << std::endl;

  // edge out
  //out.clear();
  //for(map<string,int>::iterator it = edges.begin();
  //    it != edges.end(); ++it){
  //  p.key   = it->first;
  //  p.value = it->second;
  //  out.insert(p);
  //}
  //id = 0;
  //for(set<Spair>::iterator it = out.begin();
  //    it != out.end(); ++it){
  //  ofs << "e " << id << " " << it->key << std::endl;
  //  id += 1;
  //}

  ofs.close();
}

template<typename T>
void Gspan::read(T& ins, unordered_map<string,int>& node_dict){
  Graph g;
  Triplet labels;
  Edge edge;

  int i = 0;
  int eid = 0;
  int file = 0;
  int num_v = 0;
  int num_e = 0;
  string line, name;

  while ( getline(ins,line) ){
    std::stringstream stream(line);
    if (3 < i && i <= 3 + num_v) {
      g.resize(i-4+1);
      g.label[i-4] = node_dict[line.substr(31,3)];
      //std::cout << "v " << i-4 << " " << node_dict[line.substr(31,3)] << std::endl;  
    }else if (3 + num_v < i && i <= 3 + num_v + num_e){
      int atom1  = atoi(line.substr(0,3).c_str())-1;
      int atom2  = atoi(line.substr(3,3).c_str())-1;
      int btype  = atoi(line.substr(6,3).c_str())-1;
      if (g.label[atom1] != -1 && g.label[atom2] != -1){
	labels.x = g.label[atom1];
	labels.y = btype;
	labels.z = g.label[atom2];
	edge.to = atom2; edge.labels = labels; edge.id = eid++;
	g[atom1].push_back(edge);
	edge.to = atom1; edge.labels = labels.reverse(); edge.id = eid++;
	g[atom2].push_back(edge);
      }
      //std::cout << "e " << atom1 << " " << atom2 << " " << btype << std::endl;
    }else if (i == 3){
      num_v = atoi(line.substr(0,3).c_str());
      num_e = atoi(line.substr(3,3).c_str());
    }else if (i == 0){
      g.clear();
      eid = 0;
      stream >> name;
      //std::cout << name << "->" << file << std::endl;
      mapper[name] = file;
      //std::cout << "t # " << file << std::endl;
    }
    i += 1;
    if (line.find("$$$$") < line.length()){
      file += 1;
      g.num_of_edges = eid;
      gdata.push_back(g);
      i = 0;
      //std::cout << std::endl;
    }
  }
}

#endif /*GSPAN_H_*/
