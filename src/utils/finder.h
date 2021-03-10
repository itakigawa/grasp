#ifndef FINDER_H_
#define FINDER_H_

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

#include "patricia.h"

using std::greater;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::unordered_map;
using std::vector;

struct IntPair
{
public:
  int g;
  int s;
  int rank;
};

inline bool operator<(const IntPair &left, const IntPair &right)
{
  return (left.rank < right.rank);
}

struct Counterpart
{
public:
  int target;
  int rank;
};

inline bool operator<(const Counterpart &left, const Counterpart &right)
{
  return (left.rank < right.rank);
}

struct Triplet
{
public:
  int x, y, z;
  Triplet reverse();
  explicit Triplet() {}
  explicit Triplet(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}
};

struct Pair
{
public:
  int a, b;
  void set(int, int);
};
inline void Pair::set(int _a, int _b)
{
  a = _a;
  b = _b;
}

struct VertexPair : public Pair
{
public:
  int id;
};
struct Edge
{
public:
  int to, id;
  Triplet labels;
};
struct DFSCode
{
public:
  Triplet labels;
  Pair time;
};

typedef vector<vector<Edge> > AdjacentList;

class Graph : public AdjacentList
{
public:
  vector<int> label;
  int num_of_edges;
  void resize(unsigned s)
  {
    AdjacentList::resize(s);
    label.resize(s);
  }
};

class EdgeTracer
{
public:
  VertexPair vpair;
  EdgeTracer *predec;
  explicit EdgeTracer(){};
  void set(int, int, int, EdgeTracer *);
};

inline void EdgeTracer::set(int a, int b, int id, EdgeTracer *pr)
{
  vpair.a = a;
  vpair.b = b;
  vpair.id = id;
  predec = pr;
}

typedef vector<EdgeTracer> Tracers;
typedef map<int, Tracers> GraphToTracers;
typedef map<Pair, GraphToTracers> PairSorter;

inline Triplet Triplet::reverse() { return Triplet(z, y, x); }

inline bool operator<(const Triplet &left, const Triplet &right)
{
  if (left.x != -1 && right.x != -1 && left.x != right.x)
    return (left.x < right.x);
  if (left.y != -1 && right.y != -1 && left.y != right.y)
    return (left.y < right.y);
  return (left.z < right.z);
}

inline bool operator<=(const Triplet &left, const Triplet &right)
{
  return !(right < left);
}

inline bool operator<(const Pair &left, const Pair &right)
{
  if (left.a != right.a)
    return (left.a < right.a);
  return (left.b < right.b);
}

inline bool operator==(const DFSCode &left, const DFSCode &right)
{
  if (left.time.a != right.time.a)
    return false;
  if (left.time.b != right.time.b)
    return false;
  if (left.labels.x != right.labels.x)
    return false;
  if (left.labels.y != right.labels.y)
    return false;
  return (left.labels.z == right.labels.z);
}

inline bool operator!=(const DFSCode &x, const DFSCode &y)
{
  return !(x == y);
}

inline std::ostream &operator<<(std::ostream &os, const vector<DFSCode> pattern)
{
  if (pattern.empty())
    return os;
  os << "(" << pattern[0].labels.x << ") " << pattern[0].labels.y << " (0f" << pattern[0].labels.z << ")";
  for (unsigned i = 1; i < pattern.size(); ++i)
  {
    if (pattern[i].time.a < pattern[i].time.b)
    {
      os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
    }
    else
    {
      os << " " << pattern[i].labels.y << " (b" << pattern[i].time.b << ")";
    }
  }
  return os;
}

typedef map<int, vector<int> > occurence;

struct datapack
{
  int cmp;
  int seq;
  int label;
};

class Finder
{
public:
  unordered_map<int, string> g_tname;
  unordered_map<int, string> s_tname;
  unordered_map<string, int> g_tid;
  unordered_map<string, int> s_tid;
  map<int, string> fg_name;
  map<int, string> fs_name;
  unordered_map<string, int> fg_id;
  unordered_map<string, int> fs_id;
  vector<Graph> gdata;
  vector<string> sdata;
  map<int, map<int, int> > links;
  vector<datapack> trans;
  map<Triplet, GraphToTracers> gheap;
  unsigned num_of_feat;
  patricia gpat_tree;
  patricia spat_tree;
  // member functions
  template <typename T>
  void read_sdf(T &ins, unordered_map<string, int> &node_dict);
  template <typename T>
  void read_fasta(T &is);
  template <typename T>
  bool read_pairs(T &ins);
  template <typename T>
  void read_features(T &ins);
};

template <typename T>
void Finder::read_features(T &ins)
{
  string line;
  vector<string> v, vg, vs;
  int gcount = 1, scount = 1;
  int rank = 1;
  while (getline(ins, line))
  {
    v.clear();
    vs.clear();
    vg.clear();
    boost::algorithm::split(v, line, boost::algorithm::is_any_of("|"));
    //std::cout << "feat" << count << " " << v[2] << " vs "  << v[3] << std::endl;
    if (fs_id.find(v[2]) == fs_id.end())
    {
      fs_id[v[2]] = scount;
      fs_name[scount] = v[2];
      stokenize(vs, v[2]);
      //std::cout << "   --->" << boost::algorithm::join(vs,"-") << std::endl;
      spat_tree.insert(vs, patricia::root, scount);
      scount += 1;
    }
    if (fg_id.find(v[3]) == fg_id.end())
    {
      fg_id[v[3]] = gcount;
      fg_name[gcount] = v[3];
      gtokenize(vg, v[3]);
      //std::cout << "   --->" << boost::algorithm::join(vg,"-") << std::endl;
      gpat_tree.insert(vg, patricia::root, gcount);
      gcount += 1;
    }
    links[fs_id[v[2]]][fg_id[v[3]]] = rank;
    //trans.push_back(pair<int,int>(fs_id[v[2]],fg_id[v[3]]));
    rank += 1;
  }
  std::cout << "# features: " << rank - 1 << std::endl;
  num_of_feat = rank - 1;
}

template <typename T>
bool Finder::read_pairs(T &ins)
{
  string line, cmp, seq;
  char label;

  bool first = true;
  bool classify = false;
  int count = 1;
  while (getline(ins, line))
  {
    datapack tmp;
    std::stringstream stream(line);
    stream >> seq >> cmp;
    tmp.seq = s_tid[seq];
    tmp.cmp = g_tid[cmp];
    if (first)
    {
      if (stream >> label)
      {
        classify = true;
        if (label == '+')
        {
          tmp.label = 1;
        }
        else
        {
          tmp.label = -1;
        }
      }
      first = false;
    }
    else if (classify)
    {
      if (!(stream >> label))
      {
        std::cerr << "Wrong pair file format. (with label or not?)" << std::endl;
      }
      if (label == '+')
      {
        tmp.label = 1;
      }
      else
      {
        tmp.label = -1;
      }
    }
    //trans.push_back(pair<int,int>(s_tid[seq],g_tid[cmp]));
    trans.push_back(tmp);
    count++;
  }
  std::cout << "# pairs: " << count - 1 << std::endl;

  return classify;
}

template <typename T>
void Finder::read_sdf(T &ins, unordered_map<string, int> &node_dict)
{
  Graph g;
  Triplet labels;
  Edge edge;

  int i = 0;
  int eid = 0;
  int file = 0;
  int num_v = 0;
  int num_e = 0;
  string line, name;

  EdgeTracer cursor;

  while (getline(ins, line))
  {
    std::stringstream stream(line);
    if (3 < i && i <= 3 + num_v)
    {
      g.resize(i - 4 + 1);
      if (node_dict.find(line.substr(31, 3)) != node_dict.end())
      {
        g.label[i - 4] = node_dict[line.substr(31, 3)];
      }
      else
      {
        g.label[i - 4] = -1;
      }
      //std::cout << "v " << i-4 << " " << node_dict[line.substr(31,3)] << std::endl;
    }
    else if (3 + num_v < i && i <= 3 + num_v + num_e)
    {
      int atom1 = atoi(line.substr(0, 3).c_str()) - 1;
      int atom2 = atoi(line.substr(3, 3).c_str()) - 1;
      int btype = atoi(line.substr(6, 3).c_str()) - 1;
      if (g.label[atom1] != -1 && g.label[atom2] != -1)
      {
        labels.x = g.label[atom1];
        labels.y = btype;
        labels.z = g.label[atom2];
        int etmp = eid;
        edge.to = atom2;
        edge.labels = labels;
        edge.id = eid++;
        g[atom1].push_back(edge);
        edge.to = atom1;
        edge.labels = labels.reverse();
        edge.id = eid++;
        g[atom2].push_back(edge);
        // edge 1->2
        cursor.vpair.a = atom1;
        cursor.vpair.b = atom2;
        cursor.vpair.id = etmp;
        cursor.predec = 0;
        gheap[labels][file].push_back(cursor);
        // edge 2->1
        cursor.vpair.a = atom2;
        cursor.vpair.b = atom1;
        cursor.vpair.id = etmp + 1;
        cursor.predec = 0;
        gheap[labels.reverse()][file].push_back(cursor);
      }
      //std::cout << "e " << atom1 << " " << atom2 << " " << btype << std::endl;
    }
    else if (i == 3)
    {
      num_v = atoi(line.substr(0, 3).c_str());
      num_e = atoi(line.substr(3, 3).c_str());
    }
    else if (i == 0)
    {
      g.clear();
      eid = 0;
      stream >> name;
      //std::cout << name << "->" << file << std::endl;
      g_tid[name] = file;
      g_tname[file] = name;
      //std::cout << "t # " << file << std::endl;
    }
    i += 1;
    if (line.find("$$$$") < line.length())
    {
      file += 1;
      g.num_of_edges = eid;
      gdata.push_back(g);
      i = 0;
      //std::cout << std::endl;
    }
  }
}

template <typename T>
void Finder::read_fasta(T &is)
{
  string line, name;
  string seq = "";
  bool noread = true;
  int file = 0;
  while (getline(is, line))
  {
    string str = "";
    std::istringstream stream(line);
    if (line[0] == '>')
    {
      if (noread)
      {
        noread = false;
      }
      else
      {
        //std::cout << seq << std::endl;
        sdata.push_back(seq);
        seq = "";
      }
      stream >> name;
      name.erase(0, 1);
      //std::cout << name << "->" << file << std::endl;
      s_tid[name] = file;
      s_tname[file] = name;
      file += 1;
      continue;
    }
    stream >> str;
    seq += str;
  }
  //std::cout << seq << std::endl;
  sdata.push_back(seq);
}

int buildDFSCode(const string &str, vector<DFSCode> &query, int _cur, map<int, int> &label, vector<int> &addkeys)
{

  int cur = _cur;
  DFSCode tmp;
  std::istringstream in(str);
  char c;
  int i, j, elab = -200;
  while (in.get(c))
  {
    string token(1, c);
    unsigned int n = 1;
    while (in.get(c) && c != ' ')
    {
      token += c;
      n += 1;
    }
    if (token[0] == '(' && token[n - 1] == ')')
    {
      string sstr = token.substr(1, n - 2);
      string::size_type r;
      if ((r = sstr.find("f")) != string::npos)
      {
        i = atoi(sstr.substr(0, r).c_str());
        j = atoi(sstr.substr(r + 1, string::npos).c_str());
        addkeys.push_back(cur);
        label[cur] = j;
        tmp.time.a = i;
        tmp.time.b = cur;
        tmp.labels.x = label[i];
        tmp.labels.y = elab;
        tmp.labels.z = j;
        query.push_back(tmp);
        cur += 1;
      }
      else if ((r = sstr.find("b")) != string::npos)
      {
        i = atoi(sstr.substr(r + 1, string::npos).c_str());
        tmp.time.a = cur - 1;
        tmp.time.b = i;
        tmp.labels.x = label[cur - 1];
        tmp.labels.y = elab;
        tmp.labels.z = label[i];
        query.push_back(tmp);
      }
      else
      {
        i = atoi(sstr.c_str());
        label[cur] = i;
        addkeys.push_back(cur);
        cur += 1;
      }
    }
    else
    {
      elab = atoi(token.c_str());
    }
  }

  return cur;
}

#endif /*FINDER_H_*/
