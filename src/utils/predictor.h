#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "finder.h"

#include <sstream>
#include <boost/dynamic_bitset.hpp>

using boost::dynamic_bitset;
using std::istringstream;
using std::multiset;

struct PairPattern
{
public:
  string g;
  string s;
  double score;
  int ptype;
  int count1;
  int count2;
};

inline bool operator<(const PairPattern &left, const PairPattern &right)
{
  return (left.score > right.score);
}

struct Container
{
public:
  dynamic_bitset<> key;
  int count;
  int index;
  Container(){};
  Container(dynamic_bitset<> &_key, int _index)
  {
    key = _key;
    index = _index;
    count = key.count();
  };
};

class Predictor : public Finder
{
public:
  multiset<PairPattern> result;
  vector<Container> queue;
  template <typename T>
  void read_table(T &, unsigned);
};

struct larger_than : public std::binary_function<Container, Container, bool>
{
public:
  bool operator()(const Container &lh, const Container &rh)
  {
    return lh.count > rh.count;
  }
};

template <typename T>
void Predictor::read_table(T &ins, unsigned ulimit)
{

  std::string line;
  std::multiset<Container, larger_than> b;

  int nline = 0;
  while (getline(ins, line))
  {
    std::stringstream stream(line);
    int ith;
    dynamic_bitset<> a(10000);
    nline += 1;
    while (stream >> ith)
    {
      a.set(ith - 1, true);
    }
    if (a.count() == 0)
      continue;
    Container bag(a, nline - 1);
    b.insert(bag);
  }

  for (std::multiset<Container>::iterator u = b.begin(); u != b.end(); ++u)
  {
    queue.push_back(*u);
  }

  std::cout << "# training samples: " << queue.size() << std::endl;
  return;
}

#endif /*PREDICTOR_H_*/
