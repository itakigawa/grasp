#ifndef TREE_H_
#define TREE_H_

#include <list>
#include <vector>
//#include <tr1/unordered_set>

#include<set>
using std::set;

using std::list;
using std::string;
using std::vector;
//using std::tr1::unordered_set;

template <class T>
class tree {
  int newest;
 public:
  const static int root;
  vector<T> node;
  vector< list<int> > children;
  int add_leaf(T&, int);
 tree(): newest(1) {};
};

template <class T>
const int tree<T>::root = 0;

template <class T>
inline int tree<T>::add_leaf(T& c, int parent_node){
  int current = newest++;
  node.resize(newest);
  children.resize(newest);
  node[current] = c;
  children[parent_node].push_back(current);
  return current;
}

struct vertex {
  string name;
  int num;
  int support;
  int neg_support;
  set<int> exist;
  //unordered_set<int> exist;
};

#endif /* TREE_H_*/
