#ifndef PATRICIA_H_
#define PATRICIA_H_

#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <iostream>

#include <boost/xpressive/xpressive.hpp>

using std::list;
using std::map;
using std::set;
using std::string;
using std::vector;

typedef vector<string> strSeq;

class patricia
{
  int newest;
  int add(strSeq &, int);

public:
  const static int root;
  vector<strSeq> node;
  vector<int> outQ;
  map<int, set<int> > instance;
  vector<list<int> > children;
  int insert(strSeq &, int, int);
  patricia() : newest(1){};
};

const int patricia::root = 0;

inline int patricia::add(strSeq &c, int parent_node)
{
  int current = newest++;
  node.resize(newest);
  outQ.resize(newest);
  children.resize(newest);
  node[current] = c;
  children[parent_node].push_back(current);
  return current;
}

inline int patricia::insert(strSeq &query_str, int node_id, int tag)
{
  if (newest > 1)
  {
    for (list<int>::const_iterator it = children[node_id].begin();
         it != children[node_id].end(); ++it)
    {
      int child_id = *it;
      strSeq &node_str = node[*it];
      strSeq prefix, query_postfix, node_postfix;
      unsigned matched = 0;
      while (node_str[matched] == query_str[matched])
      {
        prefix.push_back(query_str[matched]);
        matched += 1;
        if (matched == node_str.size())
          break;
        if (matched == query_str.size())
          break;
      }

      for (unsigned j = matched; j < node_str.size(); ++j)
      {
        node_postfix.push_back(node_str[j]);
      }
      for (unsigned j = matched; j < query_str.size(); ++j)
      {
        query_postfix.push_back(query_str[j]);
      }

      if (matched == 0)
      {
        // no match
        continue;
      }
      else if (matched == node_str.size() && matched == query_str.size())
      {
        outQ[child_id] = tag;
        return child_id;
      }
      else if (matched == node_str.size())
      {
        // extend
        return insert(query_postfix, child_id, tag);
      }
      else if (matched == query_str.size())
      {
        // split and extend
        node_str = node_postfix;
        int cur = add(prefix, node_id);
        children[node_id].remove(child_id);
        children[cur].push_back(child_id);
        outQ[cur] = tag;
        return cur;
      }
      else
      {
        // split into two
        node_str = node_postfix;
        int cur = add(prefix, node_id);
        outQ[cur] = -1;
        children[node_id].remove(child_id);
        children[cur].push_back(child_id);
        int cur2 = add(query_postfix, cur);
        outQ[cur2] = tag;
        return cur2;
      }
    }
  }
  int cur = add(query_str, node_id);
  outQ[cur] = tag;
  return cur;
}

void gtokenize(vector<string> &v, string s)
{
  using namespace boost::xpressive;

  sregex re = sregex::compile("(\\(\\d+\\) )*\\d \\([^\\)]+\\)");
  sregex_iterator i(s.begin(), s.end(), re);
  sregex_iterator end;

  while (i != end)
  {
    v.push_back(i->str());
    i++;
  }
}

void stokenize(vector<string> &v, string s)
{
  using namespace boost::xpressive;

  sregex re = sregex::compile(".");
  sregex_iterator i(s.begin(), s.end(), re);
  sregex_iterator end;

  while (i != end)
  {
    v.push_back(i->str());
    i++;
  }
}

#endif /*PATRICIA_H_*/
