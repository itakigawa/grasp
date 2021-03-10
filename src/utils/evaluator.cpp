#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <boost/xpressive/xpressive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "predictor.h"

using std::deque;
using std::string;
using std::vector;

using std::istringstream;

using std::unordered_map;

using boost::posix_time::microsec_clock;
using boost::posix_time::ptime;
using boost::posix_time::time_duration;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using boost::iostreams::filtering_istream;
using boost::iostreams::gzip_decompressor;

#define USAGE " features atom.dict graphs sequences pairs train.table"

void straverse(int node_id, string &buf, unsigned pnt,
               occurence &oc, vector<string> &D, unsigned delta, patricia &P,
               map<int, set<int> > &saver)
{
  if (buf.length() > pnt)
  {
    occurence oc_new;
    for (occurence::iterator p = oc.begin(); p != oc.end(); ++p)
    {
      int id = p->first;
      for (vector<int>::iterator q = p->second.begin(); q != p->second.end(); ++q)
      {
        unsigned ulimit = D[id].length();
        if (*q + delta + 1 < ulimit)
          ulimit = *q + delta + 1;
        for (unsigned pos = *q + 1; pos < ulimit; ++pos)
        {
          char key = D[id][pos];
          if (key == buf[pnt])
          {
            oc_new[id].push_back(pos);
          }
        }
      }
    }
    straverse(node_id, buf, pnt + 1, oc_new, D, delta, P, saver);
  }
  else
  {
    if (P.outQ[node_id] > 0)
    {
      //std::cout << node_id << " " << buf << " [" << oc.size() << "]" << std::endl;
      for (occurence::iterator p = oc.begin(); p != oc.end(); ++p)
      {
        //saver[node_id].insert(p->first);
        //saver[P.outQ[node_id]].insert(p->first);
        saver[p->first].insert(P.outQ[node_id]);
      }
    }
    list<int> &child = P.children[node_id];
    for (list<int>::iterator it = child.begin(); it != child.end(); ++it)
    {
      string newstr = boost::algorithm::join(P.node[*it], "");
      string newbuf(buf + newstr);
      straverse(*it, newbuf, pnt, oc, D, delta, P, saver);
    }
  }
}

void gtraverse(int node_id, vector<DFSCode> &query, unsigned pnt, int cur,
               GraphToTracers &g2tracers, vector<Graph> &gdata,
               const patricia &P, map<int, set<int> > &saver, int treelevel,
               map<int, int> &dfslabel)
{

  if (g2tracers.size() == 0)
    return;

  if (query.size() > pnt)
  {
    DFSCode &dfscur = query[pnt];
    int s = dfscur.time.a;
    int t = dfscur.time.b;
    EdgeTracer *tracer;
    EdgeTracer cursor;
    GraphToTracers newtracer;

    for (GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x)
    {
      int gid = x->first;
      Graph &g = gdata[gid];
      for (vector<EdgeTracer>::iterator it = x->second.begin(); it != x->second.end(); ++it)
      {

        tracer = &(*it);

        vector<unsigned> discovered(g.size(), 0);
        int growpoint = -999;
        int backpoint = -999;
        for (int j = pnt - 1; j >= 0; --j, tracer = tracer->predec)
        {
          discovered[tracer->vpair.a] = 1;
          discovered[tracer->vpair.b] = 1;
          if (query[j].time.a == s)
            growpoint = tracer->vpair.a;
          if (query[j].time.b == s)
            growpoint = tracer->vpair.b;
          if (s > t)
          {
            if (query[j].time.a == t)
              backpoint = tracer->vpair.a;
            if (query[j].time.b == t)
              backpoint = tracer->vpair.b;
          }
        }

        for (unsigned int i = 0; i < g[growpoint].size(); ++i)
        {
          const Edge &e = g[growpoint][i];
          if (e.labels.y == dfscur.labels.y && e.labels.z == dfscur.labels.z)
          {
            if (s < t)
            {
              if (discovered[e.to] != 1)
              {
                cursor.vpair.a = growpoint;
                cursor.vpair.b = e.to;
                cursor.vpair.id = e.id;
                cursor.predec = &(*it);
                newtracer[gid].push_back(cursor);
              }
            }
            else
            {
              if (backpoint == e.to)
              {
                cursor.vpair.a = growpoint;
                cursor.vpair.b = e.to;
                cursor.vpair.id = e.id;
                cursor.predec = &(*it);
                newtracer[gid].push_back(cursor);
              }
            }
          }
        }
      }
    }
    gtraverse(node_id, query, pnt + 1, cur, newtracer, gdata, P, saver, treelevel, dfslabel);
  }
  else
  {

    if (P.outQ[node_id] > 0)
    {
      for (GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x)
      {
        saver[x->first].insert(P.outQ[node_id]);
      }
    }

    const list<int> &child = P.children[node_id];
    for (list<int>::const_iterator it = child.begin(); it != child.end(); ++it)
    {

      const vector<string> &vec = P.node[*it];
      int addcount = 0;
      int c2 = cur;
      vector<int> addkeys;
      for (vector<string>::const_iterator it2 = vec.begin(); it2 != vec.end(); ++it2)
      {
        c2 = buildDFSCode(*it2, query, c2, dfslabel, addkeys);
        addcount += 1;
      }
      gtraverse(*it, query, pnt, c2, g2tracers, gdata, P, saver, treelevel + 1, dfslabel);
      for (vector<int>::iterator lit = addkeys.begin(); lit != addkeys.end(); ++lit)
      {
        dfslabel.erase(*lit);
      }
      for (int j = 0; j < addcount; ++j)
      {
        query.pop_back();
      }
    }
  }
}

int main(int argc, char **argv)
{

  int opt;
  unsigned ulimit = 100;
  unsigned maxnum = 10000;
  unsigned delta = 1;
  string outfile = "prediction";
  while ((opt = getopt(argc, argv, "d:o:n:u:")) != -1)
  {
    switch (opt)
    {
    case 'd':
      delta = atoi(optarg);
      break;
    case 'o':
      outfile = string(optarg);
      break;
    case 'n':
      maxnum = atoi(optarg);
      break;
    case 'u':
      ulimit = atoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0] << USAGE << std::endl;
      return -1;
    }
  }

  // *********** i/o check ***********

  if (argc - optind != 6)
  {
    std::cerr << "Usage: " << argv[0] << USAGE << std::endl;
    return -1;
  }

  // *********** mode check ***********

  Predictor f;

  string filefeat(argv[optind]);
  string fileatom(argv[optind + 1]);
  string filegraph(argv[optind + 2]);
  string filestring(argv[optind + 3]);
  string filepair(argv[optind + 4]);
  string filetable(argv[optind + 5]);

  // *********** file check ***********

  std::ifstream feat_file(filefeat.c_str());
  if (feat_file.fail())
  {
    std::cerr << "File not found: " << filefeat << std::endl;
    return -1;
  }
  std::ifstream graph_file(filegraph.c_str());
  if (graph_file.fail())
  {
    std::cerr << "File not found: " << filegraph << std::endl;
    return -1;
  }
  std::ifstream atom_file(fileatom.c_str());
  if (atom_file.fail())
  {
    std::cerr << "File not found: " << fileatom << std::endl;
    return -1;
  }
  std::ifstream sequence_file(filestring.c_str());
  if (sequence_file.fail())
  {
    std::cerr << "File not found: " << filestring << std::endl;
    return -1;
  }
  std::ifstream pair_file;
  pair_file.open(filepair.c_str());
  if (pair_file.fail())
  {
    std::cerr << "File not found: " << filepair << std::endl;
    return -1;
  }
  std::ifstream table_file(filetable.c_str());
  if (table_file.fail())
  {
    std::cerr << "File not found: " << filetable << std::endl;
    return -1;
  }

  ptime now = microsec_clock::local_time();

  std::cout << "Evaluator version 0.1" << std::endl;
  std::cout << "      by Ichigaku Takigawa" << std::endl;
  std::cout << ">starting at " << now << std::endl;
  std::cout << ">settings are:" << std::endl;
  std::cout << "  Feature => " << filefeat << std::endl;
  std::cout << "  Graph   => " << filegraph << std::endl;
  std::cout << "  String  => " << filestring << std::endl;
  std::cout << "  Pair    => " << filepair << std::endl;
  std::cout << "  Table   => " << filetable << std::endl;

  // *********** main ***********

  unordered_map<string, int> node_dict;

  boost::xpressive::sregex re = boost::xpressive::sregex::compile("^v (\\d+) (.+)$");
  boost::xpressive::smatch m;

  string line;
  while (getline(atom_file, line))
  {
    if (boost::xpressive::regex_match(line, m, re))
    {
      string i = m[1];
      string nm = m[2];
      node_dict[nm] = atoi(i.c_str());
      //std::cout << "[" << atoi(i.c_str()) << "]..[" << nm << "]" << std::endl;
    }
  }

  if (filefeat.substr(filefeat.size() - 3, 3) == ".gz")
  {
    // fs is a stream
    filtering_istream fs;
    fs.push(gzip_decompressor());
    fs.push(feat_file);
    f.read_features(fs);
  }
  else
  {
    // feat_file is a stream
    f.read_features(feat_file);
  }

  if (filegraph.substr(filegraph.size() - 3, 3) == ".gz")
  {
    // fs is a stream
    filtering_istream fs;
    fs.push(gzip_decompressor());
    fs.push(graph_file);
    f.read_sdf(fs, node_dict);
  }
  else
  {
    // graph_file is a stream
    f.read_sdf(graph_file, node_dict);
  }

  if (filestring.substr(filestring.size() - 3, 3) == ".gz")
  {
    // fs is a stream
    filtering_istream fs;
    fs.push(gzip_decompressor());
    fs.push(sequence_file);
    f.read_fasta(fs);
  }
  else
  {
    // sequence_file is a stream
    f.read_fasta(sequence_file);
  }

  if (filepair.substr(filepair.size() - 3, 3) == ".gz")
  {
    filtering_istream fs;
    fs.push(gzip_decompressor());
    fs.push(pair_file);
    f.read_pairs(fs);
  }
  else
  {
    f.read_pairs(pair_file);
  }

  if (filetable.substr(filetable.size() - 3, 3) == ".gz")
  {
    // fs is a stream
    filtering_istream fs;
    fs.push(gzip_decompressor());
    fs.push(table_file);
    f.read_table(table_file, ulimit);
  }
  else
  {
    // table_file is a stream
    f.read_table(table_file, ulimit);
  }

  std::cout << "# graphs: " << f.gdata.size() << std::endl;
  std::cout << "# sequences: " << f.sdata.size() << std::endl;

  //std::ofstream o_file("hoge1.out");
  //print(f.gdata,o_file);

  time_duration elapsed = microsec_clock::local_time() - now;
  std::cout << "done. " << elapsed << std::endl;

  map<char, occurence> octracer;
  for (unsigned int id = 0; id < f.sdata.size(); ++id)
  {
    for (unsigned int pos = 0; pos < f.sdata[id].length(); ++pos)
    {
      octracer[f.sdata[id][pos]][id].push_back(pos);
    }
  }

  map<int, set<int> > ssave;

  list<int> &childs = f.spat_tree.children[patricia::root];
  for (list<int>::iterator it = childs.begin(); it != childs.end(); ++it)
  {
    string str = boost::algorithm::join(f.spat_tree.node[*it], "");
    occurence &oc_new = octracer[str[0]];
    straverse(*it, str, 1, oc_new, f.sdata, delta, f.spat_tree, ssave);
  }

  //   for(map<int,set<int> >::iterator ait = ssave.begin(); ait != ssave.end(); ++ait){
  //     std::cout << ait->first << " " << ait->second.size() << std::endl;
  //   }

  elapsed = microsec_clock::local_time() - now;
  std::cout << "smatch done. " << elapsed << std::endl;

  map<int, set<int> > gsave;

  vector<DFSCode> gquery;
  map<int, int> dfslabel;

  list<int> &childg = f.gpat_tree.children[patricia::root];
  for (list<int>::iterator it = childg.begin(); it != childg.end(); ++it)
  {
    vector<string> &vec = f.gpat_tree.node[*it];
    int cur = 0;
    int addcount = 0;
    vector<int> addkeys;
    for (vector<string>::iterator it2 = vec.begin(); it2 != vec.end(); ++it2)
    {
      cur = buildDFSCode(*it2, gquery, cur, dfslabel, addkeys);
      addcount += 1;
    }
    GraphToTracers &newtracer = f.gheap[gquery[0].labels];
    gtraverse(*it, gquery, 1, cur, newtracer, f.gdata, f.gpat_tree, gsave, 1, dfslabel);
    for (vector<int>::iterator lit = addkeys.begin(); lit != addkeys.end(); ++lit)
    {
      dfslabel.erase(*lit);
    }
    for (int j = 0; j < addcount; ++j)
    {
      gquery.pop_back();
    }
  }

  elapsed = microsec_clock::local_time() - now;
  std::cout << "gmatch done. " << elapsed << std::endl;

  ptime main_start = microsec_clock::local_time();

  int count = 1;
  int allsum = (f.gdata.size() * f.sdata.size() / 100);

  double max = 0.0, max_prev = 0.0;

  for (vector<datapack>::iterator itr = f.trans.begin(); itr != f.trans.end(); ++itr)
  {
    int id_s = itr->seq;
    int id_g = itr->cmp;
    set<int> &tmpset_s = ssave[id_s];
    set<int> &tmpset_g = gsave[id_g];

    if (allsum != 0 and count % allsum == 0)
    {
      std::cout << count / allsum << "% " << count << " " << f.s_tname[id_s] << "<->" << f.g_tname[id_g] << std::endl;
    }
    count += 1;

    dynamic_bitset<> row(f.num_of_feat);

    for (set<int>::iterator it1 = tmpset_s.begin(); it1 != tmpset_s.end(); ++it1)
    {
      map<int, map<int, int> >::iterator mapper = f.links.find(*it1);
      if (mapper == f.links.end())
        continue;
      for (set<int>::iterator it2 = tmpset_g.begin(); it2 != tmpset_g.end(); ++it2)
      {
        map<int, int>::iterator mapper2 = mapper->second.find(*it2);
        if (mapper2 != mapper->second.end())
        {
          row.set(mapper2->second - 1, true);
        }
      }
    }

    // model
    Container q(row, 0);
    //if(q.count==0) continue;

    vector<Container>::iterator forward = std::lower_bound(f.queue.begin(), f.queue.end(), q, larger_than());
    vector<Container>::reverse_iterator back(forward);

    int target = -1;
    int c = 0;

    bool bend = true, fend = true;
    while (bend || fend)
    {
      if (max == 1.0)
        break;
      // bounding
      if (back == f.queue.rend() || max > (double)q.count / (double)back->count)
      {
        bend = false;
      }
      if (forward == f.queue.end() || max > (double)forward->count / (double)q.count)
      {
        fend = false;
      }

      if (bend)
      {
        std::size_t size = std::max(row.size(), back->key.size());
        row.resize(size);
        back->key.resize(size);
        double val = (double)(row & back->key).count() / (double)(row | back->key).count();
        if (val > max)
        {
          max = val;
          target = back->index;
          c = back->count;
        }
        back++;
      }
      if (fend)
      {
        std::size_t size = std::max(row.size(), forward->key.size());
        row.resize(size);
        forward->key.resize(size);
        double val = (double)(row & forward->key).count() / (double)(row | forward->key).count();
        if (val > max)
        {
          max = val;
          target = forward->index;
          c = forward->count;
        }
        forward++;
      }
    }
    //if(max == max_prev) continue;

    PairPattern pp;
    pp.count1 = q.count;
    pp.count2 = c;

    pp.s = f.s_tname[id_s];
    pp.g = f.g_tname[id_g];
    pp.score = max;
    pp.ptype = target;
    f.result.insert(pp);

    set<PairPattern>::iterator rit = f.result.end();
    //if(f.result.size()>maxnum){
    //	f.result.erase(--rit);
    //}
    max = rit->score;
    max_prev = max;
  }

  elapsed = microsec_clock::local_time() - main_start;
  std::cout << "main: " << elapsed << std::endl;

  string filename = outfile;
  std::ofstream o_file(filename.c_str());
  for (multiset<PairPattern>::iterator it = f.result.begin(); it != f.result.end(); ++it)
  {
    o_file << it->s << " " << it->g << " " << it->score << " " << it->ptype << " " << it->count1 << " " << it->count2 << std::endl;
  }
  o_file.close();

  elapsed = microsec_clock::local_time() - now;
  std::cout << "all done. " << elapsed << std::endl;

  return 0;
}
