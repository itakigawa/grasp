#include "gspan.h"
#include "prefixspan.h"
#include "tree.h"

#include <boost/date_time/posix_time/posix_time.hpp>

using std::istringstream;

using boost::posix_time::microsec_clock;
using boost::posix_time::ptime;
using boost::posix_time::time_duration;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using boost::iostreams::filtering_istream;
using boost::iostreams::gzip_decompressor;

#include "boost/filesystem/operations.hpp"

#define USAGE " [-m minsup] [-d delta] [-n numout] [-o outdir] [-v] graphs sequences pairs"

struct datapack {
public:
  int seq;
  int cmp;
  string seqname;
  string cmpname;
  int label;
};

template<typename T> void read_pairs(T& ins,
				     unordered_map<int,vector<Destpair> >& links,
				     Gspan& gs, Prefixspan& ps, vector<datapack>& input){
  //std::ofstream cmpdict("ctrans.dict");
  //std::ofstream prtdict("ptrans.dict");
  string line, cmp, seq, label;
  bool first = true;
  ps.classify = false;
  gs.classify = false;
  int count = 1;
  while ( getline(ins, line) ) {
    std::stringstream stream(line);
    stream >> seq >> cmp;
    if (first) {
      if (stream >> label){
	std::cout << ">interaction datasets with label" << std::endl;
	gs.classify = true;
	ps.classify = true;
      }
      first = false;
    }else if (ps.classify){
      if (!(stream >> label)){
	std::cerr << "Wrong pair file format. (with label or not?)" << std::endl;
      }
    }
    if (gs.mapper.find(cmp) == gs.mapper.end() || ps.mapper.find(seq) == ps.mapper.end()){
      continue;
    }
    if(ps.verbose){
      std::cout << cmp << "-->" << gs.mapper[cmp] << std::endl;
      std::cout << seq << "-->" << ps.mapper[seq] << std::endl;
      //cmpdict << cmp << "-->" << gs.mapper[cmp] << std::endl;
      //prtdict << seq << "-->" << ps.mapper[seq] << std::endl;
    }
    Destpair dest;
    int gi = gs.mapper[cmp];
    int si = ps.mapper[seq];
    datapack tmp;
    tmp.cmp = gi;
    tmp.cmpname = cmp;
    tmp.seq = si;
    tmp.seqname = seq;

    dest.graph = gi;
    dest.id = count;
    count++;

    if (ps.classify){
      if (label == "+"){
	dest.is_positive = true;
	tmp.label = 1;
	ps.num_links++;
	gs.degree[gi]++;
	ps.degree[si]++;
      } else{
	dest.is_positive = false;
	tmp.label = -1;
	ps.neg_num_links++;
	gs.neg_degree[gi]++;
	ps.neg_degree[si]++;
      }
      links[si].push_back( dest );
    } else{
      dest.is_positive = true;
      gs.degree[gi]++;
      ps.degree[si]++;
      ps.num_links++;
      links[si].push_back( dest );
    }
    input.push_back(tmp);
  }
  //cmpdict.close();
  //prtdict.close();
}

int main(int argc, char **argv) {
  unsigned int maxout = 100;
  unsigned int minsup = 0;
  unsigned int delta = 0;
  string outdir;
  bool verbose = false;

  // *********** opt check ***********

  int opt;

  while ((opt = getopt(argc, argv, "m:d:n:o:v")) != -1) {
    switch (opt) {
    case 'm':
      minsup = atoi (optarg);
      break;
    case 'd':
      delta = atoi (optarg);
      break;
    case 'n':
      maxout = atoi (optarg);
      break;
    case 'v':
      verbose = true;
      break;
    case 'o':
      outdir = string(optarg);
      if( !boost::filesystem::exists( outdir ) ){
	std::cerr << "Directory does not exist: "<< outdir << std::endl;
	return -1;
      }
      break;
    default:
      std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
      return -1;
    }
  }

  // *********** i/o check ***********

  if(argc-optind != 3){
    std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
    return -1;
  }
  
  // *********** file check ***********

  std::ifstream graph_file(argv[optind]);
  if(graph_file.fail()){
    std::cerr << "File not found: " << argv[optind] << std::endl;
    return -1;
  }
  std::ifstream sequence_file(argv[optind+1]);
  if(sequence_file.fail()){
    std::cerr << "File not found: " << argv[optind+1] << std::endl;
    return -1;
  }
  std::ifstream pair_file(argv[optind+2]);
  if(pair_file.fail()){
    std::cerr << "File not found: " << argv[optind+2] << std::endl;
    return -1;
  }


  string filegraph  (argv[optind]);
  string filestring (argv[optind+1]);
  string filepair   (argv[optind+2]);

  ptime now = microsec_clock::local_time();

  std::cout << "GRASP version 0.2"          << std::endl;
  std::cout << "      by Ichigaku Takigawa" << std::endl;
  std::cout << ">starting at " << now       << std::endl;
  std::cout << ">settings are:"             << std::endl;
  std::cout << "  Graph  => " << filegraph  << std::endl;
  std::cout << "  String => " << filestring << std::endl;
  std::cout << "  Pair   => " << filepair   << std::endl;
  std::cout << "  minsup => " << minsup     << std::endl;
  std::cout << "  delta  => " << delta      << std::endl;
  std::cout << "  maxout => " << maxout     << std::endl;

  delta += 1;

  // *********** main ***********

  Gspan gspan;
  gspan.minpat = 0;
  gspan.minsup = minsup;
  unordered_map<string,int> node_dict;

  string outf("atom.dict");
  if(!outdir.empty()) outf = outdir+"/"+outf;
  std::cout << ">output to: " << outf << std::endl;

  if (filegraph.substr(filegraph.size()-3,3) == ".gz"){
    filtering_istream fs;
    fs.push( gzip_decompressor() );
    fs.push( graph_file );
    make_dict(fs,minsup,node_dict,outf);

    graph_file.clear();
    graph_file.seekg(0,std::ios::beg);
    fs.reset();
    fs.push( gzip_decompressor() );
    fs.push( graph_file );
    gspan.read(fs,node_dict);
  }else{
    make_dict(graph_file,minsup,node_dict,outf);
    graph_file.clear();
    graph_file.seekg(0,std::ios::beg);
    gspan.read(graph_file,node_dict);
  }

  Prefixspan pspan;
  pspan.delta  = delta;
  pspan.minsup = minsup;
  pspan.maxout = maxout;
  pspan.num_out = 0;
  //pspan.verbose = verbose;
  pspan.verbose = false;
  pspan.num_links = 0;
  pspan.neg_num_links = 0;
  if (filestring.substr(filestring.size()-3,3) == ".gz"){
    filtering_istream fs;
    fs.push( gzip_decompressor() );
    fs.push( sequence_file );
    pspan.read(fs);
  }else{
    pspan.read(sequence_file);
  }

  std::cout << gspan.gdata.size() << " graphs, ";
  std::cout << pspan.sdata.size() << " sequences" << std::endl;

  unordered_map< int, vector<Destpair> > links;
  vector<datapack> input;

  if (filepair.substr(filepair.size()-3,3) == ".gz"){
    filtering_istream fs;
    fs.push( gzip_decompressor() );
    fs.push( pair_file );
    read_pairs(fs,links,gspan,pspan,input);
  }else{
    read_pairs(pair_file,links,gspan,pspan,input);
  }
  if (pspan.classify){
    std::cout << pspan.num_links << " positive interactions read." << std::endl;
    std::cout << pspan.neg_num_links << " negative interactions read." << std::endl;
  }else{
    std::cout << pspan.num_links << " interactions read." << std::endl;
  }

  std::cout << ">invoke gspan.." << std::endl;

  gspan.build_tree();

  time_duration elapsed = microsec_clock::local_time() - now; 
  std::cout << " gspan done. " << elapsed << std::endl;

  std::cout << ">co-mining by pspan with gspan tree" << std::endl;

  pspan.co_mine(links,gspan);

  elapsed = microsec_clock::local_time() - now; 
  std::cout << " pair mining done. " << elapsed << std::endl;

  unsigned rnk = 0;
  int feat_id;
  //map<int,vector<int> > attr;

  outf = "pair.patterns";
  if(!outdir.empty()) outf = outdir+"/"+outf;
  std::ofstream pair_out(outf.c_str());
  std::cout << ">output to: " << outf << std::endl;

  string fe = "feat.instances";
  if(!outdir.empty()) fe = outdir+"/"+fe;
  std::cout << ">output to: " << fe << std::endl;
  std::ofstream flist_out(fe.c_str());

  for(multimap<double,OutObject>::iterator it = pspan.ranking.begin();
      it != pspan.ranking.end(); ++it){
    feat_id = ++rnk;
    pair_out << feat_id << "|" << it->first << "|";
    pair_out << it->second.substr << "|" << it->second.subgraph << "|";
    pair_out << "gNsN" << it->second.pos(0) << ",gYsN" << it->second.pos(1) << ",";
    pair_out << "gNsY" << it->second.pos(2) << ",gYsY" << it->second.pos(3) << "|";
    pair_out << "gNsN" << it->second.neg(0) << ",gYsN" << it->second.neg(1) << ",";
    pair_out << "gNsY" << it->second.neg(2) << ",gYsY" << it->second.neg(3) << "|";
    pair_out << "pair" << it->second.sup << "/" << pspan.num_links << "(g" << it->second.pos(1)+it->second.pos(3) << ",s" << it->second.pos(2)+it->second.pos(3) << ")";

    if (!pspan.classify){
      pair_out << "|g" << it->second.hit_g <<  "/" << it->second.num_g << ",s" << it->second.hit_s <<  "/" << it->second.num_s;
    }

    pair_out << std::endl;

    if(verbose){
      std::cout << "feat" << feat_id << " " << it->second.substr << " "
		<< it->second.subgraph << std::endl;
      std::cout << feat_id << " " << it->second.notes << std::endl;
      std::cout << "[" << it->second.sup << "]" << std::endl;
    }
    flist_out << feat_id << "|";
    bool start = true;
    for(set<int>::iterator it2=it->second.instances.begin();
	it2!=it->second.instances.end(); ++it2){
      if(start){
	start = false;
      }else{
	flist_out << ",";
      }
      flist_out << *it2;
    }
    flist_out << std::endl;
  }
  flist_out.close();
  pair_out.close();
  
  std::cout << "# frequent pairs: " << pspan.num_out << std::endl;
  
  elapsed = microsec_clock::local_time() - now; 
  std::cout << "all done. " << elapsed << std::endl;

  return 0;
}
