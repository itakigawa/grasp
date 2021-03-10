#include "gspan.h"
#include <sstream>
#include <deque>

using std::map; using std::vector; using std::deque;

inline void Gspan::scan_rm(vector<DFSCode>& pat, vector<int>& idx){
  int prev = -1;
  for(int i=pat.size()-1; i>=0; --i){
    if(pat[i].time.a < pat[i].time.b){ // forward edge
      if(prev == pat[i].time.b || prev == -1){
	idx.push_back(i);
	prev = pat[i].time.a;
      }
    }
  }
}

bool Gspan::min_checker(vector<DFSCode>& comp, Graph& g, Tracers& _tracers){
  // build right most path
  vector<int> rm_path_index;
  scan_rm(comp,rm_path_index);

  int minlabel = comp[0].labels.x;
  int maxtoc = comp[rm_path_index[0]].time.b;
  Tracers& tracers =_tracers;
	
  vector<VertexPair> vpairs(comp.size());
  map<Pair,Tracers>  b_heap, f_heap;
  EdgeTracer *tracer;
  EdgeTracer cursor;
  Pair pkey;
	
  for(Tracers::iterator it = tracers.begin(); it != tracers.end(); ++it){
    // an instance (a sequence of vertex pairs) as vector "vpair"
    tracer = &(*it);
    //vector<bool> discovered(g.size());
    //vector<bool> tested(g.num_of_edges);
    deque<bool> discovered(g.size());
    deque<bool> tested(g.num_of_edges);

    for(int i = vpairs.size()-1; i >= 0; --i, tracer = tracer->predec){
      vpairs[i] = tracer->vpair;
      tested[vpairs[i].id] = true;
      discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
    }
		
    Pair& rm_vpair = vpairs[rm_path_index[0]]; 
		
    // grow from the right most vertex
    for(unsigned i=0; i<g[rm_vpair.b].size(); ++i){
      Edge& added_edge = g[rm_vpair.b][i];
      //backward from the right most vertex
      for(unsigned j=1; j<rm_path_index.size(); ++j){
	int idx = rm_path_index[j];
	if(tested[added_edge.id]) continue;
	if(vpairs[idx].a != added_edge.to) continue;
	if(comp[idx].labels <= added_edge.labels.reverse()){
	  pkey.set(comp[idx].time.a,added_edge.labels.y);
	  cursor.set(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
	  b_heap[pkey].push_back(cursor);
	}
      }
      // forward from the right most vertex
      if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
      pkey.set(added_edge.labels.y,added_edge.labels.z);
      cursor.set(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
      f_heap[pkey].push_back(cursor);
    }
  }
	
  DFSCode dcode;
  map<Pair,Tracers>::iterator b_it=b_heap.begin();
  if(b_it != b_heap.end()){
    dcode.labels = Triplet(-1,b_it->first.b,-1);
    dcode.time.set(maxtoc,b_it->first.a);
    if(dcode != pattern[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,b_it->second);
  }

  map<Pair,Tracers>::iterator f_it=f_heap.begin();
  if(f_it != f_heap.end()){
    dcode.labels = Triplet(-1,f_it->first.a,f_it->first.b);
    dcode.time.set(maxtoc,maxtoc+1);
    if(dcode != pattern[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,f_it->second);
  }
	
  map<Pair,Tracers> ff_heap;
  int from = -1;
  // forward from the other nodes on the right most path
  for(unsigned j = 0; j < rm_path_index.size(); ++j){
    int i = rm_path_index[j];
    for(Tracers::iterator it =tracers.begin(); it != tracers.end(); ++it){
      tracer = &(*it);
      //vector<bool> discovered(g.size());
      //vector<bool> tested(g.num_of_edges);
      deque<bool> discovered(g.size());
      deque<bool> tested(g.num_of_edges);

      for(int k = vpairs.size()-1; k >= 0; --k, tracer = tracer->predec){
	vpairs[k] = tracer->vpair;
	tested[vpairs[k].id] = true;
	discovered[vpairs[k].a] = discovered[vpairs[k].b] = true;
      }
			
      Pair& from_vpair = vpairs[i];
      for(unsigned k=0; k<g[from_vpair.a].size(); ++k){
	Edge& added_edge = g[from_vpair.a][k];
	if (minlabel > added_edge.labels.z || discovered[added_edge.to]) continue; 
	if(comp[i].labels <= added_edge.labels){
	  from = comp[i].time.a;
	  pkey.set(added_edge.labels.y,added_edge.labels.z);
	  cursor.set(from_vpair.a,added_edge.to,added_edge.id,&(*it));
	  ff_heap[pkey].push_back(cursor);
	}
      }
    }
    if(from != -1) break;
  }
	
  map<Pair,Tracers>::iterator ff_it=ff_heap.begin();
  if(ff_it != ff_heap.end()){
    dcode.labels = Triplet(-1,ff_it->first.a,ff_it->first.b);
    dcode.time.set(from,maxtoc+1);
    if(dcode != pattern[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,ff_it->second);
  }
  return true;
}

bool Gspan::is_min(){
  if(pattern.size() == 1) return true;
	
  Graph g = toGraph(pattern);	
  map<Triplet,Tracers> heap;
  EdgeTracer cursor;
  for(unsigned v=0; v<g.size(); ++v){
    for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
      if (e->labels.x <= e->labels.z){
	cursor.set(v,e->to,e->id,0);
	heap[e->labels].push_back(cursor);
      }
    }
  }
	
  vector<DFSCode> comp;
  comp.resize(1);
  map<Triplet,Tracers>::iterator it = heap.begin();

  comp[0].labels = it->first;
  comp[0].time.set(0,1);

  return min_checker(comp, g,it->second);
}

int Gspan::scan_gspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int,PairSorter,greater<int> >& f_heap){
  // build right most path
  vector<int> rm_path_index;
  scan_rm(pattern,rm_path_index);
	
  int maxtoc = pattern[rm_path_index[0]].time.b;
	
  vector<VertexPair> vpairs(pattern.size());
  int minlabel = pattern[0].labels.x;
  EdgeTracer *tracer;

  Pair pkey;
  EdgeTracer cursor;

  for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
    int gid = x->first;
    for(vector<EdgeTracer>::iterator it = x->second.begin(); it != x->second.end(); ++it){
      // an instance (a sequence of vertex pairs) as vector "vpair"
      tracer = &(*it);
      Graph& g = gdata[gid];
			
      //vector<bool> discovered(g.size());
      //vector<bool> tested(g.num_of_edges);
      deque<bool> discovered(g.size());
      deque<bool> tested(g.num_of_edges);

      for(int i = vpairs.size()-1; i >= 0; --i, tracer = tracer->predec){
	vpairs[i] = tracer->vpair;
	tested[vpairs[i].id] = true;
	discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
      }
			
      Pair& rm_vpair = vpairs[rm_path_index[0]];
		
      for(unsigned i=0; i<g[rm_vpair.b].size(); ++i){
	Edge& added_edge = g[rm_vpair.b][i];
	//backward from the right most vertex
	for(unsigned j=1; j<rm_path_index.size(); ++j){
	  int idx = rm_path_index[j];
	  if(tested[added_edge.id]) continue;
	  if(vpairs[idx].a != added_edge.to) continue;
	  if(pattern[idx].labels <= added_edge.labels.reverse()){
	    pkey.set(pattern[idx].time.a,added_edge.labels.y);
	    cursor.set(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
	    b_heap[pkey][gid].push_back(cursor);
	  }
	}
	// forward from the right most vertex
	if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
	pkey.set(added_edge.labels.y,added_edge.labels.z);
	cursor.set(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
	f_heap[maxtoc][pkey][gid].push_back(cursor);
      }
			
      // forward from the other nodes on the right most path
      for(unsigned j=0; j<rm_path_index.size(); ++j){
	int i = rm_path_index[j];
	Pair& from_vpair = vpairs[i];
	for(unsigned k=0; k<g[from_vpair.a].size(); ++k){
	  Edge& added_edge = g[from_vpair.a][k];
	  if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
					
	  if(pattern[i].labels <= added_edge.labels){
	    pkey.set(added_edge.labels.y,added_edge.labels.z);
	    cursor.set(from_vpair.a,added_edge.to,added_edge.id,&(*it));
	    f_heap[pattern[i].time.a][pkey][gid].push_back(cursor);
	  }
	}
      }
    }
  }
  return maxtoc;

}

unsigned Gspan::support(GraphToTracers& g2tracers){
  int sup = 0;
  for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
    sup += degree[x->first];
  }
  return sup;
}

void Gspan::edge_grow(GraphToTracers& g2tracers, int parent){

  if (support(g2tracers) < minsup) return;
  if (!is_min()){ return; }
	
  int current = report(g2tracers,parent);
	
  PairSorter b_heap;
  map<int,PairSorter,greater<int> > f_heap;
	
  int maxtoc = scan_gspan(g2tracers,b_heap,f_heap);
		
  // projecting...
  DFSCode  dcode;
  for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
    dcode.labels = Triplet(-1,it->first.b,-1);
    dcode.time.set(maxtoc, it->first.a);
    pattern.push_back(dcode);
    edge_grow(it->second,current);
    pattern.pop_back();
  }
	
  for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
    for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
      dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
      dcode.time.set(it->first,maxtoc+1);
      pattern.push_back(dcode);
      edge_grow(it2->second,current);
      pattern.pop_back();
    }
  }
}


int Gspan::report(GraphToTracers& tracers, int parent){
  int psup = 0;
  int nsup = 0;
  vertex ret;
  for(GraphToTracers::iterator x = tracers.begin(); x != tracers.end(); ++x){
    if(classify){
      psup += degree[x->first];
      nsup += neg_degree[x->first];
    }else{
      psup += degree[x->first];
    }
    ret.exist.insert(x->first);
  }

  //std::cout << tracers.size() << " : " << pattern << std::endl;

  std::ostringstream sout;
  sout << pattern;
  ret.name    = sout.str();
  if(classify){
    ret.support = psup;
    ret.neg_support = nsup;
  }else{
    ret.support = psup;
  }
  ret.num     = tracers.size();
  return search_tree.add_leaf(ret,parent);
}

const vector<Graph> readGraphs(std::istream &is) {
  vector<Graph> graphs;
  Graph g;
  Triplet labels;
  Edge edge;
	
  char c; unsigned i, j;
  std::string line;
  int eid=0;

  while (getline(is, line)) {
    if (line.empty()) {
      g.num_of_edges = eid;
      graphs.push_back(g);
    }
    std::stringstream stream(line);
    if (line[0] == 't') {
      g.clear();
      eid = 0;
      stream >> c >> c >> i;
    } else if (line[0] == 'v') {
      stream >> c >> i >> j;
      g.resize(i+1);
      g.label[i] = j;
    } else if (line[0] == 'e') {
      stream >> c >> i >> j >> labels.y;
      labels.x = g.label[i];
      labels.z = g.label[j];
      edge.to = j; edge.labels = labels; edge.id = eid++;
      g[i].push_back(edge);
      edge.to = i; edge.labels = labels.reverse(); edge.id = eid++;
      g[j].push_back(edge);
    }
  }

  return graphs;
}

Graph Gspan::toGraph(vector<DFSCode>& codes){
  Graph g;
  Edge edge;
  int eid = 0;
  for(vector<DFSCode>::iterator p=codes.begin(); p!=codes.end(); ++p){
    g.resize(std::max(p->time.a,p->time.b)+1);
    if(p->labels.x != -1) g.label[p->time.a] = p->labels.x;
    if(p->labels.z != -1) g.label[p->time.b] = p->labels.z;
    edge.to = p->time.b;
    edge.labels.x = g.label[p->time.a];
    edge.labels.y = p->labels.y;
    edge.labels.z = g.label[p->time.b];
    edge.id = eid++;
    g[p->time.a].push_back(edge);
    edge.to = p->time.a; edge.labels = edge.labels.reverse(); edge.id = eid++;
    g[p->time.b].push_back(edge);
  }
  g.num_of_edges = eid;
  return g;
}

void Gspan::build_tree(){
  map<Triplet,GraphToTracers> heap;	

  for(unsigned gid = 0; gid < gdata.size(); ++gid){
    EdgeTracer cursor;
    Graph& g = gdata[gid];
    for(unsigned v=0; v<g.size(); ++v){
      for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
	if (e->labels.x <= e->labels.z){
	  cursor.set(v,e->to,e->id,0);
	  heap[e->labels][gid].push_back(cursor);
	}
      }
    }
  }

  pattern.resize(1);
  for(map<Triplet,GraphToTracers>::iterator it = heap.begin(); it != heap.end(); ++it){		
    pattern[0].labels = it->first;
    pattern[0].time.set(0,1);
    edge_grow(it->second,search_tree.root);
  }
}
