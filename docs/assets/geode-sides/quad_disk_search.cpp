// Enumerate regular quadrilateral disks with a prescribed left/right reflection.
// Independent small-disk enumerator; no prescribed interior vertex valences.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using Quad = std::array<int,4>;
using Edge = std::pair<int,int>;
Edge edge(int a,int b) { return std::minmax(a,b); }
Quad cycle(Quad q) {
  Quad best=q;
  for(int i=0;i<3;++i) { std::rotate(q.begin(),q.begin()+1,q.end()); best=std::min(best,q); }
  return best;
}
bool isedge(Quad q,int a,int b) {
  for(int i=0;i<4;++i) if(edge(q[i],q[(i+1)%4])==edge(a,b)) return true;
  return false;
}
struct State {
  std::vector<Quad> quads;
  std::vector<int> mirror,color;
  std::map<Edge,Edge> front;
};
struct Search {
  int r,cap,boundary,limit,bottom,top;
  bool symmetry;
  std::vector<int> side_mask;
  std::set<Edge> rim;
  std::set<std::vector<Quad>> emitted;
  std::ostream& output;
  unsigned long long visits=0,solutions=0;
  Search(int subdivisions,int max_quads,std::ostream& out,bool use_orbits=true,int b=2,int t=2)
    :r(subdivisions),cap(max_quads),boundary(b+t+2+2*r),limit(cap+boundary/2+1),bottom(b),top(t),symmetry(use_orbits),output(out) {}

  State initial() {
    State s;
    // CCW boundary: bottom to the right, right upward, top to the left,
    // then left downward. The reflection sends vertex i to bottom-i.
    side_mask.resize(boundary);
    for(int i=0;i<boundary;++i) {
      int mask=0;
      if(i<=bottom) mask|=1;
      if(i>=bottom&&i<=bottom+r+1) mask|=2;
      if(i>=bottom+r+1&&i<=bottom+r+1+top) mask|=4;
      if(i>=bottom+r+1+top||i==0) mask|=8;
      side_mask[i]=mask;
      s.mirror.push_back((bottom-i+boundary)%boundary);
      s.color.push_back(i%2);
      Edge e{i,(i+1)%boundary};rim.insert(edge(e.first,e.second));s.front[edge(e.first,e.second)]=e;
    }
    return s;
  }
  bool compatible(Quad a,Quad b) const {
    std::vector<int> common;
    for(int v:a) if(std::find(b.begin(),b.end(),v)!=b.end()) common.push_back(v);
    if(common.size()>2) return false;
    return common.size()!=2||(isedge(a,common[0],common[1])&&isedge(b,common[0],common[1]));
  }
  bool allowed(Quad q) const {
    auto sorted=q;std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end()) return false;
    for(int i=0;i<4;++i) {
      int a=q[i],b=q[(i+1)%4],c=q[(i+2)%4];
      if(a<boundary&&b<boundary) {
        if((side_mask[a]&side_mask[b])&&!rim.count(edge(a,b))) return false;
      }
      // A diagonal lying on a box side cannot have its other endpoints on
      // opposite sides of its supporting line inside the square.
      if(a<boundary&&c<boundary&&(side_mask[a]&side_mask[c])) return false;
      if(a<boundary&&b<boundary&&c<boundary) {
        if(side_mask[a]&side_mask[b]&side_mask[c]) return false;
        // Three boundary vertices of a CCW convex quad must be cyclically CCW.
        if((b-a+boundary)%boundary >= (c-a+boundary)%boundary) return false;
      }
    }
    return true;
  }
  bool links(const State& s,bool final) const {
    std::set<int> open;
    for(auto [key,e]:s.front) { (void)key;open.insert(e.first);open.insert(e.second); }
    int deficits=0;
    for(int v=0;v<int(s.mirror.size());++v) {
      std::map<int,std::vector<int>> graph;int incidence=0;
      for(auto q:s.quads) for(int j=0;j<4;++j) if(q[j]==v) {
        ++incidence;int a=q[(j+3)%4],b=q[(j+1)%4];graph[a].push_back(b);graph[b].push_back(a);
      }
      bool corner=v<boundary&&__builtin_popcount(unsigned(side_mask[v]))==2;
      int minimum=v>=boundary?3:corner?1:2;
      deficits+=std::max(0,minimum-incidence);
      for(const auto& [a,ns]:graph) { (void)a;if(ns.size()>2) return false; }
      if(!final&&open.count(v)) continue;
      if(incidence<minimum) return false;
      std::set<int> seen;std::vector<int> stack{graph.begin()->first};
      while(!stack.empty()) { int a=stack.back();stack.pop_back();if(!seen.insert(a).second) continue;
        for(int b:graph[a]) stack.push_back(b); }
      if(seen.size()!=graph.size()) return false;
      int ends=0;for(auto [a,ns]:graph) { (void)a;ends+=ns.size()==1; }
      if(ends!=(v<boundary?2:0)) return false;
      if(v<boundary) for(int a:{(v+boundary-1)%boundary,(v+1)%boundary})
        if(graph[a].size()!=1) return false;
    }
    return deficits<=4*(cap-int(s.quads.size()));
  }
  bool insert(State& s,Quad q) const {
    Quad reflected{s.mirror[q[0]],s.mirror[q[3]],s.mirror[q[2]],s.mirror[q[1]]};
    std::vector<Quad> orbit{q};
    if(symmetry&&cycle(q)!=cycle(reflected)) orbit.push_back(reflected);
    if(s.quads.size()+orbit.size()>size_t(cap)) return false;
    for(auto cell:orbit) {
      if(!allowed(cell)) return false;
      for(auto old:s.quads) if(!compatible(cell,old)) return false;
      std::map<Edge,int> incidence;
      for(auto old:s.quads) for(int j=0;j<4;++j) ++incidence[edge(old[j],old[(j+1)%4])];
      for(int j=0;j<4;++j) {
        Edge e{cell[j],cell[(j+1)%4]},k=edge(e.first,e.second);
        if(incidence[k]>=(rim.count(k)?1:2)) return false;
        auto found=s.front.find(k);
        if(found!=s.front.end()) { if(found->second!=e) return false;s.front.erase(found); }
        else s.front[k]={e.second,e.first};
      }
      s.quads.push_back(cell);
    }
    if(s.front.size()>4*(cap-s.quads.size())) return false;
    return links(s,false);
  }
  std::vector<Quad> canonical(const State& s,std::vector<int>& labels) const {
    // Unique traversal of the oriented disk from its fixed, labeled rim.
    std::map<Edge,std::vector<int>> owners;
    for(int i=0;i<int(s.quads.size());++i) for(int j=0;j<4;++j)
      owners[edge(s.quads[i][j],s.quads[i][(j+1)%4])].push_back(i);
    labels.assign(s.mirror.size(),-1);for(int v=0;v<boundary;++v) labels[v]=v;
    int next=boundary;std::vector<bool> used(s.quads.size());
    std::vector<Edge> queue{{0,1}};std::vector<Quad> result;
    for(size_t i=0;i<queue.size();++i) for(int cell:owners[edge(queue[i].first,queue[i].second)]) {
      if(used[cell]) continue;
      auto q=s.quads[cell];int start=-1;
      for(int j=0;j<4;++j) if(q[j]==queue[i].first&&q[(j+1)%4]==queue[i].second) start=j;
      if(start<0) continue;
      std::rotate(q.begin(),q.begin()+start,q.end());used[cell]=true;
      Quad encoded;
      for(int j=0;j<4;++j) { if(labels[q[j]]<0) labels[q[j]]=next++;encoded[j]=labels[q[j]]; }
      result.push_back(encoded);
      for(int j=1;j<4;++j) queue.push_back({q[(j+1)%4],q[j]});
    }
    if(result.size()!=s.quads.size()) throw std::runtime_error("Disconnected dual");
    return result;
  }
  std::vector<int> reflection(const State& s) const {
    std::map<Edge,Quad> directed;
    for(auto q:s.quads) for(int j=0;j<4;++j) {
      Quad rotated=q;std::rotate(rotated.begin(),rotated.begin()+j,rotated.end());
      directed[{rotated[0],rotated[1]}]=rotated;
    }
    std::vector<int> action(s.mirror.size(),-1),inverse(s.mirror.size(),-1);
    for(int v=0;v<boundary;++v) { action[v]=(bottom-v+boundary)%boundary;inverse[action[v]]=v; }
    bool changed=true;
    while(changed) {
      changed=false;
      for(auto [e,q]:directed) if(action[e.first]>=0&&action[e.second]>=0) {
        auto target=directed.find({action[e.second],action[e.first]});
        if(target==directed.end()) return {};
        for(int j=2;j<4;++j) {
          int a=q[j],b=target->second[5-j];
          if(action[a]>=0&&action[a]!=b) return {};
          if(inverse[b]>=0&&inverse[b]!=a) return {};
          if(action[a]<0) { action[a]=b;inverse[b]=a;changed=true; }
        }
      }
    }
    for(int v=0;v<int(action.size());++v) if(action[v]<0||action[action[v]]!=v) return {};
    return action;
  }
  void emit(const State& s) {
    int faces=int(s.quads.size()),vertices=int(s.mirror.size());
    if(vertices!=faces+boundary/2+1||!links(s,true)) return;
    auto reflected=symmetry?s.mirror:reflection(s);
    if(reflected.empty()) return;
    std::vector<int> labels;auto cells=canonical(s,labels);
    if(!emitted.insert(cells).second) return;
    std::vector<int> action(vertices);for(int v=0;v<vertices;++v) action[labels[v]]=labels[reflected[v]];
    output<<"{\"quads\":[";
    for(int i=0;i<faces;++i) { if(i) output<<',';output<<'[';
      for(int j=0;j<4;++j) { if(j) output<<',';output<<cells[i][j]; }output<<']'; }
    output<<"],\"mirror\":[";
    for(int i=0;i<vertices;++i) { if(i) output<<',';output<<action[i]; }
    output<<"],\"vertices\":"<<vertices<<",\"boundary_vertices\":"<<boundary
          <<",\"vertical_interior_vertices\":"<<r<<"}\n";
    ++solutions;
  }
  void complete_quad(State s,Quad q,int slot) {
    if(slot==4) { if(insert(s,q)) dfs(s);return; }
    int wanted=s.color[q[slot%2]],existing=int(s.mirror.size());
    for(int v=0;v<existing;++v) if(s.color[v]==wanted&&std::find(q.begin(),q.begin()+slot,v)==q.begin()+slot) {
      q[slot]=v;complete_quad(s,q,slot+1);
    }
    // A new reflection orbit has either one axis vertex or two exchanged vertices.
    for(int size:{1,2}) if((symmetry||size==1)&&existing+size<=limit) {
      // An odd bottom subdivision makes reflection exchange bipartite colors.
      // Such an action cannot fix a vertex.
      if(symmetry&&size==1&&bottom%2) continue;
      State next=s;
      for(int j=0;j<size;++j) { next.mirror.push_back(existing+(size==1?0:1-j));next.color.push_back(wanted^((bottom%2)&&j)); }
      q[slot]=existing;complete_quad(std::move(next),q,slot+1);
    }
  }
  void dfs(const State& s) {
    ++visits;
    if(s.front.empty()) { emit(s);return; }
    if(s.quads.size()>=size_t(cap)) return;
    auto e=s.front.begin()->second;
    complete_quad(s,Quad{e.first,e.second,-1,-1},2);
  }
};
int main(int argc,char** argv) {
  try {
    if(argc<4) throw std::runtime_error("Usage: quad_disk_search MAX_QUADS VERTICAL_INTERIOR_VERTICES OUTPUT.jsonl [--no-orbits] [--bottom-segments N] [--top-segments N]");
    int cap=std::stoi(argv[1]),r=std::stoi(argv[2]),bottom=2,top=2;
    bool orbits=true;
    for(int i=4;i<argc;++i) {
      std::string option=argv[i];
      if(option=="--no-orbits") orbits=false;
      else if((option=="--bottom-segments"||option=="--top-segments")&&i+1<argc) {
        int value=std::stoi(argv[++i]);
        if(value<1||value>12) throw std::runtime_error("Invalid edge subdivision");
        (option=="--bottom-segments"?bottom:top)=value;
      } else throw std::runtime_error("Unknown or incomplete option");
    }
    if(cap<1||cap>12||r<0||(bottom+top)%2||r>cap-(bottom+top)/2) throw std::runtime_error("Invalid bounds");
    std::ofstream out(argv[3]);if(!out) throw std::runtime_error("Cannot open output");
    Search search(r,cap,out,orbits,bottom,top);auto state=search.initial();search.dfs(state);
    std::cerr<<"QUAD_SEARCH_FINISHED max_quads="<<cap<<" vertical_interior_vertices="<<r
             <<" symmetry_orbits="<<search.symmetry<<" visits="<<search.visits<<" candidates="<<search.solutions<<'\n';
    if(!out) throw std::runtime_error("Output failed");
  } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 2; }
}
