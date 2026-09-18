// Independent checks for the compact completion calculations.
#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main
#include <random>
int main() {
  uint64_t complexes=0;
  for(int mask=0;mask<4096;mask++) {
    std::vector<std::vector<int>> edges;
    for(int e=0;e<12;e++)if(mask&(1<<e))edges.push_back({EI[e][0],EI[e][1]});
    int cycle=int(edges.size())-rank_gf2(edges);
    if(edge_cycles.rank[mask]!=cycle)throw std::runtime_error("Edge cycle lookup mismatch");
    for(int fs=0;fs<64;fs++) {
      std::vector<std::vector<int>> faces;bool valid=true;
      for(int f=0;f<6;f++)if(fs&(1<<f)) {
        std::vector<int> column;
        for(int j=0;j<4;j++)for(int e=0;e<12;e++)if(edgekey(FI[f][j],FI[f][(j+1)%4])==edgekey(EI[e][0],EI[e][1])) {
          if(!(mask&(1<<e)))valid=false;
          column.push_back(e);
        }
        faces.push_back(column);
      }
      if(!valid)continue;
      int rank=rank_gf2(faces),beta1=cycle-rank,beta2=int(faces.size())-rank;
      if(beta2!=(fs==63)||beta1!=cycle-int(faces.size())+(fs==63)||beta1<0)
        throw std::runtime_error("Completion cancellation identity failed");
      ++complexes;
    }
  }
  std::mt19937 rng(9249);
  for(int sample=0;sample<512;sample++) {
    State s;int n=sample%65;
    while(int(s.front.size())<n) {
      Quad q;for(int& v:q)v=int(rng()%32);
      auto vertices=sorted(q);if(std::adjacent_find(vertices.begin(),vertices.end())!=vertices.end())continue;
      s.front[key(q)]=q;
    }
    int fast=front_component_upper(s),slow=front_component_upper_reference(s);
    std::map<std::pair<int,int>,int> edges;std::vector<std::vector<int>> faces;
    for(auto [k,q]:s.front) {
      (void)k;std::vector<int> column;
      for(int j=0;j<4;j++){auto e=edgekey(q[j],q[(j+1)%4]);if(!edges.count(e))edges[e]=int(edges.size());column.push_back(edges[e]);}
      faces.push_back(column);
    }
    if(fast!=slow||fast!=n-rank_gf2(faces))throw std::runtime_error("Compact frontier rank mismatch");
  }
  std::cout<<"SPEEDUP_TABLES_OK edge_masks=4096 boundary_subcomplexes="<<complexes<<" random_fronts=512\n";
}
