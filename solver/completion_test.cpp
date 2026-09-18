// Check the component bound against independent matrix rank and actual mesh
// components, including disconnected cavities and cells touching at an edge.
#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main
#include <random>

uint64_t checked=0,disconnected=0;
void verify(const std::vector<Hex>& cells) {
  State state;
  std::map<Quad,std::vector<int>> owners;
  std::map<Quad,Quad> oriented;
  for(int i=0;i<int(cells.size());i++)for(int f=0;f<6;f++) {
    auto q=face(cells[i],f);auto k=key(q);owners[k].push_back(i);oriented[k]=q;
  }
  std::vector<int> parent(cells.size());std::iota(parent.begin(),parent.end(),0);
  auto root=[&](int x){while(parent[x]!=x)x=parent[x];return x;};
  int actual=int(cells.size());
  for(auto& [k,ids]:owners) {
    if(ids.size()==1)state.front[k]=oriented[k];
    else if(ids.size()==2){int a=root(ids[0]),b=root(ids[1]);if(a!=b){parent[a]=b;--actual;}}
    else throw std::runtime_error("Invalid fixture incidence");
  }
  int upper=front_component_upper(state);
  std::map<std::pair<int,int>,int> edges;std::vector<std::vector<int>> columns;
  for(auto [k,q]:state.front) {
    (void)k;std::vector<int> c;
    for(int j=0;j<4;j++){auto e=edgekey(q[j],q[(j+1)%4]);if(!edges.count(e))edges[e]=int(edges.size());c.push_back(edges[e]);}
    columns.push_back(c);
  }
  int independent=int(columns.size())-rank_gf2(columns);
  if(upper!=independent||upper<actual||int(state.front.size())>4*int(cells.size())+2*upper)
    throw std::runtime_error("Completion bound lost a known filling");
  ++checked;disconnected+=actual>1;
}
void verify_orbit_sequences(const std::string& path,bool axial,int mirrors=2) {
  auto mesh=readmesh(path);std::set<int> boundary_vertices;
  for(auto q:mesh.boundary)for(int v:q)boundary_vertices.insert(v);
  int bv=int(boundary_vertices.size());
  for(int v=0;v<bv;v++)if(!boundary_vertices.count(v))throw std::runtime_error("Fixture boundary not first");
  auto action=actions(mesh.points,axial,1<<mirrors);
  std::map<Hex,int> ids;for(int i=0;i<int(mesh.hexes.size());i++)ids[sorted(mesh.hexes[i])]=i;
  std::set<int> unseen;for(int i=0;i<int(mesh.hexes.size());i++)unseen.insert(i);
  std::vector<int> representatives;
  while(!unseen.empty()) {
    int i=*unseen.begin();representatives.push_back(i);
    for(int g=0;g<(1<<mirrors);g++)unseen.erase(ids.at(sorted(reflect(mesh.hexes[i],action,g))));
  }
  std::mt19937 rng(3389);
  for(int sample=0;sample<32;sample++) {
    std::shuffle(representatives.begin(),representatives.end(),rng);
    Search search;search.cap=int(mesh.hexes.size());search.axial=axial;search.mirrors=mirrors;search.component_bound=true;
    search.counts=std::make_unique<Counts[]>(1);search.boundary=mesh;
    search.boundary.points.resize(bv);search.boundary.hexes.clear();
    State state=search.initial();state.a=action;
    for(int i:representatives) {
      Budget budget(state,search.cap);auto placed=state.relation;
      if(!search.completion_feasible(state,mesh.hexes[i],placed,budget)||!search.insert(state,mesh.hexes[i]))
        throw std::runtime_error("Known orbit sequence rejected");
    }
    if(!state.front.empty()||state.cells.size()!=mesh.hexes.size())throw std::runtime_error("Incomplete orbit sequence");
    ++checked;
  }
}
int main(int argc,char** argv) {
  if(argc!=2)throw std::runtime_error("Pass input fixture directory");
  std::string path=argv[1];std::mt19937 rng(9072);
  for(auto name:{"grid2-seed.mesh","slab4-seed.mesh","grid3-seed.mesh","published36-symmetric.mesh","published44-symmetric.mesh"}) {
    auto mesh=readmesh(path+"/"+name);int n=int(mesh.hexes.size());
    int samples=n<=8?(1<<n):256;
    for(int sample=0;sample<samples;sample++) {
      std::vector<Hex> cells;
      for(int i=0;i<n;i++)if(n<=8?bool(sample&(1<<i)):bool(rng()&1))cells.push_back(mesh.hexes[i]);
      verify(cells);
    }
    verify(mesh.hexes);
  }
  // Two cubes touching only on one edge have two unfilled components, despite
  // a connected front. The four-face edge must not collapse those components.
  verify({Hex{0,1,2,3,4,5,6,7},Hex{0,8,9,10,4,11,12,13}});
  verify({Hex{0,1,2,3,4,5,6,7},Hex{0,8,9,10,11,12,13,14}});
  verify({Hex{0,1,2,3,4,5,6,7},Hex{8,9,10,11,12,13,14,15}});
  verify_orbit_sequences(path+"/published36-symmetric.mesh",false);
  verify_orbit_sequences(path+"/published44-symmetric.mesh",true);
  verify_orbit_sequences(path+"/grid2-seed.mesh",false);
  verify_orbit_sequences(path+"/grid2-seed.mesh",true);
  verify_orbit_sequences(path+"/grid3-seed.mesh",false);
  verify_orbit_sequences(path+"/grid3-seed.mesh",true);
  verify_orbit_sequences(path+"/published36-symmetric.mesh",false,1);
  verify_orbit_sequences(path+"/published44-symmetric.mesh",true,1);
  verify_orbit_sequences(path+"/grid3-seed.mesh",false,1);
  verify_orbit_sequences(path+"/grid3-seed.mesh",true,1);
  verify_orbit_sequences(path+"/one-axial-seed.mesh",true,1);
  verify_orbit_sequences(path+"/one-diagonal-seed.mesh",false,1);
  std::cout<<"COMPLETION_BOUND_OK fillings="<<checked<<" disconnected="<<disconnected<<'\n';
}
