#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main
#include <random>
uint64_t table_cases=0,snapshots=0,groups=0,cache_checks=0,disconnected=0;
void check_tables() {
  auto verify=[](const auto& full,const auto& tables) {
    for(int mask=0;mask<256;mask++) {
      std::map<uint32_t,uint8_t> expected;
      for(auto [code,flags]:full) {
        uint32_t projected=0;
        for(int i=0;i<8;i++) {
          int target=int((code>>(4*i))&15)-1;
          if(target>=0&&(mask&(1<<i))&&(mask&(1<<target)))projected|=uint32_t(target+1)<<(4*i);
        }
        expected[projected]|=flags;
      }
      if(expected.size()!=tables[mask].size())throw std::runtime_error("Masked table size mismatch");
      for(auto [code,flags]:expected){if(tables[mask].at(code)!=flags)throw std::runtime_error("Masked table mismatch");++table_cases;}
    }
  };
  for(int odd=0;odd<2;odd++)verify(overlap_patterns.prefix[odd][8],overlap_patterns.masked[odd]);
  verify(overlap_patterns.reflection_prefix[8],overlap_patterns.reflection_masked);
  verify(overlap_patterns.fixed_prefix[8],overlap_patterns.fixed_masked);
  std::unordered_map<uint32_t,uint8_t> exchanged;
  for(auto [code,flags]:overlap_patterns.reflection_prefix[8])if(!(flags&2))exchanged[code]=flags;
  verify(exchanged,overlap_patterns.exchanged_masked);
}
void fixture(const std::string& path,bool axial,int mirrors) {
  Mesh mesh=readmesh(path);std::set<int> boundary_vertices;
  for(auto q:mesh.boundary)for(int v:q)boundary_vertices.insert(v);
  int bv=int(boundary_vertices.size());for(int v=0;v<bv;v++)if(!boundary_vertices.count(v))throw std::runtime_error("Boundary not first");
  std::map<Quad,int> owners;std::vector<std::pair<int,int>> dual_edges;
  for(int i=0;i<int(mesh.hexes.size());i++)for(int f=0;f<6;f++){auto q=key(face(mesh.hexes[i],f));auto it=owners.find(q);if(it==owners.end())owners[q]=i;else dual_edges.push_back({it->second,i});}
  auto action=actions(mesh.points,axial,1<<mirrors);
  std::map<Hex,int> ids;for(int i=0;i<int(mesh.hexes.size());i++)ids[sorted(mesh.hexes[i])]=i;
  std::set<int> ungrouped;for(int i=0;i<int(mesh.hexes.size());i++)ungrouped.insert(i);
  std::vector<std::vector<int>> orbits;
  while(!ungrouped.empty()) {
    int i=*ungrouped.begin();std::set<int> group;
    for(int g=0;g<(1<<mirrors);g++)group.insert(ids.at(sorted(reflect(mesh.hexes[i],action,g))));
    orbits.emplace_back(group.begin(),group.end());for(int id:group)ungrouped.erase(id);
  }
  Search search;search.cap=int(mesh.hexes.size());search.mirrors=mirrors;search.axial=axial;
  search.boundary=mesh;search.boundary.points.resize(bv);search.boundary.hexes.clear();search.counts=std::make_unique<Counts[]>(1);
  search.geometry.initialize(search.boundary.points,axial,1<<mirrors,true,true);
  std::mt19937 rng(7193);
  for(int sample=0;sample<8;sample++) {
    State state=search.initial();std::vector<int> label(mesh.points.size(),-1);for(int v=0;v<bv;v++)label[v]=v;
    std::set<int> remaining;for(int i=0;i<int(orbits.size());i++)remaining.insert(i);
    while(!remaining.empty()) {
      std::vector<bool> active(mesh.hexes.size());std::vector<int> parent(mesh.hexes.size());std::iota(parent.begin(),parent.end(),0);
      int components=0;for(int orbit:remaining)for(int id:orbits[orbit]){active[id]=true;++components;}
      auto root=[&](int v){while(parent[v]!=v)v=parent[v];return v;};
      for(auto [a,b]:dual_edges)if(active[a]&&active[b]){a=root(a);b=root(b);if(a!=b){parent[a]=b;--components;}}
      disconnected+=components>1;
      if(!search.cover_feasible(state))throw std::runtime_error("Cover bound rejected known filling: "+path);
      ++snapshots;
      Parity parity;parity.initialize(state);auto [anchor,p]=parity.root(0);
      for(int v=0;v<int(state.a.size());v++){auto [root,q]=parity.root(v);if(root!=anchor)throw std::runtime_error("Disconnected fixture prefix");state.cover_color[v]=p^q;}
      auto diagonals=known_diagonals(state,search.boundary);
      std::vector<int> choices;
      for(int orbit:remaining) {
        bool touches=false;
        for(int id:orbits[orbit]) {
          std::vector<Quad> open;
          for(int f=0;f<6;f++) {
            auto q=face(mesh.hexes[id],f);bool known=true;
            for(int& v:q){v=label[v];known&=v>=0;}
            if(known&&state.front.count(key(q))) {
              if(!same_orientation(q,state.front.at(key(q))))throw std::runtime_error("Fixture face orientation");
              open.push_back(q);touches=true;
            }
          }
          std::vector<Hex> witness;
          for(auto q:open) {
            witness=search.extend_cube_group(state,witness,q,diagonals);
            if(witness.empty())throw std::runtime_error("Known cell lost its cover witness: "+path);
          }
          if(!open.empty())++groups;
          for(int i=0;i<int(open.size());i++)for(int j=0;j<i;j++) {
            search.cover_cache=0;auto uncached=search.face_templates(state,open[i],open[j]);
            search.cover_cache=1;auto cached=search.face_templates(state,open[i],open[j]);
            if(!cached||cached!=uncached)throw std::runtime_error("Known pair/cache mismatch");
            if(!search.faces_can_share(state,open[j],open[i],diagonals))throw std::runtime_error("Asymmetric known pair");
            ++cache_checks;
          }
        }
        if(touches)choices.push_back(orbit);
      }
      if(choices.empty())throw std::runtime_error("No attachable known orbit");
      int selected=choices[rng()%choices.size()];Hex h=mesh.hexes[orbits[selected][0]];
      for(int& v:h) {
        if(label[v]<0) {
          bool a=action[v][1]==v,b=state.group_size==4&&action[v][2]==v;
          int type=state.group_size==2?(a?1:0):(a?(b?3:1):(b?2:0));int base=allocate(state,type);
          for(int g=0;g<state.group_size;g++){int old=action[v][g],now=state.a[base][g];if(label[old]>=0&&label[old]!=now)throw std::runtime_error("Label conflict");label[old]=now;}
        }
        v=label[v];
      }
      if(!search.insert(state,h))throw std::runtime_error("Known orbit insertion failed: "+path);
      remaining.erase(selected);
    }
    if(!state.front.empty()||state.cells.size()!=mesh.hexes.size())throw std::runtime_error("Incomplete known fixture");
  }
}
int main(int argc,char** argv) {
  if(argc!=2)throw std::runtime_error("Pass input fixture directory");
  check_tables();std::string p=argv[1];
  for(int mirrors:{1,2}) {
    fixture(p+"/published36-symmetric.mesh",false,mirrors);
    fixture(p+"/published44-symmetric.mesh",true,mirrors);
    for(bool axial:{false,true}) {
      fixture(p+"/grid2-seed.mesh",axial,mirrors);
      fixture(p+"/grid3-seed.mesh",axial,mirrors);
      fixture(p+"/slab4-seed.mesh",axial,mirrors);
    }
  }
  std::cout<<"COVER_TEST_OK masked_cases="<<table_cases<<" known_prefixes="<<snapshots<<" known_cell_groups="<<groups<<" cache_pairs="<<cache_checks<<" disconnected_remainders="<<disconnected<<'\n';
  return 0;
}
