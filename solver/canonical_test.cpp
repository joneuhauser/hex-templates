#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main
#include <random>
uint64_t partial_cases=0,prefixes=0,graph_pairs=0,intersection_cases=0,table_cases=0,snapshots=0,groups=0,cache_checks=0,disconnected=0;
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
// Check every tested optimistic successor against an actual remaining cell.
// Labels are allocated in the same order as choose(), including new vertices.
void check_partial(Search& search,const State& state,const Mesh& mesh,const std::vector<int>& labels,Hex cell,Quad q) {
  State prefix=state;auto label=labels;
  std::vector<int> inverse(state.a.size(),-1);
  for(int v=0;v<int(label.size());v++)if(label[v]>=0)inverse[label[v]]=v;
  Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1},source{};
  for(int j=0;j<4;j++)source[j]=inverse[h[j]];
  for(int j=0;j<4;j++) {
    int top=-1;
    for(auto e:EI) {
      int a=cell[e[0]],b=cell[e[1]];if(b==source[j])std::swap(a,b);
      if(a==source[j]&&std::find(source.begin(),source.begin()+4,b)==source.begin()+4)top=b;
    }
    if(top<0)throw std::runtime_error("Missing transverse fixture edge");
    source[j+4]=top;
  }
  auto placed=prefix.relation;auto diagonals=known_diagonals(prefix,search.boundary);
  Budget budget(prefix,search.cap,true);auto domains=budget.initial;
  for(int j=0;j<4;j++)budget.extend(h,j,domains);
  for(int j=4;j<8;j++) {
    if(label[source[j]]<0)label[source[j]]=allocate(prefix,0);
    h[j]=label[source[j]];budget.extend(h,j,domains);
    if(j==7)continue;
    for(int work:{0,1000}) {
      search.partial_work=work;
      if(!search.partial_cover_feasible(prefix,h,j+1,budget,domains,diagonals,placed))
        throw std::runtime_error("Optimistic successor rejected a known remaining cell");
      ++partial_cases;
    }
  }
  search.partial_work=1000;
  (void)mesh;
}

uint64_t canonical_subsets=0,canonical_orbits=0;
State owner_state(Search& search,const std::vector<Hex>& cells) {
  State s=search.initial();s.cells=cells;s.boundary_changed=true;
  if(s.boundary_owner.empty())return s;
  for(int i=0;i<int(cells.size());i++)for(int f=0;f<6;f++) {
    auto it=search.boundary_ids.find(key(face(cells[i],f)));
    if(it!=search.boundary_ids.end())s.boundary_owner[it->second]=i;
  }
  search.update_root_floor(s);return s;
}
void canonicalize_fixture(Search& search,Mesh& mesh) {
  search.face_order=5;search.prepare_boundary_symmetry();
  if(search.boundary_symmetries.size()<2)return;
  const int bv=int(search.boundary.points.size());
  std::vector<std::vector<Hex>> representatives;std::mt19937 rng(15382);
  for(const auto& g:search.boundary_symmetries) {
    auto cells=mesh.hexes;
    for(auto& h:cells){for(int& v:h)if(v<bv)v=g.permutation[v];if(g.odd){std::swap(h[1],h[3]);std::swap(h[5],h[7]);}}
    bool feasible=search.boundary_canonical_feasible(owner_state(search,cells));
    auto relabeled=cells;std::vector<int> labels(mesh.points.size());std::iota(labels.begin(),labels.end(),0);
    std::shuffle(labels.begin()+bv,labels.end(),rng);std::shuffle(relabeled.begin(),relabeled.end(),rng);
    for(auto& h:relabeled)for(int& v:h)v=labels[v];
    if(feasible!=search.boundary_canonical_feasible(owner_state(search,relabeled)))throw std::runtime_error("Canonicalization depends on interior labels or cell order");
    if(feasible)representatives.push_back(cells);
  }
  if(representatives.empty())throw std::runtime_error("Canonicalization discarded an entire boundary-isometry orbit");
  for(const auto& cells:representatives)for(int trial=0;trial<128;trial++) {
    std::vector<Hex> subset;for(auto h:cells)if(rng()%2)subset.push_back(h);
    auto state=owner_state(search,subset);
    if(!search.boundary_canonical_feasible(state))throw std::runtime_error("Unknown-owner comparison rejected a canonical completion");
    for(auto h:cells)if(!search.leader_feasible(state,h))throw std::runtime_error("Leader bound rejected canonical completion");
    ++canonical_subsets;
  }
  mesh.hexes=representatives.front();++canonical_orbits;
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
  canonicalize_fixture(search,mesh);
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
      if(state.cover_graph&&snapshots%4==0) {
        search.cover_cache=0;
        const auto& graph=*state.cover_graph;
        for(int i=0;i<int(graph.faces.size());i++)for(int j=0;j<i;j++) {
          bool expected=search.faces_can_share(state,graph.faces[i],graph.faces[j],diagonals);
          if(expected!=graph.compatible[i][j])throw std::runtime_error("Incremental/reflected cover graph differs from recomputation");
          ++graph_pairs;
        }
        search.cover_cache=1;
      }
      bool partial_checked=false;
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
          if(!open.empty()&&!partial_checked) {
            check_partial(search,state,mesh,label,mesh.hexes[id],open[0]);partial_checked=true;
          }
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
      Budget budget(state,search.cap,true);VertexDomains domains(state);auto placed=state.relation;
      for(int& v:h) {
        if(label[v]<0) {
          bool a=action[v][1]==v,b=state.group_size==4&&action[v][2]==v;
          int type=state.group_size==1?0:state.group_size==2?(a?1:0):(a?(b?3:1):(b?2:0));int base=allocate(state,type);
          for(int g=0;g<state.group_size;g++){int old=action[v][g],now=state.a[base][g];if(label[old]>=0&&label[old]!=now)throw std::runtime_error("Label conflict");label[old]=now;}
        }
        v=label[v];
      }
      if(h[0]<budget.vertices)for(int assigned=1;assigned<=8;assigned++)for(bool use_domains:{false,true}) {
        if(!edge_completion_feasible(state,h,assigned,placed,budget,Overlaps{},domains,use_domains))throw std::runtime_error("Zero-mirror prefix bound rejected a known future cell");
        ++prefixes;
      }
      if(!search.insert(state,h))throw std::runtime_error("Known orbit insertion failed: "+path);
      remaining.erase(selected);
    }
    if(!state.front.empty()||state.cells.size()!=mesh.hexes.size())throw std::runtime_error("Incomplete known fixture");
  }
}
void check_intersections() {
  std::vector<uint32_t> full{0};
  for(int i=0;i<8;i++)for(int j=0;j<8;j++)full.push_back(uint32_t(j+1)<<(4*i));
  for(const auto& a:EI)for(const auto& b:EI)for(int flip=0;flip<2;flip++)
    full.push_back((uint32_t(b[flip]+1)<<(4*a[0]))|(uint32_t(b[1-flip]+1)<<(4*a[1])));
  for(const auto& a:FI)for(const auto& b:FI)for(int r=0;r<4;r++) {
    uint32_t code=0;for(int k=0;k<4;k++)code|=uint32_t(b[(r-k+4)%4]+1)<<(4*a[k]);full.push_back(code);
  }
  std::mt19937 rng(92755);std::array<int,16> vertices;std::iota(vertices.begin(),vertices.end(),0);
  for(int trial=0;trial<20000;trial++) {
    std::shuffle(vertices.begin(),vertices.end(),rng);
    Hex a{0,1,2,3,4,5,6,7},b;
    unsigned ma=trial<1000?255:rng()%256,mb=trial<1000?255:rng()%256;uint32_t actual=0;
    for(int i=0;i<8;i++){if(!(ma&(1u<<i)))a[i]=-1;b[i]=(mb&(1u<<i))?vertices[i]:-1;}
    for(int i=0;i<8;i++)if(a[i]>=0)for(int j=0;j<8;j++)if(a[i]==b[j])actual|=uint32_t(j+1)<<(4*i);
    bool expected=false;
    for(auto code:full) {
      uint32_t projected=0;
      for(int i=0;i<8;i++)if(ma&(1u<<i)){int j=int((code>>(4*i))&15)-1;if(j>=0&&(mb&(1u<<j)))projected|=uint32_t(j+1)<<(4*i);}
      if(projected==actual){expected=true;break;}
    }
    if(partial_cells_compatible(a,b)!=expected)throw std::runtime_error("Partial intersection differs from exhaustive legal-map projection");
    ++intersection_cases;
  }
}

int main(int argc,char** argv) {
  if(argc!=2)throw std::runtime_error("Pass input fixture directory");
  State single;single.group_size=1;
  for(int i=0;i<20;i++){int v=allocate(single,0);if(v!=i||single.a.size()!=size_t(i+1)||single.a[i]!=Action{i,-1,-1,-1})throw std::runtime_error("Identity allocation changed more than one vertex");}
  check_intersections();check_tables();std::string p=argv[1];
  for(int mirrors:{0}) {
    fixture(p+"/published36-symmetric.mesh",false,mirrors);
    fixture(p+"/published44-symmetric.mesh",true,mirrors);
    for(bool axial:{false,true}) {
      fixture(p+"/grid2-seed.mesh",axial,mirrors);
      fixture(p+"/grid3-seed.mesh",axial,mirrors);
      fixture(p+"/slab4-seed.mesh",axial,mirrors);
    }
  }
  std::cout<<"CANONICAL_TEST_OK partial_cases="<<partial_cases<<" canonical_subsets="<<canonical_subsets<<" canonical_orbits="<<canonical_orbits<<" prefixes="<<prefixes<<" graph_pairs="<<graph_pairs<<" intersection_cases="<<intersection_cases<<" masked_cases="<<table_cases<<" known_prefixes="<<snapshots<<" known_cell_groups="<<groups<<" cache_pairs="<<cache_checks<<" disconnected_remainders="<<disconnected<<'\n';
  return 0;
}
