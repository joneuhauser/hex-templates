// Compare rollback parity constraints with an independent graph coloring check.
#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main

struct Coloring {
  bool bipartite=true;
  std::array<int,6> component,color;
  explicit Coloring(const std::vector<std::pair<int,int>>& edges) {
    component.fill(-1);color.fill(0);
    for(int root=0;root<6;root++)if(component[root]<0) {
      std::vector<int> todo{root};component[root]=root;
      while(!todo.empty()) {
        int v=todo.back();todo.pop_back();
        for(auto [a,b]:edges) {
          int w=a==v?b:b==v?a:-1;if(w<0)continue;
          if(component[w]<0){component[w]=root;color[w]=color[v]^1;todo.push_back(w);}
          else if(color[w]==color[v])bipartite=false;
        }
      }
    }
  }
};
int main() {
  Parity parity;std::vector<std::pair<int,int>> edges,choices;
  for(int a=0;a<6;a++)for(int b=0;b<a;b++)choices.push_back({a,b});
  uint64_t attempts=0,contradictions=0;
  auto check=[&] {
    Coloring expected(edges);if(!expected.bipartite)throw std::runtime_error("Invalid test prefix");
    for(int a=0;a<6;a++)for(int b=0;b<6;b++) {
      auto [ra,pa]=parity.root(a);auto [rb,pb]=parity.root(b);
      bool connected=expected.component[a]==expected.component[b];
      if((ra==rb)!=connected||(connected&&(pa^pb)!=(expected.color[a]^expected.color[b])))
        throw std::runtime_error("Rollback connectivity/parity mismatch");
    }
  };
  std::function<void(size_t)> visit=[&](size_t position) {
    check();if(position==choices.size())return;
    visit(position+1);
    auto checkpoint=parity.checkpoint();auto [a,b]=choices[position];edges.push_back({a,b});
    bool expected=Coloring(edges).bipartite,actual=parity.edge(a,b);++attempts;
    if(actual!=expected)throw std::runtime_error("Odd-cycle result mismatch");
    if(actual) {
      if(!parity.edge(a,b))throw std::runtime_error("Repeated edge rejected");
      visit(position+1);
    } else ++contradictions;
    parity.rollback(checkpoint);edges.pop_back();check();
  };
  visit(0);
  for(int v=0;v<6;v++)if(parity.edge(v,v))throw std::runtime_error("Self-loop accepted");
  State state;state.a.resize(6);
  for(auto [a,b]:std::vector<std::pair<int,int>>{{0,1},{1,2},{3,4}}) {
    state.relation[a*MAXV+b]=state.relation[b*MAXV+a]=1;edges.push_back({a,b});
  }
  parity.initialize(state);check();
  if(parity.checkpoint()!=0)throw std::runtime_error("Initial graph not committed");
  auto checkpoint=parity.checkpoint();
  if(!parity.edge(2,3)||parity.edge(0,4))throw std::runtime_error("Nested parity constraint failed");
  parity.rollback(checkpoint);check();
  std::cout<<"PARITY_OK edge_attempts="<<attempts<<" contradictions="<<contradictions<<'\n';
}
