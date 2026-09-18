// Reflection-orbit hexahedral enumeration. SPDX-License-Identifier: MIT
#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <csignal>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>

using Quad=std::array<int,4>;
using Hex=std::array<int,8>;
using Vec=std::array<double,3>;
using Action=std::array<int,4>;
using Diagonals=std::map<std::pair<int,int>,Quad>;
using Rules=std::vector<std::array<int,3>>; // source top index, target top index, group element
#ifndef HEX_NATIVE_ARCH
#define HEX_NATIVE_ARCH 0
#endif
constexpr int MAXV=112;
constexpr int FI[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
constexpr int EI[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
constexpr int BITS[8]={0,1,3,2,4,5,7,6};
std::atomic<bool> interrupted{false};
static_assert(std::atomic<bool>::is_always_lock_free,"Signal handling requires lock-free atomic bool");
extern "C" void on_signal(int) { interrupted.store(true,std::memory_order_relaxed); }

Quad face(const Hex& h,int f) { return {h[FI[f][0]],h[FI[f][1]],h[FI[f][2]],h[FI[f][3]]}; }
Quad reverse(Quad q) { std::swap(q[1],q[3]);return q; }
template<size_t N> std::array<int,N> sorted(std::array<int,N> a) { std::sort(a.begin(),a.end());return a; }
Quad oriented_key(Quad q) {
  Quad best=q;
  for(int k=1;k<4;k++) { Quad r;for(int j=0;j<4;j++)r[j]=q[(j+k)%4];best=std::min(best,r); }
  return best;
}
Quad key(Quad q) { return std::min(oriented_key(q),oriented_key(reverse(q))); }
bool same_orientation(Quad a,Quad b) { return oriented_key(a)==oriented_key(b); }
std::pair<int,int> edgekey(int a,int b) { return std::minmax(a,b); }
bool isedge(const Hex& h,int a,int b) {
  for(auto& e:EI)if(edgekey(h[e[0]],h[e[1]])==edgekey(a,b))return true;
  return false;
}
bool isedge(const Quad& q,int a,int b) {
  for(int i=0;i<4;i++)if(edgekey(q[i],q[(i+1)%4])==edgekey(a,b))return true;
  return false;
}
std::vector<int> intersection(const Hex& a,const Hex& b) {
  std::vector<int> r;for(int x:a)if(std::find(b.begin(),b.end(),x)!=b.end())r.push_back(x);return r;
}
bool compatible(const Hex& a,const Hex& b) {
  auto v=intersection(a,b);
  if(v.size()<2)return true;
  if(v.size()==2)return isedge(a,v[0],v[1])&&isedge(b,v[0],v[1]);
  if(v.size()!=4)return false;
  Quad target=sorted(Quad{v[0],v[1],v[2],v[3]});
  for(int i=0;i<6;i++)if(sorted(face(a,i))==target)
    for(int j=0;j<6;j++)if(sorted(face(b,j))==target)
      return same_orientation(face(a,i),reverse(face(b,j)));
  return false;
}
bool compatible_boundary(const Hex& h,const Quad& q) {
  std::vector<int> v;for(int x:q)if(std::find(h.begin(),h.end(),x)!=h.end())v.push_back(x);
  if(v.size()<2)return true;
  if(v.size()==2)return isedge(h,v[0],v[1])&&isedge(q,v[0],v[1]);
  if(v.size()!=4)return false;
  for(int i=0;i<6;i++)if(sorted(face(h,i))==sorted(q))return same_orientation(face(h,i),q);
  return false;
}

// For each involution, the common vertices of a cell and its image must form
// an empty set, vertex, edge, face, or the entire coincident cell. Enumerate
// every legal involutive mapping on these sets, then index its prefix views.
// A zero nibble means that the image has not appeared in the assigned prefix.
struct Overlaps {
  std::array<uint32_t,3> code{};
  std::array<bool,3> fixed{};
  std::array<uint8_t,3> possible{7,7,7}; // legal, coincident cell, shared face
};
struct OverlapPatterns {
  std::array<std::array<std::unordered_map<uint32_t,uint8_t>,9>,2> prefix;
  std::array<std::unordered_map<uint32_t,uint8_t>,9> reflection_prefix,fixed_prefix;
  std::array<std::array<std::unordered_map<uint32_t,uint8_t>,256>,2> masked;
  std::array<std::unordered_map<uint32_t,uint8_t>,256> reflection_masked,fixed_masked,exchanged_masked;
  std::array<std::unordered_map<uint32_t,uint8_t>,256> rotation_fixed_masked;
  OverlapPatterns() {
    std::vector<std::vector<int>> domains(1);
    for(int i=0;i<8;i++)domains.push_back({i});
    for(auto& e:EI)domains.push_back({e[0],e[1]});
    for(auto& f:FI)domains.push_back({f[0],f[1],f[2],f[3]});
    domains.push_back({0,1,2,3,4,5,6,7});
    Hex original{0,1,2,3,4,5,6,7};
    for(auto domain:domains) {
      std::sort(domain.begin(),domain.end());auto permutation=domain;
      do {
        std::array<int,8> mapping;mapping.fill(-1);
        for(size_t i=0;i<domain.size();i++)mapping[domain[i]]=permutation[i];
        bool involution=true;
        for(int i:domain)if(mapping[mapping[i]]!=i)involution=false;
        if(!involution)continue;
        for(int odd=0;odd<2;odd++) {
          Hex other;for(int i=0;i<8;i++)other[i]=mapping[i]<0?8+i:mapping[i];
          if(odd){std::swap(other[1],other[3]);std::swap(other[5],other[7]);}
          if(domain.size()==8) {
            // Compare oriented faces, including the cell's cube connectivity.
            std::array<Quad,6> a,b;
            for(int f=0;f<6;f++){a[f]=oriented_key(face(original,f));b[f]=oriented_key(face(other,f));}
            std::sort(a.begin(),a.end());std::sort(b.begin(),b.end());if(a!=b)continue;
          } else if(!compatible(original,other))continue;
          for(int n=0;n<=8;n++) {
            uint32_t code=0;
            for(int i=0;i<n;i++)if(mapping[i]>=0&&mapping[i]<n)code|=uint32_t(mapping[i]+1)<<(4*i);
            uint8_t flags=uint8_t(1|(domain.size()==8?2:0)|(domain.size()==4?4:0));
            prefix[odd][n][code]|=flags;
            if(odd) {
              bool pointwise=true;for(int i:domain)pointwise&=mapping[i]==i;
              if(domain.size()==8||pointwise)reflection_prefix[n][code]|=flags;
              if(domain.size()==8)fixed_prefix[n][code]|=flags;
            }
          }
        }
      } while(std::next_permutation(permutation.begin(),permutation.end()));
    }
    auto project=[](const auto& full,auto& tables) {
      for(auto [code,flags]:full)for(int mask=0;mask<256;mask++) {
        uint32_t partial=0;
        for(int i=0;i<8;i++)if(mask&(1<<i)){int j=int((code>>(4*i))&15)-1;if(j>=0&&(mask&(1<<j)))partial|=uint32_t(j+1)<<(4*i);}
        tables[mask][partial]|=flags;
      }
    };
    for(int odd=0;odd<2;odd++)project(prefix[odd][8],masked[odd]);
    project(reflection_prefix[8],reflection_masked);project(fixed_prefix[8],fixed_masked);
    std::unordered_map<uint32_t,uint8_t> rotation_fixed;
    for(auto [code,flags]:prefix[0][8])if(flags&2)rotation_fixed[code]=flags;
    project(rotation_fixed,rotation_fixed_masked);
    std::unordered_map<uint32_t,uint8_t> exchanged;for(auto [code,flags]:reflection_prefix[8])if(!(flags&2))exchanged[code]=flags;project(exchanged,exchanged_masked);
  }
};
const OverlapPatterns overlap_patterns;

struct Mesh { std::vector<Vec> points;std::vector<Quad> boundary;std::vector<Hex> hexes; };
Mesh readmesh(const std::string& path) {
  std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read "+path);
  Mesh m;std::string word;
  while(f>>word) {
    int n,tag;
    if(word=="MeshVersionFormatted"||word=="Dimension") { f>>n; }
    else if(word=="Vertices") { f>>n;for(int i=0;i<n;i++){Vec p;f>>p[0]>>p[1]>>p[2]>>tag;m.points.push_back(p);} }
    else if(word=="Quadrilaterals") { f>>n;for(int i=0;i<n;i++){Quad q;for(int& x:q){f>>x;--x;}f>>tag;m.boundary.push_back(q);} }
    else if(word=="Hexahedra") { f>>n;for(int i=0;i<n;i++){Hex h;for(int& x:h){f>>x;--x;}f>>tag;m.hexes.push_back(h);} }
    else if(word=="End")break;
    else throw std::runtime_error("Unsupported mesh keyword: "+word);
    if(!f)throw std::runtime_error("Malformed mesh: "+path);
  }
  if(m.points.empty()||m.boundary.empty())throw std::runtime_error("Empty mesh boundary");
  return m;
}
Vec transform(Vec p,int g,bool axial,bool half_turn=false) {
  if(half_turn) { if(g&1){p[0]=-p[0];p[1]=-p[1];} }
  else if(axial) { if(g&1)p[0]=-p[0];if(g&2)p[1]=-p[1]; }
  else { if(g&1)std::swap(p[0],p[1]);if(g&2){double x=p[0];p[0]=-p[1];p[1]=-x;} }
  return p;
}
std::vector<Action> actions(const std::vector<Vec>& p,bool axial,int group_size=4,bool half_turn=false) {
  if(half_turn&&group_size!=2)throw std::runtime_error("Half-turn action requires a group of order two");
  std::vector<Action> a(p.size(),Action{-1,-1,-1,-1});
  for(size_t i=0;i<p.size();i++)for(int g=0;g<group_size;g++) {
    Vec target=transform(p[i],g,axial,half_turn);int found=-1;
    for(size_t j=0;j<p.size();j++) {
      double d=0;for(int k=0;k<3;k++)d=std::max(d,std::abs(target[k]-p[j][k]));
      if(d<1e-12){if(found>=0)throw std::runtime_error("Ambiguous reflected vertex");found=int(j);}
    }
    if(found<0)throw std::runtime_error("Input coordinates do not have the requested reflections");
    a[i][g]=found;
  }
  for(size_t i=0;i<p.size();i++)for(int g=0;g<group_size;g++)for(int k=0;k<group_size;k++)
    if(a[a[i][g]][k]!=a[i][g^k])throw std::runtime_error("Invalid group action");
  return a;
}
Hex reflect(Hex h,const std::vector<Action>& a,int g,bool half_turn=false) {
  for(int& x:h)x=a[x][g];
  if(!half_turn&&(g==1||g==2)){std::swap(h[1],h[3]);std::swap(h[5],h[7]);}
  return h;
}
std::array<Quad,6> hexkey(const Hex& h) {
  std::array<Quad,6> r;for(int i=0;i<6;i++)r[i]=oriented_key(face(h,i));std::sort(r.begin(),r.end());return r;
}

int reflection_count(int group_size,bool half_turn=false){return half_turn?0:group_size==4?2:group_size==2?1:0;}

struct CoverGraph {
  uint64_t canonical_floor=0;
  std::vector<Quad> faces,keys;
  std::array<std::bitset<256>,256> compatible;
  std::array<uint8_t,MAXV*MAXV> relation;
};
enum AuditCounter { AC_BRANCHES, AC_CACHE_HITS, AC_CACHE_MISSES, AC_CANDIDATES, AC_CANONICAL_REJECTIONS, AC_COMPONENT_REJECTIONS, AC_COVER_REJECTIONS, AC_EARLY_BUDGET, AC_EARLY_COMPLETION, AC_EARLY_DOMAIN, AC_EARLY_GEOMETRY, AC_EARLY_OVERLAP, AC_EARLY_PARITY, AC_PAIR_REUSES, AC_PARTIAL_COMPLETION, AC_PARTIAL_COVER_CALLS, AC_PARTIAL_COVER_REJECTIONS, AC_QUICK_PACKING_REJECTIONS, AC_REJECT_LOCATION, AC_REJECT_SIDE, AC_REJECTED, AC_STATES, AC_TRIALS, AC_COVER_WORK, AC_COVER_UNKNOWN, AC_INSERT_ORBIT_SHAPE, AC_INSERT_CELL_CAP, AC_INSERT_BOUNDARY, AC_INSERT_PLACED_CELLS, AC_INSERT_ORBIT_CELLS, AC_INSERT_FACE_INCIDENCE, AC_INSERT_FACE_BUDGET, AC_INSERT_ODD_CYCLE, AC_RULE_CONFLICT, AC_FORCED_CHOICE, AC_VERTEX_CAP, AC_RELATIONS, AC_FACE_CONSISTENCY, AC_INITIAL_OVERLAP, AC_INITIAL_BUDGET, AC_CHILDREN, AC_COUNT };
using AuditCounters=std::array<uint64_t,AC_COUNT>;
struct State {
  AuditCounters* audit=nullptr;uint64_t audit_id=0;uint8_t audit_cover_reason=0;
  uint64_t canonical_floor=0;
  int group_size=4;
  bool half_turn=false;
  bool odd(int g)const{return !half_turn&&(g==1||g==2);}
  bool symmetry_parity=false,mirror_separation=false;
  std::array<std::array<int8_t,MAXV>,2> boundary_side;
  std::array<int,4> symmetry_flip{};
  std::vector<Action> a;
  std::vector<Hex> cells;
  std::vector<int> boundary_owner;bool boundary_changed=false;
  std::vector<std::pair<int,int>> locations;
  std::map<Quad,Quad> front; // canonical cyclic key -> required outward orientation
  std::map<Quad,int> incidence;
  std::array<uint8_t,MAXV*MAXV> relation{};
  std::bitset<256> cover_packing;
  std::shared_ptr<const CoverGraph> cover_graph;
  Quad cover_choice{}; bool cover_choice_valid=false;
  std::array<uint8_t,MAXV> cover_color{};
  std::array<std::pair<int,int>,MAXV> cover_side;
  int orbit_count=0,component_upper=1;
};
// Each remaining face-connected cell component with h cells has at most
// 4h+2 boundary faces. Component boundary chains are independent cycles of
// the front's edge/face incidence matrix. Its nullity therefore bounds the
// number of remaining components from above, even for a nonmanifold front.
int front_component_upper_reference(const State& s) {
  int n=int(s.front.size());if(n==0)return 0;
  if(n>256)throw std::runtime_error("Front component capacity exceeded");
  std::map<std::pair<int,int>,std::vector<int>> incidence;
  int index=0;for(auto [key,q]:s.front){(void)key;for(int j=0;j<4;j++)incidence[edgekey(q[j],q[(j+1)%4])].push_back(index);++index;}
  std::array<int,256> parent;std::iota(parent.begin(),parent.end(),0);
  auto root=[&](int v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
  int components=n;
  for(auto& [edge,faces]:incidence)if(faces.size()==2) {
    (void)edge;int a=root(faces[0]),b=root(faces[1]);if(a!=b){parent[a]=b;--components;}
  }
  std::array<std::bitset<256>,256> pivots{};int rank=0;
  for(auto& [edge,faces]:incidence)if(faces.size()!=2) {
    (void)edge;std::bitset<256> row;for(int f:faces)row.flip(root(f));
    for(int j=n-1;j>=0;j--)if(row[j]) {
      if(pivots[j].none()){pivots[j]=row;++rank;break;}row^=pivots[j];
    }
  }
  return components-rank;
}
// Compact edge incidence lists avoid allocating a map node and vector
// for every edge of every accepted frontier. The GF(2) calculation is the
// same as the reference implementation, including edges with >2 faces.
int front_component_upper(const State& s) {
  int n=int(s.front.size());if(n==0)return 0;
  if(n>256)throw std::runtime_error("Front component capacity exceeded");
  std::array<int16_t,MAXV*MAXV> head;head.fill(-1);
  std::array<int16_t,1024> next,edge_ids;
  int index=0,used=0;
  for(const auto& entry:s.front) {
    const auto& q=entry.second;
    for(int j=0;j<4;j++) {
      auto [a,b]=edgekey(q[j],q[(j+1)%4]);int e=a*MAXV+b;
      if(head[e]<0)edge_ids[used++]=e;
      next[index]=head[e];head[e]=index++;
    }
  }
  std::array<int,256> parent;std::iota(parent.begin(),parent.end(),0);
  auto root=[&](int v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
  int components=n;
  for(int k=0;k<used;k++) {
    int a=head[edge_ids[k]],b=next[a];
    if(b>=0&&next[b]<0){a=root(a/4);b=root(b/4);if(a!=b){parent[a]=b;--components;}}
  }
  std::array<std::bitset<256>,256> pivots{};int rank=0;
  for(int k=0;k<used;k++) {
    int a=head[edge_ids[k]],b=next[a];if(b>=0&&next[b]<0)continue;
    std::bitset<256> row;for(int p=a;p>=0;p=next[p])row.flip(root(p/4));
    for(int j=n-1;j>=0;j--)if(row[j]) {
      if(pivots[j].none()){pivots[j]=row;++rank;break;}row^=pivots[j];
    }
  }
  return components-rank;
}
void audit_add(const State& s,int counter,uint64_t amount=1){if(s.audit)(*s.audit)[counter]+=amount;}

struct ChoiceMask {
  std::array<uint64_t,2> word{~uint64_t(0),~uint64_t(0)};
  void reset(int v){word[v/64]&=~(uint64_t(1)<<(v%64));}
  void set(int v){word[v/64]|=uint64_t(1)<<(v%64);}
  void intersect(const ChoiceMask& other){word[0]&=other.word[0];word[1]&=other.word[1];}
  void limit(int n){if(n<64){word[0]&=(uint64_t(1)<<n)-1;word[1]=0;}else word[1]&=(uint64_t(1)<<(n-64))-1;}
  int pop(){int w=word[0]?0:1;if(!word[w])return -1;int bit=__builtin_ctzll(word[w]);word[w]&=word[w]-1;return w*64+bit;}
  int count()const{return __builtin_popcountll(word[0])+__builtin_popcountll(word[1]);}
};

// Exact signs for fixed boundary corners, with one unknown point relaxed to
// the convex hull of all boundary vertices. Doubles are treated as exact dyadics.
struct CornerGeometry {
  using Integer=boost::multiprecision::cpp_int;
  int boundary_count=0,side=0;
  std::vector<uint8_t> allowed;
  std::vector<ChoiceMask> domains;
  std::vector<int> triple_plane;
  std::vector<std::vector<Integer>> plane_values;
  mutable std::unique_ptr<std::atomic<uint8_t>[]> pair_cache;
  std::array<std::array<int,4>,8> corners{};
  size_t index(const std::array<int,4>& q) const {
    return ((size_t(q[0])*side+q[1])*side+q[2])*side+q[3];
  }
  void initialize(const std::vector<Vec>& points,bool axial=false,int group_size=4,bool build_domains=true,bool build_locations=false,bool half_turn=false) {
    boundary_count=int(points.size());side=boundary_count+1;
    allowed.assign(size_t(side)*side*side*side,1);
    int exponent=0;
    for(auto p:points)for(double x:p)if(x!=0){int e;std::frexp(x,&e);exponent=std::min(exponent,e-53);}
    std::vector<std::array<Integer,3>> exact(points.size());
    for(int i=0;i<boundary_count;i++)for(int k=0;k<3;k++) {
      double x=points[i][k];if(x==0)continue;
      int e;double m=std::frexp(x,&e);
      exact[i][k]=Integer(int64_t(std::ldexp(m,53)))<<(e-53-exponent);
    }
    auto orientation=[&](int a,int b,int c,int d) {
      std::array<Integer,3> u,v,w;
      for(int k=0;k<3;k++){u[k]=exact[b][k]-exact[a][k];v[k]=exact[c][k]-exact[a][k];w[k]=exact[d][k]-exact[a][k];}
      Integer det=u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0]);
      return det>0?1:det<0?-1:0;
    };
    for(int a=0;a<boundary_count;a++)for(int b=a+1;b<boundary_count;b++)
      for(int c=b+1;c<boundary_count;c++)for(int d=c+1;d<=boundary_count;d++) {
        bool positive=false,negative=false;
        if(d==boundary_count) {
          for(int v=0;v<boundary_count&&!(positive&&negative);v++) {
            int sign=orientation(a,b,c,v);positive|=sign>0;negative|=sign<0;
          }
        } else {int sign=orientation(a,b,c,d);positive=sign>0;negative=sign<0;}
        std::array<int,4> ids{a,b,c,d},perm{0,1,2,3};
        do {
          int inversions=0;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)inversions+=perm[i]>perm[j];
          std::array<int,4> q;for(int i=0;i<4;i++)q[i]=ids[perm[i]];
          allowed[index(q)]=(inversions&1)?negative:positive;
        } while(std::next_permutation(perm.begin(),perm.end()));
      }
    for(int i=0;i<8;i++) {
      corners[i][0]=i;
      for(int axis=0;axis<3;axis++)for(int j=0;j<8;j++)if(BITS[j]==(BITS[i]^(1<<axis)))corners[i][axis+1]=j;
      if(__builtin_popcount(unsigned(BITS[i]))&1)std::swap(corners[i][1],corners[i][2]);
    }
    if(build_domains)this->build_domains(points,axial,group_size,half_turn);
    if(build_locations)initialize_locations(exact);
  }
  void initialize_locations(const std::vector<std::array<Integer,3>>& points) {
    if(boundary_count>24)return; // Conservative fallback for large generic inputs.
    triple_plane.assign(size_t(boundary_count)*boundary_count*boundary_count,-1);
    std::map<std::array<Integer,4>,int> ids;
    auto gcd=[](Integer a,Integer b){if(a<0)a=-a;if(b<0)b=-b;while(b!=0){Integer r=a%b;a=b;b=r;}return a;};
    for(int a=0;a<boundary_count;a++)for(int b=a+1;b<boundary_count;b++)for(int c=b+1;c<boundary_count;c++) {
      std::array<Integer,3> u,v;for(int k=0;k<3;k++){u[k]=points[b][k]-points[a][k];v[k]=points[c][k]-points[a][k];}
      std::array<Integer,4> plane{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0],0};
      for(int k=0;k<3;k++)plane[3]-=plane[k]*points[a][k];
      Integer divisor=0;for(auto x:plane)divisor=gcd(divisor,x);if(divisor==0)continue;
      for(auto& x:plane)x/=divisor;
      bool reversed=false;for(auto x:plane)if(x!=0){reversed=x<0;break;}
      if(reversed)for(auto& x:plane)x=-x;
      auto found=ids.find(plane);int id;
      if(found!=ids.end())id=found->second;
      else {
        std::vector<Integer> values;bool positive=false,negative=false;
        for(auto point:points){Integer value=plane[3];for(int k=0;k<3;k++)value+=plane[k]*point[k];positive|=value>0;negative|=value<0;values.push_back(value);}
        if(!positive||!negative){ids[plane]=-1;continue;}
        id=int(plane_values.size());ids[plane]=id;plane_values.push_back(values);
        for(auto& value:values)value=-value;
        plane_values.push_back(values);
      }
      if(id>=0)triple_plane[(size_t(a)*boundary_count+b)*boundary_count+c]=id+int(reversed);
    }
    size_t n=plane_values.size();pair_cache=std::make_unique<std::atomic<uint8_t>[]>(n*n);
    for(size_t i=0;i<n*n;i++)pair_cache[i].store(0,std::memory_order_relaxed);
  }
  bool locations_compatible(int a,int b) const {
    if(a==b)return true;
    size_t n=plane_values.size(),key=size_t(a)*n+b;
    int cached=pair_cache[key].load(std::memory_order_relaxed);if(cached)return cached==1;
    const auto& x=plane_values[a];const auto& y=plane_values[b];bool possible=false;
    for(int v=0;v<boundary_count&&!possible;v++)possible=x[v]>0&&y[v]>0;
    // If no hull vertex works, a segment can enter both strict halfspaces iff
    // its two positive/negative endpoint ratios overlap. All products are exact.
    for(int v=0;v<boundary_count&&!possible;v++)if(x[v]>0)
      for(int w=0;w<boundary_count&&!possible;w++)if(y[w]>0)
        possible=x[v]*y[w]>x[w]*y[v];
    pair_cache[key].store(possible?1:2,std::memory_order_relaxed);
    pair_cache[size_t(b)*n+a].store(possible?1:2,std::memory_order_relaxed);
    return possible;
  }
  bool impose_locations(State& s,const std::vector<Hex>& cells)const {
    if(triple_plane.empty())return true;
    for(auto h:cells)for(auto positions:corners) {
      int missing=-1,vertex=-1,count=0,k=0;std::array<int,3> triple;
      for(int j=0;j<4;j++)if(h[positions[j]]>=boundary_count){missing=j;vertex=h[positions[j]];++count;}
      if(count!=1)continue;
      for(int j=0;j<4;j++)if(j!=missing)triple[k++]=h[positions[j]];
      int parity=(3-missing)&1;
      for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)parity^=triple[i]>triple[j];
      std::sort(triple.begin(),triple.end());
      int code=triple_plane[(size_t(triple[0])*boundary_count+triple[1])*boundary_count+triple[2]];
      if(code<0)continue;
      code^=parity;bool duplicate=false;
      for(auto [v,previous]:s.locations)if(v==vertex) {
        if(previous==code){duplicate=true;break;}
        if(!locations_compatible(code,previous))return false;
      }
      if(!duplicate)s.locations.push_back({vertex,code});
    }
    return true;
  }
  void build_domains(const std::vector<Vec>& points,bool axial,int group_size,bool half_turn=false) {
    // Larger generic fixtures use the exact scalar checker to bound table memory.
    if(side>33)return;
    auto action=actions(points,axial,group_size,half_turn);
    domains.resize(size_t(4)*side*side*side);
    for(int slot=0;slot<4;slot++)for(int a=0;a<side;a++)for(int b=0;b<side;b++)for(int c=0;c<side;c++) {
      auto& mask=domains[((size_t(slot)*side+a)*side+b)*side+c];mask.word={0,0};
      for(int candidate=0;candidate<side;candidate++) {
        std::array<int,3> rest{a,b,c};std::array<int,4> q;int k=0;
        for(int j=0;j<4;j++)q[j]=j==slot?candidate:rest[k++];
        bool good=true;
        for(int g=0;g<group_size&&good;g++) {
          auto image=q;for(int& v:image)if(v<boundary_count)v=action[v][g];
          if(!half_turn&&(g==1||g==2))std::swap(image[1],image[2]);
          good=allowed[index(image)];
        }
        if(good)mask.set(candidate);
      }
      if(mask.word[boundary_count/64]&(uint64_t(1)<<(boundary_count%64)))
        for(int v=boundary_count;v<128;v++)mask.set(v);
    }
  }
  ChoiceMask choices(const Hex& h,int position)const {
    ChoiceMask result;
    for(auto positions:corners) {
      if(std::find(positions.begin(),positions.end(),position)==positions.end())continue;
      int slot=-1,known=0,k=0;std::array<int,3> q;
      for(int j=0;j<4;j++) {
        int p=positions[j];if(p==position){slot=j;continue;}
        int v=p<position?h[p]:boundary_count;
        if(v>=boundary_count)v=boundary_count;
        q[k++]=v;known+=v<boundary_count;
      }
      if(slot<0||known<2)continue;
      result.intersect(domains[((size_t(slot)*side+q[0])*side+q[1])*side+q[2]]);
    }
    return result;
  }
  bool feasible(const State& s,const Hex& h,int assigned,int changed=-1) const {
    for(auto positions:corners) {
      if(changed>=0&&std::find(positions.begin(),positions.end(),changed)==positions.end())continue;
      for(int g=0;g<s.group_size;g++) {
        std::array<int,4> q;int known=0;
        for(int j=0;j<4;j++) {
          int pos=positions[j];int v=pos<assigned?h[pos]:-1;
          q[j]=v>=0&&v<boundary_count?s.a[v][g]:boundary_count;known+=q[j]<boundary_count;
        }
        if(known<3)continue;
        // Reflection reverses the determinant; reflect() compensates by a cube
        // ordering reversal. Equivalently swap two tetrahedron slots here.
        if(s.odd(g))std::swap(q[1],q[2]);
        if(!allowed[index(q)])return false;
      }
    }
    return true;
  }
};

// Rollback bipartiteness constraints. Each edge requires color(a)^color(b)=1.
// No path compression: every union can be undone without touching older unions.
struct Parity {
  std::array<int,MAXV> parent,size,difference{};
  std::vector<std::pair<int,int>> trail;
  Parity() {std::iota(parent.begin(),parent.end(),0);size.fill(1);trail.reserve(MAXV);}
  std::pair<int,int> root(int v) const {
    int parity=0;while(parent[v]!=v){parity^=difference[v];v=parent[v];}return {v,parity};
  }
  bool relate(int a,int b,int delta) {
    auto [ra,pa]=root(a);auto [rb,pb]=root(b);
    if(ra==rb)return (pa^pb)==delta;
    if(size[ra]<size[rb])std::swap(ra,rb);
    trail.push_back({rb,ra});parent[rb]=ra;difference[rb]=pa^pb^delta;size[ra]+=size[rb];return true;
  }
  bool edge(int a,int b){return relate(a,b,1);}
  bool symmetry(const State& s) {
    for(int v=0;v<int(s.a.size());v++)for(int g=1;g<s.group_size;g++)
      if(!relate(v,s.a[v][g],s.symmetry_flip[g]))return false;
    return true;
  }
  size_t checkpoint() const {return trail.size();}
  void rollback(size_t checkpoint) {
    while(trail.size()>checkpoint){auto [child,p]=trail.back();trail.pop_back();size[p]-=size[child];parent[child]=child;difference[child]=0;}
  }
  void initialize(const State& s) {
    for(int v=0;v<int(s.a.size());v++)for(int w=0;w<v;w++)if(s.relation[v*MAXV+w]&1)
      if(!edge(v,w))throw std::runtime_error("Accepted state has an odd cycle");
    if(s.symmetry_parity&&!symmetry(s))throw std::runtime_error("Accepted state violates symmetry colors");
    trail.clear();
  }
  bool extend(const State& s,const Hex& h,int position) {
    if(s.symmetry_parity)for(int g=1;g<s.group_size;g++)if(!relate(h[position],s.a[h[position]][g],s.symmetry_flip[g]))return false;
    for(auto& e:EI) {
      int other=-1;
      if(e[0]==position&&e[1]<position)other=e[1];
      if(e[1]==position&&e[0]<position)other=e[0];
      if(other<0)continue;
      for(int g=0;g<s.group_size;g++)if(!edge(s.a[h[position]][g],s.a[h[other]][g]))return false;
    }
    return true;
  }
};
// A cell exchanged with a distinct reflected cell lies wholly on one side
// of the mirror: otherwise its connected interior meets the fixed plane and
// overlaps the reflected cell's interior. Fixed cells may cross the mirror.
bool side_parity(const State& s,int axis,Parity& parity) {
  const int g=1<<axis,anchor=MAXV-1;
  for(int v=0;v<int(s.a.size());v++)if(s.a[v][g]!=v) {
    if(!parity.relate(v,s.a[v][g],1))return false;
    if(s.boundary_side[axis][v]>=0&&!parity.relate(v,anchor,s.boundary_side[axis][v]))return false;
  }
  for(auto h:s.cells) {
    bool fixed=true;for(int v:h)fixed&=std::find(h.begin(),h.end(),s.a[v][g])!=h.end();
    if(fixed)continue;
    int first=-1;for(int v:h)if(s.a[v][g]!=v) {
      if(first<0)first=v;else if(!parity.relate(first,v,0))return false;
    }
  }
  return true;
}

constexpr int MAX_FRONT=256;
using FaceDomains=std::array<std::bitset<MAX_FRONT>,6>;
struct Budget {
  std::array<std::bitset<MAX_FRONT>,MAXV> vertex_faces{};
  FaceDomains initial;
  int open,remaining,vertices,group_size;
  bool exact_domains;
  int packing_size=0;
  std::array<std::bitset<MAX_FRONT>,4> packing;
  Budget(const State& s,int cap,bool exact=false):open(int(s.front.size())),remaining(cap-int(s.cells.size())),vertices(int(s.a.size())),group_size(s.group_size),exact_domains(exact) {
    if(open>MAX_FRONT)throw std::runtime_error("Front capacity exceeded");
    int index=0;for(auto [key,q]:s.front){(void)key;for(int v:q)vertex_faces[v].set(index);++index;}
    for(auto& d:initial)for(int i=0;i<open;i++)d.set(i);
    packing_size=int(s.cover_packing.count());
    if(packing_size) {
      std::map<Quad,int> ids;int i=0;for(auto [k,q]:s.front){(void)q;ids[k]=i++;}
      i=0;for(auto [k,q]:s.front) {
        (void)k;for(int g=0;g<group_size;g++){auto r=q;for(int& v:r)v=s.a[v][g];if(s.cover_packing[ids.at(key(r))])packing[g].set(i);}
        ++i;
      }
    }
  }
  void extend(const Hex& h,int position,FaceDomains& domains) const {
    for(int f=0;f<6;f++) {
      bool contains=false;for(int j:FI[f])contains|=j==position;
      if(contains)domains[f]&=vertex_faces[h[position]];
      else if(exact_domains)domains[f]&=~vertex_faces[h[position]];
    }
  }
  bool feasible(const Overlaps& overlaps,const FaceDomains& domains) const {
    int can_cover=0,can_fix=0,can_share=0;
    for(auto& d:domains)can_cover+=d.any();
    for(int g=1;g<group_size;g++){auto p=overlaps.possible[g-1];can_fix+=bool(p&2);can_share+=bool(p&4);}
    int packing_images=0;
    if(packing_size) {
      std::bitset<MAX_FRONT> possible;for(const auto& d:domains)possible|=d;
      for(int g=0;g<group_size;g++)packing_images+=(possible&packing[g]).any();
    }
    for(int n=1;n<=group_size;n*=2) {
      if(packing_size-std::min(n,packing_images)>remaining-n)continue;
      if(n>remaining||can_fix<group_size/n-1)continue;
      // m <= n*can_cover old front faces can disappear. Each nontrivial
      // group element can pair at most n/2 cells along internal faces;
      // duplicate images only make this upper bound more optimistic.
      int internal=std::min(n*(n-1)/2,(n/2)*can_share);
      if(open+12*n-2*n*can_cover-2*internal<=6*remaining)return true;
    }
    return false;
  }
};

// Cycle ranks of all edge subgraphs of a cube. Isolated vertices do not
// affect E - V + components, so no vertex bookkeeping is needed at runtime.
struct EdgeCycles {
  std::array<uint8_t,4096> rank{};
  EdgeCycles() {
    for(int mask=0;mask<4096;mask++) {
      std::array<int,8> parent;std::iota(parent.begin(),parent.end(),0);
      auto root=[&](int v){while(parent[v]!=v)v=parent[v];return v;};
      for(int e=0;e<12;e++)if(mask&(1<<e)) {
        int a=root(EI[e][0]),b=root(EI[e][1]);
        if(a==b)++rank[mask];else parent[a]=b;
      }
    }
  }
  int upper(int mask)const{return rank[mask]+(mask==4095);}
};
const EdgeCycles edge_cycles;
// In the component bound the shared-face terms cancel the faces consumed
// from the front. For a valid insertion this leaves
//   F + 10*n <= 4*R + 2*C + 2*sum(cycle_rank + beta_2).
// Assume all unassigned edges are shared, giving a monotone upper bound.
// For a distinct reflected pair, their common edges are pointwise fixed.

using CoverKey=std::array<uint64_t,6>;
struct CoverHash {size_t operator()(const CoverKey& key)const {
  uint64_t h=0x9e3779b97f4a7c15ULL;
  for(auto x:key){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;h^=x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);}
  return h;
}};
struct alignas(64) CoverEntry {
  CoverKey key;
  uint32_t value=UINT32_MAX;
};
struct StaticPairEntry {uint64_t boundary=0;uint32_t shape=UINT32_MAX,value=0;};
struct alignas(64) Counts {
  uint64_t audit_serial=0;std::ofstream audit_states,audit_expansions;
  std::unique_ptr<StaticPairEntry[]> static_pair_cache;
  std::unordered_map<CoverKey,uint32_t,CoverHash> cover_cache;
  std::unique_ptr<CoverEntry[]> flat_cache;
  std::atomic<uint64_t> canonical_rejections{0};
  std::atomic<uint64_t> cover_fixed_checks{0},cover_fixed_rejections{0};
  std::atomic<uint64_t> twopack_calls{0},twopack_rejections{0},twopack_triples{0},twopack_valid{0};
  std::atomic<uint64_t> partial_cover_calls{0},partial_cover_rejections{0},quick_packing_rejections{0};
  std::atomic<uint64_t> cache_hits{0},cache_misses{0},pair_reuses{0};
  std::atomic<uint64_t> branches{0},trials{0},states{0},candidates{0},rejected{0};
  std::array<std::atomic<uint64_t>,8> insertion_rejections{};
  std::atomic<uint64_t> early_overlap{0},early_budget{0},early_parity{0},component_rejections{0},early_completion{0},early_geometry{0},early_domain{0},reject_location{0},reject_side{0},partial_completion{0},cover_rejections{0},cover_work{0},cover_unknown{0};
};

struct VertexDomains;
#include "boundary_automorphisms.hpp"

struct Search {
  // Optional single-thread proposal adapter. Empty hooks preserve exhaustive DFS.
  std::function<void(State&&)> proposal_sink;
  std::function<bool()> proposal_stop;
  std::function<void(std::vector<int>&)> proposal_order;
  int root_face=-1;
  int cover_fixed=1;
  bool supports_fixed_mesh(const State& s,const Hex& h)const;
  bool supports_fixed_mesh_uncached(const State& s,const Hex& h)const;
  std::vector<Hex> extend_cube_group_raw(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const;
  std::string tree_trace;
  Mesh boundary;
  int cover_closed=1,cover_memo_enabled=1,two_face_packing=1;
  int fast_pairs=1;
  int leader_bound=1;
  int boundary_canonical=1,topology_symmetry=0;
  struct BoundarySymmetry {std::vector<int> permutation,source;bool odd=false;};
  std::vector<BoundarySymmetry> boundary_symmetries;
  std::map<Quad,int> boundary_ids;std::vector<int> boundary_order;
  int partial_cover=6,partial_work=1000,quick_packing=1;
  int cap=35,threads=64,vertex_cap=0,progress=5,face_order=3,task_depth=0,mirrors=2,face_cover=1,cover_work_limit=1000,cover_cubes=1,cover_cache=1,mirror_cover=1,cover_global=1,cover_symmetry=1,cover_incremental=1,cover_seed=1,cover_flat=1,cover_cache_bits=18;
  double seconds=0;
  bool half_turn=false;
  int group_size()const{return half_turn?2:1<<mirrors;}
  bool axial=false,early_overlap=true,early_budget=true,early_parity=true,component_bound=true,early_completion=true,corner_geometry=true,vertex_domains=true,face_domains=true,symmetry_parity=true,color_domains=true,geometry_domains=true,location_domains=true,mirror_separation=true,edge_completion=true,incremental_faces=true,overlap_domains=true,component_index=true,completion_domains=true;
  CornerGeometry geometry;
  std::string output="candidates.jsonl",state_trace,replay_actions;
  std::unique_ptr<Counts[]> counts;
  std::atomic<int> tasks{0};
  std::atomic<bool> failed{false},finished{false};
  std::mutex io,monitor_mutex;
  std::condition_variable monitor_wake;
  std::ofstream out,trace;
  std::chrono::steady_clock::time_point start;
  std::thread monitor;
  Counts& local()const{return counts[omp_get_thread_num()];}
  bool stopping() const { return interrupted.load(std::memory_order_relaxed)||failed.load(std::memory_order_relaxed); }
  void error(const std::exception& e) { failed=true;std::lock_guard<std::mutex> l(io);std::cerr<<"ERROR "<<e.what()<<'\n'; }
  double elapsed()const{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
  std::array<uint64_t,5> totals()const {
    std::array<uint64_t,5> r{};
    for(int i=0;i<threads;i++){r[0]+=counts[i].branches.load();r[1]+=counts[i].trials.load();r[2]+=counts[i].states.load();r[3]+=counts[i].candidates.load();r[4]+=counts[i].rejected.load();}
    return r;
  }
  void report() {
    auto n=totals();std::lock_guard<std::mutex> l(io);
    std::cerr<<"PROGRESS elapsed="<<elapsed()<<" branches="<<n[0]<<" orbit_trials="<<n[1]<<" states="<<n[2]
             <<" candidates="<<n[3]<<" topology_rejections="<<n[4]<<" tasks="<<tasks.load();
    const char* names[]={"orbit_shape","cell_cap","boundary","placed_cells","orbit_cells","face_incidence","face_budget","odd_cycle"};
    for(int k=0;k<8;k++){uint64_t n=0;for(int i=0;i<threads;i++)n+=counts[i].insertion_rejections[k].load();std::cerr<<" reject_"<<names[k]<<'='<<n;}
    uint64_t early=0;for(int i=0;i<threads;i++)early+=counts[i].early_overlap.load();
    std::cerr<<" early_overlap="<<early;
    uint64_t budget=0;for(int i=0;i<threads;i++)budget+=counts[i].early_budget.load();
    std::cerr<<" early_budget="<<budget;
    uint64_t parity=0;for(int i=0;i<threads;i++)parity+=counts[i].early_parity.load();
    std::cerr<<" early_parity="<<parity;
    uint64_t component=0;for(int i=0;i<threads;i++)component+=counts[i].component_rejections.load();
    std::cerr<<" reject_components="<<component;
    uint64_t completion=0;for(int i=0;i<threads;i++)completion+=counts[i].early_completion.load();
    std::cerr<<" early_completion="<<completion;
    uint64_t geo=0;for(int i=0;i<threads;i++)geo+=counts[i].early_geometry.load();
    std::cerr<<" early_geometry="<<geo;
    uint64_t domains=0;for(int i=0;i<threads;i++)domains+=counts[i].early_domain.load();
    std::cerr<<" early_domain="<<domains;
    uint64_t locations=0;for(int i=0;i<threads;i++)locations+=counts[i].reject_location.load();
    std::cerr<<" reject_location="<<locations;
    uint64_t sides=0;for(int i=0;i<threads;i++)sides+=counts[i].reject_side.load();
    std::cerr<<" reject_side="<<sides;
    uint64_t partial=0;for(int i=0;i<threads;i++)partial+=counts[i].partial_completion.load();
    std::cerr<<" partial_completion="<<partial;
    uint64_t cr=0;for(int i=0;i<threads;i++)cr+=counts[i].cover_rejections.load();
    std::cerr<<" cover_rejections="<<cr;
    uint64_t ch=0,cm=0;for(int i=0;i<threads;i++){ch+=counts[i].cache_hits.load();cm+=counts[i].cache_misses.load();}
    std::cerr<<" cover_cache_hits="<<ch<<" cover_cache_misses="<<cm;
    uint64_t work=0,unknown=0;for(int i=0;i<threads;i++){work+=counts[i].cover_work.load();unknown+=counts[i].cover_unknown.load();}
    uint64_t reuses=0;for(int i=0;i<threads;i++)reuses+=counts[i].pair_reuses.load();
    std::cerr<<" cover_pair_reuses="<<reuses;
    uint64_t pc=0,pr=0,qr=0;for(int i=0;i<threads;i++){pc+=counts[i].partial_cover_calls.load();pr+=counts[i].partial_cover_rejections.load();qr+=counts[i].quick_packing_rejections.load();}
    uint64_t canonical=0;for(int i=0;i<threads;i++)canonical+=counts[i].canonical_rejections.load();
    std::cerr<<" canonical_rejections="<<canonical;
    uint64_t fc=0,fr=0;for(int i=0;i<threads;i++){fc+=counts[i].cover_fixed_checks;fr+=counts[i].cover_fixed_rejections;}
    std::cerr<<" cover_fixed_checks="<<fc<<" cover_fixed_rejections="<<fr;
    uint64_t tc=0,tr=0,tt=0,tv=0;for(int i=0;i<threads;i++){tc+=counts[i].twopack_calls;tr+=counts[i].twopack_rejections;tt+=counts[i].twopack_triples;tv+=counts[i].twopack_valid;}
    std::cerr<<" twopack_calls="<<tc<<" twopack_rejections="<<tr<<" twopack_triples="<<tt<<" twopack_valid="<<tv;
    std::cerr<<" partial_cover_calls="<<pc<<" partial_cover_rejections="<<pr<<" quick_packing_rejections="<<qr;
    std::cerr<<" cover_work="<<work<<" cover_unknown="<<unknown<<'\n';
  }
  void begin() {
    counts=std::make_unique<Counts[]>(threads);start=std::chrono::steady_clock::now();
    if(!tree_trace.empty()) {
      uint16_t endian=1;if(*reinterpret_cast<const unsigned char*>(&endian)!=1)throw std::runtime_error("Tree trace requires a little-endian host");
      std::ofstream format(tree_trace+"/format.json");format<<R"TRACE({"version":2,"byte_order":"little","state_words":["id","parent_id","metadata","added_hex_vertices_plus_one","local_cover_work","local_choose_calls"],"metadata":{"placed_cells":[0,8],"vertex_count":[8,8],"front_size":[16,16],"outcome":[32,8],"cover_reason":[40,8]},"outcomes":{"1":"cover_infeasible","2":"cell_or_face_budget","3":"exhausted_expansion","4":"candidate","5":"topology_rejected","6":"interrupted"},"cover_reasons":{"0":"none","1":"inherited_quick_packing","2":"greedy_incompatible_packing","3":"bounded_cover_proof","4":"mirror_cover","5":"two_face_packing"},"expansion_words":["id","branches","cache_hits","cache_misses","candidates","canonical_rejections","component_rejections","cover_rejections","early_budget","early_completion","early_domain","early_geometry","early_overlap","early_parity","pair_reuses","partial_completion","partial_cover_calls","partial_cover_rejections","quick_packing_rejections","reject_location","reject_side","rejected","states","trials","cover_work","cover_unknown","insert_orbit_shape","insert_cell_cap","insert_boundary","insert_placed_cells","insert_orbit_cells","insert_face_incidence","insert_face_budget","insert_odd_cycle","rule_conflict","forced_choice","vertex_cap","relations","face_consistency","initial_overlap","initial_budget","children"],"id":"(worker_local_serial << 8) | worker; serial starts at 1; parent 0 denotes root"})TRACE"<<'\n';
      if(!format)throw std::runtime_error("Cannot write tree trace schema");
    }
    if(!tree_trace.empty())for(int i=0;i<threads;i++) {
      counts[i].audit_states.open(tree_trace+"/states-"+std::to_string(i)+".bin",std::ios::binary);
      counts[i].audit_expansions.open(tree_trace+"/expansions-"+std::to_string(i)+".bin",std::ios::binary);
      if(!counts[i].audit_states||!counts[i].audit_expansions)throw std::runtime_error("Cannot open tree trace directory");
    }
    out.open(output);if(!out)throw std::runtime_error("Cannot open candidate output");
    if(!state_trace.empty()){trace.open(state_trace);if(!trace)throw std::runtime_error("Cannot open state trace");}
    monitor=std::thread([this]{
      std::unique_lock<std::mutex> l(monitor_mutex);
      while(!finished.load()) {
        monitor_wake.wait_for(l,std::chrono::seconds(1),[this]{return finished.load();});
        if(finished.load())break;
        if(seconds>0&&elapsed()>=seconds)interrupted=true;
        if(progress>0&&int(elapsed())%progress==0)report();
      }
    });
  }
  void end(){
    // All search tasks have joined. Release each worker's bounded cache in
    // parallel, avoiding a serial walk over millions of nodes at shutdown.
    #pragma omp parallel for num_threads(threads) schedule(static)
    for(int i=0;i<threads;i++) {
      std::unordered_map<CoverKey,uint32_t,CoverHash> empty;
      empty.swap(counts[i].cover_cache);counts[i].flat_cache.reset();counts[i].static_pair_cache.reset();
    }
    if(!tree_trace.empty())for(int i=0;i<threads;i++) {
      counts[i].audit_states.flush();counts[i].audit_expansions.flush();
      if(!counts[i].audit_states||!counts[i].audit_expansions){failed=true;std::cerr<<"ERROR tree trace flush failed\n";}
    }
    finished=true;monitor_wake.notify_all();if(monitor.joinable())monitor.join();report();
  }
  static Hex orient_partial(const Hex& c,Quad q);
  uint64_t root_code(const Hex& h,bool upper)const;
  bool leader_feasible(const State& s,const Hex& h)const;
  void update_root_floor(State& s)const;
  void prepare_boundary_symmetry();
  bool boundary_canonical_feasible(const State& s)const;
  Quad select_face(const State& s) const;
  State initial();
  bool completion_feasible(const State& s,const Hex& h,const std::array<uint8_t,MAXV*MAXV>& placed,const Budget& budget) const;
  bool insert(State& s,const Hex& h);
  bool topology(const State& s)const;
  void emit(const State& s);
  void dfs(State s,int depth);
  void expand_face(State& s,Quad q,int depth);
  uint32_t face_templates(const State& s,Quad q,Quad r)const;
  uint32_t face_templates_reference(const State& s,Quad q,Quad r)const;
  bool faces_can_share(const State& s,Quad q,Quad r,const Diagonals& diagonals)const;
  std::vector<Hex> extend_cube_group_uncached(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const;
  std::vector<Hex> extend_cube_group(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const;
  bool cover_feasible(State& s,int override_remaining=-1,const Diagonals* override_diagonals=nullptr)const;
  bool partial_cover_feasible(const State& s,const Hex& h,int assigned,const Budget& budget,const FaceDomains& domains,const Diagonals& diagonals,const std::array<uint8_t,MAXV*MAXV>& placed)const;
  bool orbit_cover_feasible(const State& s,const std::vector<Quad>& faces,const Diagonals& diagonals)const;
  void choose(State& s,Hex& h,int position,int depth,const std::array<uint8_t,MAXV*MAXV>& placed,const Diagonals& diagonals,const Rules& rules,Overlaps overlaps,const Budget& budget,FaceDomains domains,Parity& parity,const VertexDomains& vertices,std::array<ChoiceMask,2> colors);
  bool replay(const Mesh& seed);
};

bool extend_overlaps(const State& s,const Hex& h,int position,Overlaps& overlaps) {
  for(int g=1;g<s.group_size;g++) {
    int image=s.a[h[position]][g];
    for(int j=0;j<=position;j++)if(h[j]==image) {
      if(s.symmetry_parity&&(__builtin_popcount(unsigned(BITS[position]^BITS[j]))&1)!=s.symmetry_flip[g])return false;
      overlaps.code[g-1]|=uint32_t(j+1)<<(4*position);
      overlaps.code[g-1]|=uint32_t(position+1)<<(4*j);
      break;
    }
    const auto& allowed=s.mirror_separation&&g!=3?
      (overlaps.fixed[g-1]?overlap_patterns.fixed_prefix[position+1]:overlap_patterns.reflection_prefix[position+1]):
      overlap_patterns.prefix[s.odd(g)][position+1];
    auto it=allowed.find(overlaps.code[g-1]);if(it==allowed.end())return false;
    overlaps.possible[g-1]=it->second;
  }
  return true;
}

int relation_type(int i,int j){int n=__builtin_popcount(unsigned(BITS[i]^BITS[j]));return 1<<(n-1);}
bool conflict(int old,int type,bool placed) {
  if(type==1)return old&(2|4);
  if(type==2)return old&(1|4);
  return old&(placed?7:3);
}

struct PairFits {
  struct Fit {int bit;Hex source;};
  struct Pattern {std::vector<Fit> fits;uint32_t shape=0;std::array<uint32_t,2> parity{};std::array<std::array<uint32_t,8>,16> relation{};};
  std::array<Pattern,625> patterns;
  PairFits() {
    for(int code=0;code<625;code++) {
      auto& pattern=patterns[code];int digits=code;std::array<int,8> vertices{0,1,2,3,4,5,6,7};
      for(int j=0;j<4;j++){int digit=digits%5;digits/=5;if(digit<4)vertices[4+j]=digit;}
      for(int f=1;f<6;f++)for(int rotation=0;rotation<4;rotation++) {
        Hex h{0,3,2,1,-1,-1,-1,-1};bool valid=true;
        for(int j=0;j<4;j++){int pos=FI[f][j],v=vertices[4+(j+rotation)%4];if(h[pos]>=0&&h[pos]!=v){valid=false;break;}h[pos]=v;}
        for(int i=0;i<8;i++)if(h[i]>=0)for(int j=0;j<i;j++)if(h[i]==h[j])valid=false;
        if(!valid)continue;
        int bit=(f-1)*4+rotation;uint32_t flag=uint32_t(1)<<bit;pattern.fits.push_back({bit,h});pattern.shape|=flag;
        std::array<int,8> positions;positions.fill(-1);
        for(int i=0;i<8;i++)if(h[i]>=0)positions[h[i]]=i;
        int p=__builtin_popcount(unsigned(BITS[positions[0]]^BITS[positions[vertices[4]]]))&1;pattern.parity[p]|=flag;
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int old=0;old<8;old++) {
          int a=positions[i],b=positions[vertices[4+j]];
          if(a==b||!conflict(old,relation_type(a,b),true))pattern.relation[4*i+j][old]|=flag;
        }
      }
    }
  }
};
const PairFits pair_fits;
// Static necessary pair relations against the already placed complex. New
// vertices remain unrestricted here; reflected-orbit checks still run afterward.
struct VertexDomains {
  std::array<std::array<ChoiceMask,MAXV>,3> pair;
  std::array<ChoiceMask,2> colors;
  std::array<ChoiceMask,MAXV> neighbors;
  std::array<int,MAXV> color{};
  int initial;
  std::array<ChoiceMask,3> fixed_vertices;
  std::array<std::array<ChoiceMask,MAXV>,2> opposite;
  VertexDomains(const State& s):initial(int(s.a.size())) {
    for(int g=1;g<s.group_size;g++) {
      fixed_vertices[g-1].word={0,0};
      for(int v=0;v<initial;v++)if(s.a[v][g]==v)fixed_vertices[g-1].set(v);
    }
    for(auto& axis:opposite)for(auto& mask:axis)mask.word={0,0};
    if(s.mirror_separation)for(int axis=0;axis<reflection_count(s.group_size,s.half_turn);axis++) {
      Parity sides;if(!side_parity(s,axis,sides))throw std::runtime_error("Accepted state violates mirror separation");
      int g=1<<axis;std::array<std::pair<int,int>,MAXV> component;
      for(int v=0;v<initial;v++)component[v]=sides.root(v);
      for(int v=0;v<initial;v++)if(s.a[v][g]!=v)for(int w=0;w<initial;w++)if(s.a[w][g]!=w) {
        auto [a,p]=component[v];auto [b,q]=component[w];
        if(a==b&&p!=q)opposite[axis][v].set(w);
      }
    }
    for(auto& mask:colors)mask.word={0,0};
    for(auto& mask:neighbors)mask.word={0,0};
    Parity parity;parity.initialize(s);auto [anchor,p]=parity.root(0);
    for(int v=0;v<int(s.a.size());v++) {
      auto [root,q]=parity.root(v);
      if(root!=anchor)throw std::runtime_error("Placed vertex disconnected from boundary");
      color[v]=p^q;colors[color[v]].set(v);
    }
    for(int type=0;type<3;type++)for(int a=0;a<MAXV;a++)pair[type][a].reset(a);
    for(int a=0;a<initial;a++)for(int b=0;b<initial;b++) {
      int old=s.relation[a*MAXV+b];
      if(old&1)neighbors[a].set(b);
      if(old&6)pair[0][a].reset(b);
      if(old&5)pair[1][a].reset(b);
      if(old&7)pair[2][a].reset(b);
    }
  }
  ChoiceMask choices(const Hex& h,int position,int existing,int group_size,const Diagonals& diagonals)const {
    ChoiceMask result;bool allow_new=true;
    for(int j=0;j<position;j++)result.intersect(pair[__builtin_popcount(unsigned(BITS[position]^BITS[j]))-1][h[j]]);
    for(auto& f:FI)for(int k=0;k<2;k++) {
      if(position!=f[(k+1)%4]&&position!=f[(k+3)%4])continue;
      int i=f[k],j=f[k+2];if(i>=position||j>=position)continue;
      auto it=diagonals.find(edgekey(h[i],h[j]));if(it==diagonals.end())continue;
      ChoiceMask required;required.word={0,0};for(int v:it->second)required.set(v);
      result.intersect(required);allow_new=false;
    }
    result.limit(existing);
    if(allow_new)for(int type=0;type<group_size;type++)result.set(existing+type);
    return result;
  }
};
bool edge_completion_feasible(const State& s,const Hex& h,int assigned,
    const std::array<uint8_t,MAXV*MAXV>& placed,const Budget& budget,const Overlaps& overlaps,const VertexDomains& vertices,bool use_domains) {
  if(s.group_size>2)return true;
  int old=0,fixed=0;
  for(int e=0;e<12;e++) {
    int a=EI[e][0],b=EI[e][1];
    if((a>=assigned||h[a]<budget.vertices)&&(b>=assigned||h[b]<budget.vertices)&&
       (a>=assigned||b>=assigned||(placed[h[a]*MAXV+h[b]]&1)))old|=1<<e;
    if(s.group_size==2&&(a>=assigned||s.a[h[a]][1]==h[a])&&(b>=assigned||s.a[h[b]][1]==h[b]))fixed|=1<<e;
  }
  if(!s.mirror_separation)fixed=4095;
  int rhs=4*budget.remaining+2*s.component_upper;
  auto feasible=[&](int mask) {
    if(s.group_size==1)return budget.remaining>=1&&budget.open+10<=rhs+2*edge_cycles.upper(mask);
    if((overlaps.possible[0]&2)&&budget.remaining>=1&&budget.open+10<=rhs+2*edge_cycles.upper(mask))return true;
    return budget.remaining>=2&&budget.open+20<=rhs+2*(edge_cycles.upper(mask)+edge_cycles.upper(mask|fixed));
  };
  if(!feasible(old))return false;
  if(!use_domains||assigned==8)return true;
  std::array<ChoiceMask,8> domains;
  for(int i=0;i<8;i++) {
    if(i<assigned) {domains[i].word={0,0};if(h[i]<budget.vertices)domains[i].set(h[i]);}
    else {
      domains[i]=vertices.colors[vertices.color[h[0]]^(__builtin_popcount(unsigned(BITS[i]))&1)];
      for(int j=0;j<assigned;j++)domains[i].intersect(vertices.pair[__builtin_popcount(unsigned(BITS[i]^BITS[j]))-1][h[j]]);
    }
  }
  for(int e=0;e<12;e++)if(old&(1<<e)) {
    int a=EI[e][0],b=EI[e][1];if(a<assigned&&b<assigned)continue;
    bool possible=false;auto choices=domains[a];
    for(int v=choices.pop();v>=0&&!possible;v=choices.pop())
      possible=(vertices.neighbors[v].word[0]&domains[b].word[0])||(vertices.neighbors[v].word[1]&domains[b].word[1]);
    if(!possible)old&=~(1<<e);
  }
  return feasible(old);
}

// A new vertex's image is either itself, an earlier cube vertex, or
// outside the assigned prefix. Query these categories once, then intersect
// their bitsets instead of trying every vertex and rejecting afterward.
ChoiceMask overlap_choices(const State& s,const Hex& h,int position,
                           const Overlaps& overlaps,const VertexDomains& vertices) {
  const int existing=int(s.a.size());ChoiceMask result;result.limit(existing+s.group_size);
  for(int g=1;g<s.group_size;g++) {
    auto fixed=vertices.fixed_vertices[g-1];
    for(int v=vertices.initial;v<existing;v++)if(s.a[v][g]==v)fixed.set(v);
    for(int type=0;type<s.group_size;type++)if((type&g)==g)fixed.set(existing+type);
    const auto& allowed=s.mirror_separation&&g!=3?
      (overlaps.fixed[g-1]?overlap_patterns.fixed_prefix[position+1]:overlap_patterns.reflection_prefix[position+1]):
      overlap_patterns.prefix[s.odd(g)][position+1];
    auto code=overlaps.code[g-1];
    auto flags=[&](uint32_t next){auto it=allowed.find(next);return it==allowed.end()?0:int(it->second);};
    int outside=flags(code),self=flags(code|(uint32_t(position+1)<<(4*position)));
    ChoiceMask valid;valid.word={0,0};
    for(int w=0;w<2;w++)valid.word[w]=(outside?~fixed.word[w]:0)|(self?fixed.word[w]:0);
    // Opposite known sides force a coincident cell, whose flag is bit 1.
    ChoiceMask opposite;opposite.word={0,0};
    if(s.mirror_separation&&g!=3) {
      int axis=g==1?0:1;
      for(int j=0;j<position;j++)if(h[j]<vertices.initial)
        for(int w=0;w<2;w++)opposite.word[w]|=vertices.opposite[axis][h[j]].word[w];
      if(!(outside&2))for(int w=0;w<2;w++)valid.word[w]&=~(opposite.word[w]&~fixed.word[w]);
      if(!(self&2))for(int w=0;w<2;w++)valid.word[w]&=~(opposite.word[w]&fixed.word[w]);
    }
    for(int j=0;j<position;j++) {
      int v=s.a[h[j]][g];valid.reset(v);
      auto next=code|(uint32_t(j+1)<<(4*position))|(uint32_t(position+1)<<(4*j));
      int f=flags(next);
      if(f&&(!(opposite.word[v/64]&(uint64_t(1)<<(v%64)))||(f&2)))valid.set(v);
    }
    // All cube vertices must be distinct, including previously fixed ones.
    for(int j=0;j<position;j++)valid.reset(h[j]);
    result.intersect(valid);
  }
  return result;
}
struct Trail {
  struct Entry {int first;uint8_t second;};
  std::array<Entry,56> entries;
  size_t used=0;
  void push_back(Entry value){if(used==entries.size())throw std::runtime_error("Relation trail capacity exceeded");entries[used++]=value;}
  auto rbegin(){return std::make_reverse_iterator(entries.begin()+used);}
  auto rend(){return std::make_reverse_iterator(entries.begin());}
};
bool impose(State& s,const Hex& h,int position,const std::array<uint8_t,MAXV*MAXV>& placed,Trail& trail) {
  int v=h[position];
  for(int j=0;j<position;j++)if(h[j]==v)return false;
  for(int j=0;j<position;j++)for(int g=0;g<s.group_size;g++) {
    int a=s.a[v][g],b=s.a[h[j]][g],type=relation_type(position,j);
    if(a==b)return false;
    int index=a*MAXV+b,other=b*MAXV+a;
    if(conflict(placed[index],type,true)||conflict(s.relation[index],type,false))return false;
    if(!(s.relation[index]&type)) {
      trail.push_back({index,s.relation[index]});trail.push_back({other,s.relation[other]});
      s.relation[index]|=type;s.relation[other]|=type;
    }
  }
  return true;
}
void undo(State& s,Trail& trail){for(auto it=trail.rbegin();it!=trail.rend();++it)s.relation[it->first]=it->second;}
Diagonals known_diagonals(const State& s,const Mesh& boundary) {
  Diagonals d;
  auto add=[&](Quad q){d[edgekey(q[0],q[2])]=q;d[edgekey(q[1],q[3])]=q;};
  for(auto q:boundary.boundary)add(q);
  for(auto c:s.cells)for(int f=0;f<6;f++)add(face(c,f));
  return d;
}
bool consistent_faces(const State& s,const Hex& h,int assigned,const Diagonals& d) {
  // Reusing an existing face diagonal obliges the candidate to share that whole
  // face. Apply the obligation to all reflected images while it is still partial.
  for(auto& f:FI)for(int k=0;k<2;k++) {
    int i=f[k],j=f[k+2];if(i>=assigned||j>=assigned)continue;
    for(int g=0;g<s.group_size;g++) {
      int a=s.a[h[i]][g],b=s.a[h[j]][g];auto it=d.find(edgekey(a,b));if(it==d.end())continue;
      const Quad& q=it->second;
      for(int l:{f[(k+1)%4],f[(k+3)%4]})if(l<assigned) {
        int v=s.a[h[l]][g];if(std::find(q.begin(),q.end(),v)==q.end()||v==a||v==b)return false;
      }
    }
  }
  return true;
}
// All earlier obligations have already passed. Only test a diagonal
// completed by this vertex, or this vertex completing a diagonal's face.
// The placed complex and diagonal map are equivariant, so testing the
// representative also tests every reflected image.
bool consistent_new_faces(const Hex& h,int position,const Diagonals& d,
                          const std::array<uint8_t,MAXV*MAXV>& placed) {
  for(auto& f:FI)for(int k=0;k<2;k++) {
    int i=f[k],j=f[k+2],a=f[(k+1)%4],b=f[(k+3)%4];
    if(i>position||j>position)continue;
    if(i!=position&&j!=position&&a!=position&&b!=position)continue;
    if(!(placed[h[i]*MAXV+h[j]]&2))continue;
    auto it=d.find(edgekey(h[i],h[j]));
    if(it==d.end())throw std::runtime_error("Missing placed face diagonal");
    for(int v:{a,b})if(v<=position&&(i==position||j==position||v==position))
      if(std::find(it->second.begin(),it->second.end(),h[v])==it->second.end()||h[v]==h[i]||h[v]==h[j])return false;
  }
  return true;
}
Rules stabilizer_rules(const State& s,const Hex& h) {
  Rules rules;
  for(int g=1;g<s.group_size;g++) {
    std::array<int,4> targets{};bool fixes=true;
    for(int i=0;i<4;i++) {
      auto it=std::find(h.begin(),h.begin()+4,s.a[h[i]][g]);
      if(it==h.begin()+4){fixes=false;break;}targets[i]=int(it-h.begin());
    }
    if(fixes)for(int i=0;i<4;i++)rules.push_back({i+4,targets[i]+4,g});
  }
  return rules;
}
bool consistent_rules(const State& s,const Hex& h,int assigned,const Rules& rules) {
  for(auto r:rules)if(r[0]<assigned&&r[1]<assigned&&s.a[h[r[0]]][r[2]]!=h[r[1]])return false;
  return true;
}
int allocate(State& s,int stabilizer) {
  // stabilizer: 0=off both planes; 1=on A; 2=on B; 3=on their axis.
  int base=int(s.a.size());std::vector<int> representatives;
  auto coset=[stabilizer](int g){
    if(stabilizer==3)return 0;
    if(stabilizer==1)return g&2;
    if(stabilizer==2)return g&1;
    return g;
  };
  for(int g=0;g<s.group_size;g++)if(coset(g)==g)representatives.push_back(g);
  for(int r:representatives){Action a{-1,-1,-1,-1};for(int g=0;g<s.group_size;g++)a[g]=base+int(std::find(representatives.begin(),representatives.end(),coset(r^g))-representatives.begin());s.a.push_back(a);}
  return base;
}
Hex Search::orient_partial(const Hex& c,Quad q) {
  std::array<int,4> positions{};Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};
  for(int j=0;j<4;j++) {
    auto it=std::find(c.begin(),c.end(),h[j]);
    if(it==c.end())throw std::runtime_error("Missing boundary face corner");
    positions[j]=int(it-c.begin());
  }
  for(int j=0;j<4;j++)for(auto e:EI) {
    int a=e[0],b=e[1];if(b==positions[j])std::swap(a,b);
    if(a==positions[j]&&std::find(positions.begin(),positions.end(),b)==positions.end())h[j+4]=c[b];
  }
  return h;
}
uint64_t Search::root_code(const Hex& h,bool upper)const {
  int next=int(boundary.points.size());uint64_t code=uint64_t(1)<<63;
  for(int j=4;j<8;j++) {
    int v=h[j];
    if(v<0)v=upper?next++:0;
    else if(v>=int(boundary.points.size()))v=next++;
    code|=uint64_t(v)<<(8*(7-j));
  }
  return code;
}
void Search::update_root_floor(State& s)const {
  s.canonical_floor=0;
  if(!leader_bound||boundary_order.empty()||s.boundary_owner.empty())return;
  int root=boundary_order.front(),owner=s.boundary_owner[root];
  if(owner>=0)s.canonical_floor=root_code(orient_partial(s.cells[owner],boundary.boundary[root]),false);
}
bool Search::leader_feasible(const State& s,const Hex& h)const {
  if(!s.canonical_floor)return true;
  int root=boundary_order.front(),n=int(boundary.points.size());
  for(int f=0;f<6;f++) {
    bool complete=true;for(int j:FI[f])if(h[j]<0||h[j]>=n){complete=false;break;}
    if(!complete)continue;
    auto found=boundary_ids.find(key(face(h,f)));if(found==boundary_ids.end())continue;
    for(const auto& g:boundary_symmetries)if(g.source[root]==found->second) {
      Hex transformed=h;for(int& v:transformed)if(v>=0&&v<n)v=g.permutation[v];
      auto oriented=orient_partial(transformed,boundary.boundary[root]);
      if(root_code(oriented,true)<s.canonical_floor)return false;
    }
  }
  return true;
}

// Boundary automorphisms quotient equivalent meshes without imposing symmetry
// on the mesh. General topological maps are restricted to geometry-free inputs.
void Search::prepare_boundary_symmetry() {
  boundary_symmetries.clear();boundary_ids.clear();boundary_order.clear();
  if(!boundary_canonical||mirrors>1)return;
  int n=int(boundary.points.size());
  for(int i=0;i<int(boundary.boundary.size());i++)boundary_ids[key(boundary.boundary[i])]=i;
  if(topology_symmetry) {
    for(auto a:surface_automorphisms(boundary.boundary,n)) {
      BoundarySymmetry g;g.permutation=std::move(a.permutation);g.odd=a.odd;
      g.source.resize(boundary.boundary.size(),-1);
      for(int i=0;i<int(boundary.boundary.size());i++) {
        auto q=boundary.boundary[i];for(int& v:q)v=g.permutation[v];if(g.odd)q=reverse(q);
        auto it=boundary_ids.find(key(q));
        if(it==boundary_ids.end()||!same_orientation(q,boundary.boundary[it->second]))throw std::runtime_error("Invalid surface automorphism");
        g.source[it->second]=i;
      }
      boundary_symmetries.push_back(std::move(g));
    }
  } else
  for(bool swap:{false,true})for(int sx:{1,-1})for(int sy:{1,-1}) {
    BoundarySymmetry g;g.permutation.resize(n,-1);g.odd=(swap?-1:1)*sx*sy<0;
    std::vector<bool> used(n);bool valid=true;
    for(int v=0;v<n;v++) {
      auto q=boundary.points[v];Vec image{sx*q[swap?1:0],sy*q[swap?0:1],q[2]};
      for(int w=0;w<n;w++)if(image==boundary.points[w]){g.permutation[v]=w;break;}
      if(g.permutation[v]<0||used[g.permutation[v]]){valid=false;break;}used[g.permutation[v]]=true;
    }
    if(!valid)continue;
    // Quotient only isometries preserving the requested reflection subgroup.
    if(mirrors==1||half_turn) {
      auto action=actions(boundary.points,axial,2,half_turn);
      for(int v=0;v<n;v++)if(g.permutation[action[v][1]]!=action[g.permutation[v]][1]){valid=false;break;}
      if(!valid)continue;
    }
    g.source.resize(boundary.boundary.size(),-1);
    for(int i=0;i<int(boundary.boundary.size());i++) {
      auto q=boundary.boundary[i];for(int& v:q)v=g.permutation[v];if(g.odd)q=reverse(q);
      auto it=boundary_ids.find(key(q));
      if(it==boundary_ids.end()||!same_orientation(q,boundary.boundary[it->second])){valid=false;break;}
      g.source[it->second]=i;
    }
    if(valid)boundary_symmetries.push_back(std::move(g));
  }
  if(boundary_symmetries.size()<2)return;
  auto state=initial();auto root=key(select_face(state));
  boundary_order.push_back(boundary_ids.at(root));
  for(auto [q,i]:boundary_ids)if(q!=root)boundary_order.push_back(i);
}

bool Search::boundary_canonical_feasible(const State& s)const {
  if(boundary_symmetries.size()<2||!s.boundary_changed)return true;
  const int boundary_count=int(boundary.points.size());
  auto orient=[](const Hex& c,Quad q) {
    Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};
    for(int j=0;j<4;j++)for(auto e:EI) {
      int a=c[e[0]],b=c[e[1]];if(b==h[j])std::swap(a,b);
      if(a==h[j]&&std::find(h.begin(),h.begin()+4,b)==h.begin()+4)h[j+4]=b;
    }
    if(std::find(h.begin(),h.end(),-1)!=h.end())throw std::runtime_error("Cannot orient boundary owner");
    return h;
  };
  for(const auto& g:boundary_symmetries) {
    std::array<int,MAXV> left,right;left.fill(-1);right.fill(-1);
    int next_left=boundary_count,next_right=boundary_count;bool decided=false;
    for(int i:boundary_order) {
      int owner=s.boundary_owner[i],other=s.boundary_owner[g.source[i]];
      // Unknown earlier owners make later lexicographic comparisons unsafe.
      if(owner<0||other<0)break;
      Hex a=orient(s.cells[owner],boundary.boundary[i]),b=s.cells[other];
      for(int& v:b)if(v<boundary_count)v=g.permutation[v];
      if(g.odd){std::swap(b[1],b[3]);std::swap(b[5],b[7]);}
      b=orient(b,boundary.boundary[i]);
      for(int j=0;j<8;j++) {
        int x=a[j],y=b[j];
        if(x>=boundary_count){if(left[x]<0)left[x]=next_left++;x=left[x];}
        if(y>=boundary_count){if(right[y]<0)right[y]=next_right++;y=right[y];}
        if(y<x)return false;
        if(y>x){decided=true;break;}
      }
      if(decided)break;
    }
  }
  return true;
}

Quad Search::select_face(const State& s) const {
  // Every completion owns this prescribed face. Choosing it first changes
  // traversal only; prepare_boundary_symmetry uses the same face to order
  // its lexicographic comparisons, without imposing symmetry on the mesh.
  if(root_face>=0&&root_face<int(boundary.boundary.size())&&s.cells.empty())return boundary.boundary[root_face];
  if(face_order==4&&s.cover_choice_valid)return s.cover_choice;
  if(face_order==0)return s.front.begin()->second;
  std::array<int,MAXV> degree{};
  if(face_order>=2)for(auto c:s.cells)for(int v:c)++degree[v];
  // Estimate the four transverse vertex domains. This only chooses which
  // front face to fill next; it never removes a face or a vertex choice.
  if(face_order==5) {
    VertexDomains vertices(s);auto diagonals=known_diagonals(s,boundary);
    Quad best=s.front.begin()->second;std::array<long long,3> best_score{LLONG_MIN,LLONG_MIN,LLONG_MIN};
    for(auto [k,q]:s.front) {
      (void)k;long long product=1,sum=0;int d=0;for(int v:q)d+=degree[v];
      for(int r=0;r<4;r++) {
        Hex h{q[r],q[(r+3)%4],q[(r+2)%4],q[(r+1)%4],-1,-1,-1,-1};
        auto options=vertices.choices(h,4,int(s.a.size()),s.group_size,diagonals);
        auto color=vertices.colors[vertices.color[h[0]]^1];for(int g=0;g<s.group_size;g++)color.set(int(s.a.size())+g);options.intersect(color);
        if(corner_geometry&&!geometry.domains.empty())options.intersect(geometry.choices(h,4));
        int count=options.count();product*=count;sum+=count;
      }
      std::array<long long,3> score{-product,d,-sum};
      if(score>best_score){best_score=score;best=q;}
    }
    return best;
  }
  Quad best=s.front.begin()->second;std::pair<int,int> best_score{-1,-1};
  for(auto [key,q]:s.front) {
    (void)key;auto vertices=sorted(q);int fixes=0,d=0;
    for(int g=1;g<s.group_size;g++){Quad image=q;for(int& v:image)v=s.a[v][g];fixes+=sorted(image)==vertices;}
    for(int v:q)d+=degree[v];
    auto score=face_order>=3?std::make_pair(d,fixes):std::make_pair(fixes,face_order==1?0:d);
    if(score>best_score){best_score=score;best=q;}
  }
  return best;
}
State Search::initial() {
  State s;if(boundary_symmetries.size()>1)s.boundary_owner.assign(boundary.boundary.size(),-1);s.half_turn=half_turn;s.group_size=group_size();s.a=actions(boundary.points,axial,s.group_size,half_turn);
  for(auto q:boundary.boundary) {
    if(!s.front.emplace(key(q),q).second)throw std::runtime_error("Duplicate boundary quad");
    for(int i=0;i<4;i++)for(int j=i+1;j<4;j++) {
      int t=((j-i)==2)?2:1;
      s.relation[q[i]*MAXV+q[j]]|=t;s.relation[q[j]*MAXV+q[i]]|=t;
    }
  }
  for(auto q:boundary.boundary)for(int g=0;g<s.group_size;g++) {
    Quad r=q;for(int& x:r)x=s.a[x][g];if(s.odd(g))r=reverse(r);
    auto it=s.front.find(key(r));if(it==s.front.end()||!same_orientation(it->second,r))throw std::runtime_error("Boundary orientation/reflection mismatch");
  }
  s.component_upper=component_index?front_component_upper(s):front_component_upper_reference(s);
  Parity colors;colors.initialize(s);
  for(int g=0;g<s.group_size;g++) {
    auto [root,p]=colors.root(0);auto [image_root,q]=colors.root(s.a[0][g]);
    if(root!=image_root)throw std::runtime_error("Boundary edge graph must be connected");
    s.symmetry_flip[g]=p^q;
    for(int v=0;v<int(s.a.size());v++) {
      auto [r,a]=colors.root(v);auto [t,b]=colors.root(s.a[v][g]);
      if(r!=root||t!=root||(a^b)!=s.symmetry_flip[g])throw std::runtime_error("Inconsistent boundary symmetry colors");
    }
  }
  s.symmetry_parity=symmetry_parity&&s.group_size>1;s.mirror_separation=mirror_separation&&s.group_size>1&&!half_turn;
  for(auto& side:s.boundary_side)side.fill(-1);
  for(int axis=0;axis<reflection_count(s.group_size,s.half_turn);axis++)for(int v=0;v<int(s.a.size());v++)if(s.a[v][1<<axis]!=v) {
    auto p=boundary.points[v];double coordinate=axial?p[axis]:(axis==0?p[0]-p[1]:p[0]+p[1]);
    if(coordinate==0)throw std::runtime_error("Nonfixed boundary vertex lies on mirror");
    s.boundary_side[axis][v]=coordinate>0?1:0;
  }
  return s;
}
// Mayer-Vietoris: adding a contractible cell can increase beta_2 of
// (prescribed boundary + placed cells) by at most beta_1 of its intersection.
// Component boundary cycles give the corresponding complement-component count.
bool Search::completion_feasible(const State& s,const Hex& h,const std::array<uint8_t,MAXV*MAXV>& placed,const Budget& budget) const {
  std::vector<Hex> added;
  for(int g=0;g<s.group_size;g++) {
    Hex r=reflect(h,s.a,g,s.half_turn);bool duplicate=false;
    for(auto c:added)if(sorted(c)==sorted(r)){if(hexkey(c)!=hexkey(r))return true;duplicate=true;break;}
    if(!duplicate)added.push_back(r);
  }
  int n=int(added.size());if(n>budget.remaining)return false;
  std::array<std::array<Quad,6>,4> keys;
  for(int i=0;i<n;i++)for(int f=0;f<6;f++)keys[i][f]=key(face(added[i],f));
  int components=s.component_upper,covered=0,internal=0;
  for(int i=0;i<n;i++) {
    auto c=added[i];std::array<int,8> parent;std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int x){while(parent[x]!=x)x=parent[x];return x;};
    int vertices=0,edges=0,faces=0,connected=0;
    for(int v:c) {
      bool shared=v<budget.vertices;
      for(int j=0;j<i&&!shared;j++)shared=std::find(added[j].begin(),added[j].end(),v)!=added[j].end();
      vertices+=shared;
    }
    connected=vertices;
    for(auto& e:EI) {
      bool shared=placed[c[e[0]]*MAXV+c[e[1]]]&1;
      for(int j=0;j<i&&!shared;j++)shared=isedge(added[j],c[e[0]],c[e[1]]);
      if(shared){++edges;int a=root(e[0]),b=root(e[1]);if(a!=b){parent[a]=b;--connected;}}
    }
    for(int f=0;f<6;f++) {
      auto k=keys[i][f];bool old_front=s.front.count(k),old_face=old_front||s.incidence.count(k),previous=false;
      for(int j=0;j<i&&!previous;j++)for(int g=0;g<6;g++)if(keys[j][g]==k){previous=true;break;}
      faces+=old_face||previous;covered+=old_front;internal+=previous;
    }
    // A subcomplex of the cube boundary has beta_2=1 only if all six
    // faces are present. Clamping is conservative for incompatible proposals.
    components+=std::max(0,edges-vertices+connected-faces+(faces==6));
  }
  int front_after=budget.open+6*n-2*covered-2*internal;
  return front_after<=4*(budget.remaining-n)+2*components;
}
bool Search::insert(State& s,const Hex& h) {
  auto reject=[&](int reason){audit_add(s,AC_INSERT_ORBIT_SHAPE+reason);local().insertion_rejections[reason].fetch_add(1,std::memory_order_relaxed);return false;};
  std::vector<Hex> added;
  for(int g=0;g<s.group_size;g++) {
    Hex r=reflect(h,s.a,g,s.half_turn);bool duplicate=false;
    for(auto& c:added)if(sorted(c)==sorted(r)) {if(hexkey(c)!=hexkey(r))return reject(0);duplicate=true;break;}
    if(!duplicate)added.push_back(r);
  }
  if(s.cells.size()+added.size()>size_t(cap))return reject(1);
  if(location_domains&&!geometry.impose_locations(s,added)){(audit_add(s,AC_REJECT_LOCATION,1),local().reject_location.fetch_add(1,std::memory_order_relaxed));return false;}
  for(size_t i=0;i<added.size();i++) {
    for(auto q:boundary.boundary)if(!compatible_boundary(added[i],q))return reject(2);
    for(auto& c:s.cells)if(!compatible(added[i],c))return reject(3);
    for(size_t j=0;j<i;j++)if(!compatible(added[i],added[j]))return reject(4);
  }
  s.boundary_changed=false;
  for(auto& c:added)for(int f=0;f<6;f++) {
    Quad q=face(c,f),k=key(q);auto it=s.front.find(k);
    if(!s.boundary_owner.empty()) {
      auto boundary_face=boundary_ids.find(k);
      if(boundary_face!=boundary_ids.end()){s.boundary_owner[boundary_face->second]=int(s.cells.size()+(&c-added.data()));s.boundary_changed=true;}
    }
    int& occurrences=s.incidence[k];
    if(it!=s.front.end()) {
      if(!same_orientation(q,it->second))return reject(5);
      s.front.erase(it);
    } else {
      if(occurrences)return reject(5); // a closed interior or covered boundary face
      s.front.emplace(k,reverse(q));
    }
    ++occurrences;
  }
  for(auto& c:added) {
    for(int i=0;i<8;i++)for(int j=i+1;j<8;j++) {
      int t=relation_type(i,j);s.relation[c[i]*MAXV+c[j]]|=t;s.relation[c[j]*MAXV+c[i]]|=t;
    }
    s.cells.push_back(c);
  }
  ++s.orbit_count;update_root_floor(s);
  if(!boundary_canonical_feasible(s)){(audit_add(s,AC_CANONICAL_REJECTIONS,1),local().canonical_rejections.fetch_add(1,std::memory_order_relaxed));return false;}
  // Every remaining cell can consume at most six open faces. This deliberately
  // weak bound remains valid for disconnected or non-shellable remainders.
  if(s.front.size()>6*(size_t(cap)-s.cells.size()))return reject(6);
  if(component_bound) {
    s.component_upper=component_index?front_component_upper(s):front_component_upper_reference(s);
    if(int(s.front.size())>4*(cap-int(s.cells.size()))+2*s.component_upper) {
      (audit_add(s,AC_COMPONENT_REJECTIONS,1),local().component_rejections.fetch_add(1,std::memory_order_relaxed));return false;
    }
  }
  // A cubical ball has a bipartite edge graph. Detect incompatible odd cycles.
  std::vector<int> colors(s.a.size(),-1);
  for(size_t root=0;root<s.a.size();root++)if(colors[root]<0) {
    std::vector<int> todo{int(root)};colors[root]=0;
    while(!todo.empty()) {int v=todo.back();todo.pop_back();for(size_t w=0;w<s.a.size();w++)if(s.relation[v*MAXV+w]&1) {
      if(colors[w]<0){colors[w]=colors[v]^1;todo.push_back(int(w));}
      else if(colors[w]==colors[v])return reject(7);
    }}
  }
  if(s.symmetry_parity) {
    Parity parity;
    for(int v=0;v<int(s.a.size());v++)for(int w=0;w<v;w++)if(s.relation[v*MAXV+w]&1)
      if(!parity.edge(v,w))return reject(7);
    if(!parity.symmetry(s))return reject(7);
  }
  if(s.mirror_separation)for(int axis=0;axis<reflection_count(s.group_size,s.half_turn);axis++) {
    Parity sides;if(!side_parity(s,axis,sides)){(audit_add(s,AC_REJECT_SIDE,1),local().reject_side.fetch_add(1,std::memory_order_relaxed));return false;}
  }
  return true;
}

int rank_gf2(const std::vector<std::vector<int>>& columns) {
  std::array<std::bitset<512>,512> pivots{};int rank=0;
  for(auto& c:columns) {
    std::bitset<512> x;for(int r:c){if(r>=512)throw std::runtime_error("Homology capacity exceeded");x.flip(r);}
    for(int r=511;r>=0;r--)if(x[r]) {
      if(pivots[r].none()){pivots[r]=x;++rank;break;}x^=pivots[r];
    }
  }
  return rank;
}
bool Search::topology(const State& s)const {
  std::map<std::pair<int,int>,int> edges;
  std::map<Quad,int> faces;
  std::vector<Quad> face_list;
  std::vector<std::vector<int>> owners;
  for(size_t i=0;i<s.cells.size();i++) {
    auto h=s.cells[i];
    for(auto& e:EI){auto k=edgekey(h[e[0]],h[e[1]]);if(!edges.count(k))edges[k]=int(edges.size());}
    for(int f=0;f<6;f++){auto k=key(face(h,f));if(!faces.count(k)){faces[k]=int(faces.size());face_list.push_back(k);owners.emplace_back();}owners[faces[k]].push_back(int(i));}
  }
  const int V=int(s.a.size()),E=int(edges.size()),F=int(faces.size()),H=int(s.cells.size());
  if(V-E+F-H!=1)return false;
  std::vector<std::vector<int>> d1,d2,d3;
  for(auto [e,index]:edges){(void)index;d1.push_back({e.first,e.second});}
  for(auto q:face_list){std::vector<int> c;for(int i=0;i<4;i++)c.push_back(edges.at(edgekey(q[i],q[(i+1)%4])));d2.push_back(c);}
  for(auto h:s.cells){std::vector<int> c;for(int f=0;f<6;f++)c.push_back(faces.at(key(face(h,f))));d3.push_back(c);}
  int r1=rank_gf2(d1),r2=rank_gf2(d2),r3=rank_gf2(d3);
  if(V-r1!=1||E-r1-r2!=0||F-r2-r3!=0||H-r3!=0)return false;
  // These are necessary conditions, not a full 3-ball recognition algorithm.
  return true;
}
void Search::emit(const State& s) {
  if(!topology(s)){(audit_add(s,AC_REJECTED,1),local().rejected.fetch_add(1,std::memory_order_relaxed));return;}
  (audit_add(s,AC_CANDIDATES,1),local().candidates.fetch_add(1,std::memory_order_relaxed));
  std::lock_guard<std::mutex> lock(io);
  out<<"{\"hexes\":"<<s.cells.size()<<",\"vertices\":"<<s.a.size()<<",\"cell_orbits\":"<<s.orbit_count<<",\"vertex_actions\":[";
  for(size_t i=0;i<s.a.size();i++){if(i)out<<',';out<<'[';for(int g=0;g<s.group_size;g++){if(g)out<<',';out<<s.a[i][g];}out<<']';}
  out<<"],\"cells\":[";
  for(size_t i=0;i<s.cells.size();i++){if(i)out<<',';out<<'[';for(int j=0;j<8;j++){if(j)out<<',';out<<s.cells[i][j];}out<<']';}
  out<<"]}\n";out.flush();if(!out)throw std::runtime_error("Candidate output write failed");
}

void Search::choose(State& s,Hex& h,int position,int depth,const std::array<uint8_t,MAXV*MAXV>& placed,const Diagonals& diagonals,const Rules& rules,Overlaps overlaps,const Budget& budget,FaceDomains domains,Parity& parity,const VertexDomains& vertices,std::array<ChoiceMask,2> colors) {
  if(stopping())return;
  if(proposal_stop&&proposal_stop())return;
  (audit_add(s,AC_BRANCHES,1),local().branches.fetch_add(1,std::memory_order_relaxed));
  if(edge_completion&&component_bound&&early_completion&&position>=5&&!edge_completion_feasible(s,h,position,placed,budget,overlaps,vertices,completion_domains)) {
    (audit_add(s,AC_PARTIAL_COMPLETION,1),local().partial_completion.fetch_add(1,std::memory_order_relaxed));return;
  }
  if(partial_cover&&face_cover&&s.group_size<=2&&position==partial_cover&&
     !partial_cover_feasible(s,h,position,budget,domains,diagonals,placed))return;
  if(position==8) {
    if(component_bound&&early_completion&&!completion_feasible(s,h,placed,budget)){(audit_add(s,AC_EARLY_COMPLETION,1),local().early_completion.fetch_add(1,std::memory_order_relaxed));return;}
    (audit_add(s,AC_TRIALS,1),local().trials.fetch_add(1,std::memory_order_relaxed));
    State next=s;
    if(!insert(next,h))return;
    audit_add(s,AC_CHILDREN);
    if(proposal_sink){proposal_sink(std::move(next));return;}
    // Bounded dynamic task splitting: each task owns its complete search state.
    int pending=tasks.load(std::memory_order_relaxed);
    if((task_depth==0||depth<task_depth)&&threads>1&&pending<threads*4&&tasks.compare_exchange_strong(pending,pending+1)) {
      #pragma omp task firstprivate(next,depth)
      {
        try{dfs(std::move(next),depth+1);}catch(const std::exception& e){error(e);}
        tasks.fetch_sub(1,std::memory_order_relaxed);
      }
    } else dfs(std::move(next),depth+1);
    return;
  }
  // The snapshot distinguishes an existing cell's body diagonal from repeated
  // occurrences within the candidate orbit (possibly the same cell).
  const int existing=int(s.a.size());
  int forced=-1;
  for(auto r:rules)if(r[1]==position&&r[0]<position) {
    int v=s.a[h[r[0]]][r[2]];if(forced>=0&&forced!=v){audit_add(s,AC_RULE_CONFLICT);return;}forced=v;
  }
  ChoiceMask options;options.limit(existing+s.group_size);
  if(vertex_domains)options=vertices.choices(h,position,existing,s.group_size,diagonals);
  int required_color=vertices.color[h[0]]^(__builtin_popcount(unsigned(BITS[position]))&1);
  if(color_domains) {
    auto valid=colors[required_color];for(int type=0;type<s.group_size;type++)valid.set(existing+type);
    options.intersect(valid);
  }
  if(vertex_domains||color_domains)(audit_add(s,AC_EARLY_DOMAIN,existing+s.group_size-options.count()),local().early_domain.fetch_add(existing+s.group_size-options.count(),std::memory_order_relaxed));
  if(corner_geometry&&!geometry.domains.empty()) {
    int before=options.count();options.intersect(geometry.choices(h,position));
    (audit_add(s,AC_EARLY_GEOMETRY,before-options.count()),local().early_geometry.fetch_add(before-options.count(),std::memory_order_relaxed));
  }
  if(overlap_domains&&early_overlap) {
    int before=options.count();options.intersect(overlap_choices(s,h,position,overlaps,vertices));
    (audit_add(s,AC_EARLY_OVERLAP,before-options.count()),local().early_overlap.fetch_add(before-options.count(),std::memory_order_relaxed));
  }
  std::vector<int> ordered;size_t cursor=0;
  if(proposal_order){for(int option=options.pop();option>=0;option=options.pop())ordered.push_back(option);proposal_order(ordered);}
  auto pop=[&](){return proposal_order?(cursor<ordered.size()?ordered[cursor++]:-1):options.pop();};
  for(int option=pop();option>=0&&!stopping();option=pop()) {
    if(proposal_stop&&proposal_stop())break;
    if(forced>=0&&option!=forced){audit_add(s,AC_FORCED_CHOICE);continue;}
    int candidate=option;auto next_colors=colors;
    if(option>=existing) {
      int type=option-existing,size=s.group_size==1?1:s.group_size==2?(type==0?2:1):(type==0?4:type==3?1:2);
      if(existing+size>vertex_cap){audit_add(s,AC_VERTEX_CAP);continue;}
      candidate=allocate(s,type);
      for(int g=0;g<s.group_size;g++)next_colors[required_color^s.symmetry_flip[g]].set(s.a[candidate][g]);
    }
    h[position]=candidate;Trail trail;
    if(corner_geometry&&geometry.domains.empty()&&!geometry.feasible(s,h,position+1,position)) {(audit_add(s,AC_EARLY_GEOMETRY,1),local().early_geometry.fetch_add(1,std::memory_order_relaxed));s.a.resize(existing);continue;}
    auto next_overlaps=overlaps;
    if(s.mirror_separation&&candidate<vertices.initial)for(int axis=0;axis<reflection_count(s.group_size,s.half_turn);axis++) {
      int g=1<<axis;
      for(int j=0;j<position;j++)if(h[j]<vertices.initial&&((vertices.opposite[axis][h[j]].word[candidate/64]>>(candidate%64))&1))next_overlaps.fixed[g-1]=true;
    }
    auto next_domains=domains;
    if(early_overlap&&!extend_overlaps(s,h,position,next_overlaps)) {
      (audit_add(s,AC_EARLY_OVERLAP,1),local().early_overlap.fetch_add(1,std::memory_order_relaxed));s.a.resize(existing);continue;
    }
    if(early_budget) {
      budget.extend(h,position,next_domains);
      if(!budget.feasible(next_overlaps,next_domains)){(audit_add(s,AC_EARLY_BUDGET,1),local().early_budget.fetch_add(1,std::memory_order_relaxed));s.a.resize(existing);continue;}
    }
    auto checkpoint=parity.checkpoint();
    if(early_parity&&!parity.extend(s,h,position)) {
      (audit_add(s,AC_EARLY_PARITY,1),local().early_parity.fetch_add(1,std::memory_order_relaxed));parity.rollback(checkpoint);s.a.resize(existing);continue;
    }
    if(!consistent_rules(s,h,position+1,rules))audit_add(s,AC_RULE_CONFLICT);
    else if(!impose(s,h,position,placed,trail))audit_add(s,AC_RELATIONS);
    else if(!(incremental_faces?consistent_new_faces(h,position,diagonals,placed):consistent_faces(s,h,position+1,diagonals)))audit_add(s,AC_FACE_CONSISTENCY);
    else choose(s,h,position+1,depth,placed,diagonals,rules,next_overlaps,budget,next_domains,parity,vertices,next_colors);
    parity.rollback(checkpoint);undo(s,trail);s.a.resize(existing);
  }
}
void Search::dfs(State s,int depth) {
  if(stopping())return;
  AuditCounters audit{};uint64_t parent=s.audit_id;
  if(!tree_trace.empty()){s.audit=&audit;s.audit_id=(++local().audit_serial<<8)|uint64_t(omp_get_thread_num());}
  auto record=[&](unsigned outcome) {
    if(tree_trace.empty())return;
    uint64_t hex=0;if(!s.cells.empty())for(int j=0;j<8;j++)hex|=uint64_t(s.cells.back()[j]+1)<<(8*j);
    uint64_t meta=uint64_t(s.cells.size())|(uint64_t(s.a.size())<<8)|(uint64_t(s.front.size())<<16)|(uint64_t(outcome)<<32)|(uint64_t(s.audit_cover_reason)<<40);
    std::array<uint64_t,6> row{s.audit_id,parent,meta,hex,audit[AC_COVER_WORK],audit[AC_BRANCHES]};
    local().audit_states.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    if(outcome==3||outcome==6){local().audit_expansions.write(reinterpret_cast<const char*>(&s.audit_id),8);local().audit_expansions.write(reinterpret_cast<const char*>(audit.data()),sizeof(audit));}
    if(!local().audit_states||!local().audit_expansions)throw std::runtime_error("Tree trace write failed");
  };
  (audit_add(s,AC_STATES,1),local().states.fetch_add(1,std::memory_order_relaxed));
  if(!state_trace.empty()) {
    std::ostringstream line;line<<"a";for(auto a:s.a)for(int v:a)line<<' '<<v;
    line<<" h";for(auto h:s.cells)for(int v:h)line<<' '<<v;line<<'\n';
    std::lock_guard<std::mutex> lock(io);trace<<line.str();if(!trace)throw std::runtime_error("State trace write failed");
  }
  if(s.front.empty()){emit(s);record(audit[AC_REJECTED]?5:4);return;}
  if(s.cells.size()>=size_t(cap)||s.front.size()>6*(size_t(cap)-s.cells.size())){record(2);return;}
  if(face_cover&&!cover_feasible(s)){(audit_add(s,AC_COVER_REJECTIONS,1),local().cover_rejections.fetch_add(1,std::memory_order_relaxed));record(1);return;}
  Quad q=select_face(s);
  expand_face(s,q,depth);record(stopping()?6:3);
}
void Search::expand_face(State& s,Quad q,int depth) {
  Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};
  if(corner_geometry&&!geometry.feasible(s,h,4)){(audit_add(s,AC_EARLY_GEOMETRY,1),local().early_geometry.fetch_add(1,std::memory_order_relaxed));return;}
  // Existing front-face relations are already known. No cell is inserted until
  // all four images have passed compatibility checks as a single transaction.
  auto placed=s.relation;
  auto diagonals=known_diagonals(s,boundary);
  auto rules=stabilizer_rules(s,h);
  Overlaps overlaps{};
  if(early_overlap)for(int j=0;j<4;j++)if(!extend_overlaps(s,h,j,overlaps)){audit_add(s,AC_INITIAL_OVERLAP);return;}
  Budget budget(s,cap,face_domains);auto domains=budget.initial;
  if(early_budget){for(int j=0;j<4;j++)budget.extend(h,j,domains);if(!budget.feasible(overlaps,domains)){audit_add(s,AC_INITIAL_BUDGET);return;}}
  Parity parity;if(early_parity)parity.initialize(s);
  VertexDomains vertices(s);
  choose(s,h,4,depth,placed,diagonals,rules,overlaps,budget,domains,parity,vertices,vertices.colors);
}

// Necessary regular intersection test for two incompletely specified cubes.
// Unknown slots can acquire any distinct vertex; only known matches constrain it.

struct FaceMasks {
  std::array<unsigned,256> containing{};std::array<bool,256> edge{};
  FaceMasks(){for(int m=0;m<256;m++)for(int f=0;f<6;f++){unsigned mask=0;for(int v:FI[f])mask|=1u<<v;if(!(m&~mask))containing[m]|=1u<<f;}
    for(auto e:EI)edge[(1u<<e[0])|(1u<<e[1])]=true;}
};
const FaceMasks face_masks;
bool partial_cells_compatible_uncached(const Hex& a,const Hex& b) {
  unsigned ma=0,mb=0;std::array<int,8> mapping;mapping.fill(-1);int count=0;
  for(int i=0;i<8;i++)if(a[i]>=0)for(int j=0;j<8;j++)if(a[i]==b[j]){ma|=1u<<i;mb|=1u<<j;mapping[i]=j;++count;break;}
  if(count<=1)return true;
  if(count>4)return false;
  if(count==2&&face_masks.edge[ma]&&face_masks.edge[mb])return true;
  unsigned fa=face_masks.containing[ma],fb=face_masks.containing[mb];
  if(!fa||!fb)return false;
  int first=__builtin_ctz(ma);
  while(fa) {
    int f=__builtin_ctz(fa);fa&=fa-1;int pos_a=0;while(FI[f][pos_a]!=first)++pos_a;
    auto available=fb;
    while(available) {
      int g=__builtin_ctz(available);available&=available-1;int pos_b=0;while(FI[g][pos_b]!=mapping[first])++pos_b;
      int rotation=(pos_a+pos_b)%4;bool valid=true;
      for(int k=0;k<4;k++) {
        int ia=FI[f][k],ib=FI[g][(rotation-k+4)%4];
        if((mapping[ia]>=0&&mapping[ia]!=ib)||(a[ia]>=0&&b[ib]>=0&&a[ia]!=b[ib])){valid=false;break;}
      }
      if(valid)return true;
    }
  }
  return false;
}
// This predicate depends only on the two partial cubes, so the exact-key
// cache remains valid across states. Negative slots all mean unspecified.
bool partial_cells_compatible(const Hex& a,const Hex& b) {
  uint64_t x=0,y=0;
  for(int i=0;i<8;i++) {
    if(a[i]>=255||b[i]>=255)return partial_cells_compatible_uncached(a,b);
    x|=uint64_t(std::max(0,a[i]+1))<<(8*i);y|=uint64_t(std::max(0,b[i]+1))<<(8*i);
  }
  if(x>y)std::swap(x,y);
  struct Entry {uint64_t a=0,b=0;uint8_t result=0;};
  thread_local std::array<Entry,16384> cache;
  uint64_t hash=x^(y*0x9e3779b97f4a7c15ULL);hash^=hash>>31;hash*=0xbf58476d1ce4e5b9ULL;hash^=hash>>29;
  auto& entry=cache[hash&16383];
  if(entry.result&&entry.a==x&&entry.b==y)return entry.result==1;
  bool result=partial_cells_compatible_uncached(a,b);entry={x,y,uint8_t(result?1:2)};return result;
}
bool partial_cells_compatible_reference(const Hex& a,const Hex& b) {
  unsigned ma=0,mb=0;int count=0;
  for(int i=0;i<8;i++)if(a[i]>=0)for(int j=0;j<8;j++)if(a[i]==b[j]){ma|=1u<<i;mb|=1u<<j;++count;}
  if(count<=1)return true;
  if(count>4)return false;
  if(count==2) {
    bool ea=false,eb=false;
    for(auto e:EI){unsigned m=(1u<<e[0])|(1u<<e[1]);ea|=m==ma;eb|=m==mb;}
    if(ea&&eb)return true;
  }
  for(const auto& fa:FI) {
    unsigned maska=0;for(int i:fa)maska|=1u<<i;if(ma&~maska)continue;
    for(const auto& fb:FI) {
      unsigned maskb=0;for(int i:fb)maskb|=1u<<i;if(mb&~maskb)continue;
      for(int rotation=0;rotation<4;rotation++) {
        bool ok=true;
        for(int k=0;k<4;k++) {
          int ia=fa[k],ib=fb[(rotation-k+4)%4];
          if(a[ia]>=0&&b[ib]>=0&&a[ia]!=b[ib]){ok=false;break;}
          if(a[ia]>=0)for(int j=0;j<8;j++)if(b[j]==a[ia]&&j!=ib)ok=false;
          if(b[ib]>=0)for(int j=0;j<8;j++)if(a[j]==b[ib]&&j!=ia)ok=false;
        }
        if(ok)return true;
      }
    }
  }
  return false;
}

bool partial_orbit_support(const State& s,const Hex& h) {
  int mask=0;for(int i=0;i<8;i++)if(h[i]>=0) {
    mask|=1<<i;
    if((s.cover_color[h[i]]^s.cover_color[h[0]])!=(__builtin_popcount(unsigned(BITS[i]))&1))return false;
  }
  for(int g=1;g<s.group_size;g++) {
    uint32_t code=0;
    for(int i=0;i<8;i++)if(h[i]>=0)for(int j=0;j<8;j++)if(h[j]>=0&&s.a[h[i]][g]==h[j]){code|=uint32_t(j+1)<<(4*i);break;}
    const auto& allowed=s.mirror_separation&&g!=3?overlap_patterns.reflection_masked[mask]:overlap_patterns.masked[s.odd(g)][mask];
    if(!allowed.count(code))return false;
  }
  return true;
}


struct CoverMemo {
  struct FixedEntry {uint64_t generation=0;Hex h{};bool result=false;};
  std::array<FixedEntry,4096> fixed;
  uint64_t fixed_generation=0;
  std::array<uint64_t,MAXV> fixed_incidence{};
  struct FaceEntry {uint64_t generation=0;Hex h{};bool result=false;};
  struct GroupEntry {uint64_t generation=0;Quad q{};std::vector<Hex> group,result;};
  std::array<FaceEntry,4096> face;
  std::array<GroupEntry,1024> group;
  uint64_t serial=0,generation=0;
  const State* state=nullptr;const Diagonals* diagonals=nullptr;
};
thread_local CoverMemo cover_memo;
uint64_t hex_hash(const Hex& h) {
  uint64_t a=0,b=0;for(int i=0;i<8;i++)(i<4?a:b)|=uint64_t(h[i]+1)<<(16*(i%4));
  a^=b*0x9e3779b97f4a7c15ULL;a^=a>>31;a*=0xbf58476d1ce4e5b9ULL;return a^(a>>29);
}
struct CoverMemoScope {
  const State* old_state;const Diagonals* old_diagonals;uint64_t old_generation;
  CoverMemoScope(const State& s,const Diagonals& d):old_state(cover_memo.state),old_diagonals(cover_memo.diagonals),old_generation(cover_memo.generation) {
    cover_memo.state=&s;cover_memo.diagonals=&d;cover_memo.generation=++cover_memo.serial;
  }
  ~CoverMemoScope(){cover_memo.state=old_state;cover_memo.diagonals=old_diagonals;cover_memo.generation=old_generation;}
};
bool partial_face_support_uncached(const State& s,const Hex& h,const Diagonals& diagonals) {
  ChoiceMask available;for(int v:h)if(v>=0)available.reset(v);
  std::array<ChoiceMask,8> domains;for(int i=0;i<8;i++)if(h[i]<0)domains[i]=available;
  for(auto f:FI)for(int k=0;k<2;k++) {
    int a=h[f[k]],b=h[f[k+2]];
    if(a<0||b<0||!(s.relation[a*MAXV+b]&2))continue;
    auto it=diagonals.find(edgekey(a,b));if(it==diagonals.end())throw std::runtime_error("Missing cover diagonal");
    std::array<int,2> other{};int count=0;for(int v:it->second)if(v!=a&&v!=b)other[count++]=v;
    for(int j:{f[(k+1)%4],f[(k+3)%4]}) {
      if(h[j]>=0){if(h[j]!=other[0]&&h[j]!=other[1])return false;}
      else {ChoiceMask required;required.word={0,0};for(int v:other)required.set(v);domains[j].intersect(required);}
    }
  }
  int previous=-1;
  for(int i=0;i<8;i++)if(h[i]<0) {
    int count=domains[i].count();if(!count)return false;
    if(previous>=0&&count==1&&domains[previous].count()==1&&domains[previous].word==domains[i].word)return false;
    previous=i;
  }
  return true;
}
bool partial_face_support(const State& s,const Hex& h,const Diagonals& diagonals) {
  auto& cache=cover_memo;
  if(cache.state!=&s||cache.diagonals!=&diagonals)return partial_face_support_uncached(s,h,diagonals);
  auto& entry=cache.face[hex_hash(h)&4095];
  if(entry.generation==cache.generation&&entry.h==h)return entry.result;
  bool result=partial_face_support_uncached(s,h,diagonals);
  entry={cache.generation,h,result};return result;
}
bool partial_face_support_reference(const State& s,const Hex& h,const Diagonals& diagonals) {
  std::array<ChoiceMask,8> domain;
  for(int i=0;i<8;i++) {
    if(h[i]>=0){domain[i].word={0,0};domain[i].set(h[i]);}
    else for(int v:h)if(v>=0)domain[i].reset(v);
  }
  for(auto f:FI)for(int k=0;k<2;k++) {
    int a=h[f[k]],b=h[f[k+2]];
    if(a<0||b<0||!(s.relation[a*MAXV+b]&2))continue;
    auto it=diagonals.find(edgekey(a,b));if(it==diagonals.end())throw std::runtime_error("Missing cover diagonal");
    ChoiceMask required;required.word={0,0};for(int v:it->second)if(v!=a&&v!=b)required.set(v);
    domain[f[(k+1)%4]].intersect(required);domain[f[(k+3)%4]].intersect(required);
  }
  int unknown=-1;
  for(int i=0;i<8;i++) {
    if(!domain[i].count())return false;
    if(h[i]<0) {
      if(unknown>=0&&domain[i].count()==1&&domain[unknown].count()==1&&domain[i].word==domain[unknown].word)return false;
      unknown=i;
    }
  }
  return true;
}

// Any two old front faces assigned to one future cube must fit an
// adjacent or opposite pair of cube faces. Enumerate all oriented fits.

uint32_t Search::face_templates(const State& s,Quad q,Quad r)const {
  if(!fast_pairs||s.group_size>2)return face_templates_reference(s,q,r);
  int shape=0,multiplier=1,common=0;
  std::array<int,8> vertices;for(int j=0;j<4;j++){vertices[j]=q[j];vertices[j+4]=r[j];}
  for(int j=0;j<4;j++) {
    int digit=4;for(int i=0;i<4;i++)if(q[i]==r[j]){digit=i;++common;break;}
    shape+=multiplier*digit;multiplier*=5;
  }
  if(common!=0&&common!=2)return 0;
  const auto& pattern=pair_fits.patterns[shape];
  uint32_t result=pattern.parity[s.cover_color[q[0]]^s.cover_color[r[0]]];
  if(!result)return 0;
  // Input faces already have their prescribed edge/face-diagonal relations
  // and alternating colors. Only cross-face constraints can change the fits.
  for(int i=0;i<4&&result;i++)for(int j=0;j<4;j++) {
    int old=s.relation[q[i]*MAXV+r[j]];
    if(old)result&=pattern.relation[4*i+j][old];
  }
  if(!result)return 0;
  if(corner_geometry) {
    uint64_t boundary_key=0;for(int i=0;i<8;i++)boundary_key|=uint64_t(std::min(vertices[i],geometry.boundary_count))<<(8*i);
    uint64_t hash=boundary_key^(uint64_t(shape)*0x9e3779b97f4a7c15ULL);hash^=hash>>33;hash*=0xff51afd7ed558ccdULL;hash^=hash>>33;
    auto& cache=local().static_pair_cache;if(!cache)cache=std::make_unique<StaticPairEntry[]>(65536);
    auto& entry=cache[hash&65535];uint32_t allowed=0;
    if(cover_cache&&entry.shape==uint32_t(shape)&&entry.boundary==boundary_key)allowed=entry.value;
    else {
      for(const auto& fit:pattern.fits) {
        std::array<int,8> h;for(int i=0;i<8;i++)h[i]=fit.source[i]<0?geometry.boundary_count:std::min(vertices[fit.source[i]],geometry.boundary_count);
        bool valid=true;
        for(auto corner:geometry.corners) {
          std::array<int,4> ids;for(int j=0;j<4;j++)ids[j]=h[corner[j]];
          if(!geometry.allowed[geometry.index(ids)]){valid=false;break;}
        }
        if(valid)allowed|=uint32_t(1)<<fit.bit;
      }
      if(cover_cache)entry={boundary_key,uint32_t(shape),allowed};
    }
    result&=allowed;
  }
  if(s.canonical_floor||s.group_size>1)for(const auto& fit:pattern.fits)if(result&(uint32_t(1)<<fit.bit)) {
    Hex h;for(int i=0;i<8;i++)h[i]=fit.source[i]<0?-1:vertices[fit.source[i]];
    bool valid=partial_orbit_support(s,h)&&leader_feasible(s,h);
    if(corner_geometry)for(int g=1;g<s.group_size&&valid;g++)for(auto corner:geometry.corners) {
      std::array<int,4> ids;
      for(int j=0;j<4;j++){int v=h[corner[j]];ids[j]=v<0||v>=geometry.boundary_count?geometry.boundary_count:s.a[v][g];}
      if(s.odd(g))std::swap(ids[1],ids[2]);
      if(!geometry.allowed[geometry.index(ids)]){valid=false;break;}
    }
    if(!valid)result&=~(uint32_t(1)<<fit.bit);
  }
  return result;
}
uint32_t Search::face_templates_reference(const State& s,Quad q,Quad r)const {
  int common=0;for(int v:q)common+=std::find(r.begin(),r.end(),v)!=r.end();
  if(common!=0&&common!=2)return 0;
  CoverKey cache_key{};if(s.group_size<=2)cache_key[5]=s.canonical_floor;std::array<int,8> vertices;
  for(int i=0;i<4;i++){vertices[i]=q[i];vertices[i+4]=r[i];}
  if(cover_cache) {
    for(int i=0;i<8;i++) {
      cache_key[0]|=uint64_t(vertices[i]<int(boundary.points.size())?vertices[i]:255)<<(8*i);
      int label=i;for(int j=0;j<i;j++)if(vertices[j]==vertices[i]){label=j;break;}
      cache_key[1]|=uint64_t(label)<<(4*i);
    }
    for(int i=0;i<8;i++)cache_key[1]|=uint64_t(s.cover_color[vertices[i]]^s.cover_color[vertices[0]])<<(32+i);
    for(int g=1;g<s.group_size;g++)for(int i=0;i<8;i++) {
      int image=s.a[vertices[i]][g],target=0;for(int j=0;j<8;j++)if(vertices[j]==image){target=j+1;break;}
      cache_key[4+(g-1)/2]|=uint64_t(target)<<(((g-1)%2)*32+4*i);
    }
    int k=0;for(int i=0;i<8;i++)for(int j=0;j<i;j++) {
      cache_key[2+k/21]|=uint64_t(s.relation[vertices[i]*MAXV+vertices[j]])<<(3*(k%21));++k;
    }
    if(cover_flat) {
      auto& cache=local().flat_cache;
      if(!cache)cache=std::make_unique<CoverEntry[]>(size_t(1)<<cover_cache_bits);
      const auto& entry=cache[CoverHash{}(cache_key)&((size_t(1)<<cover_cache_bits)-1)];
      if(entry.value!=UINT32_MAX&&entry.key==cache_key){(audit_add(s,AC_CACHE_HITS,1),local().cache_hits.fetch_add(1,std::memory_order_relaxed));return entry.value;}
    } else {
      auto& cache=local().cover_cache;auto found=cache.find(cache_key);
      if(found!=cache.end()){(audit_add(s,AC_CACHE_HITS,1),local().cache_hits.fetch_add(1,std::memory_order_relaxed));return found->second;}
    }
    (audit_add(s,AC_CACHE_MISSES,1),local().cache_misses.fetch_add(1,std::memory_order_relaxed));
  }
  uint32_t result=0;
  auto save=[&](){if(cover_cache){
    if(cover_flat) {
      auto& entry=local().flat_cache[CoverHash{}(cache_key)&((size_t(1)<<cover_cache_bits)-1)];entry.key=cache_key;entry.value=result;
    } else {
      auto& cache=local().cover_cache;if(cache.size()>=131072)cache.clear();cache.emplace(cache_key,result);
    }
  }return result;};

  for(int f=1;f<6;f++) {
    if((f==1)!=(common==0))continue;
    for(int rotation=0;rotation<4;rotation++) {
      Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};bool valid=true;
      for(int j=0;j<4;j++) {
        int pos=FI[f][j],v=r[(j+rotation)%4];
        if(h[pos]>=0&&h[pos]!=v){valid=false;break;}h[pos]=v;
      }
      if(!valid)continue;
      for(int i=0;i<8&&valid;i++)if(h[i]>=0)for(int j=0;j<i;j++)if(h[j]>=0)
        if(h[i]==h[j]||conflict(s.relation[h[i]*MAXV+h[j]],relation_type(i,j),true)){valid=false;break;}
      if(!valid||!partial_orbit_support(s,h))continue;
      if(corner_geometry)for(int g=0;g<s.group_size&&valid;g++)for(auto corner:geometry.corners) {
        std::array<int,4> ids;
        for(int j=0;j<4;j++){int v=h[corner[j]];ids[j]=v<0||v>=geometry.boundary_count?geometry.boundary_count:s.a[v][g];}
        if(s.odd(g))std::swap(ids[1],ids[2]);
        if(!geometry.allowed[geometry.index(ids)]){valid=false;break;}
      }
      if(valid&&leader_feasible(s,h))result|=uint32_t(1)<<((f-1)*4+rotation);
    }
  }
  return save();
}
bool Search::faces_can_share(const State& s,Quad q,Quad r,const Diagonals& diagonals)const {
  auto templates=face_templates(s,q,r);
  while(templates) {
    int bit=__builtin_ctz(templates);templates&=templates-1;
    int f=bit/4+1,rotation=bit%4;Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};
    for(int j=0;j<4;j++)h[FI[f][j]]=r[(j+rotation)%4];
    if(partial_face_support(s,h,diagonals))return true;
  }
  return false;
}


// A future cell must have a regular, consistently oriented intersection
// with every placed cell. Pair relations and known face diagonals alone do
// not enforce this jointly. Unknown slots stay existentially unspecified.
// This predicate is memoized only within an immutable cover scope. It does
// not alter the face-pair graph, whose positive entries can be inherited.
bool Search::supports_fixed_mesh(const State& s,const Hex& h)const {
  auto& cache=cover_memo;
  if(cache.state!=&s)return supports_fixed_mesh_uncached(s,h);
  auto& entry=cache.fixed[hex_hash(h)&4095];
  if(entry.generation==cache.generation&&entry.h==h)return entry.result;
  bool result=supports_fixed_mesh_uncached(s,h);entry={cache.generation,h,result};return result;
}
bool Search::supports_fixed_mesh_uncached(const State& s,const Hex& h)const {
  auto test=[&](const Hex& k) {
    local().cover_fixed_checks.fetch_add(1,std::memory_order_relaxed);
    bool ok=partial_cells_compatible(h,k);
    if(!ok)local().cover_fixed_rejections.fetch_add(1,std::memory_order_relaxed);
    return ok;
  };
  auto& cache=cover_memo;
  if(cache.state==&s&&s.cells.size()<=64) {
    if(cache.fixed_generation!=cache.generation) {
      cache.fixed_incidence.fill(0);
      for(int i=0;i<int(s.cells.size());i++)for(int v:s.cells[i])if(v>=0)cache.fixed_incidence[v]|=1ULL<<i;
      cache.fixed_generation=cache.generation;
    }
    uint64_t seen=0,repeated=0;
    for(int v:h)if(v>=0){repeated|=seen&cache.fixed_incidence[v];seen|=cache.fixed_incidence[v];}
    while(repeated){int i=__builtin_ctzll(repeated);repeated&=repeated-1;if(!test(s.cells[i]))return false;}
    return true;
  }
  for(const auto& k:s.cells)if(!test(k))return false;
  return true;
}

std::vector<Hex> Search::extend_cube_group_raw(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const {
  if(group.empty())return {Hex{q[0],q[3],q[2],q[1],-1,-1,-1,-1}};
  std::vector<Hex> result;
  if(group.size()==1&&std::all_of(group[0].begin()+4,group[0].end(),[](int v){return v<0;})) {
    auto old=group[0];Quad first{old[0],old[3],old[2],old[1]};auto templates=face_templates(s,first,q);
    while(templates){int bit=__builtin_ctz(templates);templates&=templates-1;auto h=old;
      for(int j=0;j<4;j++)h[FI[bit/4+1][j]]=q[(j+bit%4)%4];
      if(partial_face_support(s,h,diagonals))result.push_back(h);
    }
    return result;
  }
  for(auto old:group) {
    if(std::find(old.begin(),old.end(),-1)==old.end()) {
      for(int f=0;f<6;f++)if(same_orientation(face(old,f),q)){result.push_back(old);break;}
      continue;
    }
    for(int f=0;f<6;f++)for(int rotation=0;rotation<4;rotation++) {
    Hex h=old;bool valid=true;
    for(int j=0;j<4;j++) {
      int pos=FI[f][j],v=q[(j+rotation)%4];
      if(h[pos]>=0&&h[pos]!=v){valid=false;break;}h[pos]=v;
    }
    if(!valid)continue;
    for(int i=0;i<8&&valid;i++)if(h[i]>=0)for(int j=0;j<i;j++)if(h[j]>=0)
      if(h[i]==h[j]||conflict(s.relation[h[i]*MAXV+h[j]],relation_type(i,j),true)){valid=false;break;}
    if(!valid||!partial_orbit_support(s,h)||!partial_face_support(s,h,diagonals))continue;
    if(corner_geometry)for(int g=0;g<s.group_size&&valid;g++)for(auto corner:geometry.corners) {
      std::array<int,4> ids;
      for(int j=0;j<4;j++){int v=h[corner[j]];ids[j]=v<0||v>=geometry.boundary_count?geometry.boundary_count:s.a[v][g];}
      if(s.odd(g))std::swap(ids[1],ids[2]);
      if(!geometry.allowed[geometry.index(ids)]){valid=false;break;}
    }
    if(valid&&leader_feasible(s,h)&&std::find(result.begin(),result.end(),h)==result.end())result.push_back(h);
  }
  }
  return result;
}

std::vector<Hex> Search::extend_cube_group_uncached(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const {
  auto result=extend_cube_group_raw(s,group,q,diagonals);
  if(cover_fixed&&s.group_size==1)result.erase(std::remove_if(result.begin(),result.end(),[&](const Hex& h){return !supports_fixed_mesh(s,h);}),result.end());
  return result;
}
std::vector<Hex> Search::extend_cube_group(const State& s,const std::vector<Hex>& group,Quad q,const Diagonals& diagonals)const {
  auto& cache=cover_memo;
  if(cache.state!=&s||cache.diagonals!=&diagonals||group.empty())return extend_cube_group_uncached(s,group,q,diagonals);
  uint64_t hash=0;for(int v:q)hash=hash*131+uint64_t(v+1);
  for(const auto& h:group)hash=(hash*0x9e3779b97f4a7c15ULL)^hex_hash(h);
  auto& entry=cache.group[hash&1023];
  if(entry.generation==cache.generation&&entry.q==q&&entry.group==group)return entry.result;
  auto result=extend_cube_group_uncached(s,group,q,diagonals);
  entry.generation=cache.generation;entry.q=q;entry.group=group;entry.result=result;return result;
}

struct CoverWork {
  Counts& count;int& nodes;int limit;const State& state;
  ~CoverWork(){audit_add(state,AC_COVER_WORK,nodes);if(nodes>limit)audit_add(state,AC_COVER_UNKNOWN);count.cover_work.fetch_add(nodes,std::memory_order_relaxed);if(nodes>limit)count.cover_unknown.fetch_add(1,std::memory_order_relaxed);}
};
bool Search::orbit_cover_feasible(const State& s,const std::vector<Quad>& faces,const Diagonals& diagonals)const {
  struct FaceOrbit {Quad first,image;int size;};
  std::vector<FaceOrbit> orbits;std::set<Quad> seen;
  for(auto q:faces)if(!seen.count(key(q))) {
    auto image=q;for(int& v:image)v=s.a[v][1];if(s.odd(1))image=reverse(image);
    seen.insert(key(q));seen.insert(key(image));orbits.push_back({q,image,key(q)==key(image)?1:2});
  }
  auto legal=[&](const Hex& h,bool fixed) {
    int mask=0;uint32_t code=0;
    for(int i=0;i<8;i++)if(h[i]>=0) {
      mask|=1<<i;for(int j=0;j<8;j++)if(h[j]>=0&&s.a[h[i]][1]==h[j]){code|=uint32_t(j+1)<<(4*i);break;}
      if(!fixed&&s.mirror_separation&&s.a[h[i]][1]!=h[i]) {
        // Compare every inferred side component, not just the anchored one.
        for(int j=0;j<i;j++)if(h[j]>=0&&s.a[h[j]][1]!=h[j]) {
          auto [a,p]=s.cover_side[h[i]];auto [b,q]=s.cover_side[h[j]];if(a==b&&p!=q)return false;
        }
      }
    }
    const auto& fixed_table=s.half_turn?overlap_patterns.rotation_fixed_masked[mask]:overlap_patterns.fixed_masked[mask];
    const auto& table=fixed?fixed_table:overlap_patterns.exchanged_masked[mask];
    // Without geometric separation, the exchanged-cell table must also
    // admit non-pointwise proper intersections. Use the unrestricted table
    // conservatively; full search still enforces the actual orbit type.
    if(!fixed&&!s.mirror_separation)return overlap_patterns.masked[s.odd(1)][mask].count(code)!=0;
    return table.count(code)!=0;
  };
  auto extend=[&](const std::vector<Hex>& old,int index,bool fixed) {
    std::vector<Hex> result;const auto& orbit=orbits[index];
    if(!fixed&&orbit.size==1)return result;
    if(fixed) {
      auto first=extend_cube_group(s,old,orbit.first,diagonals);
      if(first.empty())return result;
      if(orbit.size==2)first=extend_cube_group(s,first,orbit.image,diagonals);
      for(auto h:first)if(legal(h,true))result.push_back(h);
    } else {
      // The first face chooses the representative of an exchanged pair.
      for(int flip=0;flip<(old.empty()?1:2);flip++) {
        auto possible=extend_cube_group(s,old,flip?orbit.image:orbit.first,diagonals);
        for(auto h:possible)if(legal(h,false)&&std::find(result.begin(),result.end(),h)==result.end())result.push_back(h);
      }
    }
    return result;
  };
  int n=int(orbits.size()),remaining=cap-int(s.cells.size());
  std::array<std::array<std::vector<Hex>,2>,256> initial;
  for(int v=0;v<n;v++) {
    initial[v][0]=extend({},v,true);initial[v][1]=extend({},v,false);
    if(initial[v][0].empty()&&initial[v][1].empty())return false;
  }
  std::array<int,256> assigned;assigned.fill(-1);
  std::array<int,64> cost{},face_count{};
  std::array<std::vector<Hex>,64> witness;
  int nodes=0;CoverWork work{local(),nodes,cover_work_limit,s};
  std::function<bool(int,int,int,std::array<uint64_t,256>)> solve;
  solve=[&](int left,int used,int spent,std::array<uint64_t,256> forbidden) {
    if(!left)return true;
    if(++nodes>cover_work_limit||stopping())return true;
    int best=-1,best_options=1000,best_size=0;
    for(int v=0;v<n;v++)if(assigned[v]<0) {
      int options=used-__builtin_popcountll(forbidden[v]);
      options+=spent<remaining&&!initial[v][0].empty();
      options+=spent+2<=remaining&&!initial[v][1].empty();
      if(options==0)return false;
      if(options<best_options||(options==best_options&&orbits[v].size<best_size)){best=v;best_options=options;best_size=orbits[v].size;}
    }
    auto attempt=[&](int c,int c_cost,std::vector<Hex> possible) {
      if(possible.empty())return false;
      int increment=c_cost==1?orbits[best].size:1;if(face_count[c]+increment>6)return false;
      auto previous=std::move(witness[c]);witness[c]=std::move(possible);
      int old_cost=cost[c];cost[c]=c_cost;face_count[c]+=increment;assigned[best]=c;
      auto next=forbidden;int new_used=std::max(used,c+1);
      for(int v=0;v<n;v++)if(assigned[v]<0) {
        int add=c_cost==1?orbits[v].size:1;
        if(face_count[c]+add>6||extend(witness[c],v,c_cost==1).empty())next[v]|=uint64_t(1)<<c;
      }
      bool ok=solve(left-1,new_used,spent+(c==used?c_cost:0),next);
      assigned[best]=-1;face_count[c]-=increment;cost[c]=old_cost;witness[c]=std::move(previous);return ok;
    };
    for(int c=0;c<used;c++)if(!(forbidden[best]&(uint64_t(1)<<c)))
      if(attempt(c,cost[c],extend(witness[c],best,cost[c]==1)))return true;
    if(spent<remaining&&attempt(used,1,initial[best][0]))return true;
    if(spent+2<=remaining&&attempt(used,2,initial[best][1]))return true;
    return false;
  };
  return solve(n,0,0,{});
}

bool Search::partial_cover_feasible(const State& s,const Hex& h,int assigned,const Budget& budget,const FaceDomains& domains,const Diagonals& diagonals,const std::array<uint8_t,MAXV*MAXV>& placed)const {
  (audit_add(s,AC_PARTIAL_COVER_CALLS,1),local().partial_cover_calls.fetch_add(1,std::memory_order_relaxed));
  State optimistic=s;
  if(leader_bound&&!s.canonical_floor&&!boundary_order.empty()) {
    auto q=face(h,0);int root=boundary_order.front();
    if(key(q)==key(boundary.boundary[root])) {
      Hex partial=h;for(int j=assigned;j<8;j++)partial[j]=-1;
      optimistic.canonical_floor=root_code(orient_partial(partial,boundary.boundary[root]),false);
    }
  }
  optimistic.relation=placed;
  auto next_diagonals=diagonals;
  // Remove the union of faces consumable by either orbit member. This
  // relaxation can omit real successor faces, but cannot invent a demand.
  std::bitset<256> can_consume;
  for(const auto& d:domains)can_consume|=d;
  std::set<Quad> removed;int index=0;
  for(auto [k,q]:s.front)if(can_consume[index++])for(int g=0;g<s.group_size;g++) {
    for(int& v:q)v=s.a[v][g];
    removed.insert(key(q));
  }
  for(auto k:removed)optimistic.front.erase(k);
  for(int a=0;a<assigned;a++)for(int b=0;b<a;b++) {
    int type=relation_type(a,b);
    if(type!=2)for(int g=0;g<s.group_size;g++) {
      int u=s.a[h[a]][g],v=s.a[h[b]][g];
      optimistic.relation[u*MAXV+v]|=type;optimistic.relation[v*MAXV+u]|=type;
    }
  }
  Hex partial=h;for(int j=assigned;j<8;j++)partial[j]=-1;
  auto partner_can_contain=[&](Quad q) {
    if(s.group_size==1)return false;
    // A face potentially internal to the current orbit is not mandatory.
    for(auto positions:FI) {
      bool possible=true;
      for(int j=0;j<8;j++)if(partial[j]>=0) {
        bool on_face=std::find(positions,positions+4,j)!=positions+4;
        bool in_quad=std::find(q.begin(),q.end(),s.a[partial[j]][1])!=q.end();
        if(on_face!=in_quad){possible=false;break;}
      }
      if(possible)return true;
    }
    return false;
  };
  for(int f=0;f<6;f++) {
    bool complete=true;for(int j:FI[f])complete&=j<assigned;
    if(!complete)continue;
    auto original=face(h,f);bool mandatory=!partner_can_contain(original);
    for(int g=0;g<s.group_size;g++) {
      auto q=original;for(int& v:q)v=s.a[v][g];if(s.odd(g))q=reverse(q);
      auto k=key(q);
      if(mandatory&&!s.front.count(k))optimistic.front.emplace(k,reverse(q));
      for(int j=0;j<2;j++) {
        int a=q[j],b=q[j+2];next_diagonals[edgekey(a,b)]=q;
        optimistic.relation[a*MAXV+b]|=2;optimistic.relation[b*MAXV+a]|=2;
      }
    }
  }
  int cost=1;
  if(s.group_size==2) {
    uint32_t code=0;int mask=(1<<assigned)-1;
    for(int i=0;i<assigned;i++)for(int j=0;j<assigned;j++)if(s.a[h[i]][1]==h[j])code|=uint32_t(j+1)<<(4*i);
    const auto& fixed_table=s.half_turn?overlap_patterns.rotation_fixed_masked[mask]:overlap_patterns.fixed_masked[mask];
    if(!fixed_table.count(code))cost=2;
  }
  bool possible=budget.remaining>=cost&&cover_feasible(optimistic,budget.remaining-cost,&next_diagonals);
  if(!possible)(audit_add(s,AC_PARTIAL_COVER_REJECTIONS,1),local().partial_cover_rejections.fetch_add(1,std::memory_order_relaxed));
  return possible;
}

bool Search::cover_feasible(State& s,int override_remaining,const Diagonals* override_diagonals)const {
  s.audit_cover_reason=0;
  int work_limit=override_remaining>=0?partial_work:cover_work_limit;
  auto previous_packing=s.cover_packing;
  s.cover_packing.reset();s.cover_choice_valid=false;
  int remaining=override_remaining>=0?override_remaining:cap-int(s.cells.size()),n=int(s.front.size());if(n<=remaining){s.cover_graph.reset();return true;}
  Parity colors;colors.initialize(s);auto [anchor,base]=colors.root(0);
  for(int v=0;v<int(s.a.size());v++){auto [root,p]=colors.root(v);if(root!=anchor)throw std::runtime_error("Disconnected cover vertex");s.cover_color[v]=base^p;}
  if(s.group_size==2&&s.mirror_separation){Parity sides;if(!side_parity(s,0,sides))throw std::runtime_error("Invalid cover side state");for(int v=0;v<int(s.a.size());v++)s.cover_side[v]=sides.root(v);}
  auto diagonals=override_diagonals?*override_diagonals:known_diagonals(s,boundary);
  std::vector<Quad> faces,keys;for(auto [k,q]:s.front){keys.push_back(k);faces.push_back(q);}

  std::array<std::bitset<256>,256> compatible{};
  std::array<std::array<int,256>,4> images{};

  for(int g=0;g<(cover_symmetry?s.group_size:1);g++)for(int i=0;i<n;i++) {
    auto q=faces[i];for(int& v:q)v=s.a[v][g];images[g][i]=int(std::lower_bound(keys.begin(),keys.end(),key(q))-keys.begin());
  }
  std::array<int,256> previous;previous.fill(-1);
  if(cover_incremental&&s.cover_graph) {
    const auto& old=s.cover_graph->keys;int j=0;
    for(int i=0;i<n;i++) {
      while(j<int(old.size())&&old[j]<keys[i])++j;
      if(j<int(old.size())&&old[j]==keys[i]&&same_orientation(faces[i],s.cover_graph->faces[j]))previous[i]=j;
    }
  }
  // Incompatible old pairs stay incompatible as relations are added.
  // Try extending the parent's packing with newly exposed faces before
  // paying to construct the entire compatibility graph.
  if(quick_packing&&s.group_size<=2&&s.cover_graph) {
    std::vector<int> packing;
    for(int i=0;i<n;i++)if(previous[i]>=0&&previous_packing[previous[i]])packing.push_back(i);
    if(int(packing.size())>remaining){(audit_add(s,AC_QUICK_PACKING_REJECTIONS,1),local().quick_packing_rejections.fetch_add(1,std::memory_order_relaxed));s.audit_cover_reason=1;return false;}
    for(int pass=0;pass<2;pass++)for(int i=0;i<n;i++)if((previous[i]<0)==(pass==0)) {
      if(previous[i]>=0&&previous_packing[previous[i]])continue;
      bool possible=true;
      for(int j:packing) {
        bool share=previous[i]>=0&&previous[j]>=0?s.cover_graph->compatible[previous[i]][previous[j]]:
          faces_can_share(s,faces[i],faces[j],diagonals);
        if(share){possible=false;break;}
      }
      if(possible){packing.push_back(i);if(int(packing.size())>remaining){(audit_add(s,AC_QUICK_PACKING_REJECTIONS,1),local().quick_packing_rejections.fetch_add(1,std::memory_order_relaxed));s.audit_cover_reason=1;return false;}}
    }
    for(int i:packing)s.cover_packing.set(i);
  }
  std::array<std::bitset<MAXV>,MAXV> changed{};
  std::array<std::bitset<MAXV>,256> face_vertices{},dirty_neighbors{};
  std::array<bool,256> dirty_self{};
  if(cover_incremental&&s.cover_graph) {
    for(int a=0;a<int(s.a.size());a++)for(int b=0;b<a;b++)if(s.relation[a*MAXV+b]!=s.cover_graph->relation[a*MAXV+b]){changed[a].set(b);changed[b].set(a);}
    for(int i=0;i<n;i++)if(previous[i]>=0) {
      for(int v:faces[i]){face_vertices[i].set(v);dirty_neighbors[i]|=changed[v];}
      dirty_self[i]=(face_vertices[i]&dirty_neighbors[i]).any();
    }
  }
  uint64_t reused_pairs=0;
  std::array<std::bitset<256>,256> evaluated{};
  for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(!evaluated[i][j]) {
    bool reuse=previous[i]>=0&&previous[j]>=0;
    if(reuse&&s.canonical_floor!=s.cover_graph->canonical_floor&&s.cover_graph->compatible[previous[i]][previous[j]])reuse=false;
    if(reuse)reuse=!dirty_self[i]&&!dirty_self[j]&&!(dirty_neighbors[i]&face_vertices[j]).any();
    reused_pairs+=reuse;
    bool possible=reuse?s.cover_graph->compatible[previous[i]][previous[j]]:faces_can_share(s,faces[i],faces[j],diagonals);
    for(int g=0;g<(cover_symmetry?s.group_size:1);g++) {
      int a=images[g][i],b=images[g][j];evaluated[a].set(b);evaluated[b].set(a);
      if(possible){compatible[a].set(b);compatible[b].set(a);}
    }
  }
  (audit_add(s,AC_PAIR_REUSES,reused_pairs),local().pair_reuses.fetch_add(reused_pairs,std::memory_order_relaxed));
  if(face_order==4&&override_remaining<0) {
    std::array<int,MAXV> degrees{};for(auto c:s.cells)for(int v:c)++degrees[v];
    std::array<int,3> best_score{-100000,-100000,-100000};
    for(int i=0;i<n;i++) {
      int degree=0;for(int v:faces[i])degree+=degrees[v];
      int share=int(compatible[i].count());
      std::array<int,3> score{degree,-share,0};
      if(score>best_score){best_score=score;s.cover_choice=faces[i];s.cover_choice_valid=true;}
    }
  }
  // Every pairwise incompatible face needs its own future cell. Greedy
  // packings are lower bounds regardless of whether they are maximum.
  for(int start=-1;start<0;start++) {
    std::bitset<256> available;for(int i=0;i<n;i++)available.set(i);
    int bound=0,forced=start;std::bitset<256> picked;
    while(available.any()) {
      int best=forced,degree=257;forced=-1;
      if(best<0)for(int i=0;i<n;i++)if(available[i]) {
        int d=int((compatible[i]&available).count());
        if(d<degree){best=i;degree=d;if(d==0)break;}
      }
      picked.set(best);if(++bound>remaining){s.audit_cover_reason=2;return false;}
      available.reset(best);available&=~compatible[best];
    }
    if(picked.count()>s.cover_packing.count())s.cover_packing=picked;
  }
  auto save_graph=[&](){
    if(cover_incremental&&override_remaining<0) {
      auto graph=std::make_shared<CoverGraph>();graph->faces=faces;graph->keys=keys;graph->compatible=compatible;graph->relation=s.relation;graph->canonical_floor=s.canonical_floor;s.cover_graph=std::move(graph);
    }
    return true;
  };
  if(!work_limit)return save_graph();
  std::unique_ptr<CoverMemoScope> memo_scope;
  if(cover_memo_enabled)memo_scope=std::make_unique<CoverMemoScope>(s,diagonals);
  // If S contains no feasible triple of front faces, each future cube
  // contains at most two members of S. Thus |S| > 2R proves infeasibility.
  if(two_face_packing&&cover_cubes&&s.group_size==1&&override_remaining<0&&n<=64&&n>2*remaining) {
    std::vector<uint64_t> triples;uint64_t trials=0;
    std::array<std::bitset<MAXV>,64> fv;for(int i=0;i<n;i++)for(int v:faces[i])fv[i].set(v);
    for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(compatible[i][j]) {
      auto common=compatible[i]&compatible[j];bool any=false;for(int k=0;k<j;k++)any|=common[k];if(!any)continue;
      for(int k=0;k<j;k++)if(common[k]) {
        ++trials;
        int count=(fv[i]|fv[j]|fv[k]).count(),common3=(fv[i]&fv[j]&fv[k]).count();
        int opposite=!(fv[i]&fv[j]).any();opposite+=!(fv[i]&fv[k]).any();opposite+=!(fv[j]&fv[k]).any();
        // Three cube faces either meet at one corner (seven vertices), or
        // contain one opposite pair (eight vertices). Retain every triple
        // with either pattern; extra geometric impossibilities only weaken S.
        if(!((opposite==0&&count==7&&common3==1)||(opposite==1&&count==8&&common3==0)))continue;
        triples.push_back((uint64_t(1)<<i)|(uint64_t(1)<<j)|(uint64_t(1)<<k));
      }
    }
    local().twopack_calls.fetch_add(1,std::memory_order_relaxed);local().twopack_triples.fetch_add(trials,std::memory_order_relaxed);local().twopack_valid.fetch_add(triples.size(),std::memory_order_relaxed);
    int removed=0;
    while(!triples.empty()&&n-removed>2*remaining) {
      std::array<int,64> degree{};
      for(auto mask:triples)while(mask){int v=__builtin_ctzll(mask);mask&=mask-1;++degree[v];}
      int v=int(std::max_element(degree.begin(),degree.end())-degree.begin());++removed;
      triples.erase(std::remove_if(triples.begin(),triples.end(),[&](uint64_t t){return t&(uint64_t(1)<<v);}),triples.end());
    }
    if(triples.empty()&&n-removed>2*remaining) {
      local().twopack_rejections.fetch_add(1,std::memory_order_relaxed);s.audit_cover_reason=5;return false;
    }
  }
  // A completion partitions the front into <=remaining face cliques, each
  // of size <=6. Try coloring the incompatibility graph. Reaching the work
  // limit returns "possibly feasible", never a rejection.
  std::array<int,256> color; color.fill(-1);
  std::array<int,64> sizes{};
  std::array<std::vector<Hex>,64> witnesses;
  int nodes=0,color_limit=remaining;CoverWork work{local(),nodes,work_limit,s};
  int component_size=0;
  auto future_faces_supported=[&](const Hex& h,int cell,int used) {
    if(!cover_closed||override_remaining>=0||s.group_size>2||used!=remaining||component_size!=n)return true;
    for(int f=0;f<6;f++) {
      auto q=face(h,f);if(std::find(q.begin(),q.end(),-1)!=q.end()||s.front.count(key(q)))continue;
      bool supported=false;
      for(int other=0;other<used&&!supported;other++)if(other!=cell) {
        auto possible=extend_cube_group(s,witnesses[other],reverse(q),diagonals);
        for(const auto& k:possible)if(partial_cells_compatible(h,k)){supported=true;break;}
      }
      if(!supported)return false;
    }
    return true;
  };

  std::array<uint64_t,256> forbidden{};
  std::array<int,256> incompatibility_degree{};
  for(int v=0;v<n;v++)incompatibility_degree[v]=n-1-int(compatible[v].count());
  auto colorable=[&](auto&& self,int left,int used)->bool {
    if(!left) {
      for(int c=0;c<used;c++)if(cover_cubes) {
        bool supported=false;
        for(const auto& h:witnesses[c])if(future_faces_supported(h,c,used)){supported=true;break;}
        if(!supported)return false;
      }
      return true;
    }
    if(++nodes>work_limit||stopping())return true;
    int best=-1,sat=-1,degree=-1;
    for(int v=0;v<n;v++)if(color[v]==-1) {
      int count=__builtin_popcountll(forbidden[v]),d=incompatibility_degree[v];
      if(count>sat||(count==sat&&d>degree)){best=v;sat=count;degree=d;}
    }
    for(int c=0;c<std::min(color_limit,used+1);c++)if(!(forbidden[best]&(uint64_t(1)<<c))&&sizes[c]<6) {
      std::vector<Hex> previous;
      if(cover_cubes) {
        auto possible=extend_cube_group(s,witnesses[c],faces[best],diagonals);
        if(cover_global)possible.erase(std::remove_if(possible.begin(),possible.end(),[&](const Hex& h){
          if(!future_faces_supported(h,c,used))return true;
          for(int other=0;other<used;other++)if(other!=c) {
            bool supported=false;for(const auto& k:witnesses[other])if(partial_cells_compatible(h,k)){supported=true;break;}
            if(!supported)return true;
          }
          return false;
        }),possible.end());
        if(possible.empty())continue;
        previous=std::move(witnesses[c]);witnesses[c]=std::move(possible);
      }
      color[best]=c;++sizes[c];
      std::array<uint8_t,256> changed;int changed_count=0;uint64_t bit=uint64_t(1)<<c;
      for(int v=0;v<n;v++)if(color[v]==-1&&!(forbidden[v]&bit)) {
        if(!compatible[best][v]||(cover_cubes&&sizes[c]>=2&&extend_cube_group(s,witnesses[c],faces[v],diagonals).empty())) {
          forbidden[v]|=bit;changed[changed_count++]=uint8_t(v);
        }
      }
      bool ok=self(self,left-1,std::max(used,c+1));
      for(int i=0;i<changed_count;i++)forbidden[changed[i]]&=~bit;
      color[best]=-1;--sizes[c];if(cover_cubes)witnesses[c]=std::move(previous);if(ok)return true;
    }
    return false;
  };
  // Disconnected components of the compatibility graph cannot share a
  // future cell. Prove their minimum cover costs independently, avoiding
  // backtracking over unrelated components.
  std::bitset<256> unseen;for(int i=0;i<n;i++)unseen.set(i);
  int total=int(s.cover_packing.count());
  while(unseen.any()) {
    int first=0;while(!unseen[first])++first;
    std::bitset<256> component,pending;pending.set(first);unseen.reset(first);
    while(pending.any()) {
      int v=0;while(!pending[v])++v;pending.reset(v);component.set(v);
      auto more=compatible[v]&unseen;unseen&=~more;pending|=more;
    }
    component_size=int(component.count());
    int size=int(component.count()),lower=int((component&s.cover_packing).count());
    if(lower==size)continue;
    for(int i=0;i<n;i++)color[i]=component[i]?-1:-2;
    while(lower<size) {
      color_limit=lower;
      int seeded=0;forbidden.fill(0);
      if(cover_seed)for(int v=0;v<n;v++)if(component[v]&&s.cover_packing[v]) {
        int c=seeded++;color[v]=c;sizes[c]=1;
        if(cover_cubes)witnesses[c]=extend_cube_group(s,{},faces[v],diagonals);
        for(int w=0;w<n;w++)if(component[w]&&!compatible[v][w])forbidden[w]|=uint64_t(1)<<c;
      }
      bool possible=colorable(colorable,size-seeded,seeded);
      for(int v=0;v<n;v++)if(component[v])color[v]=-1;
      for(int c=0;c<seeded;c++){sizes[c]=0;witnesses[c].clear();}
      if(possible)break;
      ++lower;if(++total>remaining){s.audit_cover_reason=3;return false;}
    }
  }
  if(mirror_cover&&s.group_size==2&&!orbit_cover_feasible(s,faces,diagonals)){s.audit_cover_reason=4;return false;}
  return save_graph();
}

bool Search::replay(const Mesh& seed) {
  State s=initial();
  if(seed.points.size()<boundary.points.size())throw std::runtime_error("Replay has too few boundary vertices");
  std::vector<Action> seed_actions;
  if(replay_actions.empty())seed_actions=actions(seed.points,axial,s.group_size,half_turn);
  else {
    std::ifstream f(replay_actions);if(!f)throw std::runtime_error("Cannot open replay actions");
    seed_actions.resize(seed.points.size(),Action{-1,-1,-1,-1});
    for(auto& a:seed_actions)for(int g=0;g<s.group_size;g++)
      if(!(f>>a[g])||a[g]<0||a[g]>=int(seed.points.size()))throw std::runtime_error("Invalid replay action vertex");
    std::string extra;if(f>>extra)throw std::runtime_error("Extra replay action data");
    for(int i=0;i<int(seed.points.size());i++) {
      if(seed_actions[i][0]!=i)throw std::runtime_error("Replay identity action mismatch");
      for(int g=0;g<s.group_size;g++)for(int k=0;k<s.group_size;k++)
        if(seed_actions[seed_actions[i][g]][k]!=seed_actions[i][g^k])throw std::runtime_error("Invalid replay group action");
      if(s.group_size==4&&seed_actions[i][3]==i&&(seed_actions[i][1]!=i||seed_actions[i][2]!=i))
        throw std::runtime_error("Replay vertex has a geometrically impossible stabilizer");
    }
    std::set<std::array<Quad,6>> cells;
    for(auto h:seed.hexes)cells.insert(hexkey(h));
    for(auto h:seed.hexes)for(int g=0;g<s.group_size;g++)if(!cells.count(hexkey(reflect(h,seed_actions,g,s.half_turn))))
      throw std::runtime_error("Replay action does not preserve oriented cells");
  }
  for(size_t v=0;v<boundary.points.size();v++) {
    if(seed_actions[v]!=s.a[v])throw std::runtime_error("Replay boundary action mismatch");
    for(int k=0;k<3;k++)if(std::abs(seed.points[v][k]-boundary.points[v][k])>1e-12)
      throw std::runtime_error("Replay boundary coordinate mismatch");
  }
  std::vector<int> to_search(seed.points.size(),-1),to_seed;
  for(size_t i=0;i<boundary.points.size();i++){to_search[i]=int(i);to_seed.push_back(int(i));}
  std::set<Hex> used;
  while(!s.front.empty()) {
    if(face_cover&&!cover_feasible(s))return false;
    Quad q=select_face(s);Hex h{q[0],q[3],q[2],q[1],-1,-1,-1,-1};
    Hex source{};bool found=false;
    // Recover an oriented cube ordering with the required bottom face and
    // uniquely determined opposite endpoints along its four transverse edges.
    for(auto c:seed.hexes)if(!used.count(sorted(c)))for(int f=0;f<6&&!found;f++) {
      Quad b=face(c,f),wanted;for(int j=0;j<4;j++)wanted[j]=to_seed[q[j]];
      if(!same_orientation(b,wanted))continue;
      for(int j=0;j<4;j++)source[j]=to_seed[h[j]];
      for(int j=0;j<4;j++) {
        int top=-1;
        for(auto& e:EI){int a=c[e[0]],v=c[e[1]];if(v==source[j])std::swap(a,v);
          if(a==source[j]&&std::find(source.begin(),source.begin()+4,v)==source.begin()+4)top=v;}
        if(top<0)throw std::runtime_error("Replay could not orient seed cell");
        source[j+4]=top;
      }
      found=true;
    }
    if(!found)throw std::runtime_error("Replay seed does not fill the selected face");
    auto placed=s.relation;auto diagonals=known_diagonals(s,boundary);auto rules=stabilizer_rules(s,h);
    Overlaps overlaps{};
    if(early_overlap)for(int j=0;j<4;j++)if(!extend_overlaps(s,h,j,overlaps))throw std::runtime_error("Replay bottom overlap rejected");
    Budget budget(s,cap,face_domains);auto domains=budget.initial;
    if(early_budget)for(int j=0;j<4;j++)budget.extend(h,j,domains);
    Parity parity;if(early_parity)parity.initialize(s);
    VertexDomains vertices(s);
    for(int j=4;j<8;j++) {
      auto options=overlap_choices(s,h,j,overlaps,vertices);
      int option=to_search[source[j]];
      int v=source[j];
      if(to_search[v]<0) {
        bool fa=seed_actions[v][1]==v,fb=seed_actions[v][2]==v;
        int type=s.group_size==1?0:s.group_size==2?(fa?1:0):(fa?(fb?3:1):(fb?2:0));
        option=int(s.a.size())+type;
        int base=allocate(s,type);if(int(s.a.size())>vertex_cap)return false;
        to_seed.resize(s.a.size(),-1);
        for(int g=0;g<s.group_size;g++){int old=seed_actions[v][g],now=s.a[base][g];
          if(to_search[old]>=0&&to_search[old]!=now)throw std::runtime_error("Replay orbit-label conflict");
          to_search[old]=now;to_seed[now]=old;}
      }
      h[j]=to_search[v];Trail trail;
      if(overlap_domains&&early_overlap&&!(options.word[option/64]&(uint64_t(1)<<(option%64))))throw std::runtime_error("Replay overlap domain rejected");
      if(corner_geometry&&!geometry.feasible(s,h,j+1))throw std::runtime_error("Replay corner geometry rejected");
      if(early_overlap&&!extend_overlaps(s,h,j,overlaps))throw std::runtime_error("Replay overlap rejected");
      if(early_budget){budget.extend(h,j,domains);if(!budget.feasible(overlaps,domains))return false;}
      if(early_parity&&!parity.extend(s,h,j))throw std::runtime_error("Replay parity rejected");
      if(edge_completion&&component_bound&&early_completion&&!edge_completion_feasible(s,h,j+1,placed,budget,overlaps,vertices,completion_domains))return false;
      if(partial_cover&&face_cover&&s.group_size<=2&&j+1==partial_cover&&
        !partial_cover_feasible(s,h,j+1,budget,domains,diagonals,placed))throw std::runtime_error("Replay partial cover rejected");
      if(incremental_faces&&!consistent_new_faces(h,j,diagonals,placed))throw std::runtime_error("Replay incremental faces rejected");
      if(!consistent_rules(s,h,j+1,rules)||!impose(s,h,j,placed,trail)||!consistent_faces(s,h,j+1,diagonals))throw std::runtime_error("Replay partial orbit rejected");
    }
    if(component_bound&&early_completion&&!completion_feasible(s,h,placed,budget))return false;
    size_t before=s.cells.size();if(!insert(s,h))return false;
    for(size_t i=before;i<s.cells.size();i++) {
      Hex old;for(int j=0;j<8;j++)old[j]=to_seed[s.cells[i][j]];
      if(std::none_of(seed.hexes.begin(),seed.hexes.end(),[&](Hex c){return hexkey(c)==hexkey(old);}))throw std::runtime_error("Replay generated a cell absent from seed");
      used.insert(sorted(old));
    }
  }
  if(used.size()!=seed.hexes.size())throw std::runtime_error("Replay did not recover all seed cells");
  if(!topology(s))throw std::runtime_error("Replay seed failed necessary topology checks");
  emit(s);
  std::ofstream meshout(output+".replay.mesh");
  if(!meshout)throw std::runtime_error("Cannot write replay mesh");
  meshout<<std::setprecision(17)<<"MeshVersionFormatted 2\nDimension 3\nVertices\n"<<to_seed.size()<<'\n';
  for(int v:to_seed){auto p=seed.points[v];meshout<<p[0]<<' '<<p[1]<<' '<<p[2]<<" 0\n";}
  meshout<<"Hexahedra\n"<<s.cells.size()<<'\n';
  for(auto c:s.cells){for(int v:c)meshout<<v+1<<' ';meshout<<"0\n";}
  meshout<<"Quadrilaterals\n"<<boundary.boundary.size()<<'\n';
  for(auto q:boundary.boundary){for(int v:q)meshout<<v+1<<' ';meshout<<"0\n";}
  meshout<<"End\n";if(!meshout)throw std::runtime_error("Replay mesh write failed");
  std::lock_guard<std::mutex> l(io);
  std::cerr<<"REPLAY_OK cells="<<s.cells.size()<<" cell_orbits="<<s.orbit_count<<" vertices="<<s.a.size()<<'\n';
  return true;
}

int main(int argc,char** argv) {
  Search search;bool cover_work_given=false;std::string input="input/pyramid.mesh",replay;int interior=-1;bool face_order_given=false,planes_none=false,mirrors_given=false;
  try {
    for(int i=1;i<argc;i++) {
      std::string k=argv[i];if(k=="--help"){std::cout<<"--input FILE --cap N --threads N --interior N --planes none|diagonal|axial --mirrors 0|1|2 --half-turn 0|1 --seconds N --progress N --output FILE --replay SEED --early-overlap 0|1 --early-budget 0|1 --early-parity 0|1 --state-trace FILE --tree-trace DIR --replay-actions FILE --face-order 0|1|2|3|4|5 --root-face N --partial-cover 0|5|6|7 --partial-work N --quick-packing 0|1 --boundary-canonical 0|1 --topology-symmetry 0|1 --leader-bound 0|1 --fast-pairs 0|1 --cover-closed 0|1 --cover-fixed 0|1 --cover-memo 0|1 --two-face-packing 0|1 --component-bound 0|1 --early-completion 0|1 --edge-completion 0|1 --incremental-faces 0|1 --overlap-domains 0|1 --component-index 0|1 --completion-domains 0|1 --task-depth N --corner-geometry 0|1 --vertex-domains 0|1 --face-domains 0|1 --symmetry-parity 0|1 --color-domains 0|1 --geometry-domains 0|1 --location-domains 0|1 --mirror-separation 0|1 --face-cover 0|1 --cover-work N --cover-cubes 0|1 --cover-cache 0|1 --mirror-cover 0|1 --cover-global 0|1 --cover-symmetry 0|1 --cover-incremental 0|1 --cover-seed 0|1 --cover-flat 0|1 --cover-cache-bits N\n";return 0;}
      if(i+1==argc)throw std::runtime_error("Missing option value");
      std::string v=argv[++i];
      if(k=="--input")input=v;else if(k=="--cap")search.cap=std::stoi(v);
      else if(k=="--threads")search.threads=std::stoi(v);else if(k=="--interior")interior=std::stoi(v);
      else if(k=="--seconds")search.seconds=std::stod(v);else if(k=="--progress")search.progress=std::stoi(v);
      else if(k=="--task-depth"){search.task_depth=std::stoi(v);if(search.task_depth<0||search.task_depth>48)throw std::runtime_error("Invalid task depth");}
      else if(k=="--early-completion"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.early_completion=v=="1";}
      else if(k=="--edge-completion"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.edge_completion=v=="1";}
      else if(k=="--face-cover"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.face_cover=std::stoi(v);}
      else if(k=="--cover-work"){cover_work_given=true;search.cover_work_limit=std::stoi(v);if(search.cover_work_limit<0||search.cover_work_limit>1000000)throw std::runtime_error("Cover work must be 0..1000000");}
      else if(k=="--cover-cache-bits"){search.cover_cache_bits=std::stoi(v);if(search.cover_cache_bits<10||search.cover_cache_bits>22)throw std::runtime_error("Cover cache bits must be 10..22");}
      else if(k=="--cover-flat"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_flat=std::stoi(v);}
      else if(k=="--cover-seed"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_seed=std::stoi(v);}
      else if(k=="--cover-incremental"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_incremental=std::stoi(v);}
      else if(k=="--cover-symmetry"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_symmetry=std::stoi(v);}
      else if(k=="--cover-global"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_global=std::stoi(v);}
      else if(k=="--cover-cubes"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_cubes=std::stoi(v);}
      else if(k=="--cover-cache"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_cache=std::stoi(v);}
      else if(k=="--mirror-cover"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.mirror_cover=std::stoi(v);}
      else if(k=="--completion-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.completion_domains=v=="1";}
      else if(k=="--component-index"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.component_index=v=="1";}
      else if(k=="--overlap-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.overlap_domains=v=="1";}
      else if(k=="--incremental-faces"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.incremental_faces=v=="1";}
      else if(k=="--component-bound"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.component_bound=v=="1";}
      else if(k=="--cover-closed"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_closed=std::stoi(v);}
      else if(k=="--fast-pairs"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.fast_pairs=std::stoi(v);}
      else if(k=="--two-face-packing"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.two_face_packing=std::stoi(v);}
      else if(k=="--cover-fixed"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_fixed=std::stoi(v);}
      else if(k=="--root-face")search.root_face=std::stoi(v);
      else if(k=="--leader-bound"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.leader_bound=std::stoi(v);}
      else if(k=="--topology-symmetry"){search.topology_symmetry=std::stoi(v);if(search.topology_symmetry<0||search.topology_symmetry>1)throw std::runtime_error("Topology symmetry must be 0 or 1");}
      else if(k=="--boundary-canonical"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.boundary_canonical=std::stoi(v);}
      else if(k=="--partial-work"){search.partial_work=std::stoi(v);if(search.partial_work<0||search.partial_work>1000000)throw std::runtime_error("Partial work must be 0..1000000");}
      else if(k=="--cover-memo"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.cover_memo_enabled=std::stoi(v);}
      else if(k=="--quick-packing"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.quick_packing=std::stoi(v);}
      else if(k=="--partial-cover"){search.partial_cover=std::stoi(v);if(search.partial_cover!=0&&(search.partial_cover<5||search.partial_cover>7))throw std::runtime_error("Partial cover must be 0, 5, 6, or 7");}
      else if(k=="--face-order"){face_order_given=true;search.face_order=std::stoi(v);if(search.face_order<0||search.face_order>5)throw std::runtime_error("Invalid face order");}
      else if(k=="--early-parity"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.early_parity=v=="1";}
      else if(k=="--early-budget"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.early_budget=v=="1";}
      else if(k=="--early-overlap"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.early_overlap=v=="1";}
      else if(k=="--replay-actions")search.replay_actions=v;
      else if(k=="--tree-trace")search.tree_trace=v;
      else if(k=="--state-trace")search.state_trace=v;
      else if(k=="--output")search.output=v;else if(k=="--replay")replay=v;
      else if(k=="--mirror-separation"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.mirror_separation=v=="1";}
      else if(k=="--location-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.location_domains=v=="1";}
      else if(k=="--geometry-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.geometry_domains=v=="1";}
      else if(k=="--color-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.color_domains=v=="1";}
      else if(k=="--symmetry-parity"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.symmetry_parity=v=="1";}
      else if(k=="--face-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.face_domains=v=="1";}
      else if(k=="--vertex-domains"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.vertex_domains=v=="1";}
      else if(k=="--corner-geometry"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.corner_geometry=v=="1";}
      else if(k=="--mirrors"){if(v!="0"&&v!="1"&&v!="2")throw std::runtime_error("Expected zero, one, or two mirrors");search.mirrors=std::stoi(v);mirrors_given=true;}
      else if(k=="--half-turn"){if(v!="0"&&v!="1")throw std::runtime_error("Expected 0 or 1");search.half_turn=v=="1";}
      else if(k=="--planes"){if(v!="none"&&v!="diagonal"&&v!="axial")throw std::runtime_error("Unknown mirror planes");search.axial=v=="axial";planes_none=v=="none";}
      else throw std::runtime_error("Unknown option "+k);
    }
    if(search.half_turn) {
      if(mirrors_given&&search.mirrors!=0)throw std::runtime_error("Half-turn mode cannot be combined with mirrors");
      search.mirrors=0;search.axial=false;search.mirror_separation=false;planes_none=true;
    }
    if(planes_none&&search.mirrors!=0)throw std::runtime_error("--planes none requires --mirrors 0");
    if(search.mirrors==0&&!search.half_turn){search.axial=false;search.mirror_separation=false;search.symmetry_parity=false;search.overlap_domains=false;search.cover_symmetry=false;search.mirror_cover=false;}
    if(search.group_size()==1&&!cover_work_given)search.cover_work_limit=100000;
    if(!face_order_given){if(search.mirrors==0)search.face_order=5;else if(search.axial&&search.mirrors==1)search.face_order=4;}
    if(search.cap<1||search.cap>48||search.threads<1||search.threads>256||search.progress<0||search.seconds<0||interior< -1)throw std::runtime_error("Invalid search bounds");
    if(replay.empty()&&!search.replay_actions.empty())throw std::runtime_error("Replay actions require --replay");
    if(!search.early_overlap){search.early_budget=false;search.early_parity=false;}
    if(!search.tree_trace.empty()&&search.group_size()!=1)throw std::runtime_error("Tree trace currently requires zero mirrors");
    search.boundary=readmesh(input);
    if(search.root_face< -1||search.root_face>=int(search.boundary.boundary.size()))throw std::runtime_error("Root face must be -1 or a zero-based boundary-face index");
    std::set<int> bv;for(auto q:search.boundary.boundary)for(int v:q)bv.insert(v);
    if(bv.size()!=search.boundary.points.size()||*bv.begin()!=0||*bv.rbegin()!=int(bv.size())-1)throw std::runtime_error("Input must contain only contiguous boundary vertices");
    int full=2*search.cap+1-int(search.boundary.boundary.size())/2;
    if(full<0)throw std::runtime_error("Cell cap is below the incidence/degree bound");
    if(interior<0)interior=full;
    if(interior>full)throw std::runtime_error("Interior allowance exceeds proved bound");
    search.vertex_cap=int(bv.size())+interior;
    if(search.vertex_cap>=MAXV)throw std::runtime_error("Vertex capacity exceeded");
    std::signal(SIGTERM,on_signal);std::signal(SIGINT,on_signal);
    omp_set_dynamic(0);omp_set_num_threads(search.threads);
    if(search.corner_geometry||search.location_domains)search.geometry.initialize(search.boundary.points,search.axial,search.group_size(),search.geometry_domains,search.location_domains,search.half_turn);
    if(search.topology_symmetry&&(search.mirrors||search.half_turn||search.corner_geometry||search.location_domains||
        search.geometry_domains||!search.boundary_canonical))
      throw std::runtime_error("Topological boundary automorphisms require zero mirrors, disabled geometry, and boundary canonicalization");
    search.prepare_boundary_symmetry();
    search.begin();
    std::cerr<<"SEARCH cap="<<search.cap<<" interior_cap="<<interior<<" full_interior_bound="<<full<<" planes="<<(search.mirrors==0?"none":search.axial?"axial":"diagonal")<<" mirrors="<<search.mirrors<<" half_turn="<<search.half_turn<<" enforcement="<<(search.group_size()==1?"individual_cells":"cell_orbits")<<" early_overlap="<<search.early_overlap<<" early_budget="<<search.early_budget<<" early_parity="<<search.early_parity<<" two_face_packing="<<search.two_face_packing<<" cover_memo="<<search.cover_memo_enabled<<" cover_fixed="<<search.cover_fixed<<" cover_closed="<<search.cover_closed<<" fast_pairs="<<search.fast_pairs<<" leader_bound="<<search.leader_bound<<" boundary_canonical="<<search.boundary_canonical<<" topology_symmetry="<<search.topology_symmetry<<" boundary_group="<<search.boundary_symmetries.size()<<" native_arch="<<HEX_NATIVE_ARCH<<" partial_work="<<search.partial_work<<" quick_packing="<<search.quick_packing<<" partial_cover="<<search.partial_cover<<" root_face="<<search.root_face<<" face_order="<<search.face_order<<" component_bound="<<search.component_bound<<" early_completion="<<search.early_completion<<" face_cover="<<search.face_cover<<" cover_work_limit="<<search.cover_work_limit<<" cover_cubes="<<search.cover_cubes<<" cover_cache="<<search.cover_cache<<" cover_global="<<search.cover_global<<" cover_symmetry="<<search.cover_symmetry<<" cover_incremental="<<search.cover_incremental<<" cover_seed="<<search.cover_seed<<" cover_flat="<<search.cover_flat<<" cover_cache_bits="<<search.cover_cache_bits<<" mirror_cover="<<search.mirror_cover<<" completion_domains="<<search.completion_domains<<" component_index="<<search.component_index<<" overlap_domains="<<search.overlap_domains<<" incremental_faces="<<search.incremental_faces<<" edge_completion="<<search.edge_completion<<" mirror_separation="<<search.mirror_separation<<" location_domains="<<search.location_domains<<" geometry_domains="<<search.geometry_domains<<" color_domains="<<search.color_domains<<" symmetry_parity="<<search.symmetry_parity<<" face_domains="<<search.face_domains<<" vertex_domains="<<search.vertex_domains<<" corner_geometry="<<search.corner_geometry<<" task_depth="<<search.task_depth<<'\n';
    if(!replay.empty()) {
      bool ok=search.replay(readmesh(replay));search.end();if(!ok){std::cerr<<"REPLAY_REJECTED\n";return 3;}
    } else {
      State root=search.initial();
      #pragma omp parallel
      {
        #pragma omp single
        {
          std::cerr<<"SEARCH_WORKERS count="<<omp_get_num_threads()<<'\n';
          try{search.dfs(std::move(root),0);}catch(const std::exception& e){search.error(e);}
        }
      }
      search.end();
    }
    if(search.failed)return 2;
    if(interrupted){std::cerr<<"SEARCH_INCOMPLETE\n";return 124;}
    std::cerr<<"SEARCH_FINISHED candidates="<<search.totals()[3]<<'\n';return 0;
  } catch(const std::exception& e) {
    if(search.monitor.joinable()){search.failed=true;search.end();}
    std::cerr<<"ERROR "<<e.what()<<'\n';return 2;
  }
}
