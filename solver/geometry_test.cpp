// Independent exact 4x4 determinant checks and valid geometric fixture checks.
#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main
#include <random>
using Rational=boost::multiprecision::cpp_rational;
Rational rational(double x) {
  if(x==0)return 0;
  int e;double m=std::frexp(x,&e);Rational r=int64_t(std::ldexp(m,53));
  if(e>=53)r*=CornerGeometry::Integer(1)<<(e-53);
  else r/=CornerGeometry::Integer(1)<<(53-e);
  return r;
}
int reference(const std::vector<Vec>& p,std::array<int,4> q) {
  std::array<std::array<Rational,4>,4> a;
  for(int i=0;i<4;i++){a[i][0]=1;for(int j=0;j<3;j++)a[i][j+1]=rational(p[q[i]][j]);}
  Rational det=0;std::array<int,4> perm{0,1,2,3};
  do {
    Rational term=1;int inversions=0;
    for(int i=0;i<4;i++){term*=a[i][perm[i]];for(int j=i+1;j<4;j++)inversions+=perm[i]>perm[j];}
    if(inversions&1)det-=term;else det+=term;
  }while(std::next_permutation(perm.begin(),perm.end()));
  return det>0?1:det<0?-1:0;
}
int main(int argc,char** argv) {
  if(argc!=2)throw std::runtime_error("Pass input directory");
  std::string path=argv[1];std::mt19937 random(129417);int signs=0,prefixes=0,masks=0,locations=0;
  for(auto name:{"pyramid.mesh","cube.mesh"}) {
    auto mesh=readmesh(path+"/"+name);int n=int(mesh.points.size());
    for(bool axial:{false,true})for(int mirrors:{1,2}) {
      CornerGeometry geometry;geometry.initialize(mesh.points,axial,1<<mirrors,true,true);
      for(int sample=0;sample<256;sample++) {
        std::array<int,4> q;
        do{for(int& v:q)v=int(random()%(n+1));}while(std::set<int>(q.begin(),q.end()).size()!=4);
        bool positive=false;auto it=std::find(q.begin(),q.end(),n);
        if(it==q.end())positive=reference(mesh.points,q)>0;
        else {int slot=int(it-q.begin());for(int v=0;v<n;v++){auto trial=q;trial[slot]=v;positive|=reference(mesh.points,trial)>0;}}
        if(bool(geometry.allowed[geometry.index(q)])!=positive)throw std::runtime_error("Exact determinant mismatch");
        ++signs;
      }
      for(int sample=0;sample<512&&!geometry.plane_values.empty();sample++) {
        int a=random()%geometry.plane_values.size(),b=random()%geometry.plane_values.size();
        const auto& x=geometry.plane_values[a];const auto& y=geometry.plane_values[b];
        // Independent dual test: search for a nonnegative combination of the
        // inequalities that is nonpositive throughout the boundary convex hull.
        auto separates=[&](CornerGeometry::Integer alpha,CornerGeometry::Integer beta) {
          if(alpha==0&&beta==0)return false;
          for(int v=0;v<n;v++)if(alpha*x[v]+beta*y[v]>0)return false;
          return true;
        };
        bool impossible=separates(1,0)||separates(0,1);
        for(int v=0;v<n&&!impossible;v++)if(x[v]*y[v]<=0) {
          CornerGeometry::Integer alpha=y[v]<0?-y[v]:y[v],beta=x[v]<0?-x[v]:x[v];
          impossible=separates(alpha,beta);
        }
        if(geometry.locations_compatible(a,b)==impossible)throw std::runtime_error("Halfspace primal/dual mismatch");
        ++locations;
      }
      State state;state.group_size=1<<mirrors;state.a=actions(mesh.points,axial,state.group_size);
      // Add one generic interior orbit and compare mask filtering with scalar
      // corner checks for every candidate, including new-vertex placeholders.
      allocate(state,0);
      for(int sample=0;sample<512;sample++) {
        Hex h;std::set<int> used;
        for(int& v:h){do{v=int(random()%state.a.size());}while(used.count(v));used.insert(v);}
        int position=4+int(random()%4);auto mask=geometry.choices(h,position);
        for(int v=0;v<int(state.a.size());v++) {
          h[position]=v;
          bool scalar=geometry.feasible(state,h,position+1,position);
          bool bit=(mask.word[v/64]>>(v%64))&1;
          if(scalar!=bit)throw std::runtime_error("Geometry mask mismatch");
          ++masks;
        }
      }
    }
  }
  for(auto entry:{std::pair<const char*,bool>{"published36-symmetric.mesh",false},{"published44-symmetric.mesh",true},{"grid3-seed.mesh",false},{"grid3-seed.mesh",true},{"one-axial-seed.mesh",true},{"one-diagonal-seed.mesh",false}}) {
    auto mesh=readmesh(path+"/"+entry.first);std::set<int> boundary;
    for(auto q:mesh.boundary)for(int v:q)boundary.insert(v);
    auto points=mesh.points;points.resize(boundary.size());
    CornerGeometry geometry;geometry.initialize(points,entry.second,2,true,true);
    State state;state.group_size=2;state.a=actions(mesh.points,entry.second,2);
    if(!geometry.impose_locations(state,mesh.hexes))throw std::runtime_error("Valid fixture locations rejected");
    for(auto h:mesh.hexes)for(int n=0;n<=8;n++) {
      if(!geometry.feasible(state,h,n))throw std::runtime_error("Valid fixture corner rejected");
      ++prefixes;
    }
  }
  {
    State paired;paired.group_size=2;paired.a.resize(16);
    for(auto& side:paired.boundary_side)side.fill(-1);
    for(int v=0;v<16;v++)paired.a[v]={v,v^8,-1,-1};
    paired.boundary_side[0][0]=1;paired.boundary_side[0][8]=0;
    paired.cells={Hex{0,1,2,3,4,5,6,7},Hex{8,9,10,11,12,13,14,15}};
    Parity sides;if(!side_parity(paired,0,sides))throw std::runtime_error("Separated cells rejected");
    paired.cells.push_back(Hex{0,1,2,3,4,5,6,15});
    Parity contradiction;if(side_parity(paired,0,contradiction))throw std::runtime_error("Crossing nonfixed cell accepted");
    auto cube=readmesh(path+"/cube.mesh");State fixed;fixed.group_size=2;fixed.a=actions(cube.points,true,2);
    for(auto& side:fixed.boundary_side)side.fill(-1);
    for(int v=0;v<8;v++)fixed.boundary_side[0][v]=cube.points[v][0]>0;
    fixed.cells={Hex{0,1,2,3,4,5,6,7}};
    Parity crossing;if(!side_parity(fixed,0,crossing))throw std::runtime_error("Fixed mirror-crossing cell rejected");
  }
  // Exact zero, subnormal input, and very different coordinate exponents.
  std::vector<Vec> probes{{0,0,0},{1,0,0},{0,1,0},{0.25,0.25,0},{0,0,std::ldexp(1.,-1074)},{0,0,std::ldexp(1.,900)}};
  CornerGeometry geometry;geometry.initialize(probes,false,2,false);
  for(auto q:{std::array<int,4>{0,1,2,3},{0,1,2,4},{0,1,2,5},{0,2,1,4}}) {
    if(bool(geometry.allowed[geometry.index(q)])!=(reference(probes,q)>0))throw std::runtime_error("Extreme determinant mismatch");
    ++signs;
  }
  std::cout<<"GEOMETRY_OK signs="<<signs<<" mask_comparisons="<<masks<<" halfspace_pairs="<<locations<<" valid_prefixes="<<prefixes<<'\n';
}
