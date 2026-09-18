// Fixed-boundary, two-reflection adapter for the authors' HexOpt sJGrad kernel.
// Compile with upstream meshQuality.cpp and geometry.cpp; see HEXOPT_NOTES.md.
#include <cmath>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <tuple>
#include "meshQuality.h"
using V=std::array<double,3>;
V reflect(V v,int g,int axial) {
 if(axial){if(g&1)v[0]=-v[0];if(g&2)v[1]=-v[1];}
 else {if(g&1)std::swap(v[0],v[1]);if(g&2){std::swap(v[0],v[1]);v[0]=-v[0];v[1]=-v[1];}}
 return v;
}
int main(int argc,char**argv) {
 if(argc!=3&&argc!=4)return 2;
 std::ifstream in(argv[1]);int n,ne,fixed,axial;in>>n>>ne>>fixed>>axial;
 if(!in||n<8||ne<1||fixed<0||fixed>n)return 2;
 std::vector<V>x(n);std::vector<std::array<int,4>>actions(n);std::vector<std::array<int,8>>hex(ne);
 for(auto& v:x)for(auto& c:v)in>>c;
 for(auto& a:actions)for(auto& i:a){in>>i;if(i<0||i>=n)return 2;}
 for(auto& h:hex)for(auto& i:h){in>>i;if(i<0||i>=n)return 2;}
 if(!in)return 2;
 auto limit=[](const char*name,int fallback) {
  const char* value=std::getenv(name);if(!value)return fallback;
  char* end=nullptr;long parsed=std::strtol(value,&end,10);
  return end==value||*end||parsed<1||parsed>10000000?-1:int(parsed);
 };
 const int max_iterations=limit("HEXOPT_MAX_ITERATIONS",300000);
 const int stall_limit=limit("HEXOPT_STALL_LIMIT",5000);
 if(max_iterations<1||stall_limit<1)return 2;
 const auto original=x;auto best=x;
 // Optional orthogonal affine columns permit matched side vertices to move.
 // Each row: initial_value lower_bound upper_bound number_of_terms, followed
 // by (vertex coordinate coefficient) terms. Coordinates are zero based.
 struct Column { double value,lower,upper,norm=0;std::vector<std::tuple<int,int,double>> terms; };
 std::vector<Column> columns;auto offset=original;
 const bool linear=argc==4;
 if(linear) {
  std::ifstream constraint(argv[3]);int count;constraint>>count;
  if(!constraint||count<1||count>3*n)return 2;
  columns.resize(count);std::vector<int> owner(3*n,-1);
  for(int j=0;j<count;j++) {
   auto& c=columns[j];int terms;constraint>>c.value>>c.lower>>c.upper>>terms;
   if(!constraint||!std::isfinite(c.value)||!std::isfinite(c.lower)||!std::isfinite(c.upper)||
      c.lower>c.value||c.value>c.upper||terms<1||terms>3*n)return 2;
   for(int t=0;t<terms;t++) {
    int v,k;double a;constraint>>v>>k>>a;
    if(!constraint||v<0||v>=n||k<0||k>=3||!std::isfinite(a)||a==0||owner[3*v+k]>=0)return 2;
    owner[3*v+k]=j;c.terms.emplace_back(v,k,a);c.norm+=a*a;offset[v][k]-=a*c.value;
   }
  }
  std::string extra;if(constraint>>extra)return 2;
 }
 auto affine_project=[&](std::vector<V>&y,bool coordinates) {
  const auto old=y;y=coordinates?offset:std::vector<V>(n,V{0,0,0});
  for(const auto& c:columns) {
   double value=0;for(auto [v,k,a]:c.terms)value+=a*(old[v][k]-(coordinates?offset[v][k]:0));
   value/=c.norm;if(coordinates)value=std::clamp(value,c.lower,c.upper);
   for(auto [v,k,a]:c.terms)y[v][k]+=a*value;
  }
 };
 auto project=[&](std::vector<V>&y){auto old=y;for(int i=0;i<n;i++){
  y[i]={0,0,0};for(int g=0;g<4;g++){auto v=reflect(old[actions[i][g]],g,axial);for(int k=0;k<3;k++)y[i][k]+=v[k]/4;}
 }};
 omp_set_dynamic(0);omp_set_num_threads(1);
 int last=0,it=0;double threshold=.01,achieved=0;
 auto start=std::chrono::steady_clock::now();
 for(;it<max_iterations&&it-last<stall_limit&&threshold<.91;it++) {
  std::vector<V> grad(n,V{0,0,0});bool all=true;
  for(auto h:hex){std::array<V,8>v,d{};for(int j=0;j<8;j++)v[j]=x[h[j]];bool good=true;
   sJGrad(v[0][0],v[0][1],v[0][2],v[1][0],v[1][1],v[1][2],v[2][0],v[2][1],v[2][2],v[3][0],v[3][1],v[3][2],v[4][0],v[4][1],v[4][2],v[5][0],v[5][1],v[5][2],v[6][0],v[6][1],v[6][2],v[7][0],v[7][1],v[7][2],
   d[0][0],d[0][1],d[0][2],d[1][0],d[1][1],d[1][2],d[2][0],d[2][1],d[2][2],d[3][0],d[3][1],d[3][2],d[4][0],d[4][1],d[4][2],d[5][0],d[5][1],d[5][2],d[6][0],d[6][1],d[6][2],d[7][0],d[7][1],d[7][2],-threshold,good);
   all&=good;
   for(int j=0;j<8;j++)for(int k=0;k<3;k++)grad[h[j]][k]+=d[j][k];
  }
  if(all){best=x;achieved=threshold;threshold+=.01;last=it;std::cout<<"THRESHOLD "<<achieved<<" iteration="<<it<<'\n';continue;}
  if(linear)affine_project(grad,false);else project(grad);
  for(int i=linear?0:fixed;i<n;i++){
   double norm=0;for(double v:grad[i])norm+=v*v;
   if(!std::isfinite(norm)){std::cerr<<"NONFINITE_GRADIENT\n";return 3;}
   double step=1e-3/std::max(1.,std::sqrt(norm));
   for(int k=0;k<3;k++)x[i][k]-=step*grad[i][k];
  }
  if(linear)affine_project(x,true);
  else {project(x);for(int i=0;i<fixed;i++)x[i]=original[i];}
 }
 std::ofstream out(argv[2]);out<<std::setprecision(17);for(auto v:best)out<<v[0]<<' '<<v[1]<<' '<<v[2]<<'\n';
 std::cout<<"FINISHED achieved="<<achieved<<" iterations="<<it<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
 return out?0:2;
}
