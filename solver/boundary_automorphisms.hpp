// Exact automorphisms of a spherical quadrangulation. SPDX-License-Identifier: MIT
#pragma once

struct SurfaceAutomorphism {std::vector<int> permutation;bool odd;};
std::vector<SurfaceAutomorphism> surface_automorphisms(const std::vector<Quad>& faces,int vertices) {
  const int n=int(faces.size());
  if(n<1||n>256||vertices!=n+2||vertices>MAXV)
    throw std::runtime_error("Automorphisms require a spherical quad surface within vertex capacity");
  std::array<int16_t,MAXV*MAXV> half;half.fill(-1);
  for(int i=0;i<n;i++) {
    auto q=faces[i];auto sorted=q;std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Repeated quad vertex");
    for(int j=0;j<4;j++) {
      int a=q[j],b=q[(j+1)%4];
      if(a<0||a>=vertices||b<0||b>=vertices||half[a*MAXV+b]>=0)
        throw std::runtime_error("Invalid oriented surface edge");
      half[a*MAXV+b]=int16_t(i);
    }
  }
  for(auto q:faces)for(int j=0;j<4;j++)if(half[q[(j+1)%4]*MAXV+q[j]]<0)
    throw std::runtime_error("Surface is not closed and consistently oriented");
  std::vector<SurfaceAutomorphism> result;
  for(bool odd:{false,true})for(int target=0;target<n;target++)for(int rotation=0;rotation<4;rotation++) {
    std::vector<int> map(vertices,-1),inverse(vertices,-1),face_map(n,-1),face_inverse(n,-1);
    bool valid=true;
    auto set_vertex=[&](int a,int b){
      if((map[a]>=0&&map[a]!=b)||(inverse[b]>=0&&inverse[b]!=a)){valid=false;return;}
      map[a]=b;inverse[b]=a;
    };
    for(int j=0;j<4;j++)set_vertex(faces[0][j],faces[target][(rotation+(odd?-j:j)+4)%4]);
    std::vector<std::pair<int,int>> queue{{0,target}};face_map[0]=target;face_inverse[target]=0;
    for(size_t head=0;head<queue.size()&&valid;head++) {
      int source=queue[head].first;auto q=faces[source];
      for(int j=0;j<4&&valid;j++) {
        int a=q[j],b=q[(j+1)%4],ma=map[a],mb=map[b];
        if(ma<0||mb<0){valid=false;break;}
        int next_source=half[b*MAXV+a],next_target=half[odd?ma*MAXV+mb:mb*MAXV+ma];
        if(next_target<0){valid=false;break;}
        if(face_map[next_source]>=0){if(face_map[next_source]!=next_target)valid=false;continue;}
        if(face_inverse[next_target]>=0){valid=false;break;}
        auto u=faces[next_source],v=faces[next_target];
        int si=0,ti=0;while(u[si]!=a)++si;while(v[ti]!=ma)++ti;
        for(int k=0;k<4;k++)set_vertex(u[(si+k)%4],v[(ti+(odd?-k:k)+4)%4]);
        face_map[next_source]=next_target;face_inverse[next_target]=next_source;
        queue.push_back({next_source,next_target});
      }
    }
    if(valid&&queue.size()==size_t(n)&&std::find(map.begin(),map.end(),-1)==map.end())
      result.push_back({std::move(map),odd});
  }
  // A successful identity traversal also establishes dual connectivity and
  // that every boundary vertex is used. With closed oriented edge incidence
  // and Euler characteristic two, this is a spherical surface.
  if(result.empty())throw std::runtime_error("Surface is disconnected or has unused vertices");
  return result;
}
