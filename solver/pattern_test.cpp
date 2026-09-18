// Exhaustive cross-check of the early overlap lookup against the original
// complete-cell compatibility test, for every partial involution on 8 slots.
#define main mirror_search_program_main
#include "mirror_search.cpp"
#undef main

int main() {
  std::array<int,8> mapping;mapping.fill(-2);
  uint64_t tested=0,legal=0,geometric=0;
  std::function<void()> visit=[&] {
    int i=0;while(i<8&&mapping[i]!=-2)++i;
    if(i<8) {
      mapping[i]=-1;visit();mapping[i]=i;visit();
      for(int j=i+1;j<8;j++)if(mapping[j]==-2) {
        mapping[i]=j;mapping[j]=i;visit();mapping[j]=-2;
      }
      mapping[i]=-2;return;
    }
    Hex original{0,1,2,3,4,5,6,7};int shared=0;
    for(int v:mapping)shared+=v>=0;
    for(int odd=0;odd<2;odd++) {
      Hex other;for(int j=0;j<8;j++)other[j]=mapping[j]<0?8+j:mapping[j];
      if(odd){std::swap(other[1],other[3]);std::swap(other[5],other[7]);}
      bool expected=shared==8?hexkey(original)==hexkey(other):compatible(original,other);
      uint32_t code=0;for(int j=0;j<8;j++)if(mapping[j]>=0)code|=uint32_t(mapping[j]+1)<<(4*j);
      bool actual=overlap_patterns.prefix[odd][8].count(code);
      if(actual!=expected)throw std::runtime_error("Overlap table mismatch");
      if(odd) {
        bool pointwise=true;for(int j=0;j<8;j++)if(mapping[j]>=0)pointwise&=mapping[j]==j;
        bool valid=expected&&(shared==8||pointwise);
        if(bool(overlap_patterns.reflection_prefix[8].count(code))!=valid)throw std::runtime_error("Reflection intersection mismatch");
        if(bool(overlap_patterns.fixed_prefix[8].count(code))!=(expected&&shared==8))throw std::runtime_error("Fixed-cell table mismatch");
        if(valid) {
          ++geometric;
          for(int n=0;n<=8;n++) {
            uint32_t partial=0;for(int j=0;j<n;j++)if(mapping[j]>=0&&mapping[j]<n)partial|=uint32_t(mapping[j]+1)<<(4*j);
            if(!overlap_patterns.reflection_prefix[n].count(partial))throw std::runtime_error("Geometric prefix lost");
            if(shared==8&&!overlap_patterns.fixed_prefix[n].count(partial))throw std::runtime_error("Fixed-cell prefix lost");
          }
        }
      }
      ++tested;if(!expected)continue;++legal;
      for(int n=0;n<=8;n++) {
        code=0;for(int j=0;j<n;j++)if(mapping[j]>=0&&mapping[j]<n)code|=uint32_t(mapping[j]+1)<<(4*j);
        auto it=overlap_patterns.prefix[odd][n].find(code);
        if(it==overlap_patterns.prefix[odd][n].end()||(shared==8&&!(it->second&2))||(shared==4&&!(it->second&4)))
          throw std::runtime_error("Legal intersection lost in a prefix");
      }
    }
  };
  visit();std::cout<<"OVERLAP_TABLES_OK mappings="<<tested<<" legal="<<legal<<" geometric_reflections="<<geometric<<'\n';
}
