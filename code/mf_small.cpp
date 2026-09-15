// Exact enumeration of alternating-form subspaces over F_2.
// Compile: g++ -O3 -std=c++17 mf_small.cpp -o mf_small
// Run: ./mf_small
// No GAP or external libraries needed. Integers encode alternating matrices.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <vector>
using namespace std;

int binary_rank(vector<unsigned> rows) {
  int r=0;
  for (int bit=0; bit<16; ++bit) {
    int j=r;
    while (j<(int)rows.size() && !(rows[j]&(1u<<bit))) ++j;
    if(j==(int)rows.size()) continue;
    swap(rows[r],rows[j]);
    for(int k=r+1;k<(int)rows.size();++k)
      if(rows[k]&(1u<<bit)) rows[k]^=rows[r];
    ++r;
  }
  return r;
}

struct FormSpace {
 int m,n,N;
 vector<int> ranks;
 vector<uint32_t> radicals;
 vector<vector<uint8_t>> compatible;
 FormSpace(int m_):m(m_),n(m*(m-1)/2),N(1<<n),ranks(N),radicals(N),
    compatible(N,vector<uint8_t>(N)) {
   for(int f=0;f<N;++f) {
     vector<unsigned> rows(m,0); int bit=0;
     for(int i=0;i<m;++i) for(int j=i+1;j<m;++j,++bit)
       if(f&(1<<bit)){ rows[i]^=1u<<j; rows[j]^=1u<<i; }
     ranks[f]=binary_rank(rows);
     for(int v=0;v<(1<<m);++v) {
       bool rad=true;
       for(auto row:rows) if(__builtin_parity(row&v)) {rad=false;break;}
       if(rad) radicals[f]|=1u<<v;
     }
   }
   for(int a=0;a<N;++a) for(int b=0;b<N;++b) {
     int count=__builtin_popcount(radicals[a]&radicals[b]);
     int dim=__builtin_ctz((unsigned)count);
     compatible[a][b]=(ranks[a]+ranks[b]+ranks[a^b]==2*(m-dim));
   }
 }
 bool extend(const vector<int>& S,int x,vector<int>& T) {
   // Old-old triples already satisfy MF. Every new triple contains an old
   // vector a and two members x+s,x+s+a of the new coset, so these tests suffice.
   for(int a:S) if(a) for(int s:S)
     if(!compatible[a][x^s]) return false;
   T=S;
   for(int s:S) T.push_back(x^s);
   return true;
 }
 int ambivalent_dimension(const vector<int>& S) {
   int r=__builtin_ctz((unsigned)S.size()),nv=m*r;
   vector<uint64_t> echelon(nv,0);
   int equation_rank=0;
   for(int lambda=1;lambda<(int)S.size();++lambda) {
     int form=S[lambda];
     for(int v=1;v<(1<<m);++v) if(radicals[form]&(1u<<v)) {
       uint64_t equation=0;
       for(int i=0;i<r;++i) if(lambda&(1<<i))
         equation|=uint64_t(v)<<(m*i);
       int rhs=0,bit=0;
       for(int i=0;i<m;++i) for(int j=i+1;j<m;++j,++bit)
         if((form&(1<<bit)) && (v&(1<<i)) && (v&(1<<j))) rhs^=1;
       if(rhs) equation|=1ull<<nv;
       for(int p=0;p<nv;++p) if(equation&(1ull<<p)) {
         if(echelon[p]) equation^=echelon[p];
         else {echelon[p]=equation;++equation_rank;equation=0;break;}
       }
       if(equation) return -1;
     }
   }
   return nv-equation_rank;
 }
 bool star_plus_symplectic(const vector<int>& S) {
   if(m!=5 || S.size()!=32) return false;
   vector<int> ranktwo;
   for(int f:S) if(ranks[f]==2) ranktwo.push_back(f);
   if(ranktwo.size()!=15) return false;
   int common=0;
   for(int a=1;a<32;++a) {
     bool good=true;
     for(int f:ranktwo) {
       int coeff[5][5]={},bit=0;
       for(int i=0;i<5;++i) for(int j=i+1;j<5;++j,++bit)
         coeff[i][j]=(f>>bit)&1;
       for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) for(int k=j+1;k<5;++k)
         if((((a>>i)&1)*coeff[j][k])^(((a>>j)&1)*coeff[i][k])^(((a>>k)&1)*coeff[i][j]))
           good=false;
     }
     if(good) ++common;
   }
   return common==1;
 }
};

void enumerate(int m,int target,bool contain_e12) {
 FormSpace F(m);
 const int offset=contain_e12?1:0, d=F.n-offset, k=target-offset;
 uint64_t visited=0, surviving=0, represented_total=0, pruned=0, ambivalent=0;
 uint64_t star_extensions=0;
 vector<uint64_t> nodes(k+1);
 vector<int> pivots;
 // Each pivot pattern specifies one RREF cell. Free entries are only in
 // nonpivot columns to the right of each row's leading 1.
 function<void(int)> choose_pivots = [&](int lo) {
   if((int)pivots.size()==k) {
     vector<vector<int>> row_options(k);
     uint64_t cell_total=1;
     for(int row=0;row<k;++row) {
       vector<int> freecols;
       for(int c=pivots[row]+1;c<d;++c)
         if(find(pivots.begin(),pivots.end(),c)==pivots.end()) freecols.push_back(c);
       for(int mask=0;mask<(1<<(int)freecols.size());++mask) {
         int x=1<<(pivots[row]+offset);
         for(int j=0;j<(int)freecols.size();++j)
           if(mask&(1<<j)) x^=1<<(freecols[j]+offset);
         row_options[row].push_back(x);
       }
       cell_total*=row_options[row].size();
     }
     represented_total+=cell_total;
     vector<int> initial=contain_e12?vector<int>{0,1}:vector<int>{0};
     function<void(int,const vector<int>&)> dfs = [&](int row,const vector<int>& S) {
       ++nodes[row];
       if(row==k) {
         ++visited; ++surviving;
         int ad=F.ambivalent_dimension(S);
         if(ad>=0) ++ambivalent;
         if(F.star_plus_symplectic(S)) ++star_extensions;
         if(surviving<=10 || ad>=0) {
           cout<<"MF survivor m="<<m<<" r="<<target<<" ambivalent_dim="<<ad<<" vectors=";
           for(int x:S) cout<<x<<",";
           cout<<"\n";
         }
         return;
       }
       for(int x:row_options[row]) {
         vector<int> T;
         if(F.extend(S,x,T)) dfs(row+1,T);
         else {
           uint64_t omitted=1;
           for(int j=row+1;j<k;++j) omitted*=row_options[j].size();
           pruned+=omitted;
         }
       }
     };
     dfs(0,initial);
     return;
   }
   for(int p=lo;p<=d-(k-(int)pivots.size());++p) {
     pivots.push_back(p); choose_pivots(p+1); pivots.pop_back();
   }
 };
 choose_pivots(0);
 cout<<"RESULT m="<<m<<" r="<<target<<" contain_e12="<<contain_e12
     <<" total_RREF_subspaces="<<represented_total<<" pruned="<<pruned
     <<" surviving_MF="<<surviving<<" surviving_MF_ambivalent="<<ambivalent
     <<" star_plus_symplectic="<<star_extensions
     <<" partition_check="<<(pruned+visited==represented_total)
     <<" nodes_by_depth=";
 for(auto x:nodes) cout<<x<<",";
 cout<<"\n";
}

int main() {
 auto start=chrono::steady_clock::now();
 enumerate(3,3,false);
 enumerate(4,4,false);
 enumerate(5,4,true);
 enumerate(5,5,true);
 cout<<"elapsed_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
}
