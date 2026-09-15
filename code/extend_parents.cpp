// Enumerate scalar-polar one-dimensional extensions of exported 512 parents.
// Compile: g++ -O3 -std=c++17 extend_parents.cpp -o extend_parents
// Run: ./extend_parents 6 [or 7] < parents512_forms.txt
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
using namespace std;
using Radical=__uint128_t;

int weight(Radical x) {
 return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64));
}

vector<int> rref(vector<int> B,int n) {
 int row=0;
 for(int p=0;p<n;++p) {
  int k=row;
  while(k<(int)B.size() && !(B[k]&(1<<p))) ++k;
  if(k==(int)B.size()) continue;
  swap(B[row],B[k]);
  for(int i=0;i<(int)B.size();++i)
   if(i!=row && (B[i]&(1<<p))) B[i]^=B[row];
  ++row;
 }
 B.resize(row); return B;
}

vector<int> span(const vector<int>& B) {
 vector<int>S{0};
 for(int a:B) {
  int len=S.size();
  for(int j=0;j<len;++j) S.push_back(a^S[j]);
 }
 return S;
}

struct Geometry {
 int m,n,N;
 vector<unsigned char> rank;
 vector<Radical> rad;
 vector<int> quadratic_monomials;
 Geometry(int mm):m(mm),n(m*(m-1)/2),N(1<<n),rank(N),rad(N),
  quadratic_monomials(1<<m) {
   for(int v=0;v<(1<<m);++v) {
    int bit=0;
    for(int i=0;i<m;++i) for(int j=i+1;j<m;++j,++bit)
     if((v&(1<<i))&&(v&(1<<j))) quadratic_monomials[v]|=1<<bit;
   }
   for(int f=0;f<N;++f) {
    unsigned rows[7]={}; int bit=0;
    for(int i=0;i<m;++i) for(int j=i+1;j<m;++j,++bit)
     if(f&(1<<bit)) {rows[i]|=1u<<j; rows[j]|=1u<<i;}
    unsigned temp[7]; copy(rows,rows+m,temp);
    int rr=0;
    for(int p=0;p<m;++p) {
     int k=rr;
     while(k<m && !(temp[k]&(1u<<p))) ++k;
     if(k==m) continue;
     swap(temp[rr],temp[k]);
     for(int i=rr+1;i<m;++i) if(temp[i]&(1u<<p)) temp[i]^=temp[rr];
     ++rr;
    }
    rank[f]=rr;
    for(int v=0;v<(1<<m);++v) {
     bool yes=true;
     for(int i=0;i<m;++i) if(__builtin_parity(rows[i]&v)) {yes=false;break;}
     if(yes) rad[f]|=Radical(1)<<v;
    }
   }
 }
 bool mf_pair(int a,int b) const {
   int count=weight(rad[a]&rad[b]);
   int dim=__builtin_ctz((unsigned)count);
   return rank[a]+rank[b]+rank[a^b]==2*(m-dim);
 }
 bool mf_extend(const vector<int>& S,int x) const {
   for(int a:S) if(a) for(int s:S) if(!mf_pair(a,x^s)) return false;
   return true;
 }
 int affine_dim(const vector<int>& B) const {
   auto S=span(B); int r=B.size(),nv=m*r,eqrank=0;
   vector<uint64_t> echelon(nv,0);
   for(int lambda=1;lambda<(int)S.size();++lambda) {
    int f=S[lambda];
    for(int v=1;v<(1<<m);++v) if(rad[f]&(Radical(1)<<v)) {
     uint64_t eq=0;
     for(int i=0;i<r;++i) if(lambda&(1<<i)) eq|=uint64_t(v)<<(m*i);
     if(__builtin_parity(f&quadratic_monomials[v])) eq|=1ull<<nv;
     for(int p=0;p<nv;++p) if(eq&(1ull<<p)) {
       if(echelon[p]) eq^=echelon[p];
       else {echelon[p]=eq;++eqrank;eq=0;break;}
     }
     if(eq) return -1;
    }
   }
   return nv-eqrank;
 }
};

int main(int argc,char**argv) {
 int m=argc>1?stoi(argv[1]):6,n=m*(m-1)/2;
 if(m!=6 && m!=7) {cerr<<"Only m=6 or 7 is supported.\n";return 2;}
 auto start=chrono::steady_clock::now();
 map<vector<int>,vector<int>> parents;
 string line; int inputcount=0;
 while(getline(cin,line)) {
  stringstream ss(line);int id,mm,r,x;
  if(!(ss>>id>>mm>>r)) continue;
  if(mm!=m || r!=9-m) continue;
  vector<int>B;
  while(ss>>x) B.push_back(x);
  if((int)B.size()!=r) {cerr<<"Malformed row\n";return 2;}
  parents[rref(B,n)].push_back(id); ++inputcount;
 }
 cout<<"m="<<m<<" parent_groups="<<inputcount<<" distinct_coordinate_parent_nets="<<parents.size()<<"\n";
 Geometry F(m);
 cout<<"Geometry initialized seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
 map<vector<int>,int> good;
 set<vector<int>> mf_nets,stem_nets;
 uint64_t total=0,mf=0,stem=0,real=0;
 int pi=0;
 for(const auto& entry:parents) {
  const auto &B=entry.first; auto S=span(B);
  for(int a:S) for(int b:S) if(!F.mf_pair(a,b)) {cerr<<"Parent MF failed\n";return 3;}
  if(F.affine_dim(B)<0) {cerr<<"Parent reality failed\n";return 3;}
  Radical parentrad=~Radical(0);
  for(int a:B) parentrad&=F.rad[a];
  vector<int> freecols;
  for(int p=0;p<n;++p) {
    bool pivot=false;
    for(int a:B) if(__builtin_ctz((unsigned)a)==p) pivot=true;
    if(!pivot) freecols.push_back(p);
  }
  uint64_t pt=0,pm=0,ps=0,pr=0;
  for(int k=1;k<(1<<(int)freecols.size());++k) {
   int x=0;
   for(int j=0;j<(int)freecols.size();++j) if(k&(1<<j)) x|=1<<freecols[j];
   ++pt;
   if(!F.mf_extend(S,x)) continue;
   ++pm;
   auto T=B;T.push_back(x);T=rref(T,n);
   mf_nets.insert(T);
   if((parentrad&F.rad[x])!=Radical(1)) continue;
   ++ps;stem_nets.insert(T);
   int ad=F.affine_dim(T);
   if(ad<0) continue;
   ++pr;good[T]=ad;
  }
  total+=pt;mf+=pm;stem+=ps;real+=pr;
  cout<<"parent_net="<<++pi<<" representative_id="<<entry.second[0]<<" input_groups="<<entry.second.size()
      <<" extensions="<<pt<<" MF="<<pm<<" stem="<<ps<<" real="<<pr<<"\n";
 }
 ofstream out("candidate_m"+to_string(m)+"_nets.txt");
 map<int,int> dims;
 for(const auto& entry:good) {
   out<<m<<" "<<entry.first.size()<<" "<<entry.second;
   for(int x:entry.first) out<<" "<<x;
   out<<"\n"; ++dims[entry.second];
 }
 cout<<"RESULT m="<<m<<" quotient_extensions="<<total<<" MF_hits="<<mf<<" stem_hits="<<stem<<" real_hits="<<real
     <<" distinct_coordinate_MF_nets="<<mf_nets.size()<<" distinct_coordinate_stem_nets="<<stem_nets.size()
     <<" distinct_coordinate_SR_nets="<<good.size()<<" affine_dimensions=";
 for(auto d:dims) cout<<d.first<<":"<<d.second<<",";
 cout<<"\nelapsed_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
 return 0;
}
