// Safe orbit reduction under explicitly verified parent-net stabilizer subgroups.
// A subgroup need not be the full stabilizer: representatives still form a
// COMPLETE COVER, but their count is NOT an isomorphism-class count.
#define main extension_probe_main
#include "extend_parents.cpp"
#undef main

vector<int> exterior_columns(const vector<int>& g,int m) {
 vector<int> result; int n=m*(m-1)/2;
 for(int a=0;a<m;++a) for(int b=a+1;b<m;++b) {
  int image=0,bit=0;
  for(int i=0;i<m;++i) for(int j=i+1;j<m;++j,++bit)
   if((((g[i]>>a)&1)*((g[j]>>b)&1))^(((g[i]>>b)&1)*((g[j]>>a)&1))) image|=1<<bit;
  result.push_back(image);
 }
 return result;
}
int linear_image(int x,const vector<int>& columns) {
 int y=0;
 while(x) {int p=__builtin_ctz((unsigned)x); y^=columns[p];x&=x-1;}
 return y;
}
int residue(int x,const vector<int>& B) {
 for(int a:B) if(x&(1<<__builtin_ctz((unsigned)a))) x^=a;
 return x;
}
vector<int> compose(const vector<int>& a,const vector<int>& b) {
 vector<int>c;
 for(int row:a) c.push_back(linear_image(row,b));
 return c;
}
vector<vector<int>> parent_actions(const vector<int>& B,int m,const vector<int>& freecols,
                                  const vector<vector<int>>& supplied={}) {
 int n=m*(m-1)/2,d=freecols.size();
 vector<int>I;for(int i=0;i<m;++i) I.push_back(1<<i);
 vector<vector<int>> elementary;
 for(int i=0;i<m;++i) for(int j=0;j<m;++j) if(i!=j) {
  auto g=I;g[i]^=1<<j;elementary.push_back(g);
 }
 for(int i=0;i<m;++i) for(int j=i+1;j<m;++j) {
  auto g=I;swap(g[i],g[j]);elementary.push_back(g);
 }
 set<vector<int>> actions;
 auto accept=[&](const vector<int>& g) {
  if((int)g.size()!=m || (int)rref(g,m).size()!=m) return false;
  for(int row:g)if(row<0 || row>=(1<<m))return false;
  auto cols=exterior_columns(g,m);
  vector<int> images;
  for(int a:B) images.push_back(linear_image(a,cols));
  if(rref(images,n)!=B) return false;
  vector<int> act;
  for(int c:freecols) {
   int y=residue(cols[c],B),z=0;
   for(int j=0;j<d;++j) if(y&(1<<freecols[j])) z|=1<<j;
   act.push_back(z);
  }
  actions.insert(act);return true;
 };
 if(!supplied.empty()) {
  for(auto g:supplied)if(!accept(g)){cerr<<"Supplied generator does not preserve parent\n";exit(4);}
 } else {
  vector<bool> already;
  for(auto g:elementary) already.push_back(accept(g));
  for(int i=0;i<(int)elementary.size();++i) for(int j=0;j<(int)elementary.size();++j)
   if(!(already[i]&&already[j])) accept(compose(elementary[i],elementary[j]));
 }
 vector<int>qI;for(int j=0;j<d;++j) qI.push_back(1<<j);actions.erase(qI);
 return vector<vector<int>>(actions.begin(),actions.end());
}

int main(int argc,char**argv) {
 int m=argc>1?stoi(argv[1]):6,n=m*(m-1)/2;
 if(m!=6 && m!=7) return 2;
 auto start=chrono::steady_clock::now();
 string prefix=argc>2?"fullparent":"subgroup";
 map<vector<int>,vector<vector<int>>> supplied;
 if(argc>2) {
  ifstream auts(argv[2]);int mm,r,ng,x;
  while(auts>>mm>>r>>ng) {
   vector<int>B;
   for(int j=0;j<r;++j){auts>>x;B.push_back(x);}
   vector<vector<int>> matrices;
   for(int k=0;k<ng;++k){vector<int>g;for(int j=0;j<mm;++j){auts>>x;g.push_back(x);}matrices.push_back(g);}
   if(mm==m)supplied[B]=matrices;
  }
 }
 map<vector<int>,vector<int>> parents;
 string line;
 while(getline(cin,line)) {
  stringstream ss(line);int id,mm,r,x;
  if(!(ss>>id>>mm>>r)||mm!=m||r!=9-m) continue;
  vector<int>B;while(ss>>x)B.push_back(x);
  parents[rref(B,n)].push_back(id);
 }
 Geometry F(m);
 map<vector<int>,int>good;
 uint64_t candidate_mass=0,total_orbits=0;int pi=0;
 ofstream masses(prefix+"_m"+to_string(m)+"_orbit_masses.txt");
 for(const auto& entry:parents) {
  auto B=entry.first;auto S=span(B);
  Radical common=~Radical(0);for(int a:B)common&=F.rad[a];
  vector<int> freecols;
  for(int p=0;p<n;++p) {
   bool isp=false;for(int a:B) if(__builtin_ctz((unsigned)a)==p) isp=true;
   if(!isp) freecols.push_back(p);
  }
  if(argc>2 && !supplied.count(B)){cerr<<"Missing parent automorphisms\n";return 4;}
  auto acts=parent_actions(B,m,freecols,supplied[B]);
  int N=1<<freecols.size();vector<bool>seen(N,false);
  uint64_t count=0,mass=0;
  auto expand=[&](int k) {
   int x=0;for(int j=0;j<(int)freecols.size();++j)if(k&(1<<j))x|=1<<freecols[j];return x;
  };
  for(int k=1;k<N;++k) {
   if(seen[k])continue;
   int x=expand(k);
   if((common&F.rad[x])!=Radical(1)||!F.mf_extend(S,x))continue;
   auto T=B;T.push_back(x);T=rref(T,n);
   int ad=F.affine_dim(T);if(ad<0)continue;
   vector<int>queue{k};seen[k]=true;
   for(size_t head=0;head<queue.size();++head)for(const auto&action:acts){
    int y=linear_image(queue[head],action);
    if(!seen[y]){seen[y]=true;queue.push_back(y);}
   }
   ++count;mass+=queue.size();good[T]=ad;
   masses<<pi+1<<" "<<k<<" "<<queue.size()<<"\n";
  }
  cout<<"parent="<<++pi<<" id="<<entry.second[0]<<" verified_subgroup_generators="<<acts.size()
      <<" SR_extension_orbits="<<count<<" SR_extension_mass="<<mass<<"\n"<<flush;
  candidate_mass+=mass;total_orbits+=count;
 }
 ofstream out(prefix+"_m"+to_string(m)+"_net_representatives.txt");
 for(const auto&entry:good){out<<m<<" "<<entry.first.size()<<" "<<entry.second;for(int b:entry.first)out<<" "<<b;out<<"\n";}
 cout<<"RESULT m="<<m<<" subgroup_orbit_representatives="<<total_orbits<<" coordinate_distinct_representatives="<<good.size()
     <<" SR_extension_mass="<<candidate_mass<<" elapsed_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
}
