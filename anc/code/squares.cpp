// (1) Does the set of square cycles generate the cycle space (over Q) of the root-poset Hasse diagram?
// (2) For every Coxeter word c, are opposite edges of every square parallel?
#include <bits/stdc++.h>
using namespace std;
int n; vector<vector<int>> A; typedef vector<int> V;
V refl(int i, V v){int c=0;for(int j=0;j<n;j++)c+=A[i][j]*v[j];v[i]-=c;return v;}
V applyw(const vector<int>&w,int upto,V v){for(int k=upto-1;k>=0;k--)v=refl(w[k],v);return v;}
bool pos(const V&v){for(int x:v)if(x<0)return false;return true;}
// Exact rank modulo the prime P=2^31-1 of an integer matrix. rank_P <= rank_Q <= (cycle rank) for the square cycles,
// so rank_P == cycle rank certifies exactly that the squares span H_1(Hasse diagram; Q).
int rankQ(vector<vector<double>> Md){ const long long P=2147483647LL; int rows=Md.size(); if(!rows) return 0; int cols=Md[0].size();
  vector<vector<long long>> M(rows,vector<long long>(cols)); for(int i=0;i<rows;i++)for(int k=0;k<cols;k++){ long long v=llround(Md[i][k]); M[i][k]=((v%P)+P)%P; }
  auto pw=[&](long long b,long long e){long long r=1;b%=P;while(e){if(e&1)r=r*b%P;b=b*b%P;e>>=1;}return r;};
  int r=0; for(int c=0;c<cols&&r<rows;c++){int p=-1;for(int i=r;i<rows;i++)if(M[i][c]){p=i;break;} if(p<0)continue; swap(M[p],M[r]);
    long long inv=pw(M[r][c],P-2); for(int k=c;k<cols;k++) M[r][k]=M[r][k]*inv%P;
    for(int i=0;i<rows;i++) if(i!=r&&M[i][c]){ long long f=M[i][c]; for(int k=c;k<cols;k++) M[i][k]=((M[i][k]-f*M[r][k])%P+P)%P; } r++; }
  return r;}
int main(int argc,char**argv){
  char t=argv[1][0]; n=atoi(argv[2]); A.assign(n,vector<int>(n,0));
  for(int i=0;i<n;i++)A[i][i]=2; for(int i=0;i+1<n;i++){A[i][i+1]=-1;A[i+1][i]=-1;}
  if(t=='B'){A[n-2][n-1]=-1;A[n-1][n-2]=-2;} if(t=='C'){A[n-2][n-1]=-2;A[n-1][n-2]=-1;}
  if(t=='D'){A[n-2][n-1]=0;A[n-1][n-2]=0;A[n-3][n-1]=-1;A[n-1][n-3]=-1;}
  if(t=='E'){A[n-2][n-1]=0;A[n-1][n-2]=0;A[2][n-1]=-1;A[n-1][2]=-1;}
  if(t=='F'){A[1][2]=-2;A[2][1]=-1;}
  if(t=='G'){A[0][1]=-1;A[1][0]=-3;}
  set<V> R; vector<V> fr; for(int i=0;i<n;i++){V e(n,0);e[i]=1;R.insert(e);fr.push_back(e);}
  while(!fr.empty()){vector<V> nf;for(auto&r:fr)for(int i=0;i<n;i++){V s=refl(i,r);bool nz=false;for(int x:s)if(x)nz=true;if(nz&&pos(s)&&!R.count(s)){R.insert(s);nf.push_back(s);}}fr=nf;}
  vector<V> roots(R.begin(),R.end());int N=roots.size();map<V,int>id;for(int i=0;i<N;i++)id[roots[i]]=i;
  vector<array<int,3>> cov; map<pair<int,int>,int> eid;
  for(int a=0;a<N;a++)for(int s=0;s<n;s++){V b=roots[a];b[s]++;if(id.count(b)){eid[{a,id[b]}]=cov.size();cov.push_back({a,id[b],s});}}
  int E=cov.size();
  // squares
  vector<array<int,4>> sq; // a,b(=a+s),d(=a+t),e
  vector<vector<double>> M;
  for(int a=0;a<N;a++)for(int s=0;s<n;s++)for(int u=s+1;u<n;u++){V b=roots[a];b[s]++;V d=roots[a];d[u]++;V e=b;e[u]++;
    if(id.count(b)&&id.count(d)&&id.count(e)){int B=id[b],D=id[d],Ee=id[e];sq.push_back({a,B,D,Ee});
      vector<double> row(E,0); row[eid[{a,B}]]+=1; row[eid[{B,Ee}]]+=1; row[eid[{D,Ee}]]-=1; row[eid[{a,D}]]-=1; M.push_back(row);}}
  int adjsq=0; for(int a=0;a<N;a++)for(int s1=0;s1<n;s1++)for(int u=s1+1;u<n;u++){V b=roots[a];b[s1]++;V d=roots[a];d[u]++;V e=b;e[u]++; if(id.count(b)&&id.count(d)&&id.count(e)&&A[s1][u]!=0) adjsq++;}
  int cyc=E-N+1, r=rankQ(M); printf("   squares with adjacent s,t: %d\n",adjsq);
  printf("%c%d: roots %d edges %d cycle rank %d squares %zu span-rank %d -> %s\n",t,n,N,E,cyc,sq.size(),r,(r==cyc?"SQUARES GENERATE":"NOT generated"));
  // parallelism for all c
  vector<int> perm(n);iota(perm.begin(),perm.end(),0);long long bad=0,cnt=0;
  do{ vector<int> w;int k=0; while((int)w.size()<N){int s=perm[k%n];k++;V e(n,0);e[s]=1;if(pos(applyw(w,w.size(),e)))w.push_back(s);}
    vector<int>P(N); for(int m=0;m<N;m++){V e(n,0);e[w[m]]=1;P[id[applyw(w,m,e)]]=m;}
    for(auto&q:sq){ bool ab=P[q[0]]<P[q[1]], de=P[q[2]]<P[q[3]], ad=P[q[0]]<P[q[2]], be=P[q[1]]<P[q[3]]; if(ab!=de||ad!=be) bad++; }
    cnt++; }while(next_permutation(perm.begin(),perm.end()));
  printf("   Coxeter words %lld, non-parallel squares total %lld\n",cnt,bad);
}
