#include <cstring>
// For crystallographic types (integer Cartan): for every Coxeter element c, check that the w0(c) inversion order
// induces an orientation of the root-poset Hasse diagram with zero circulation (existence of a potential).
#include <bits/stdc++.h>
using namespace std;
int n; vector<vector<int>> A; // A[i][j] = <alpha_i^vee, alpha_j>
typedef vector<int> V;
V refl(int i, V v){ int c=0; for(int j=0;j<n;j++) c+=A[i][j]*v[j]; v[i]-=c; return v; }
V applyw(const vector<int>& w, int upto, V v){ for(int k=upto-1;k>=0;k--) v=refl(w[k],v); return v; }
bool pos(const V&v){ for(int x:v) if(x<0) return false; return true; }
int main(int argc,char**argv){
  char t=argv[1][0]; n=atoi(argv[2]); A.assign(n,vector<int>(n,0));
  for(int i=0;i<n;i++) A[i][i]=2;
  for(int i=0;i+1<n;i++){A[i][i+1]=-1;A[i+1][i]=-1;}
  if(t=='B'){A[n-2][n-1]=-1;A[n-1][n-2]=-2;}       // alpha_{n-1} short
  if(t=='C'){A[n-2][n-1]=-2;A[n-1][n-2]=-1;}       // alpha_{n-1} long
  if(t=='E'){ // E_n: chain 0-1-...-(n-2) with node n-1 attached to node 2
    A[n-2][n-1]=0;A[n-1][n-2]=0;A[2][n-1]=-1;A[n-1][2]=-1;}
  if(t=='D'){A[n-2][n-1]=0;A[n-1][n-2]=0;A[n-3][n-1]=-1;A[n-1][n-3]=-1;}
  if(t=='F'){ if(n!=4){fprintf(stderr,"F needs n=4\n");return 1;} A[1][2]=-2;A[2][1]=-1;}
  if(t=='G'){ if(n!=2){fprintf(stderr,"G needs n=2\n");return 1;} A[0][1]=-1;A[1][0]=-3;}
  if(!strchr("ABCDEFG",t)){fprintf(stderr,"unknown type %c\n",t);return 1;}
  if(t=='E'&&(n<6||n>8)){fprintf(stderr,"E needs 6<=n<=8\n");return 1;}
  // positive roots
  set<V> R; vector<V> fr; for(int i=0;i<n;i++){V e(n,0);e[i]=1;R.insert(e);fr.push_back(e);}
  while(!fr.empty()){vector<V> nf; for(auto&r:fr) for(int i=0;i<n;i++){V s=refl(i,r); bool nz=false; for(int x:s) if(x) nz=true; if(nz&&pos(s)&&!R.count(s)){R.insert(s);nf.push_back(s);}} fr=nf;}
  vector<V> roots(R.begin(),R.end()); int N=roots.size(); map<V,int> id; for(int i=0;i<N;i++) id[roots[i]]=i;
  // covers: beta - alpha = simple root (crystallographic)
  vector<array<int,3>> cov; // (a,b,s)
  for(int a=0;a<N;a++) for(int s=0;s<n;s++){V b=roots[a]; b[s]++; if(id.count(b)) cov.push_back({a,id[b],s});}
  vector<int> perm(n); iota(perm.begin(),perm.end(),0); long long cnt=0,bad=0;
  do{
    vector<int> w; int k=0;
    while((int)w.size()<N){ int s=perm[k%n]; k++; V e(n,0); e[s]=1; if(pos(applyw(w,w.size(),e))) w.push_back(s); }
    vector<int> P(N);
    for(int m=0;m<N;m++){ V e(n,0); e[w[m]]=1; P[id[applyw(w,m,e)]]=m; }
    // potential via BFS
    vector<vector<pair<int,int>>> adj(N); for(auto&c:cov){int d=(P[c[0]]<P[c[1]])?1:-1; adj[c[0]].push_back({c[1],d}); adj[c[1]].push_back({c[0],-d});}
    vector<int> g(N,INT_MIN); g[0]=0; deque<int> q{0}; bool ok=true;
    while(!q.empty()){int u=q.front();q.pop_front(); for(auto [v,d]:adj[u]){ if(g[v]==INT_MIN){g[v]=g[u]+d;q.push_back(v);} else if(g[v]!=g[u]+d) ok=false; }}
    cnt++; if(!ok) bad++;
  }while(next_permutation(perm.begin(),perm.end()));
  printf("%c%d: roots %d covers %zu, Coxeter words checked %lld, without potential %lld\n",t,n,N,cov.size(),cnt,bad);
}
