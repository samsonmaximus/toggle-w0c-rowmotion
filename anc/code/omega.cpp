// Test: orientation of cover (g, g+a_s) is "up" (g toggled before g+a_s) iff Omega(g, a_s) > 0, where
// Omega(x,y) = <x,y> - <y,x>, <.,.> Euler form of valued quiver with orientation i->j iff s_i before s_j in c (or reverse).
#include <bits/stdc++.h>
using namespace std;
int n; vector<vector<int>> A; typedef vector<int> V;
V refl(int i, V v){int c=0;for(int j=0;j<n;j++)c+=A[i][j]*v[j];v[i]-=c;return v;}
V applyw(const vector<int>&w,int upto,V v){for(int k=upto-1;k>=0;k--)v=refl(w[k],v);return v;}
bool pos(const V&v){for(int x:v)if(x<0)return false;return true;}
int main(int argc,char**argv){
  char t=argv[1][0]; n=atoi(argv[2]); A.assign(n,vector<int>(n,0));
  for(int i=0;i<n;i++)A[i][i]=2; for(int i=0;i+1<n;i++){A[i][i+1]=-1;A[i+1][i]=-1;}
  vector<int> f(n,1);
  if(t=='B'){A[n-2][n-1]=-1;A[n-1][n-2]=-2; for(int i=0;i<n-1;i++) f[i]=2;}   // alpha_{n-1} short: f = |alpha|^2/2 -> long=2, short=1
  if(t=='C'){A[n-2][n-1]=-2;A[n-1][n-2]=-1; f[n-1]=2;}
  if(t=='D'){A[n-2][n-1]=0;A[n-1][n-2]=0;A[n-3][n-1]=-1;A[n-1][n-3]=-1;}
  if(t=='E'){A[n-2][n-1]=0;A[n-1][n-2]=0;A[2][n-1]=-1;A[n-1][2]=-1;}
  if(t=='F'){A[1][2]=-2;A[2][1]=-1; f[2]=2;f[3]=2;}   // F4: alpha2,alpha3 long (0-indexed)
  if(t=='G'){A[0][1]=-1;A[1][0]=-3; f[0]=3;}           // G2: alpha0 long (0-indexed)
  // check symmetrizability: f_i a_ij = f_j a_ji
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) if(f[i]*A[i][j]!=f[j]*A[j][i]){printf("bad symmetrizer\n");return 1;}
  set<V> R; vector<V> fr; for(int i=0;i<n;i++){V e(n,0);e[i]=1;R.insert(e);fr.push_back(e);}
  while(!fr.empty()){vector<V> nf;for(auto&r:fr)for(int i=0;i<n;i++){V s=refl(i,r);bool nz=false;for(int x:s)if(x)nz=true;if(nz&&pos(s)&&!R.count(s)){R.insert(s);nf.push_back(s);}}fr=nf;}
  vector<V> roots(R.begin(),R.end());int N=roots.size();map<V,int>id;for(int i=0;i<N;i++)id[roots[i]]=i;
  vector<array<int,3>> cov; for(int a=0;a<N;a++)for(int s=0;s<n;s++){V b=roots[a];b[s]++;if(id.count(b))cov.push_back({a,id[b],s});}
  vector<int> perm(n);iota(perm.begin(),perm.end(),0); long long cnt=0, viol[2]={0,0}, zero=0, pairviol=0;
  do{ vector<int> w;int k=0; while((int)w.size()<N){int s=perm[k%n];k++;V e(n,0);e[s]=1;if(pos(applyw(w,w.size(),e)))w.push_back(s);}
    vector<int>P(N); for(int m=0;m<N;m++){V e(n,0);e[w[m]]=1;P[id[applyw(w,m,e)]]=m;}
    vector<int> where(n); for(int i=0;i<n;i++) where[perm[i]]=i;
    for(int conv=0;conv<2;conv++){
      // Euler matrix E: E[i][i]=f_i ; arrow i->j : E[i][j]= f_i*A[i][j]
      vector<vector<long long>> E(n,vector<long long>(n,0));
      for(int i=0;i<n;i++){E[i][i]=f[i]; for(int j=0;j<n;j++) if(i!=j&&A[i][j]!=0){ bool arr = conv? (where[i]<where[j]) : (where[i]>where[j]); if(arr) E[i][j]=(long long)f[i]*A[i][j]; }}
      for(auto&cv:cov){ const V&g=roots[cv[0]]; int s=cv[2];
        long long gs=0,sg=0; for(int i=0;i<n;i++){ gs+=g[i]*E[i][s]; sg+=E[s][i]*g[i]; }
        long long Om=gs-sg; if(conv==0 && Om==0) zero++;
        bool up=P[cv[0]]<P[cv[1]]; if(Om==0 || (Om>0)!=up) viol[conv]++; }
    }
    { vector<vector<long long>> E(n,vector<long long>(n,0));
      for(int i=0;i<n;i++){E[i][i]=f[i]; for(int j=0;j<n;j++) if(i!=j&&A[i][j]!=0 && where[i]>where[j]) E[i][j]=(long long)f[i]*A[i][j];}
      for(int a=0;a<N;a++)for(int b=0;b<N;b++) if(P[a]<P[b]){ long long x=0,y=0; for(int i=0;i<n;i++)for(int j=0;j<n;j++){x+=roots[a][i]*E[i][j]*roots[b][j]; y+=roots[b][i]*E[i][j]*roots[a][j];} if(x-y<0) pairviol++; } }
    cnt++; }while(next_permutation(perm.begin(),perm.end()));
  printf("%c%d: words %lld covers %zu | Omega==0 on covers %lld | cover-rule violations %lld | pairs beta<beta' with Omega<0: %lld\n",t,n,cnt,cov.size(),zero,viol[0],pairviol);
}
