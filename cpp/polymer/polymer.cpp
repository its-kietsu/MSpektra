#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
#define API extern "C" __declspec(dllexport)
struct Summary { int used; double mn,mw,dispersity; };
struct Candidate { double mass,score; int support,single; };
static thread_local char error_text[256];
API unsigned poly_abi(){return 1;}
API unsigned poly_summary_size(){return sizeof(Summary);}
API unsigned poly_candidate_size(){return sizeof(Candidate);}
API const char* poly_error(){return error_text;}
static void check(const double* m,int n){
    if(!m||n<1||n>5000) throw std::runtime_error("Use between 1 and 5000 species");
    for(int i=0;i<n;++i) if(!std::isfinite(m[i])||m[i]<=0) throw std::runtime_error("Species masses must be finite and positive");
}
static int fail(const std::exception& e){std::strncpy(error_text,e.what(),255);error_text[255]=0;return -1;}
API int poly_distribution(const double* m,const double* amount,const double* response,int n,int basis,
                           Summary* out,double* nf,double* wf){
    error_text[0]=0;
    try{
        check(m,n);if(!amount||!response||!out||!nf||!wf||basis<0||basis>1) throw std::runtime_error("Invalid distribution input");
        std::vector<long double> a(n);long double total=0,first=0,second=0;
        int used=0;
        for(int i=0;i<n;++i){
            if(!std::isfinite(amount[i])||amount[i]<0||!std::isfinite(response[i])||response[i]<=0)
                throw std::runtime_error("Amounts must be nonnegative; response factors must be finite and positive");
            a[i]=(long double)amount[i]/response[i];if(basis) a[i]/=m[i];
            total+=a[i];first+=a[i]*m[i];second+=a[i]*m[i]*m[i];if(a[i]>0)++used;
        }
        if(total<=0) throw std::runtime_error("No positive signal remains; choose positive species areas or heights");
        Summary r{};r.used=used;r.mn=double(first/total);r.mw=double(second/first);r.dispersity=std::max(1.,r.mw/r.mn);
        if(!std::isfinite(r.mn)||!std::isfinite(r.mw)||!std::isfinite(r.dispersity)) throw std::runtime_error("Numerical scale exceeded");
        for(int i=0;i<n;++i){nf[i]=double(a[i]/total);wf[i]=double(a[i]*m[i]/first);}
        *out=r;return 0;
    }catch(const std::exception& e){return fail(e);}catch(...){std::strcpy(error_text,"Native distribution error");return -1;}
}
static void repeat_check(double repeat,double tol){
    if(!std::isfinite(repeat)||repeat<=0||!std::isfinite(tol)||tol<1e-8||tol>=repeat/2)
        throw std::runtime_error("Repeat mass must be positive; tolerance must be >= 1e-8 Da and less than half the repeat mass");
}
API int poly_assign(const double* m,int n,double repeat,double ends,double tol,
                    double* dp,double* predicted,double* delta,double* ppm,int* matched){
    error_text[0]=0;
    try{
        check(m,n);repeat_check(repeat,tol);
        if(!std::isfinite(ends)||!dp||!predicted||!delta||!ppm||!matched) throw std::runtime_error("Invalid end-group assignment input");
        for(int i=0;i<n;++i){
            double q=std::round((m[i]-ends)/repeat);
            if(!std::isfinite(q)||q<1||q>1e9){dp[i]=predicted[i]=delta[i]=ppm[i]=std::numeric_limits<double>::quiet_NaN();matched[i]=0;continue;}
            double pred=ends+q*repeat;
            if(!std::isfinite(pred)||pred<=0) throw std::runtime_error("Predicted mass exceeds valid range");
            predicted[i]=pred;delta[i]=m[i]-pred;ppm[i]=delta[i]/pred*1e6;
            matched[i]=std::abs(delta[i])<=tol;dp[i]=matched[i]?q:std::numeric_limits<double>::quiet_NaN();
        }
        return 0;
    }catch(const std::exception& e){return fail(e);}catch(...){std::strcpy(error_text,"Native assignment error");return -1;}
}
// Candidate end-group residues modulo R. This cannot identify E or DP uniquely.
API int poly_offsets(const double* m,int n,double repeat,double tol,Candidate* out,int capacity){
    error_text[0]=0;
    try{
        check(m,n);repeat_check(repeat,tol);if(!out||capacity<1||capacity>100) throw std::runtime_error("Invalid candidate buffer");
        std::vector<double> phases(n);
        for(int i=0;i<n;++i)phases[i]=std::fmod(m[i],repeat);
        std::vector<Candidate> candidates;
        for(double centre:phases){
            int count=0;long double sum=0;
            for(double phase:phases){double d=std::remainder(phase-centre,repeat);if(std::abs(d)<=tol){++count;sum+=d;}}
            double mass=double(std::fmod((long double)centre+sum/count+repeat,(long double)repeat));
            count=0;for(double phase:phases)if(std::abs(std::remainder(phase-mass,repeat))<=tol)++count;
            candidates.push_back({mass,double(count),count,0});
        }
        std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b){return a.support>b.support||(a.support==b.support&&a.mass<b.mass);});
        int k=0;for(const auto& c:candidates){
            bool duplicate=false;for(int j=0;j<k;++j)if(std::abs(std::remainder(c.mass-out[j].mass,repeat))<=2*tol)duplicate=true;
            if(!duplicate){out[k++]=c;if(k==capacity)break;}
        }return k;
    }catch(const std::exception& e){return fail(e);}catch(...){std::strcpy(error_text,"Native offset error");return -1;}
}
API int poly_repeats(const double* m,int n,double lo,double hi,double tol,Candidate* out,int capacity){
    error_text[0]=0;
    try{
        check(m,n);
        if(n<3)throw std::runtime_error("At least three species are needed to find a repeating series");
        if(!std::isfinite(lo)||!std::isfinite(hi)||lo<=0||hi<=lo||hi>1e7||!std::isfinite(tol)||tol<1e-8||tol>=lo/2)
            throw std::runtime_error("Use a positive repeat-mass range and tolerance >= 1e-8 Da, less than half its minimum");
        if(!out||capacity<1||capacity>100)throw std::runtime_error("Invalid candidate buffer");
        std::vector<double> sorted(m,m+n);std::sort(sorted.begin(),sorted.end());
        std::vector<double> gaps;
        std::map<long long,std::pair<long double,int>> bins;
        double width=tol/2;
        for(int i=0;i<n;++i)for(int j=i+1;j<std::min(n,i+5);++j){
            double d=sorted[j]-sorted[i];if(d<=0||d>hi*6+tol)continue;gaps.push_back(d);
            for(int k=1;k<=6;++k){double r=d/k;if(r<lo||r>hi)continue;auto& b=bins[std::llround(r/width)];b.first+=r;++b.second;}
        }
        std::vector<std::pair<long double,int>> ranked;for(const auto& b:bins)ranked.push_back(b.second);
        std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b){return a.second>b.second||(a.second==b.second&&a.first/a.second<b.first/b.second);});
        std::vector<Candidate> candidates;
        for(size_t i=0;i<std::min(size_t(512),ranked.size());++i){
            Candidate c{double(ranked[i].first/ranked[i].second),0,0,0};
            for(double d:gaps){double k=std::round(d/c.mass);if(k>=1&&k<=6&&std::abs(d-k*c.mass)<=tol){++c.support;if(k==1)++c.single;c.score+=1/(k*k);}}
            if(c.support>=2)candidates.push_back(c);
        }
        std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b){return a.score>b.score||(a.score==b.score&&a.mass<b.mass);});
        int k=0;for(const auto& c:candidates){
            bool duplicate=false;for(int j=0;j<k;++j)if(std::abs(c.mass-out[j].mass)<=tol)duplicate=true;
            if(!duplicate){out[k++]=c;if(k==capacity)break;}
        }return k;
    }catch(const std::exception& e){return fail(e);}catch(...){std::strcpy(error_text,"Native repeat finder error");return -1;}
}
