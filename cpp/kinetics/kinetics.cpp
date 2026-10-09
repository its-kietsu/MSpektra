// Native least-squares fitting. ABI contains no owning pointers; no external runtime needed.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
#define API extern "C" __declspec(dllexport)
// Models: 0 straight line (zero order), 1 first-order decay, 2 first-order rise, 3 second-order decay
// y = C + A / (1 + q (X - X0)), with q = k A (k: the second-order rate constant, in 1 / (signal x X)).
struct Result {
    int model, count, parameters, dof, warnings;
    double origin, values[3], errors[3], sse, rmse, r2, half_life, half_life_error;
    // ABI 2: reaction order, Akaike criteria (to compare models on the same data; AICc NaN when n <= p + 1),
    // the second-order rate constant k = q / A and its uncertainty (NaN for the other models)
    int order, reserved;
    double aic, aicc, k2, k2_error;
};
static thread_local char last_error[256];
static constexpr double missing_value = std::numeric_limits<double>::quiet_NaN();
API unsigned kin_abi() { return 2; }
API unsigned kin_result_size() { return sizeof(Result); }
API const char* kin_error() { return last_error; }
struct Trial { double c=0, a=0, sse=0; };
static double shape(int model, double q, double t) {
    if(model==3) return 1.0/(1.0+q*t);
    return model==1 ? std::exp(-q*t) : -std::expm1(-q*t);
}
// d shape / d q
static double dshape(int model, double q, double t) {
    if(model==3) { double d=1.0+q*t; return -t/(d*d); }
    return (model==1?-1:1)*t*std::exp(-q*t);
}
// Solve C and A exactly at each log-rate. Centered sums avoid intercept cancellation.
static Trial trial(const std::vector<double>& t, const std::vector<double>& y,
                   int model, double q, bool offset) {
    long double sf=0, sy=0, ff=0, fy=0;
    std::vector<double> f(t.size());
    for (size_t i=0;i<t.size();++i) { f[i]=shape(model,q,t[i]); sf+=f[i]; sy+=y[i]; }
    double mf=offset?double(sf/t.size()):0, my=offset?double(sy/t.size()):0;
    for (size_t i=0;i<t.size();++i) { ff+=(f[i]-mf)*(f[i]-mf); fy+=(f[i]-mf)*(y[i]-my); }
    Trial r;
    r.a=ff>0?std::max(0.,double(fy/ff)):0;
    r.c=offset?my-r.a*mf:0;
    long double s=0;
    for (size_t i=0;i<t.size();++i) { double d=y[i]-r.c-r.a*f[i]; s+=d*d; }
    r.sse=double(s); return r;
}
// Invert a correlation-scaled information matrix; reject near-singular uncertainty estimates.
static bool inverse(double a[3][3], int n, double out[3][3]) {
    double b[3][6]={}, scale[3]={};
    for(int i=0;i<n;++i) { scale[i]=std::sqrt(a[i][i]); if(!(scale[i]>0)) return false; }
    for(int i=0;i<n;++i) { for(int j=0;j<n;++j) b[i][j]=a[i][j]/scale[i]/scale[j]; b[i][n+i]=1; }
    for(int j=0;j<n;++j) {
        int p=j; for(int i=j+1;i<n;++i) if(std::abs(b[i][j])>std::abs(b[p][j])) p=i;
        if(std::abs(b[p][j])<1e-10) return false;
        for(int k=0;k<2*n;++k) std::swap(b[j][k],b[p][k]);
        double v=b[j][j]; for(int k=0;k<2*n;++k) b[j][k]/=v;
        for(int i=0;i<n;++i) if(i!=j) { double d=b[i][j]; for(int k=0;k<2*n;++k) b[i][k]-=d*b[j][k]; }
    }
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) out[i][j]=b[i][n+j]/scale[i]/scale[j];
    return true;
}
API int kin_fit(const double* x, const double* y, int n, int model, int offset,
                Result* result, double* fitted, double* residuals) {
    last_error[0]=0;
    try {
        if(!x||!y||!result||!fitted||!residuals) throw std::runtime_error("Missing fit buffers");
        if(n<2||n>100000) throw std::runtime_error("Provide between 2 and 100000 observations");
        if(model<0||model>3||(offset!=0&&offset!=1)) throw std::runtime_error("Invalid model or baseline option");
        int p=model==0?2:(offset?3:2);
        if(n<=p) throw std::runtime_error("More observations than fitted parameters are required");
        std::vector<double> sorted(x,x+n);
        double ymin=y[0], ymax=y[0];
        for(int i=0;i<n;++i) {
            if(!std::isfinite(x[i])||!std::isfinite(y[i])) throw std::runtime_error("X and area values must be finite");
            ymin=std::min(ymin,y[i]); ymax=std::max(ymax,y[i]);
        }
        std::sort(sorted.begin(),sorted.end());
        int unique=int(std::unique(sorted.begin(),sorted.end())-sorted.begin());
        if(unique<p) throw std::runtime_error("Too few distinct X values for this model");
        double origin=sorted.front(), span=sorted.back()-origin;
        if(!(span>0)||!std::isfinite(span)) throw std::runtime_error("X values must span a finite nonzero range");
        double ys=ymax-ymin;
        if(!(ys>0)||!std::isfinite(ys)||ys<=std::max(std::abs(ymin),std::abs(ymax))*1e-12)
            throw std::runtime_error("Areas are constant or too close to constant to determine a rate");
        std::vector<double> t(n), v(n);
        for(int i=0;i<n;++i) { t[i]=(x[i]-origin)/span; v[i]=(y[i]-ymin)/ys; }
        double c=0,a=0,q=0;
        int warnings=0;
        if(model==0) {
            long double sx=0,sy=0,xx=0,xy=0;
            for(int i=0;i<n;++i) { sx+=t[i];sy+=v[i]; }
            double mx=double(sx/n),my=double(sy/n);
            for(int i=0;i<n;++i) { xx+=(t[i]-mx)*(t[i]-mx);xy+=(t[i]-mx)*(v[i]-my); }
            a=double(xy/xx);c=my-a*mx;
        } else {
            // When baseline is fixed, scale without subtracting ymin.
            if(!offset) for(int i=0;i<n;++i) v[i]=y[i]/ys;
            constexpr int steps=480;
            const double lo=std::log(1e-6), hi=std::log(1e4), h=(hi-lo)/steps;
            int best=0; double cost=std::numeric_limits<double>::infinity();
            for(int j=0;j<=steps;++j) {
                double s=trial(t,v,model,std::exp(lo+j*h),offset).sse;
                if(s<cost) {best=j;cost=s;}
            }
            double l=lo+std::max(0,best-1)*h, r=lo+std::min(steps,best+1)*h;
            constexpr double phi=.6180339887498948482;
            double b=r-phi*(r-l),d=l+phi*(r-l);
            double fb=trial(t,v,model,std::exp(b),offset).sse,fd=trial(t,v,model,std::exp(d),offset).sse;
            for(int j=0;j<100;++j) {
                if(fb<fd) {r=d;d=b;fd=fb;b=r-phi*(r-l);fb=trial(t,v,model,std::exp(b),offset).sse;}
                else {l=b;b=d;fb=fd;d=l+phi*(r-l);fd=trial(t,v,model,std::exp(d),offset).sse;}
            }
            q=std::exp((l+r)/2); Trial z=trial(t,v,model,q,offset);c=z.c;a=z.a;
            if(a<=1e-12) throw std::runtime_error("Areas do not follow the selected growth or decay direction");
            // Sparse X spacing can make high rates unidentifiable before the nominal bound.
            if(best==0||best==steps) warnings|=1;
        }
        Result out{};
        out.model=model;out.count=n;out.parameters=p;out.dof=n-p;out.origin=origin;
        for(double& e:out.errors) e=missing_value;
        out.half_life=out.half_life_error=missing_value;
        out.values[0]=c*ys+((model==0||offset)?ymin:0);
        out.values[1]=model==0?a*ys/span:a*ys;
        out.values[2]=model==0?missing_value:q/span;
        // half of the change: ln 2 / k (first order), 1 / (k A) (second order)
        if(model!=0) out.half_life=(model==3?1.0:std::log(2.))/out.values[2];
        out.order=model==0?0:(model==3?2:1);
        out.k2=out.k2_error=missing_value;
        double info[3][3]={}, inv[3][3]={};
        long double sse=0, normalized_sse=0, mean=0;
        for(int i=0;i<n;++i) mean+=y[i];
        mean/=n; long double sst=0;
        for(int i=0;i<n;++i) {
            double f=model==0?t[i]:shape(model,q,t[i]);
            double yn=c+a*f;
            fitted[i]=yn*ys+((model==0||offset)?ymin:0);
            residuals[i]=y[i]-fitted[i]; sse+=(long double)residuals[i]*residuals[i];
            double nr=residuals[i]/ys; normalized_sse+=(long double)nr*nr;
            sst+=(y[i]-mean)*(y[i]-mean);
            double j[3]={1,f,0};
            if(model) { j[2]=a*dshape(model,q,t[i]); if(!offset) {j[0]=j[1];j[1]=j[2];} }
            for(int u=0;u<p;++u) for(int w=0;w<p;++w) info[u][w]+=j[u]*j[w];
        }
        out.sse=double(sse);out.rmse=double(std::sqrt(sse/n));out.r2=1-double(sse/sst);
        {   // Akaike: n ln(SSE/n) + 2p (+ the small sample term); a perfect fit is held at a tiny SSE
            double s=std::max(double(sse/n), 1e-300);
            out.aic=n*std::log(s)+2.0*p;
            out.aicc=(n-p-1>0)?out.aic+2.0*p*(p+1)/double(n-p-1):missing_value;
        }
        if(!std::isfinite(out.sse)||!std::isfinite(out.rmse)||!std::isfinite(out.r2))
            throw std::runtime_error("Numerical scale exceeded; rescale the input values");
        bool ok=inverse(info,p,inv);
        if(model && q*q*info[p-1][p-1] < 1e-14*n) ok=false;
        double sigma=double(normalized_sse/out.dof);
        if(ok && !(warnings&1)) {
            double e[3]={};for(int i=0;i<p;++i) e[i]=std::sqrt(std::max(0.,inv[i][i]*sigma));
            if(model==0) {out.errors[0]=e[0]*ys;out.errors[1]=e[1]*ys/span;}
            else {
                out.errors[0]=offset?e[0]*ys:0;
                out.errors[1]=e[offset?1:0]*ys;out.errors[2]=e[offset?2:1]/span;
                out.half_life_error=out.half_life*out.errors[2]/out.values[2];
                if(out.errors[2]>=out.values[2]) warnings|=4;
                if(model==3) {  // k = q / A, with the covariance of A and q
                    int ia=offset?1:0, iq=offset?2:1;
                    double ra=out.errors[1]/out.values[1], rq=out.errors[2]/out.values[2];
                    double cov=inv[ia][iq]*sigma/(a*q);  // relative covariance (scale factors cancel)
                    out.k2_error=std::abs(out.values[2]/out.values[1])*std::sqrt(std::max(0.,ra*ra+rq*rq-2*cov));
                }
            }
        } else warnings|=2;
        for(int i=0;i<3;++i) if(!std::isfinite(out.errors[i]) && (model!=0||i<2)) warnings|=2;
        if(model==3) out.k2=out.values[2]/out.values[1];
        out.warnings=warnings;
        for(int i=0;i<3;++i) if(i!=2||model!=0)
            if(!std::isfinite(out.values[i])) throw std::runtime_error("Fitted parameters exceed numerical range");
        *result=out;return 0;
    } catch(const std::exception& e) { std::strncpy(last_error,e.what(),255);last_error[255]=0;return 1; }
      catch(...) { std::strcpy(last_error,"Unexpected native fitting error");return 2; }
}
// Prediction is native too, including the dense curve shown/exported by the interface.
API int kin_predict(const Result* r, const double* x, int n, double* y) {
    last_error[0]=0;
    if(!r||!x||!y||n<1||r->model<0||r->model>3) { std::strcpy(last_error,"Invalid prediction input");return 1; }
    for(int i=0;i<n;++i) {
        double dt=x[i]-r->origin;
        y[i]=r->values[0]+r->values[1]*(r->model==0?dt:shape(r->model,r->values[2],dt));
        if(!std::isfinite(x[i])||!std::isfinite(y[i])) { std::strcpy(last_error,"Prediction is outside the finite model range");return 1; }
    }
    return 0;
}
