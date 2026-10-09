// Native MaxEnt pipeline extension. Keeps the installed solver and its ABI unchanged.
// Python supplies settings/arrays and receives results; all numerical work is native.
#define NOMINMAX
#include <windows.h>
#include "ms/core.hpp"
#include "ms/maxent_pipeline.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string error;
constexpr double missing=std::numeric_limits<double>::quiet_NaN();
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
struct Module {
    HMODULE handle;
    explicit Module(const wchar_t* path):handle(LoadLibraryExW(path,nullptr,LOAD_WITH_ALTERED_SEARCH_PATH)){
        require(handle!=nullptr,"Cannot load the installed native engine");
    }
    ~Module(){FreeLibrary(handle);}
    template<class F> F get(const char* name){auto p=GetProcAddress(handle,name);if(!p)throw std::runtime_error(std::string("Native engine function missing: ")+name);F f;static_assert(sizeof(f)==sizeof(p));std::memcpy(&f,&p,sizeof(f));return f;}
    void checked(int rc){if(rc)throw std::runtime_error(get<const char*(*)()>("ms_last_error")());}
};
struct Result {
    std::array<ms::Vector,5> x,y;
    std::vector<std::array<double,12>> peaks;
    std::array<double,9> summary{}; // base, base m/z, threshold, n, used, peaks, lo, hi, setting
    std::array<double,3> metrics{}; // R2, chi2, rounds
    std::string notes;
};
void read_array(Module& core,void* out,const char* name,ms::Vector& x,ms::Vector& y){
    double *px=nullptr,*py=nullptr;auto n=core.get<long(*)(void*,double**,double**)>(name)(out,&px,&py);
    require(n>=0&&n<=50000000&&(!n||(px&&py)),"Invalid native result array");
    if(n){x.assign(px,px+n);y.assign(py,py+n);}
}
void fractions(Result& r){double top=1e-30,total=0;for(auto& p:r.peaks){top=std::max(top,p[1]);total+=p[1];}if(total==0)total=1;for(auto& p:r.peaks){p[5]=100*p[1]/top;p[6]=100*p[1]/total;}}
double width(const ms::Vector& x,const ms::Vector& y,std::size_t i){
    double half=y[i]*.5;auto lo=i,hi=i;while(lo>0&&y[lo]>half)--lo;while(hi+1<y.size()&&y[hi]>half)++hi;
    if(lo==i||hi==i||y[lo]>half||y[hi]>half)return 0;
    double a=x[lo]+(half-y[lo])*(x[lo+1]-x[lo])/(y[lo+1]-y[lo]);
    double b=x[hi-1]+(y[hi-1]-half)*(x[hi]-x[hi-1])/(y[hi-1]-y[hi]);return b-a;
}
double resolution(const ms::Vector& x,const ms::Vector& y){
    std::vector<std::size_t> order(y.size());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return y[a]>y[b];});
    ms::Vector used,values;
    for(std::size_t k=0;k<std::min<std::size_t>(5000,order.size());++k){auto i=order[k];if(y[i]<.05*y[order[0]]||values.size()>=5)break;
        if(std::any_of(used.begin(),used.end(),[&](double v){return std::abs(x[i]-v)<.2;}))continue;
        used.push_back(x[i]);auto w=width(x,y,i);if(w>0)values.push_back(x[i]/w);
    }
    if(!values.empty())return ms::median(values);
    auto i=order[0];auto w=width(x,y,i);if(w<=0){ms::Vector d;for(std::size_t j=1;j<x.size();++j)d.push_back(x[j]-x[j-1]);w=5*ms::median(d);}return x[i]/std::max(w,1e-9);
}
double averagine(double mass){
    const double n=mass/111.1254;double var=0;
    const std::vector<std::pair<double,std::vector<ms::Range>>> atoms={
        {4.9384,{{0,.9893},{1.00336,.0107}}},{7.7583,{{0,.999885},{1.00628,.000115}}},
        {1.3577,{{0,.99636},{.99704,.00364}}},{1.4773,{{0,.99757},{1.00422,.00038},{2.00425,.00205}}},
        {.0417,{{0,.9499},{.99939,.0075},{1.99580,.0425}}}};
    for(auto& [count,iso]:atoms){double a=0,b=0;for(auto [d,w]:iso){a+=d*w;b+=d*d*w;}var+=count*n*(b-a*a);}return std::sqrt(var);
}
void isotope_list(Module& core,Result& r){
    if(r.x[0].empty())return;
    double *a=nullptr,*b=nullptr;long n=0;
    core.checked(core.get<int(*)(double*,double*,long,double,long,double**,double**,long*)>("ms_isotope_peaks")(
        r.x[0].data(),r.y[0].data(),long(r.x[0].size()),.01,400,&a,&b,&n));
    if(n){r.x[4].assign(a,a+n);r.y[4].assign(b,b+n);}
}
void grouped_peaks(Module& core,Result& r,double threshold){
    double *m=nullptr,*a=nullptr,*h=nullptr,*area=nullptr,*first=nullptr,*last=nullptr;int* ni=nullptr;long n=0;
    core.checked(core.get<int(*)(double*,double*,long,double,double,double**,double**,double**,double**,int**,double**,double**,long*)>("ms_group_isotopes2")(
        r.x[0].data(),r.y[0].data(),long(r.x[0].size()),threshold,1.00235,&m,&a,&h,&area,&ni,&first,&last,&n));
    r.peaks.clear();for(long i=0;i<n;++i)r.peaks.push_back({m[i],h[i],area[i],a[i],double(ni[i]),0,0,missing,missing,missing,0,1});
}
void refine(Module& core,Result& r,const Params& p){
    auto f=core.get<int(*)(double*,double*,long,double*,double*,long,double,int,int,double,double,double*)>("ms_refine_species");
    for(auto& q:r.peaks){if(q[4]<=0)continue;double v[6]{};int rc=f(r.x[0].data(),r.y[0].data(),long(r.x[0].size()),r.x[3].data(),r.y[3].data(),long(r.x[3].size()),p.adduct_mass,p.z_lo,p.z_hi,1.00235,q[3],v);
        core.checked(rc);if(v[0]){q[7]=q[3];if(v[5]!=0)q[3]=v[1];if(std::isfinite(v[2])){q[8]=v[2];q[9]=v[3];}q[10]=v[4];}
    }
}
void validate(const Params& p,double display){
    require(p.z_lo>=1&&p.z_hi>=p.z_lo&&p.z_hi-p.z_lo<1000,"Invalid charge range");
    require(std::isfinite(p.mass_lo)&&std::isfinite(p.mass_hi)&&p.mass_lo>0&&p.mass_hi>p.mass_lo,"Invalid mass range");
    require(std::isfinite(p.mass_step)&&p.mass_step>=.001,"Mass step must be finite and at least 0.001 Da");
    // The original solver accepts broad ranges and manages its own grids.
    // Do not reject valid settings with the former Python solver's 400k cap.
    // Only guard arithmetic/index overflow in the long-based engine ABI.
    const double intervals=(p.mass_hi-p.mass_lo)/p.mass_step;
    require(std::isfinite(intervals)&&intervals<=double(std::numeric_limits<long>::max()-1),
            "Mass range and step require too many points; use a narrower range or a larger mass step");
    require(p.rounds>=1&&p.iterations>=1&&p.min_rounds>=1,"Rounds and iterations must be positive whole numbers");
    require(std::isfinite(display)&&display>=0,"Display resolution must be finite and non-negative");
    require(std::isfinite(p.resolution)&&p.resolution>=-1&&std::isfinite(p.peak_width)&&p.peak_width>=0,"Invalid resolution or peak width");
    require(std::isfinite(p.adduct_mass)&&p.adduct_mass!=0,"Invalid charge carrier mass");
    require(std::isfinite(p.peak_threshold)&&p.peak_threshold>=0&&p.peak_threshold<=1&&std::isfinite(p.peak_window)&&p.peak_window>0,"Invalid peak threshold or window");
    require(std::isfinite(p.min_intensity)&&p.min_intensity>=0&&std::isfinite(p.baseline_width)&&p.baseline_width>0,"Invalid minimum intensity or baseline width");
    require(std::isfinite(p.chi2_target)&&p.chi2_target>=0,"Invalid chi squared target");
}
std::unique_ptr<Result> run(Module& core,const double* x,const double* y,long n,const Params& original,double display,Progress progress,void* user){
    validate(original,display);require(x&&y&&n>=10&&n<=20000000,"Too few or excessive input points");
    auto result=std::make_unique<Result>();auto& r=*result;double *rx=nullptr,*ry=nullptr;long nr=0;
    core.checked(core.get<int(*)(const double*,const double*,long,double,double,double**,double**,long*)>("ms_restrict")(
        x,y,n,original.mz_lo,original.mz_hi,&rx,&ry,&nr));require(nr>=10,"Too few data points in the m/z range");
    r.x[3].assign(rx,rx+nr);r.y[3].assign(ry,ry+nr);
    auto base=std::max_element(r.y[3].begin(),r.y[3].end());double thr=original.min_intensity;
    if(original.min_intensity_pct)thr=(thr/100.0)*(*base);
    std::vector<unsigned char> mask(nr,0);long npk=0,used=0;
    if(thr>0)core.checked(core.get<int(*)(double*,double*,long,double,unsigned char*,long*)>("ms_peak_mask")(
        r.x[3].data(),r.y[3].data(),nr,thr,mask.data(),&npk));
    else {for(long i=0;i<nr;++i)mask[i]=r.y[3][i]>0;for(long i=1;i+1<nr;++i)if(r.y[3][i]>r.y[3][i-1]&&r.y[3][i]>=r.y[3][i+1]&&r.y[3][i]>0)++npk;}
    used=std::count(mask.begin(),mask.end(),1);
    r.summary={*base,r.x[3][std::size_t(base-r.y[3].begin())],thr,double(nr),double(used),double(npk),r.x[3].front(),r.x[3].back(),original.min_intensity};
    // With a minimum intensity the solver gets the data before baseline and minimum intensity, with
    // both settings: it applies them itself, once and in the same order (peaks judged on the data as
    // shown, zeroed after the baseline subtraction), and its noise model gets the right mask. Giving
    // it the data processed below with the minimum intensity still set applied the threshold a second
    // time, on the baseline subtracted data, and weaker species were lost. Without a minimum
    // intensity it gets the processed data, as before.
    ms::Vector solver_y;
    if(thr>0)solver_y=r.y[3];
    Params p=original;
    if(p.baseline){ms::Vector out(nr);core.checked(core.get<int(*)(double*,double*,long,double,double,double*)>("ms_subtract_baseline")(
        r.x[3].data(),r.y[3].data(),nr,0,p.baseline_width,out.data()));r.y[3]=std::move(out);p.baseline=0;}
    if(thr>0){for(long i=0;i<nr;++i)if(!mask[i])r.y[3][i]=0;
        require(npk>0&&used>=5&&std::any_of(r.y[3].begin(),r.y[3].end(),[](double v){return v>0;}),"No data above the minimum intensity; lower the minimum intensity");}
    void* raw=nullptr;core.checked(core.get<int(*)(double*,double*,long,const Params*,void**,Progress,void*)>("ms_maxent")(
        r.x[3].data(),thr>0?solver_y.data():r.y[3].data(),nr,thr>0?&original:&p,&raw,progress,user));require(raw!=nullptr,"Native solver returned no result");
    auto free=core.get<void(*)(void*)>("ms_maxent_free");std::unique_ptr<void,decltype(free)> owned(raw,free);
    read_array(core,raw,"ms_maxent_mass",r.x[0],r.y[0]);read_array(core,raw,"ms_maxent_fit",r.x[1],r.y[1]);read_array(core,raw,"ms_maxent_zdist",r.x[2],r.y[2]);
    r.metrics={core.get<double(*)(void*)>("ms_maxent_r2")(raw),core.get<double(*)(void*)>("ms_maxent_chi2")(raw),double(core.get<int(*)(void*)>("ms_maxent_rounds")(raw))};
    r.notes=core.get<const char*(*)(void*)>("ms_maxent_notes")(raw);
    const bool resolved=p.resolved_isotopes&&r.notes.find("resolved isotopes")!=std::string::npos;
    double *m=nullptr,*h=nullptr,*a=nullptr,*ap=nullptr;int* ni=nullptr;
    auto count=core.get<long(*)(void*,double**,double**,double**,double**,int**)>("ms_maxent_peaks")(raw,&m,&h,&a,&ap,&ni);
    for(long i=0;i<count;++i)r.peaks.push_back({m[i],h[i],a[i],ap[i],resolved?double(ni[i]):0,0,0,missing,missing,missing,0,1});
    // Broadening changes the display spectrum, not the fit/charge distribution.
    if(display>0&&p.resolution>=0&&!r.x[0].empty()){
        double instrument=p.resolution>0?p.resolution:resolution(r.x[3],r.y[3]);
        if(display<instrument){double effective=instrument;
            if(!resolved){double centre=.5*(p.mass_lo+p.mass_hi);effective=1/std::sqrt(1/(instrument*instrument)+std::pow(2.3548*averagine(centre)/centre,2));}
            ms::Spectrum spectrum;for(std::size_t i=0;i<r.x[0].size();++i)spectrum.push_back({r.x[0][i],r.y[0][i]});
            auto shown=ms::broaden_mass(spectrum,display,effective);for(std::size_t i=0;i<shown.size();++i)r.y[0][i]=shown[i].y;
            if(resolved)grouped_peaks(core,r,p.peak_threshold);
            else {double window=std::max(p.peak_window,2*p.mass_step);window=std::max(window,.5*.5*(p.mass_lo+p.mass_hi)/effective);
                auto peaks=ms::pick_mass_peaks(shown,window,p.peak_threshold);r.peaks.clear();for(auto q:peaks)r.peaks.push_back({q.mass,q.height,q.area,q.apex,0,0,0,missing,missing,missing,0,1});}
            r.notes+="; shown at resolving power "+std::to_string(display);
        }
    }
    fractions(r);if(resolved)isotope_list(core,r);refine(core,r,p);r.notes+="; native MaxEnt pipeline";
    return result;
}
}
#define API extern "C" __declspec(dllexport)
API unsigned mx_abi(){return 1;}
API unsigned mx_params_size(){return sizeof(Params);}
API const char* mx_error(){return error.c_str();}
API int mx_run(const wchar_t* path,const double* x,const double* y,long n,const Params* params,double display,Progress progress,void* user,void** output) noexcept {
    if(output)*output=nullptr;
    try {error.clear();require(path&&params&&output,"Missing native input or output");Module core(path);auto r=run(core,x,y,n,*params,display,progress,user);*output=r.release();return 0;}
    catch(const std::exception& ex){error=ex.what();return -1;}catch(...){error="Unexpected native MaxEnt error";return -1;}
}
API void mx_free(void* handle){delete static_cast<Result*>(handle);}
API long mx_array(void* handle,int kind,double** x,double** y){if(!handle||!x||!y||kind<0||kind>=5)return -1;auto& r=*static_cast<Result*>(handle);*x=r.x[kind].data();*y=r.y[kind].data();return long(r.x[kind].size());}
API long mx_peaks(void* handle,double** rows){if(!handle||!rows)return -1;auto& r=*static_cast<Result*>(handle);*rows=r.peaks.empty()?nullptr:r.peaks[0].data();return long(r.peaks.size());}
API const double* mx_summary(void* handle){return handle?static_cast<Result*>(handle)->summary.data():nullptr;}
API const double* mx_metrics(void* handle){return handle?static_cast<Result*>(handle)->metrics.data():nullptr;}
API const char* mx_notes(void* handle){return handle?static_cast<Result*>(handle)->notes.c_str():"";}
// Direct numerical primitive for strict comparisons with the previous Python feature.
API int mx_broaden(const double* xy,long n,double display,double model,double* out) noexcept {
    try {error.clear();require(xy&&out&&n>=0,"Invalid broadening input");ms::Spectrum s;for(long i=0;i<n;++i)s.push_back({xy[2*i],xy[2*i+1]});auto b=ms::broaden_mass(s,display,model);for(long i=0;i<n;++i){out[2*i]=b[i].x;out[2*i+1]=b[i].y;}return 0;}
    catch(const std::exception& ex){error=ex.what();return -1;}
}
