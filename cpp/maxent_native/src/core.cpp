#include "ms/core.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace ms {
static void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
static bool finite(double x) { return std::isfinite(x); }
Spectrum broaden_mass(const Spectrum& s,double display_resolution,double model_resolution){
    require(finite(display_resolution)&&display_resolution>0&&finite(model_resolution)&&model_resolution>0,"Resolving powers must be positive and finite");
    if(s.size()<3)return s;
    Vector steps,y;double moment=0,weight=0;for(std::size_t i=0;i<s.size();++i){require(finite(s[i].x)&&finite(s[i].y),"Invalid mass spectrum");if(i){require(s[i].x>s[i-1].x,"Mass axis must increase");steps.push_back(s[i].x-s[i-1].x);}double w=std::max(s[i].y,0.)+1e-30;moment+=s[i].x*w;weight+=w;y.push_back(s[i].y);}
    double m=moment/weight,model=m/model_resolution/2.3548,shown=m/display_resolution/2.3548,extra=std::sqrt(std::max(shown*shown-model*model,0.));if(extra<=0)return s;
    auto blurred=gaussian(y,extra/median(steps));auto out=s;for(std::size_t i=0;i<s.size();++i)out[i].y=blurred[i]*shown/std::max(model,1e-30);return out;
}
std::vector<Peak> group_isotopes(const Spectrum& s,double threshold,double spacing){
    require(finite(threshold)&&threshold>=0&&threshold<=1&&finite(spacing)&&spacing>0,"Invalid isotope grouping parameters");
    auto iso=pick_mass_peaks(s,.35,threshold*.2);if(iso.empty())return {};
    std::sort(iso.begin(),iso.end(),[](const Peak& a,const Peak& b){return a.mass<b.mass;});Vector steps;for(std::size_t i=1;i<s.size();++i)steps.push_back(s[i].x-s[i-1].x);double dx=steps.empty()?1:median(steps);
    std::vector<Peak> out;std::size_t first=0;
    for(std::size_t end=1;end<=iso.size();++end){if(end<iso.size()){double gap=iso[end].mass-iso[end-1].mass;if(std::abs(gap-spacing)<=.12||gap<.6)continue;}
        auto top=std::max_element(iso.begin()+first,iso.begin()+end,[](const Peak& a,const Peak& b){return a.height<b.height;});double lo=iso[first].mass-.5,hi=iso[end-1].mass+.5,moment=0,sum=0;
        for(auto p:s)if(p.x>=lo&&p.x<=hi){sum+=p.y;moment+=p.x*p.y;}
        out.push_back({sum>0?moment/sum:iso[first].mass,top->mass,top->height,sum*dx,0,0,int(end-first)});first=end;
    }
    double highest=0;for(const auto& p:out)highest=std::max(highest,p.height);std::erase_if(out,[&](const Peak& p){return p.height<threshold*highest;});double total=0;for(const auto& p:out)total+=p.height;if(total==0)total=1;for(auto& p:out){p.relative=highest>0?100*p.height/highest:0;p.fraction=100*p.height/total;}return out;
}
double median(Vector v) {
    require(!v.empty(), "Median of empty array");
    const auto n=v.size(), k=n/2;
    std::nth_element(v.begin(),v.begin()+k,v.end());
    const double hi=v[k];
    return n%2 ? hi : (hi+*std::max_element(v.begin(),v.begin()+k))/2;
}
Spectrum restrict_spectrum(Spectrum s, std::optional<Range> range) {
    if(range) require(finite(range->first)&&finite(range->second)&&range->first<range->second,"Invalid m/z range");
    s.erase(std::remove_if(s.begin(),s.end(),[](Point p){return !finite(p.x)||!finite(p.y)||p.x<=0;}),s.end());
    for(auto p:s) require(p.y>=0,"Spectrum intensities must not be negative");
    std::stable_sort(s.begin(),s.end(),[](Point a,Point b){return a.x<b.x;});
    if(range) s.erase(std::remove_if(s.begin(),s.end(),[&](Point p){return p.x<range->first||p.x>range->second;}),s.end());
    return s;
}
Spectrum engine_unique(const Spectrum& s) {
    struct Bin {double x=0,y=0; std::size_t n=0;};
    std::map<float,Bin> bins;
    for(auto p:s) {
        require(finite(p.x)&&finite(p.y)&&finite(static_cast<float>(p.x)),"Nonfinite engine input");
        auto& b=bins[static_cast<float>(p.x)]; b.x+=p.x; b.y+=p.y; ++b.n;
    }
    Spectrum out; out.reserve(bins.size());
    for(const auto& [key,b]:bins) out.push_back({b.x/b.n,b.y});
    return out;
}
Spectrum read_xy(const std::filesystem::path& path) {
    std::ifstream f(path); require(f.good(),"Cannot open spectrum");
    f.imbue(std::locale::classic());
    Spectrum s; std::string line; std::size_t row=0;
    while(std::getline(f,line)) {
        ++row; if(line.empty()||line[0]=='#') continue;
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream in(line); in.imbue(std::locale::classic()); Point p; std::string extra;
        if(!(in>>p.x>>p.y)||in>>extra) throw std::invalid_argument("Expected two numeric columns at line "+std::to_string(row));
        s.push_back(p);
    }
    return restrict_spectrum(std::move(s));
}
void write_xy(const std::filesystem::path& path,const Spectrum& s) {
    std::ofstream f(path); require(f.good(),"Cannot create spectrum output");
    f.imbue(std::locale::classic()); f<<std::setprecision(17);
    for(auto p:s) f<<p.x<<'\t'<<p.y<<'\n';
    require(f.good(),"Spectrum write failed");
}
std::vector<Peak> pick_mass_peaks(const Spectrum& s,double window,double threshold) {
    require(finite(window)&&window>=0&&finite(threshold)&&threshold>=0&&threshold<=1,"Invalid peak parameters");
    if(s.size()<3) return {};
    Vector delta; double top=0;
    for(std::size_t i=0;i<s.size();++i) { require(finite(s[i].x)&&finite(s[i].y),"Nonfinite mass data"); top=std::max(top,s[i].y); if(i) delta.push_back(s[i].x-s[i-1].x); }
    if(top<=0) return {};
    double dx=median(delta); require(dx>0,"Mass coordinates must increase");
    // Half-window integration must cover the same samples as the Python
    // reference even when the peak window exceeds the displayed mass range.
    // 2*N is sufficient to cover the whole array while bounding the index.
    auto w=static_cast<std::size_t>(std::min(2.0*double(s.size()),std::max(1.0,std::nearbyint(window/std::max(dx,1e-9)))));
    std::vector<std::size_t> order(s.size()); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return s[a].y>s[b].y;});
    std::vector<bool> taken(s.size()); std::vector<Peak> out;
    for(auto i:order) {
        if(s[i].y<=0||s[i].y<threshold*top) break;
        if(taken[i]) continue;
        auto lo=i>w?i-w:0,hi=std::min(s.size(),i+w+1);
        double ymax=0; for(auto j=lo;j<hi;++j)ymax=std::max(ymax,s[j].y);
        if(s[i].y<ymax)continue;
        std::fill(taken.begin()+lo,taken.begin()+hi,true);
        lo=i>w/2?i-w/2:0; hi=std::min(s.size(),i+w/2+1);
        double area=0,sum=0,moment=0;
        for(auto j=lo;j<hi;++j) {area+=s[j].y; if(s[j].y>=.5*s[i].y){sum+=s[j].y; moment+=s[j].x*s[j].y;}}
        out.push_back({sum>0?moment/sum:s[i].x,s[i].x,s[i].y,area*dx,100*s[i].y/top,0});
    }
    double total=0; for(auto p:out)total+=p.height;
    for(auto& p:out)p.fraction=100*p.height/total;
    std::sort(out.begin(),out.end(),[](Peak a,Peak b){return a.mass<b.mass;}); return out;
}
static Vector extremum(const Vector& y,std::size_t width,bool minimum) {
    Vector out(y.size()); std::deque<std::pair<std::int64_t,double>> q;
    auto left=static_cast<std::int64_t>(width/2),right=static_cast<std::int64_t>(width)-left-1;
    for(std::int64_t j=-left;j<static_cast<std::int64_t>(y.size())+right;++j) {
        double value=y[std::clamp<std::int64_t>(j,0,static_cast<std::int64_t>(y.size())-1)];
        while(!q.empty()&&(minimum?q.back().second>=value:q.back().second<=value))q.pop_back();
        q.emplace_back(j,value);
        while(q.front().first<j-static_cast<std::int64_t>(width)+1)q.pop_front();
        auto i=j-right; if(i>=0)out[i]=q.front().second;
    } return out;
}
Spectrum baseline(const Spectrum& s,double width) {
    require(finite(width)&&width>0,"Baseline width must be positive"); if(s.size()<20)return s;
    Vector x,y,delta; for(auto p:s){x.push_back(p.x);y.push_back(p.y);}
    for(std::size_t i=1;i<x.size();++i)delta.push_back(x[i]-x[i-1]);
    const double dx=median(delta); require(dx>0,"Spectrum axis must increase");
    auto k=static_cast<std::size_t>(std::max(5.0,std::min(double(s.size()/4),std::nearbyint(width/std::max(dx,1e-12)))));
    auto base=extremum(extremum(y,k,true),k,false); auto kw=std::max<std::size_t>(3,k/2);
    // Prefix sums with nearest-edge extension reproduce scipy uniform_filter1d.
    auto left=kw/2,right=kw-left-1; Vector prefix(s.size()+kw,0);
    for(std::size_t i=0;i<s.size()+kw-1;++i) prefix[i+1]=prefix[i]+base[std::clamp<std::int64_t>(static_cast<std::int64_t>(i)-left,0,s.size()-1)];
    Spectrum out=s;
    for(std::size_t i=0;i<s.size();++i)out[i].y=std::max(0.0,y[i]-std::min(y[i],(prefix[i+left+right+1]-prefix[i])/kw));
    return out;
}
double r_squared(const Vector& y,const Vector& fit) {
    require(!y.empty()&&y.size()==fit.size(),"R squared requires equal nonempty vectors");
    double mean=std::accumulate(y.begin(),y.end(),0.0)/y.size(),ss=0,err=0;
    for(std::size_t i=0;i<y.size();++i){ss+=std::pow(y[i]-mean,2);err+=std::pow(y[i]-fit[i],2);}return ss>0?1-err/ss:0;
}
double pairwise(const double* p,std::size_t n){
    if(n<8){double s=-0.;for(std::size_t i=0;i<n;++i)s+=p[i];return s;}
    if(n<=128){double lane[8];for(int j=0;j<8;++j)lane[j]=p[j];std::size_t i=8;for(;i+7<n;i+=8)for(int j=0;j<8;++j)lane[j]+=p[i+j];double s=((lane[0]+lane[1])+(lane[2]+lane[3]))+((lane[4]+lane[5])+(lane[6]+lane[7]));for(;i<n;++i)s+=p[i];return s;}
    std::size_t mid=(n/2)/8*8;return pairwise(p,mid)+pairwise(p+mid,n-mid);
}

Vector gaussian(const Vector& v,double sigma) {
    require(finite(sigma)&&sigma>0&&sigma<100000,"Invalid Gaussian width");
    auto radius=static_cast<std::int64_t>(4*sigma+.5); Vector kernel(2*radius+1); double total=0;
    for(std::int64_t i=-radius;i<=radius;++i)kernel[i+radius]=std::exp((-.5/(sigma*sigma))*(i*i));
    total=pairwise(kernel.data(),kernel.size());for(auto& k:kernel)k/=total;
    Vector out(v.size(),0);
    // SciPy 1.18 NI_Correlate1D evaluates symmetric kernels centre first,
    // then adds paired samples from the outside towards the centre.
    for(std::size_t i=0;i<v.size();++i){out[i]=v[i]*kernel[radius];for(auto j=-radius;j<0;++j){auto left=static_cast<std::int64_t>(i)+j,right=static_cast<std::int64_t>(i)-j;double a=left>=0?v[left]:0,b=right<static_cast<std::int64_t>(v.size())?v[right]:0;out[i]+=(a+b)*kernel[j+radius];}}
    return out;
}
static double interpolate(const Vector& x,const Vector& y,double value) {
    if(value<=x.front())return y.front();
    if(value>=x.back())return y.back();
    auto i=static_cast<std::size_t>(std::upper_bound(x.begin(),x.end(),value)-x.begin()-1);
    return y[i]+((y[i+1]-y[i])/(x[i+1]-x[i]))*(value-x[i]);
}
Vector cell_average(const Vector& x,const Vector& y,const Vector& centres,double width) {
    require(x.size()==y.size()&&x.size()>=2&&finite(width)&&width>0,"Invalid cell averaging input");
    Vector cum(x.size(),0);for(std::size_t i=1;i<x.size();++i){require(x[i]>x[i-1],"Coordinates must increase");cum[i]=cum[i-1]+.5*(y[i]+y[i-1])*(x[i]-x[i-1]);}
    Vector out;out.reserve(centres.size());for(auto c:centres)out.push_back((interpolate(x,cum,c+width/2)-interpolate(x,cum,c-width/2))/width);return out;
}
Vector grid(const std::vector<Range>& windows,double d,std::size_t limit) {
    require(finite(d)&&d>0&&limit>0,"Invalid grid spacing or limit"); Vector out;
    std::size_t count=0;for(auto [a,b]:windows){require(finite(a)&&finite(b)&&b>=a,"Invalid grid interval");double n=std::ceil((b-a)/d+.5);require(n<=double(limit-count),"Grid exceeds limit");count+=static_cast<std::size_t>(n);}
    out.reserve(count);for(auto [a,b]:windows){auto n=static_cast<std::size_t>(std::ceil((b-a)/d+.5));double delta=(a+d)-a;for(std::size_t i=0;i<n;++i)out.push_back(a+i*delta);}return out;
}
Spectrum sum_profiles(const std::vector<Spectrum>& spectra) {
    const Spectrum* first=nullptr;bool same=true;
    for(const auto& s:spectra){if(s.empty())continue;if(!first){first=&s;continue;}if(s.size()!=first->size()){same=false;continue;}for(std::size_t i=0;i<s.size();++i)if(s[i].x!=(*first)[i].x){same=false;break;}}
    if(!first)return {};
    if(same){Spectrum out=*first;for(auto& p:out)p.y=0;for(const auto& s:spectra)for(std::size_t i=0;i<s.size();++i)out[i].y+=s[i].y;return out;}
    Vector diff;double last=0;for(auto p:*first)if(p.x>0){if(last>0&&p.x>last)diff.push_back(std::log(p.x)-std::log(last));last=p.x;}
    double spacing=diff.empty()?1e-6:median(diff);
    struct Bin{double moment=0,weight=0,intensity=0;};std::map<std::int64_t,Bin> bins;
    for(const auto& s:spectra)for(auto p:s)if(p.x>0){double key=std::nearbyint(std::log(p.x)/spacing);require(finite(key)&&std::abs(key)<9e18,"Profile bin overflow");auto& b=bins[static_cast<std::int64_t>(key)];double w=std::max(p.y,1e-30);b.moment+=p.x*w;b.weight+=w;b.intensity+=p.y;}
    Spectrum out;for(const auto& [key,b]:bins)out.push_back({b.moment/b.weight,b.intensity});return out;
}
double integrate(const Spectrum& s,double b0,double b1) {
    if(s.size()<2)return 0;
    require(s.back().x>s.front().x,"Time must increase");double sum=0,prev=0;
    for(std::size_t i=0;i<s.size();++i){double b=b0+(b1-b0)*(s[i].x-s.front().x)/(s.back().x-s.front().x);double y=std::max(s[i].y-b,0.0);if(i){require(s[i].x>=s[i-1].x,"Time must increase");sum+=.5*(prev+y)*(s[i].x-s[i-1].x)*60;}prev=y;}return sum;
}
Vector decode_pda(const std::vector<std::uint8_t>& raw,std::size_t records,std::size_t wavelengths) {
    require(wavelengths>0&&records<=100000000/wavelengths,"PDA matrix exceeds limit");Vector out;out.reserve(records*wavelengths);std::size_t off=0;
    auto u16=[&](std::size_t i){require(i+2<=raw.size(),"Truncated PDA block");return unsigned(raw[i])|(unsigned(raw[i+1])<<8);};
    auto u32=[&](std::size_t i){return std::uint32_t(u16(i))|(std::uint32_t(u16(i+2))<<16);};
    for(std::size_t r=0;r<records;++r){require(off+24<=raw.size()&&raw[off]=='R'&&raw[off+1]=='C',"Invalid PDA record");auto len=u32(off+12);require(len>=24&&len<=raw.size()-off,"Invalid PDA record length");auto end=off+len,pos=off+24;std::size_t count=0;
        while(pos+2<=end&&count<wavelengths){auto n=u16(pos);auto i=pos+2,stop=i+n;require(stop+2<=end,"Truncated PDA block");std::int64_t acc=0;bool first=true;
            while(i<stop){auto pre=raw[i]>>5;auto bytes=std::min<unsigned>(pre+1,4);require(i+bytes<=stop,"Truncated PDA value");std::int64_t value=raw[i]&31;for(unsigned j=1;j<bytes;++j)value=(value<<8)|raw[i+j];unsigned bits=5+8*(bytes-1);if(value&(std::int64_t(1)<<(bits-1)))value-=(std::int64_t(1)<<bits);acc=first?value:acc+value;first=false;require(count<wavelengths,"Excess PDA wavelengths");out.push_back(static_cast<float>(acc));++count;i+=bytes;}
            require(u16(stop)==n,"PDA block footer mismatch");pos=stop+2;
        }require(count==wavelengths,"Missing PDA wavelengths");off=end;
    }return out;
}
ChargeModel::ChargeModel(Vector coordinate,double d,const Vector& masses,const Vector& charges,double carrier,double fwhm,bool log):coordinate_(std::move(coordinate)),masses_(masses.size()) {
    require(coordinate_.size()>=2&&finite(d)&&d>0&&finite(fwhm)&&fwhm>0&&finite(carrier),"Invalid charge model");
    require(masses.size()<=400000&&charges.size()<=1000&&double(masses.size())*charges.size()<=60000000,"Charge model exceeds limit");
    require(std::is_sorted(coordinate_.begin(),coordinate_.end()),"Charge model grid must increase");
    if(carrier==1||carrier==-1)carrier*=proton;
    sigma_=std::max(.3,fwhm/d/2.3548);
    for(auto z:charges){require(finite(z)&&z>=1&&z==std::floor(z),"Invalid charge");index_.emplace_back();w0_.emplace_back();w1_.emplace_back();
        for(auto m:masses){double pos=(m+z*carrier)/z;if(log)pos=std::log(std::max(pos,1e-9));auto idx=std::upper_bound(coordinate_.begin(),coordinate_.end(),pos)-coordinate_.begin()-1;bool ok=idx>=0&&idx<static_cast<std::ptrdiff_t>(coordinate_.size())-1;if(!ok)idx=0;ok=ok&&(coordinate_[idx+1]-coordinate_[idx]<1.5*d);double w=ok?std::clamp((pos-coordinate_[idx])/d,0.0,1.0):0;index_.back().push_back(idx);w0_.back().push_back(ok?1-w:0);w1_.back().push_back(w);}
    }
}
Vector ChargeModel::scatter(std::size_t k,const Vector& f) const {
    require(k<index_.size()&&f.size()==masses_,"Charge scatter dimensions differ");Vector out(coordinate_.size(),0),right(coordinate_.size(),0);for(std::size_t i=0;i<f.size();++i){out[index_[k][i]]+=w0_[k][i]*f[i];right[index_[k][i]+1]+=w1_[k][i]*f[i];}for(std::size_t i=0;i<out.size();++i)out[i]+=right[i];return out;
}
Vector ChargeModel::forward(const Vector& f,const Vector& c) const {
    require(c.size()==index_.size(),"Charge envelope dimensions differ");Vector out(coordinate_.size(),0);for(std::size_t k=0;k<c.size();++k)if(c[k]>0){auto s=scatter(k,f);for(std::size_t i=0;i<out.size();++i)out[i]+=c[k]*s[i];}return gaussian(out,sigma_);
}
Vector ChargeModel::adjoint(const Vector& r,const Vector& c) const {
    require(r.size()==coordinate_.size()&&c.size()==index_.size(),"Charge adjoint dimensions differ");auto gp=gaussian(r,sigma_);Vector out(masses_,0);for(std::size_t k=0;k<c.size();++k)if(c[k]>0)for(std::size_t i=0;i<masses_;++i)out[i]+=c[k]*(w0_[k][i]*gp[index_[k][i]]+w1_[k][i]*gp[index_[k][i]+1]);return out;
}
std::vector<Vector> ChargeModel::columns(const Vector& f) const {std::vector<Vector> out;for(std::size_t k=0;k<index_.size();++k)out.push_back(gaussian(scatter(k,f),sigma_));return out;}
}
