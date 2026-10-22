// Direct +1/-1 mass conversion. No charge-state deconvolution or protein model.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <map>

struct DirectPeak {
    double mass, mz, area, height, first_mz, last_mz;
    int isotopes;
};
struct DirectParams {
    double carrier, mz_min, mz_max, threshold_percent, isotope_tolerance;
    int profile, group_isotopes, mass_mode; // 0:first observed, 1:strongest, 2:mean
};
struct DirectMember {
    int species, isotope;
    double mz, area, height;
};
static thread_local char direct_error[256];
#define DIRECT_API extern "C" __declspec(dllexport)
DIRECT_API unsigned poly_direct_peak_size() { return sizeof(DirectPeak); }
DIRECT_API unsigned poly_direct_params_size() { return sizeof(DirectParams); }
DIRECT_API const char* poly_direct_error() { return direct_error; }
DIRECT_API unsigned poly_direct_member_size() { return sizeof(DirectMember); }

namespace poly_direct {
struct Peak { double mz, area, height; };
inline double positive(double y) { return std::max(0., y); }

inline std::vector<DirectPeak> calculate(const double* x, const double* y, int n,
                                       const DirectParams& p,
                                       std::vector<DirectMember>* members=nullptr) {
    if (!x || !y || n < 1 || n > 5000000)
        throw std::runtime_error("Use a spectrum with 1 to 5 million points");
    if (!std::isfinite(p.carrier) || !std::isfinite(p.mz_min) ||
        !std::isfinite(p.mz_max) || p.mz_min < 0 || p.mz_max <= p.mz_min ||
        !std::isfinite(p.threshold_percent) || p.threshold_percent <= 0 ||
        p.threshold_percent > 100 || !std::isfinite(p.isotope_tolerance) ||
        p.isotope_tolerance < 1e-8 || p.isotope_tolerance >= .25 ||
        p.mass_mode < 0 || p.mass_mode > 2 ||
        (p.profile != 0 && p.profile != 1) ||
        (p.group_isotopes != 0 && p.group_isotopes != 1))
        throw std::runtime_error("Check the range, carrier, peak threshold and isotope tolerance");
    for (int i = 0; i < n; ++i)
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || x[i] <= 0 ||
            (i && x[i] <= x[i-1]))
            throw std::runtime_error("Spectrum m/z values must increase; all values must be finite");
    const int begin = int(std::lower_bound(x, x+n, p.mz_min)-x);
    const int end = int(std::upper_bound(x, x+n, p.mz_max)-x);
    if (begin == end) throw std::runtime_error("No spectrum points in this m/z range");
    double top = 0;
    for (int i = begin; i < end; ++i) top = std::max(top, y[i]);
    if (top <= 0) throw std::runtime_error("No positive signal in this m/z range");
    const double cutoff = top * p.threshold_percent / 100.;
    std::vector<Peak> peaks;
    for (int i = begin; i < end; ++i) {
        if (y[i] < cutoff) continue;
        if (!p.profile) {
            // Centroid intensities are discrete signals, never trapezoids across gaps.
            peaks.push_back({x[i], y[i], y[i]});
            continue;
        }
        // Require both flanks. A profile cut at a range boundary is not a full peak.
        if (i == begin || i == end-1 || y[i] <= y[i-1]) continue;
        int plateau = i;
        while (plateau+1 < end && y[plateau+1] == y[i]) ++plateau;
        if (plateau == end-1 || y[plateau+1] >= y[i]) { i = plateau; continue; }
        int left = i-1, right = plateau+1;
        while (left > begin && y[left-1] < y[left]) --left;
        while (right+1 < end && y[right+1] < y[right]) ++right;
        // Height above the local valleys rejects noise riding on a broad background.
        if (y[i]-std::max(y[left],y[right]) < cutoff) { i=plateau;continue; }
        double mz = (x[i]+x[plateau])*.5, height = y[i];
        if (plateau == i) {
            // Quadratic apex on the actual, potentially nonuniform profile axis.
            double dl=x[i]-x[i-1], dr=x[i+1]-x[i];
            double sl=(y[i]-y[i-1])/dl, sr=(y[i+1]-y[i])/dr;
            double a=(sr-sl)/(dl+dr), derivative=sl+a*dl;
            if (a < 0) {
                double dx=-derivative/(2*a);
                if (dx >= -dl && dx <= dr) {
                    mz=x[i]+dx; height=y[i]+derivative*dx+a*dx*dx;
                }
            }
        }
        double area=0, slope=(y[right]-y[left])/(x[right]-x[left]);
        for (int j=left; j<right; ++j) {
            double v0=positive(y[j]-(y[left]+slope*(x[j]-x[left])));
            double v1=positive(y[j+1]-(y[left]+slope*(x[j+1]-x[left])));
            area+=(v0+v1)*.5*(x[j+1]-x[j]);
        }
        if (area > 0) peaks.push_back({mz,area,height});
        i=plateau;
    }
    if (peaks.empty()) throw std::runtime_error("No complete peaks found; lower the threshold or widen the range");
    std::vector<char> used(peaks.size(),0);
    std::vector<DirectPeak> out;
    std::vector<std::pair<double, DirectMember>> pending;
    constexpr double isotope_spacing = 1.00335483507;
    for (size_t i=0; i<peaks.size(); ++i) {
        if (used[i]) continue;
        std::vector<size_t> group{i}; used[i]=1;
        if (p.group_isotopes) {
            double previous=peaks[i].mz;
            for (;;) {
                double target=previous+isotope_spacing;
                auto it=std::lower_bound(peaks.begin(),peaks.end(),target-p.isotope_tolerance,
                                        [](const Peak& q,double v){return q.mz<v;});
                size_t best=peaks.size(); double distance=p.isotope_tolerance;
                for (;it!=peaks.end() && it->mz<=target+p.isotope_tolerance;++it) {
                    size_t j=size_t(it-peaks.begin());
                    double d=std::abs(it->mz-target);
                    if (!used[j] && (best==peaks.size() || d<distance)) {best=j;distance=d;}
                }
                if (best==peaks.size()) break;
                used[best]=1;group.push_back(best);previous=peaks[best].mz;
            }
        }
        long double area=0, first=0;
        size_t strongest=i;
        for (auto j:group) {
            area+=peaks[j].area;first+=(long double)peaks[j].area*peaks[j].mz;
            if (peaks[j].height>peaks[strongest].height) strongest=j;
        }
        double mz=p.mass_mode==0?peaks[i].mz:p.mass_mode==1?peaks[strongest].mz:double(first/area);
        double mass=mz-p.carrier;
        if (!std::isfinite(mass) || mass<=0)
            throw std::runtime_error("Carrier subtraction gives a nonpositive mass; check the adduct and range");
        if (!std::isfinite(double(area))) throw std::runtime_error("Signal exceeds the supported numerical range");
        out.push_back({mass,mz,double(area),peaks[strongest].height,
                       peaks[i].mz,peaks[group.back()].mz,int(group.size())});
        if (members) {
            for (size_t k=0;k<group.size();++k) {
                const auto& q=peaks[group[k]];
                pending.push_back({peaks[i].mz,{0,int(k),q.mz,q.area,q.height}});
                if (pending.size()>200000)
                    throw std::runtime_error("More than 200000 isotope peaks: narrow the range or raise the threshold");
            }
        }
        if (out.size()>5000)
            throw std::runtime_error("More than 5000 species: raise the peak threshold or narrow the m/z range");
    }
    std::sort(out.begin(),out.end(),[](const DirectPeak& a,const DirectPeak& b){return a.mass<b.mass;});
    if (members) {
        std::map<double,int> indices;
        for (size_t i=0;i<out.size();++i) indices[out[i].first_mz]=int(i);
        for (auto& q:pending) {q.second.species=indices.at(q.first);members->push_back(q.second);}
    }
    return out;
}
}

DIRECT_API int poly_singly_charged(const double* x,const double* y,int n,
                                  const DirectParams* p,DirectPeak* out,int capacity) {
    direct_error[0]=0;
    try {
        if (!p || !out || capacity<1 || capacity>5000)
            throw std::runtime_error("Invalid direct polymer output buffer");
        auto result=poly_direct::calculate(x,y,n,*p);
        if (result.size()>size_t(capacity)) throw std::runtime_error("Peak output buffer is too small");
        std::copy(result.begin(),result.end(),out);
        return int(result.size());
    } catch (const std::exception& e) {
        std::strncpy(direct_error,e.what(),255);direct_error[255]=0;return -1;
    } catch (...) { std::strcpy(direct_error,"Direct polymer conversion failed");return -1; }
}
DIRECT_API int poly_singly_charged_linked(const double* x,const double* y,int n,
                                         const DirectParams* p,DirectPeak* out,int capacity,
                                         DirectMember* members,int member_capacity,int* member_count) {
    direct_error[0]=0;
    try {
        if (!p || !out || capacity<1 || capacity>5000 || !members || !member_count ||
            member_capacity<1 || member_capacity>200000)
            throw std::runtime_error("Invalid linked polymer output buffer");
        *member_count=0;
        std::vector<DirectMember> links;
        auto result=poly_direct::calculate(x,y,n,*p,&links);
        if (result.size()>size_t(capacity) || links.size()>size_t(member_capacity))
            throw std::runtime_error("Linked polymer output buffer is too small");
        std::copy(result.begin(),result.end(),out);
        std::copy(links.begin(),links.end(),members);
        *member_count=int(links.size());
        return int(result.size());
    } catch (const std::exception& e) {
        std::strncpy(direct_error,e.what(),255);direct_error[255]=0;return -1;
    } catch (...) { std::strcpy(direct_error,"Linked polymer conversion failed");return -1; }
}
#undef DIRECT_API
