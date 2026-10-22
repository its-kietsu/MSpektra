// Resolved isotope charge assignment, with explicit uncertainty and unique mass groups.
// No composition/averagine model. Each detected measured peak has one owner.
#include "singly_charged.hpp"
#include "profile_peaks.hpp"
#include <numeric>
#include <set>

struct ChargeParams {
    double carrier, mz_min, mz_max, threshold_percent, isotope_tolerance, merge_da, merge_ppm;
    int profile, charge_min, charge_max, sign, min_isotopes, mass_mode;
};
struct ChargePeak {
    double mass, mz, area, height, first_mz, last_mz, spacing_error, response;
    int charge, isotopes, state, alternate;
    // state: 0 isotope supported; 1 ambiguous; 2 unassigned; 3 manual; 4 deconvolved; 5 outside range.
};
struct NeutralPeak { double mass, area, height; int count, representative; };
static thread_local char charge_error[256];
#define API extern "C" __declspec(dllexport)
API unsigned poly_charge_abi() {return 5;}
API unsigned poly_charge_params_size() {return sizeof(ChargeParams);}
API unsigned poly_charge_peak_size() {return sizeof(ChargePeak);}
API unsigned poly_charge_neutral_size() {return sizeof(NeutralPeak);}
API const char* poly_charge_error() {return charge_error;}
namespace {
constexpr double spacing=1.00335483507;
void require(bool valid,const char* text) {if(!valid)throw std::runtime_error(text);}
void validate(const ChargeParams& p) {
    require(p.charge_min>=1 && p.charge_max>=p.charge_min && p.charge_max<=100,
            "Choose a charge range from 1 to 100");
    require(p.sign==1 || p.sign==-1,"Choose positive or negative ions");
    require(std::isfinite(p.carrier) && p.mass_mode>=0 && p.mass_mode<=2 && (p.profile==0 || p.profile==1),
            "Check the carrier, spectrum and mass basis");
    require(p.min_isotopes>=2 && p.min_isotopes<=20,"Choose 2 to 20 isotope peaks");
    require(std::isfinite(p.merge_da) && p.merge_da>=0 && p.merge_da<=1 &&
            std::isfinite(p.merge_ppm) && p.merge_ppm>=0 && p.merge_ppm<=1000,
            "Check the mass matching tolerance");
}
struct Candidate {std::vector<int> indices;int z;double error,signal;double score=0;int seeds=0;};
// Peaks of a profile need this signal to noise ratio (height above the local valleys / noise) to start a group,
// and SUPPORT_SNR to be one of its isotope peaks; noise bumps otherwise extend groups and form false ladders.
constexpr double SEED_SNR=5.,SUPPORT_SNR=4.;
// Noise of the profile at each peak of the support list (0: centroid data). Set by support_peaks.
static thread_local std::vector<double> peak_noise;
double noise_at(size_t k){return k<peak_noise.size()?peak_noise[k]:0.;}
// The isotope pattern of one species is log-concave, so its heights rise to one maximum and then fall. Two
// species side by side, or a species and noise between its peaks, break this; heights may differ from that
// shape by a quarter of the larger one plus twice the noise.
bool unimodal(const std::vector<DirectPeak>& peaks,const std::vector<int>& indices) {
    if(indices.size()<3)return true;
    size_t top=0;
    for(size_t k=1;k<indices.size();++k)if(peaks[indices[k]].height>peaks[indices[top]].height)top=k;
    auto slack=[&](size_t a,size_t b){
        return .25*std::max(peaks[indices[a]].height,peaks[indices[b]].height)+
               2.*std::max(noise_at(indices[a]),noise_at(indices[b]));};
    // before the maximum no peak lies below an earlier one, after it none above a peak between (so small
    // steps cannot add up to a second maximum)
    for(size_t k=1,high=0;k<=top;++k) {
        if(peaks[indices[k]].height<peaks[indices[high]].height-slack(k,high))return false;
        if(peaks[indices[k]].height>peaks[indices[high]].height)high=k;
    }
    for(size_t k=top+1,low=top;k<indices.size();++k) {
        if(peaks[indices[k]].height>peaks[indices[low]].height+slack(k,low))return false;
        if(peaks[indices[k]].height<peaks[indices[low]].height)low=k;
    }
    // log-concave: each height at least the geometric mean of its neighbours' (a fall and then a rise, or a
    // steep fall after a slow one, is two species); within half of the square plus twice the noise
    for(size_t k=1;k+1<indices.size();++k) {
        const double n=2.*noise_at(indices[k]);
        const double a=std::max(0.,peaks[indices[k-1]].height-n),b=std::max(0.,peaks[indices[k+1]].height-n);
        const double h=peaks[indices[k]].height+n;
        if(a*b>1.5*h*h)return false;
    }
    return true;
}
int nearest(const std::vector<DirectPeak>& peaks,double value,double tolerance) {
    auto it=std::lower_bound(peaks.begin(),peaks.end(),value-tolerance,
                            [](const DirectPeak& a,double b){return a.mz<b;});
    int best=-1;double d=tolerance;
    for(;it!=peaks.end() && it->mz<=value+tolerance;++it) {
        double delta=std::abs(it->mz-value);
        if(best<0 || delta<d){best=int(it-peaks.begin());d=delta;}
    }
    return best;
}
int overlap(const Candidate& a,const Candidate& b) {
    int n=0;size_t i=0,j=0;
    while(i<a.indices.size() && j<b.indices.size()) {
        if(a.indices[i]==b.indices[j]){++n;++i;++j;}
        else if(a.indices[i]<b.indices[j])++i;else ++j;
    }
    return n;
}
bool harmonic_subset(const Candidate& large,const Candidate& small,int shared) {
    // Every-other isotope candidates of a denser, longer envelope are harmonics.
    return shared==int(small.indices.size()) &&
        large.indices.size()>=2*small.indices.size()-1 && large.z>small.z && large.z%small.z==0;
}
double evidence_score(const std::vector<int>& indices,const std::vector<DirectPeak>& peaks,
                      const std::vector<char>& seeds,int gaps,double error,double tolerance) {
    double main_area=0;int main=0;
    for(int i:indices)if(seeds[i]){++main;main_area+=peaks[i].area;}
    double weak=0,reference=main?main_area/main:1;
    for(int i:indices)if(!seeds[i])weak+=std::min(1.,peaks[i].area/reference);
    return 4.*main+weak-1.5*gaps-.5*error/tolerance;
}
bool profile_support(const DirectPeak& peak,double reference_area,double reference_width) {
    double width=peak.last_mz-peak.first_mz;
    if(width<=0 || reference_width<=0)return true; // centroid or undersampled profile
    if(width<reference_width*.35 || width>reference_width*3.)return false;
    // Small clean peaks remain eligible. Small background-dominated bumps need
    // enough measured area relative to the seed to count as isotope evidence.
    return peak.area>=reference_area*.05 || peak.area>=.8*peak.height*width;
}
std::vector<DirectPeak> support_peaks(const double* x,const double* y,int n,const ChargeParams& p,
                                     std::vector<char>& seeds) {
    DirectParams pick{0,p.mz_min,p.mz_max,p.threshold_percent,p.isotope_tolerance,p.profile,0,p.mass_mode};
    std::vector<double> strong_noise,strong_snr,noise,snr;
    const auto strong=charge_profile::calculate(x,y,n,pick,&strong_noise,&strong_snr);
    auto peaks=strong;noise=strong_noise;snr=strong_snr;
    // A group needs one strong peak, not three peaks each above a global cutoff.
    // Weak isotope support is picked separately; excessive noise raises only this support cutoff.
    double threshold=p.threshold_percent*(p.profile?.05:.2);
    while(threshold<p.threshold_percent*.999) {
        try {pick.threshold_percent=threshold;peaks=charge_profile::calculate(x,y,n,pick,&noise,&snr);break;}
        catch(const std::runtime_error& e) {
            if(std::string(e.what()).find("More than 5000")==std::string::npos)throw;
            threshold*=1.8;
        }
    }
    std::vector<char> strong_peak(peaks.size(),0);
    for(size_t i=0;i<strong.size();++i) {
        int k=nearest(peaks,strong[i].mz,1e-8);
        if(k>=0 && strong_snr[i]>=SEED_SNR)strong_peak[k]=1;
    }
    // Support peaks above the noise only (all peaks of centroid data, whose noise cannot be told).
    std::vector<DirectPeak> kept;seeds.clear();peak_noise.clear();
    for(size_t k=0;k<peaks.size();++k) {
        if(!strong_peak[k] && snr[k]<SUPPORT_SNR)continue;
        kept.push_back(peaks[k]);seeds.push_back(strong_peak[k]);peak_noise.push_back(noise[k]);
    }
    return kept;
}
Candidate sequence(const std::vector<DirectPeak>& peaks,const std::vector<char>& seeds,
                   int anchor,int z,const ChargeParams& p,bool allow_gaps=false) {
    const double step=spacing/z,tol=std::min(p.isotope_tolerance,step*.10);
    double width=peaks[anchor].last_mz-peaks[anchor].first_mz;
    if(p.profile && width>0 && step<width*.85)return Candidate{{},z,0,0,0};
    std::map<int,int> positions{{0,anchor}};
    for(int direction:{-1,1}) {
        int missed=0;
        for(int k=1;k<=80 && missed<2;++k) {
            const int order=k*direction;
            double target=peaks[anchor].mz+order*step;
            auto start=std::lower_bound(peaks.begin(),peaks.end(),target-2*tol,
                                       [](const DirectPeak& a,double b){return a.mz<b;});
            int index=-1;double distance=2*tol;
            for(auto it=start;it!=peaks.end() && it->mz<=target+2*tol;++it) {
                int k=int(it-peaks.begin());
                if(p.profile && !seeds[k] && !profile_support(*it,peaks[anchor].area,peaks[anchor].last_mz-peaks[anchor].first_mz))continue;
                double delta=std::abs(it->mz-target);
                if(index<0 || delta<distance){index=k;distance=delta;}
            }
            if(index<0 || (direction<0?index>=anchor:index<=anchor)){++missed;continue;}
            positions[order]=index;missed=0;
        }
    }
    Candidate c{{},z,0,0,0};
    std::vector<std::pair<int,int>> ordered(positions.begin(),positions.end());
    int count=int(ordered.size()),center=0;
    std::vector<long double> w(count+1),v(count+1),v2(count+1);
    std::vector<double> signal(count+1);
    std::vector<int> strong(count+1),adjacent(count+1);
    for(int i=0;i<count;++i) {
        auto [order,index]=ordered[i];if(order==0)center=i;
        long double weight=std::sqrt(peaks[index].height);
        // Anchor-relative offsets avoid subtracting squared large m/z values.
        long double offset=peaks[index].mz-peaks[anchor].mz-order*step;
        w[i+1]=w[i]+weight;v[i+1]=v[i]+weight*offset;v2[i+1]=v2[i]+weight*offset*offset;
        signal[i+1]=signal[i]+peaks[index].height;strong[i+1]=strong[i]+(seeds[index]?1:0);
        adjacent[i+1]=adjacent[i]+(i>0 && order-ordered[i-1].first==1?1:0);
    }
    // An unrelated weak match outside a missing isotope must not invalidate a sound core.
    // Find the best supported interval containing the seed instead of accepting/rejecting
    // every possible tail match together. Missing internal isotopes remain allowed.
    for(int a=0;a<=center;++a)for(int b=std::max(center,a+p.min_isotopes-1);b<count;++b) {
        int size=b-a+1,edges=adjacent[b+1]-adjacent[a+1],run=1,longest=1;
        for(int k=a+1;k<=b;++k){run=ordered[k].first-ordered[k-1].first==1?run+1:1;longest=std::max(longest,run);}
        int gaps=ordered[b].first-ordered[a].first+1-size;
        if(!allow_gaps && gaps)continue;
        if(edges<p.min_isotopes-1 || longest<p.min_isotopes || gaps>std::max(1,size/3) ||
           ordered[a+1].first-ordered[a].first!=1 || ordered[b].first-ordered[b-1].first!=1)continue;
        // Fixed-lattice RMS alone can hide a single implausible jump, especially
        // among weak tails. Use the same local spacing bound as window joining.
        bool coherent=true;
        for(int k=a+1;k<=b;++k) {
            double delta=peaks[ordered[k].second].mz-peaks[ordered[k-1].second].mz;
            int slots=ordered[k].first-ordered[k-1].first;
            if(std::abs(delta-slots*step)>2*tol){coherent=false;break;}
        }
        if(!coherent)continue;
        long double weight=w[b+1]-w[a],mean=(v[b+1]-v[a])/weight;
        double error=std::sqrt(std::max(0.L,(v2[b+1]-v2[a])/weight-mean*mean));
        if(error>tol)continue;
        int main=strong[b+1]-strong[a];
        std::vector<int> indices;for(int k=a;k<=b;++k)indices.push_back(ordered[k].second);
        if(p.profile) {
            double area=0,width=0;int count=0;
            for(int k:indices)if(seeds[k]){area+=peaks[k].area;width+=peaks[k].last_mz-peaks[k].first_mz;++count;}
            bool valid=count>0;
            if(count)for(int k:indices)if(!seeds[k] && !profile_support(peaks[k],area/count,width/count)){valid=false;break;}
            if(!valid)continue;
        }
        double score=evidence_score(indices,peaks,seeds,gaps,error,tol);
        if(!c.indices.empty() && score<=c.score)continue;
        if(!unimodal(peaks,indices))continue;
        c.indices=std::move(indices);
        c.error=error;c.seeds=main;c.signal=signal[b+1]-signal[a];c.score=score;
    }
    return c;
}
void join_supported_windows(std::vector<Candidate>& candidates,const std::vector<DirectPeak>& peaks,
                            const std::vector<char>& seeds,const ChargeParams& p) {
    std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b){
        if(a.z!=b.z)return a.z<b.z;
        if(a.indices.front()!=b.indices.front())return a.indices.front()<b.indices.front();
        return a.score>b.score;
    });
    std::vector<Candidate> joined;
    for(auto& c:candidates) {
        bool found=false;double step=spacing/c.z,tol=std::min(p.isotope_tolerance,step*.10);
        for(auto it=joined.rbegin();it!=joined.rend();++it) {
            if(it->z!=c.z)break;
            if(it->indices.back()<c.indices.front())continue;
            if(overlap(*it,c)<std::max(2,p.min_isotopes-1))continue;
            std::vector<int> indices;
            std::set_union(it->indices.begin(),it->indices.end(),c.indices.begin(),c.indices.end(),std::back_inserter(indices));
            bool coherent=true;int gaps=0;
            for(size_t k=1;k<indices.size();++k) {
                double delta=peaks[indices[k]].mz-peaks[indices[k-1]].mz;
                int slots=int(std::round(delta/step));
                if(slots!=1 || std::abs(delta-slots*step)>2*tol){coherent=false;break;}
                gaps+=slots-1;
            }
            if(!coherent || !unimodal(peaks,indices))continue;
            it->indices=std::move(indices);it->error=std::max(it->error,c.error);it->signal=0;it->seeds=0;
            for(int k:it->indices){it->signal+=peaks[k].height;it->seeds+=seeds[k]?1:0;}
            it->score=evidence_score(it->indices,peaks,seeds,gaps,it->error,tol);
            found=true;break;
        }
        if(!found)joined.push_back(std::move(c));
    }
    candidates=std::move(joined);
}
void remove_weak_harmonics(std::vector<Candidate>& candidates,const std::vector<DirectPeak>& peaks,
                           const std::vector<char>& seeds,const ChargeParams& p) {
    // A few weak interstitial noise peaks do not establish a higher charge over
    // an otherwise supported lower-charge pattern. No peptide/composition model.
    if(!p.profile)return; // Low centroid intensities alone are not evidence of noise.
    std::vector<char> drop(candidates.size(),0);
    std::vector<std::vector<size_t>> touching(peaks.size());
    for(size_t i=0;i<candidates.size();++i)for(int k:candidates[i].indices)touching[k].push_back(i);
    for(size_t a=0;a<candidates.size();++a) {
      std::set<size_t> neighbors;
      for(int k:candidates[a].indices)for(size_t b:touching[k])if(b<a)neighbors.insert(b);
      for(size_t b:neighbors) {
        const auto& low=candidates[b];const auto& high=candidates[a];
        if(high.z<=low.z || high.z%low.z || high.indices.front()>low.indices.back() ||
           low.indices.front()>high.indices.back())continue;
        int shared=overlap(high,low);
        if(shared!=int(low.indices.size()))continue;
        double common=0,extra=0;int independent=0,clean=0;
        for(int k:high.indices) {
            if(std::binary_search(low.indices.begin(),low.indices.end(),k))common+=peaks[k].area;
            else {
                const auto& q=peaks[k];double width=q.last_mz-q.first_mz;
                extra+=q.area;independent+=seeds[k]?1:0;
                clean+=(width<=0 || q.area>=.8*q.height*width)?1:0;
            }
        }
        if(!independent && !clean && extra<.2*common)drop[a]=1;
      }
    }
    std::vector<Candidate> kept;
    for(size_t i=0;i<candidates.size();++i)if(!drop[i])kept.push_back(std::move(candidates[i]));
    candidates=std::move(kept);
}
// A strong peak half way between two strong peaks of a group: the group may be every second peak of an
// envelope of twice the charge whose shape another species hides (its peak on one of the isotope peaks).
// Returns that charge, or 0.
int half_step_charge(const std::vector<DirectPeak>& peaks,const Candidate& c,
                     const ChargeParams& p) {
    if(2*c.z>100 || c.indices.size()<2)return 0;
    const double step=spacing/c.z,tol=std::min(p.isotope_tolerance,step*.05);
    double top=0;for(int k:c.indices)top=std::max(top,peaks[k].height);
    for(size_t k=0;k+1<c.indices.size();++k) {
        const auto& a=peaks[c.indices[k]];const auto& b=peaks[c.indices[k+1]];
        if(std::min(a.height,b.height)<.1*top)continue;   // (not in the weak tail of the group)
        int mid=nearest(peaks,(a.mz+b.mz)*.5,tol);
        if(mid>=0 && peaks[mid].height>=std::max(.3*std::min(a.height,b.height),.25*top))return 2*c.z;
    }
    return 0;
}
void remove_interleaved(std::vector<Candidate>& candidates) {
    // Ladders of two lower charges side by side (e.g. a +2 envelope and a +1 envelope a quarter of a unit
    // apart) also sit on the lattice of a higher charge. A higher charge group whose peaks all belong to lower
    // charge groups that go on beyond it is such a pair, not an ion of that charge. The lower groups of a true
    // higher charge envelope (every second peak of a +4 forms a +2 ladder) stay inside it.
    std::map<int,std::vector<size_t>> touching;
    for(size_t i=0;i<candidates.size();++i)for(int k:candidates[i].indices)touching[k].push_back(i);
    std::vector<char> drop(candidates.size(),0);
    for(size_t h=0;h<candidates.size();++h) {
        const auto& high=candidates[h];
        if(high.z<2)continue;
        std::set<int> covered;
        for(int k:high.indices)for(size_t l:touching[k]) {
            const auto& low=candidates[l];
            if(low.z>=high.z)continue;
            // (a lower group that goes on for at least two peaks: one weak peak of another species after a
            // harmonic of a true envelope does not make it a separate ion)
            int outside=0;
            for(int q:low.indices)outside+=q<high.indices.front() || q>high.indices.back();
            if(outside>=2)covered.insert(k);
        }
        if(covered.size()==high.indices.size())drop[h]=1;
    }
    std::vector<Candidate> kept;
    for(size_t i=0;i<candidates.size();++i)if(!drop[i])kept.push_back(std::move(candidates[i]));
    candidates=std::move(kept);
}
ChargePeak summarize(const std::vector<DirectPeak>& peaks,const std::vector<int>& indices,
                     int z,int state,int alternate,const ChargeParams& p,double error=0,bool strict=true) {
    require(!indices.empty(),"Empty isotope group");
    int strongest=indices[0];long double area=0,weighted=0;
    for(int i:indices) {
        area+=peaks[i].area;weighted+=(long double)peaks[i].area*peaks[i].mz;
        if(peaks[i].height>peaks[strongest].height)strongest=i;
    }
    const double mz=p.mass_mode==0?peaks[indices[0]].mz:p.mass_mode==1?peaks[strongest].mz:double(weighted/area);
    double mass=z?z*(mz-p.carrier):0;
    if(!strict && z && std::isfinite(mass) && mass<=0) {
        // (a peak below the carrier mass: listed without a charge; the analysis goes on)
        mass=0;z=0;state=2;alternate=0;
    }
    require(std::isfinite(double(area)) && std::isfinite(mass) && (!z || mass>0),
            "Check the carrier mass and spectrum range");
    return {mass,mz,double(area),peaks[strongest].height,peaks[indices.front()].mz,
            peaks[indices.back()].mz,error,1.,z*p.sign,int(indices.size()),state,alternate*p.sign};
}
template<class F> int guard(F body) {
    charge_error[0]=0;
    try{return body();}
    catch(const std::exception& e){std::strncpy(charge_error,e.what(),255);charge_error[255]=0;return -1;}
    catch(...){std::strcpy(charge_error,"Polymer charge analysis failed");return -1;}
}
}
API int poly_charge_find(const double* x,const double* y,int n,const ChargeParams* params,
                         ChargePeak* output,int capacity,DirectMember* links,int link_capacity,int* link_count) {
 return guard([&] {
    require(params && output && links && link_count && capacity>=1 && capacity<=5000 &&
            link_capacity>=1 && link_capacity<=200000,"Invalid charge analysis output buffer");
    *link_count=0;const auto& p=*params;validate(p);
    std::vector<char> seeds;
    auto peaks=support_peaks(x,y,n,p,seeds);
    std::vector<Candidate> candidates;
    // Detection is independent of the output filter. Check supported competing charges.
    // Selecting +1 alone must not relabel every second isotope of a +2 envelope as +1.
    std::set<std::pair<int,std::vector<int>>> seen;
    for(int z=1;z<=100;++z) {
        for(size_t i=0;i<peaks.size();++i) {
            if(!seeds[i])continue;
            auto c=sequence(peaks,seeds,int(i),z,p);
            if(!c.indices.empty() && seen.emplace(z,c.indices).second)candidates.push_back(std::move(c));
        }
    }
    join_supported_windows(candidates,peaks,seeds,p);
    remove_weak_harmonics(candidates,peaks,seeds,p);
    remove_interleaved(candidates);
    // A charge above the range is looked for only so that its envelope is not read as a lower charge; a short
    // group of one strong peak and weak neighbours there is no such envelope and must not take the peak.
    candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[&](const Candidate& c){
        if(c.z<=p.charge_max)return false;
        int strong=0;for(int k:c.indices)strong+=seeds[k]?1:0;
        return int(c.indices.size())<std::max(p.min_isotopes,4) || strong<2;}),candidates.end());
    std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b){
        if(a.score!=b.score)return a.score>b.score;
        if(a.indices.size()!=b.indices.size())return a.indices.size()>b.indices.size();
        if(a.error!=b.error)return a.error<b.error;
        if(a.signal!=b.signal)return a.signal>b.signal;
        return a.z<b.z;
    });
    // Indexed overlap lookup avoids quadratic work over all unrelated envelopes.
    std::vector<std::vector<int>> touching(peaks.size());
    for(size_t i=0;i<candidates.size();++i)for(int k:candidates[i].indices)touching[k].push_back(int(i));
    std::vector<char> used(peaks.size(),0);
    std::vector<ChargePeak> rows;std::vector<std::vector<int>> members;
    for(size_t ci=0;ci<candidates.size();++ci) {
        Candidate c=candidates[ci];
        auto total=[&](const std::vector<int>& v){double t=0;for(int k:v)t+=peaks[k].height;return t;};
        // Removes the peaks drop from the ends of a group; false when one is inside it or too little is left.
        auto trim_ends=[&](Candidate& g,const std::set<int>& drop) {
            size_t a=0,b=g.indices.size();
            while(a<b && drop.count(g.indices[a]))++a;
            while(b>a && drop.count(g.indices[b-1]))--b;
            bool inside=false,seed=false;
            for(size_t k=a;k<b;++k){inside|=drop.count(g.indices[k])>0;seed|=seeds[g.indices[k]]!=0;}
            if(inside || !seed || int(b-a)<p.min_isotopes)return false;
            g.indices=std::vector<int>(g.indices.begin()+a,g.indices.begin()+b);
            return true;
        };
        auto shared_of=[&](const Candidate& g,const Candidate& h) {
            std::vector<int> out;
            std::set_intersection(g.indices.begin(),g.indices.end(),h.indices.begin(),h.indices.end(),std::back_inserter(out));
            return out;
        };
        {
            // Peaks taken by a stronger group can only be end peaks of this one.
            std::set<int> taken;
            for(int k:c.indices)if(used[k])taken.insert(k);
            if(!taken.empty() && !trim_ends(c,taken))continue;
        }
        {
            // Weak end peaks of this group that are a larger part of a group of another charge belong to that
            // group (the last isotope peaks of a strong +1 envelope running into a weak +2 envelope).
            std::set<int> near,yield;
            for(int k:c.indices)near.insert(touching[k].begin(),touching[k].end());
            for(int j:near) {
                const auto& other=candidates[j];
                if(size_t(j)==ci || other.z==c.z ||
                   std::any_of(other.indices.begin(),other.indices.end(),[&](int k){return used[k];}))continue;
                auto shared=shared_of(c,other);
                if(shared.empty())continue;
                const double own=total(shared)/total(c.indices),theirs=total(shared)/total(other.indices);
                if(own<.25 && theirs>own)yield.insert(shared.begin(),shared.end());
            }
            Candidate trimmed=c;
            if(!yield.empty() && trim_ends(trimmed,yield))c=std::move(trimmed);
        }
        std::set<int> contenders;
        for(int k:c.indices)contenders.insert(touching[k].begin(),touching[k].end());
        int alternative=0;std::set<int> union_indices(c.indices.begin(),c.indices.end());
        for(int j:contenders) {
            if(size_t(j)==ci)continue;
            const auto& other=candidates[j];
            if(other.z==c.z || std::any_of(other.indices.begin(),other.indices.end(),[&](int k){return used[k];}))continue;
            {
                // end peaks that the other group can give up (a small part of it) are no other reading
                auto common=shared_of(c,other);
                if(common.empty())continue;
                Candidate rest=other;
                if(total(common)<.25*total(other.indices) && trim_ends(rest,std::set<int>(common.begin(),common.end())))continue;
            }
            int shared=overlap(c,other);
            int independent=0;
            for(int k:other.indices)if(seeds[k] && !std::binary_search(c.indices.begin(),c.indices.end(),k))++independent;
            if(shared>=2 && !harmonic_subset(c,other,shared) &&
               (independent>0 || c.indices==other.indices) &&
               (other.score>=c.score*.85 || independent>=2)) {
                if(!alternative)alternative=other.z;
                union_indices.insert(other.indices.begin(),other.indices.end());
            }
        }
        if(!alternative)alternative=half_step_charge(peaks,c,p);
        std::vector<int> group(union_indices.begin(),union_indices.end());
        for(int k:group)used[k]=1;
        // Ambiguous overlaps retain measured peaks but no invented neutral mass.
        const bool outside=c.z<p.charge_min || c.z>p.charge_max;
        rows.push_back(summarize(peaks,group,alternative||outside?0:c.z,alternative?1:outside?5:0,
                                 alternative?alternative:outside?c.z:0,p,c.error,false));members.push_back(std::move(group));
    }
    for(size_t i=0;i<peaks.size();++i)if(!used[i] && seeds[i]) {
        rows.push_back(summarize(peaks,{int(i)},0,2,0,p));
        members.push_back({int(i)});
    }
    require(rows.size()<=size_t(capacity),"Too many species; narrow the range or raise the threshold");
    require(peaks.size()<=size_t(link_capacity),"Too many measured isotope peaks");
    std::vector<int> order(rows.size());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){
        if(bool(rows[a].mass)!=bool(rows[b].mass))return rows[a].mass>0;
        if(rows[a].mass!=rows[b].mass)return rows[a].mass<rows[b].mass;
        return rows[a].mz<rows[b].mz;
    });
    int total=0;
    for(size_t i=0;i<order.size();++i) {
        int old=order[i];output[i]=rows[old];
        for(size_t k=0;k<members[old].size();++k){const auto& q=peaks[members[old][k]];
            links[total++]={int(i),int(k),q.mz,q.area,q.height};}
    }
    *link_count=total;return int(rows.size());
 });
}
API int poly_charge_merge(const ChargePeak* peaks,const int* include,int n,double da,double ppm,
                          NeutralPeak* output,int capacity,int* group) {
 return guard([&]{
    require(peaks && include && output && group && n>0 && n<=5000 && capacity>=n,
            "Invalid polymer mass grouping input");
    require(std::isfinite(da) && da>=0 && da<=1 && std::isfinite(ppm) && ppm>=0 && ppm<=1000,
            "Check the mass matching tolerance");
    std::vector<int> order;
    for(int i=0;i<n;++i){group[i]=-1;if(include[i]) {
        require(peaks[i].state==0 || peaks[i].state==3 || peaks[i].state==4,"Review uncertain charges before including them");
        require(std::isfinite(peaks[i].mass) && peaks[i].mass>0 && std::isfinite(peaks[i].area) &&
                peaks[i].area>=0 && std::isfinite(peaks[i].height) && peaks[i].height>=0 &&
                std::isfinite(peaks[i].response) && peaks[i].response>0,"Invalid polymer peak or response factor");
        order.push_back(i);
    }}
    std::sort(order.begin(),order.end(),[&](int a,int b){return peaks[a].mass<peaks[b].mass;});
    int count=0;double first=0;
    for(int i:order){const auto& p=peaks[i];
        const double tol=std::max(da,p.mass*ppm*1e-6);
        if(!count || p.mass-first>tol) {first=p.mass;output[count++]={p.mass,0,0,0,i};}
        auto& q=output[count-1];group[i]=count-1;q.area+=p.area/p.response;q.height+=p.height/p.response;++q.count;
        require(std::isfinite(q.area) && std::isfinite(q.height),"Signal exceeds the supported numerical range");
        // Report an actual member mass, never an intensity weighted invented mass.
        if(p.height>peaks[q.representative].height){q.representative=i;q.mass=p.mass;}
    }
    return count;
 });
}
API int poly_charge_manual(ChargePeak* peak,DirectMember* members,int n,int charge,
                           double carrier,int mass_mode,double tolerance) {
 return guard([&]{
    require(peak && members && n>=1 && n<=200000 && charge!=0 && std::abs(charge)<=100 &&
            std::isfinite(carrier) && std::isfinite(tolerance) && tolerance>0 && tolerance<.25 &&
            mass_mode>=0 && mass_mode<=2,"Check the manual charge settings");
    std::vector<DirectPeak> picked;
    for(int i=0;i<n;++i) {
        require(std::isfinite(members[i].mz) && members[i].mz>0 && std::isfinite(members[i].area) &&
                members[i].area>0 && std::isfinite(members[i].height) && members[i].height>0 &&
                (!i || members[i].mz>members[i-1].mz),"Invalid measured peaks");
        if(i)require(std::abs(members[i].mz-members[i-1].mz-spacing/std::abs(charge))<=
                        std::min(tolerance,spacing/std::abs(charge)*.12),
                        "The isotope spacing does not match this charge");
        picked.push_back({0,members[i].mz,members[i].area,members[i].height,0,0,1});
    }
    std::vector<int> indices(n);std::iota(indices.begin(),indices.end(),0);
    ChargeParams p{};p.carrier=carrier;p.mass_mode=mass_mode;p.sign=charge<0?-1:1;
    const double response=peak->response;
    *peak=summarize(picked,indices,std::abs(charge),3,0,p);peak->response=response;return 0;
 });
}
API int poly_charge_resolve(ChargePeak* peak,const DirectMember* members,int n,int charge,
                            double carrier,int mass_mode,double tolerance,int* chosen) {
 return guard([&]{
    require(peak && members && chosen && n>=1 && n<=5000 && charge!=0 && std::abs(charge)<=100 &&
            std::isfinite(carrier) && std::isfinite(tolerance) && tolerance>0 && tolerance<.25 &&
            mass_mode>=0 && mass_mode<=2,"Check the manual charge settings");
    std::vector<DirectPeak> picked;
    for(int i=0;i<n;++i) {
        require(std::isfinite(members[i].mz) && members[i].mz>0 && std::isfinite(members[i].area) &&
                members[i].area>0 && std::isfinite(members[i].height) && members[i].height>0 &&
                (!i || members[i].mz>members[i-1].mz),"Invalid measured peaks");
        picked.push_back({0,members[i].mz,members[i].area,members[i].height,0,0,1});chosen[i]=0;
    }
    ChargeParams p{};p.carrier=carrier;p.mass_mode=mass_mode;p.sign=charge<0?-1:1;
    p.min_isotopes=2;p.isotope_tolerance=tolerance;
    peak_noise.assign(picked.size(),0.);
    std::vector<char> seeds(picked.size(),1);std::vector<Candidate> candidates;
    for(int i=0;i<n;++i) {
        auto c=sequence(picked,seeds,i,std::abs(charge),p);
        if(!c.indices.empty())candidates.push_back(std::move(c));
    }
    join_supported_windows(candidates,picked,seeds,p);
    std::vector<int> best;double best_signal=-1;
    for(const auto& c:candidates)if(c.indices.size()>best.size() ||
        (c.indices.size()==best.size() && c.signal>best_signal)){best=c.indices;best_signal=c.signal;}
    if(n==1)best={0};
    require(!best.empty(),"The isotope spacing does not match this charge");
    const double response=peak->response;
    *peak=summarize(picked,best,std::abs(charge),3,0,p);peak->response=response;
    for(int i:best)chosen[i]=1;
    return int(best.size());
 });
}
API int poly_charge_pick(const double* x,const double* y,int n,const ChargeParams* params,double target,
                         int charge,int isotopes,ChargePeak* output,DirectMember* links,int capacity,int* count) {
 return guard([&]{
    require(params && output && links && count && capacity>=1 && capacity<=5000 &&
            std::isfinite(target) && target>0 && charge!=0 && std::abs(charge)<=100,
            "Check the selected peak and charge");
    auto p=*params;validate(p);require((charge>0?1:-1)==p.sign,"Charge sign must match the spectrum");
    require(x && y && n>0,"Choose a measured spectrum");
    const double step=spacing/std::abs(charge),span=20*step;
    DirectParams pick{0,std::max(p.mz_min,target-span),std::min(p.mz_max,target+span),.2,
                      p.isotope_tolerance,p.profile,0,p.mass_mode};
    if(pick.mz_max<=pick.mz_min){pick.mz_min=std::max(0.,target-.01);pick.mz_max=target+.01;}
    // The clicked peak, rather than an unrelated high peak nearby, sets the support cutoff.
    double local_top=0,window_top=0;
    for(int i=0;i<n;++i) {
        if(x[i]>=pick.mz_min && x[i]<=pick.mz_max)window_top=std::max(window_top,y[i]);
        if(std::abs(x[i]-target)<=std::max(.02,step*.2))local_top=std::max(local_top,y[i]);
    }
    require(local_top>0 && window_top>0,"Click closer to a measured peak");
    pick.threshold_percent=.2*local_top/window_top;
    auto peaks=charge_profile::calculate(x,y,n,pick,&peak_noise);
    int anchor=nearest(peaks,target,std::max(.02,step*.2));
    require(anchor>=0,"Click closer to a measured peak");
    // A manual click confirms the charge, not every weak background bump as an
    // isotope. Apply measured profile support to the neighbors as in auto search.
    std::vector<char> seeds(peaks.size(),0);
    for(size_t i=0;i<peaks.size();++i)
        // A weak clicked peak must not lower the evidence threshold for a much
        // stronger neighboring envelope and turn background into half-step peaks.
        // Clean weak peaks remain eligible through profile_support, not raw height.
        seeds[i]=peaks[i].height>=window_top*p.threshold_percent/100.;
    seeds[anchor]=1;
    Candidate c;
    // Manual expansion stays on the contiguous ladder containing the clicked
    // peak. A missing position must not bridge into a different weak group.
    if(isotopes)c=sequence(peaks,seeds,anchor,std::abs(charge),p,false);
    if(c.indices.empty())c.indices={anchor};
    require(c.indices.size()<=size_t(capacity),"Too many isotope peaks");
    *output=summarize(peaks,c.indices,std::abs(charge),3,0,p);
    for(size_t i=0;i<c.indices.size();++i) {
        const auto& q=peaks[c.indices[i]];links[i]={0,int(i),q.mz,q.area,q.height};
    }
    *count=int(c.indices.size());return 0;
 });
}
