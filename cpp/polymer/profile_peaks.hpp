#pragma once
// Profile detection uses a short positive filter only on sufficiently sampled peaks.
// Positions, heights and areas are measured from the original profile, never the filtered signal.
namespace charge_profile {
inline std::vector<DirectPeak> calculate(const double* x,const double* y,int n,const DirectParams& p) {
    if(!p.profile)return poly_direct::calculate(x,y,n,p);
    if(!x || !y || n<1 || n>5000000 || !std::isfinite(p.mz_min) || !std::isfinite(p.mz_max) ||
       p.mz_min<0 || p.mz_max<=p.mz_min || !std::isfinite(p.threshold_percent) ||
       p.threshold_percent<=0 || p.threshold_percent>100 || !std::isfinite(p.isotope_tolerance) ||
       p.isotope_tolerance<1e-8 || p.isotope_tolerance>=.25 || p.mass_mode<0 || p.mass_mode>2)
        throw std::runtime_error("Check the range, peak threshold and isotope tolerance");
    for(int i=0;i<n;++i)if(!std::isfinite(x[i]) || !std::isfinite(y[i]) || x[i]<=0 || (i && x[i]<=x[i-1]))
        throw std::runtime_error("Spectrum m/z values must increase; all values must be finite");
    int begin=int(std::lower_bound(x,x+n,p.mz_min)-x),end=int(std::upper_bound(x,x+n,p.mz_max)-x);
    if(begin==end)throw std::runtime_error("No spectrum points in this m/z range");
    double top=*std::max_element(y+begin,y+end);
    if(top<=0)throw std::runtime_error("No positive signal in this m/z range");
    std::vector<int> apices;
    for(int i=begin+1;i<end-1;++i)if(y[i]>=top*.15 && y[i]>y[i-1] && y[i]>=y[i+1])apices.push_back(i);
    std::sort(apices.begin(),apices.end(),[&](int a,int b){return y[a]>y[b];});
    std::vector<int> widths;
    for(int i:apices) {
        int a=i,b=i;
        while(a>begin && y[a-1]>y[i]*.5)--a;
        while(b+1<end && y[b+1]>y[i]*.5)++b;
        widths.push_back(b-a+1);if(widths.size()==24)break;
    }
    if(widths.empty())return poly_direct::calculate(x,y,n,p);
    std::sort(widths.begin(),widths.end());
    if(widths[widths.size()/2]<4)return poly_direct::calculate(x,y,n,p);
    std::vector<double> smooth(y,y+n);
    constexpr double kernel[5]={1,4,6,4,1};
    for(int i=begin+2;i+2<end;++i) {
        // Do not interpolate across gaps in a sparse profile axis.
        double dx=(x[i+1]-x[i-1])*.5;
        if(x[i+2]-x[i-2]>6*dx)continue;
        double value=0;for(int k=-2;k<=2;++k)value+=kernel[k+2]*y[i+k];smooth[i]=value/16.;
    }
    const double cutoff=top*p.threshold_percent/100.;
    int resolution=widths[widths.size()/2],distance=std::max(2,int(std::round(resolution*.7)));
    int reach=std::max(6,std::min(1000,resolution*3));
    std::vector<int> maxima;
    for(int i=begin+1;i<end-1;++i) {
        if(smooth[i]<=smooth[i-1])continue;
        int plateau=i;while(plateau+1<end && smooth[plateau+1]==smooth[i])++plateau;
        if(plateau+1<end && smooth[plateau+1]<smooth[i])maxima.push_back((i+plateau)/2);
        i=plateau;
    }
    std::sort(maxima.begin(),maxima.end(),[&](int a,int b){return smooth[a]>smooth[b];});
    std::vector<char> occupied(end-begin,0);std::vector<int> kept;
    for(int i:maxima) {
        if(occupied[i-begin])continue;
        kept.push_back(i);
        for(int j=std::max(begin,i-distance);j<=std::min(end-1,i+distance);++j)occupied[j-begin]=1;
    }
    std::sort(kept.begin(),kept.end());
    std::vector<DirectPeak> result;
    for(size_t k=0;k<kept.size();++k) {
        int i=kept[k],plateau=i,left=i,right=i;
        // Topographic prominence spans small fluctuations instead of stopping at the
        // first valley on a jagged flank. The window follows observed profile width.
        for(int j=i-1;j>=std::max(begin,i-reach);--j) {
            if(smooth[j]>smooth[i])break;
            if(smooth[j]<smooth[left])left=j;
        }
        for(int j=i+1;j<=std::min(end-1,i+reach);++j) {
            if(smooth[j]>smooth[i])break;
            if(smooth[j]<smooth[right])right=j;
        }
        if(left==i || right==i || (left==begin && smooth[left]>smooth[i]*.5) ||
           (right==end-1 && smooth[right]>smooth[i]*.5) ||
           smooth[i]-std::max(smooth[left],smooth[right])<cutoff)continue;
        // Quantification regions do not overlap another resolved maximum.
        if(k>0 && kept[k-1]>left)left=int(std::min_element(smooth.begin()+kept[k-1],smooth.begin()+i+1)-smooth.begin());
        if(k+1<kept.size() && kept[k+1]<right)right=int(std::min_element(smooth.begin()+i,smooth.begin()+kept[k+1]+1)-smooth.begin());
        if(left>=i || right<=i)continue;
        int apex=int(std::max_element(y+left+1,y+right)-y);
        double mz=x[apex],height=y[apex];
        int flat=apex;while(flat+1<right && y[flat+1]==y[apex])++flat;
        if(flat!=apex)mz=(x[apex]+x[flat])*.5;
        else if(apex>begin && apex+1<end) {
            double dl=x[apex]-x[apex-1],dr=x[apex+1]-x[apex];
            double sl=(y[apex]-y[apex-1])/dl,sr=(y[apex+1]-y[apex])/dr;
            double a=(sr-sl)/(dl+dr),derivative=sl+a*dl;
            if(a<0) {
                double delta=-derivative/(2*a);
                if(delta>=-dl && delta<=dr){mz+=delta;height+=derivative*delta+a*delta*delta;}
            }
        }
        double area=0,slope=(smooth[right]-smooth[left])/(x[right]-x[left]);
        for(int j=left;j<right;++j) {
            double a=std::max(0.,y[j]-(smooth[left]+slope*(x[j]-x[left])));
            double b=std::max(0.,y[j+1]-(smooth[left]+slope*(x[j+1]-x[left])));
            area+=(a+b)*.5*(x[j+1]-x[j]);
        }
        double half=(smooth[i]+std::max(smooth[left],smooth[right]))*.5;
        int a=i,b=plateau;
        while(a>left && smooth[a-1]>half)--a;
        while(b<right && smooth[b+1]>half)++b;
        double from=x[a],to=x[b];
        if(a>left && smooth[a]!=smooth[a-1])from=x[a-1]+(half-smooth[a-1])/(smooth[a]-smooth[a-1])*(x[a]-x[a-1]);
        if(b<right && smooth[b]!=smooth[b+1])to=x[b]+(smooth[b]-half)/(smooth[b]-smooth[b+1])*(x[b+1]-x[b]);
        if(area>0 && (result.empty() || mz>result.back().mz))result.push_back({mz,mz,area,height,from,to,1});
        if(result.size()>5000)throw std::runtime_error("More than 5000 species: raise the peak threshold or narrow the m/z range");
    }
    if(result.empty())throw std::runtime_error("No complete peaks found; lower the threshold or widen the range");
    return result;
}
}
