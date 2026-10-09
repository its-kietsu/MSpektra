// Exact mass tools: port of ms_formula.py (isotope table, parse_formula,
// mono_mass, isotope_pattern, find_formulas). Element order inside a formula is
// kept as in the Python dicts (insertion order) so that sums and the isotope
// convolution run in the same order and give the same numbers.
#include "common.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <unordered_map>

namespace ms {
namespace {

const double ELECTRON = 0.00054857990946;
const long kMaxCount = 1000000;   // atoms of one element in a formula (larger counts overflowed long)

// Python's round() / std::nearbyint (half to even) without a library call: nearbyint of the
// MinGW runtime goes through the x87 control word and made the isotope pattern of a 36 kDa protein
// four times slower under Windows than under Linux. Exact for every double (|x| >= 2^52 is integral).
inline double rne(double x) {
    const double big = 4503599627370496.0;   // 2^52
    const double a = std::fabs(x);
    if (!(a < big)) return x;                 // integral already (or NaN / inf)
    const double t = a + big;                 // rounded to an integer, half to even (SSE2 doubles)
    return std::copysign(t - big, x);
}

struct Iso { double mass, p; };
struct Element { const char* sym; std::vector<Iso> iso; };

// IUPAC / NIST (AME 2012, CIAAW 2013), as listed in ms_formula.py
std::vector<Element> make_table() {
    std::vector<Element> t = {
        {"H", {{1.00782503207, 0.999885}, {2.0141017778, 0.000115}}},
        {"Li", {{6.015122795, 0.0759}, {7.01600455, 0.9241}}},
        {"B", {{10.0129370, 0.199}, {11.0093054, 0.801}}},
        {"C", {{12.0, 0.9893}, {13.0033548378, 0.0107}}},
        {"N", {{14.0030740048, 0.99636}, {15.0001088982, 0.00364}}},
        {"O", {{15.99491461956, 0.99757}, {16.99913170, 0.00038}, {17.9991610, 0.00205}}},
        {"F", {{18.99840322, 1.0}}},
        {"Na", {{22.9897692809, 1.0}}},
        {"Mg", {{23.985041700, 0.7899}, {24.98583692, 0.1000}, {25.982592929, 0.1101}}},
        {"Al", {{26.98153863, 1.0}}},
        {"Si", {{27.9769265325, 0.92223}, {28.976494700, 0.04685}, {29.97377017, 0.03092}}},
        {"P", {{30.97376163, 1.0}}},
        {"S", {{31.97207100, 0.9499}, {32.97145876, 0.0075}, {33.96786690, 0.0425}, {35.96708076, 0.0001}}},
        {"Cl", {{34.96885268, 0.7576}, {36.96590259, 0.2424}}},
        {"K", {{38.96370668, 0.932581}, {39.96399848, 0.000117}, {40.96182576, 0.067302}}},
        {"Ca", {{39.96259098, 0.96941}, {41.95861801, 0.00647}, {42.9587666, 0.00135}, {43.9554818, 0.02086}}},
        {"Mn", {{54.9380451, 1.0}}},
        {"Fe", {{53.9396105, 0.05845}, {55.9349375, 0.91754}, {56.9353940, 0.02119}, {57.9332756, 0.00282}}},
        {"Co", {{58.9331950, 1.0}}},
        {"Ni", {{57.9353429, 0.680769}, {59.9307864, 0.262231}, {60.9310560, 0.011399}, {61.9283451, 0.036345},
                {63.9279660, 0.009256}}},
        {"Cu", {{62.9295975, 0.6915}, {64.9277895, 0.3085}}},
        {"Zn", {{63.9291422, 0.4917}, {65.9260334, 0.2773}, {66.9271273, 0.0404}, {67.9248442, 0.1845},
                {69.9253193, 0.0061}}},
        {"Se", {{73.9224764, 0.0089}, {75.9192136, 0.0937}, {76.9199140, 0.0763}, {77.9173091, 0.2377},
                {79.9165213, 0.4961}, {81.9166994, 0.0873}}},
        {"Br", {{78.9183371, 0.5069}, {80.9162906, 0.4931}}},
        {"Ru", {{95.907598, 0.0554}, {97.905287, 0.0187}, {98.9059393, 0.1276}, {99.9042195, 0.1260},
                {100.9055821, 0.1706}, {101.9043493, 0.3155}, {103.905433, 0.1862}}},
        {"Rh", {{102.905504, 1.0}}},
        {"Pd", {{101.905609, 0.0102}, {103.904036, 0.1114}, {104.905085, 0.2233}, {105.903486, 0.2733},
                {107.903892, 0.2646}, {109.905153, 0.1172}}},
        {"Ag", {{106.905097, 0.51839}, {108.904752, 0.48161}}},
        {"Sn", {{111.904818, 0.0097}, {113.902779, 0.0066}, {114.903342, 0.0034}, {115.901741, 0.1454},
                {116.902952, 0.0768}, {117.901603, 0.2422}, {118.903308, 0.0859}, {119.9021947, 0.3258},
                {121.9034390, 0.0463}, {123.9052739, 0.0579}}},
        {"I", {{126.904473, 1.0}}},
        {"Cs", {{132.905451933, 1.0}}},
        {"Ir", {{190.9605940, 0.373}, {192.9629264, 0.627}}},
        {"Pt", {{189.959932, 0.00012}, {191.9610380, 0.00782}, {193.9626803, 0.3286}, {194.9647911, 0.3378},
                {195.9649515, 0.2521}, {197.9678930, 0.07356}}},
        {"Au", {{196.9665687, 1.0}}},
        {"Hg", {{195.965833, 0.0015}, {197.9667690, 0.0997}, {198.9682799, 0.1687}, {199.9683260, 0.2310},
                {200.9703023, 0.1318}, {201.9706430, 0.2986}, {203.9734939, 0.0687}}},
    };
    // most abundant first (stable, as list.sort with key -abundance)
    for (Element& e : t) std::stable_sort(e.iso.begin(), e.iso.end(), [](const Iso& a, const Iso& b) { return a.p > b.p; });
    return t;
}

const std::vector<Element>& table() {
    static const std::vector<Element> t = make_table();
    return t;
}

const Element* element(const std::string& sym) {
    for (const Element& e : table()) if (sym == e.sym) return &e;
    return nullptr;
}

double mono_of(const Element& e) { return e.iso[0].mass; }   // the most abundant isotope

// a formula: (element, count) in insertion order, as the Python dict
typedef std::vector<std::pair<const Element*, long>> Formula;

void add_to(Formula& f, const Element* e, long n) {
    for (auto& q : f) if (q.first == e) {
        if ((long long)q.second + n > kMaxCount) throw std::invalid_argument("Element count too large in the formula (at most 1000000)");
        q.second += n;
        return;
    }
    f.push_back({e, n});
}

long count_of(const Formula& f, const char* sym) {
    for (const auto& q : f) if (!strcmp(q.first->sym, sym)) return q.second;
    return 0;
}

Formula parse_one(const std::string& text);

// parse_formula: 'C6H12O6', 'Ca(OH)2', 'C2H5OH', 'CuSO4·5H2O' (a part after
// U+00B7, U+2022, U+22C5, '*' or '.' may start with a count, as ms_formula.py)
Formula parse(const std::string& text_in) {
    std::string text;
    bool parts = false;
    for (size_t i = 0; i < text_in.size(); i++) {
        const unsigned char c = (unsigned char)text_in[i];
        if (c == ' ') continue;
        if (c == '*' || c == '.') { text += '|'; parts = true; continue; }
        if (c == 0xC2 && i + 1 < text_in.size() && (unsigned char)text_in[i + 1] == 0xB7) {   // U+00B7
            text += '|'; parts = true; i++; continue;
        }
        if (c == 0xE2 && i + 2 < text_in.size() &&
            (((unsigned char)text_in[i + 1] == 0x80 && (unsigned char)text_in[i + 2] == 0xA2) ||     // U+2022
             ((unsigned char)text_in[i + 1] == 0x8B && (unsigned char)text_in[i + 2] == 0x85))) {  // U+22C5
            text += '|'; parts = true; i += 2; continue;
        }
        text += text_in[i];
    }
    if (text.empty()) throw std::invalid_argument("Empty formula");
    if (!parts) return parse_one(text);
    Formula out;
    size_t a = 0;
    for (;;) {
        const size_t b = text.find('|', a);
        const std::string part = text.substr(a, b == std::string::npos ? std::string::npos : b - a);
        size_t k = 0;
        long mult = 1;
        if (k < part.size() && isdigit((unsigned char)part[k])) {
            mult = 0;
            while (k < part.size() && isdigit((unsigned char)part[k])) mult = mult * 10 + (part[k++] - '0');
        }
        const std::string body = part.substr(k);
        if (body.empty()) throw std::invalid_argument("Cannot read the formula '" + text_in + "'");
        for (const auto& q : parse_one(body)) add_to(out, q.first, q.second * mult);
        if (b == std::string::npos) break;
        a = b + 1;
    }
    return out;
}

Formula parse_one(const std::string& text) {
    std::vector<Formula> stack(1);
    size_t pos = 0;
    while (pos < text.size()) {
        const char ch = text[pos];
        if (ch >= 'A' && ch <= 'Z') {
            std::string sym(1, ch);
            pos++;
            if (pos < text.size() && text[pos] >= 'a' && text[pos] <= 'z') sym += text[pos++];
            long n = 1;
            if (pos < text.size() && isdigit((unsigned char)text[pos])) {
                n = 0;
                while (pos < text.size() && isdigit((unsigned char)text[pos])) {
                    n = n * 10 + (text[pos++] - '0');
                    if (n > kMaxCount) throw std::invalid_argument("Element count too large in the formula (at most 1000000)");
                }
            }
            const Element* e = element(sym);
            if (!e) throw std::invalid_argument("Unknown element " + sym);
            add_to(stack.back(), e, n);
        } else if (ch == '(') {
            stack.emplace_back();
            pos++;
        } else if (ch == ')') {
            pos++;
            long k = 1;
            if (pos < text.size() && isdigit((unsigned char)text[pos])) {
                k = 0;
                while (pos < text.size() && isdigit((unsigned char)text[pos])) {
                    k = k * 10 + (text[pos++] - '0');
                    if (k > kMaxCount) throw std::invalid_argument("Element count too large in the formula (at most 1000000)");
                }
            }
            if (stack.size() < 2) throw std::invalid_argument("Unbalanced parentheses");
            Formula grp = stack.back();
            stack.pop_back();
            for (const auto& q : grp) {
                if ((long long)q.second * k > kMaxCount) throw std::invalid_argument("Element count too large in the formula (at most 1000000)");
                add_to(stack.back(), q.first, q.second * k);
            }
        } else {
            throw std::invalid_argument("Cannot read the formula at '" + text.substr(pos) + "'");
        }
    }
    if (stack.size() != 1) throw std::invalid_argument("Cannot read the formula '" + text + "'");
    return stack[0];
}

// The built-in sum() of Python: plain left to right before 3.12; from 3.12 on (the bundled Python)
// sum() of floats compensates the rounding (Neumaier, CPython's builtin_sum_impl): the first value
// starts the sum, the compensation is added at the end. ms_formula's mono_mass and the formula
// finder use sum(), so the library follows the Python that calls it (flag from msengine_py).
class PySum {
public:
    explicit PySum(bool compensated) : comp_(compensated) {}
    void add(double x) {
        if (first_) { f_ = 0.0 + x; first_ = false; return; }   // int 0 + the first float
        if (!comp_) { f_ += x; return; }
        const double t = f_ + x;
        if (std::fabs(f_) >= std::fabs(x)) c_ += (f_ - t) + x;
        else c_ += (x - t) + f_;
        f_ = t;
    }
    double value() const {
        if (comp_ && c_ != 0.0 && std::isfinite(c_)) return f_ + c_;
        return f_;
    }
private:
    bool comp_, first_ = true;
    double f_ = 0.0, c_ = 0.0;
};

double mono_mass_py(const Formula& f, bool comp) {
    PySum s(comp);
    for (const auto& q : f) s.add(mono_of(*q.first) * q.second);
    return s.value();
}

double mono_mass(const Formula& f) {
    double s = 0.0;
    for (const auto& q : f) s += mono_of(*q.first) * q.second;
    return s;
}

double average_mass(const Formula& f) {
    double s = 0.0;
    for (const auto& q : f) {
        double m = 0, p = 0;
        for (const Iso& i : q.first->iso) { m += i.mass * i.p; p += i.p; }
        s += m / p * q.second;
    }
    return s;
}

std::string format_formula(const Formula& f) {
    std::vector<std::string> order;
    const bool hasC = std::any_of(f.begin(), f.end(), [](const std::pair<const Element*, long>& q) { return !strcmp(q.first->sym, "C"); });
    std::vector<std::string> rest;
    for (const auto& q : f) {
        const std::string s = q.first->sym;
        if (hasC && (s == "C" || s == "H")) continue;
        rest.push_back(s);
    }
    std::sort(rest.begin(), rest.end());
    if (hasC) order = {"C", "H"};
    order.insert(order.end(), rest.begin(), rest.end());
    std::string out;
    for (const std::string& e : order) {
        const long n = count_of(f, e.c_str());
        if (n > 0) {
            out += e;
            if (n > 1) out += std::to_string(n);
        }
    }
    return out;
}

// Index of the 0.1 mDa keys of one convolution step: open addressing with a
// generation stamp, so that it is cleared in constant time between the atoms
// (std::unordered_map spent most of the time of a large pattern in its node
// allocations). Only the speed differs: entries keep their order of first
// occurrence (the dict insertion order of ms_formula.isotope_pattern).
class KeyIndex {
public:
    void reset(size_t expected) {
        size_t cap = 64;
        while (cap < 2 * expected + 16) cap <<= 1;
        if (cap > keys_.size()) {
            keys_.assign(cap, 0);
            val_.assign(cap, 0);
            stamp_.assign(cap, 0);
            gen_ = 0;
        }
        mask_ = keys_.size() - 1;
        bits_ = 0;
        while ((size_t(1) << bits_) < keys_.size()) bits_++;
        if (++gen_ == 0) {   // the stamps wrapped around: clear them once
            std::fill(stamp_.begin(), stamp_.end(), 0u);
            gen_ = 1;
        }
    }
    // slot of key (found: its value is valid); insert with set()
    size_t find(long long key, bool& found) const {
        size_t h = (size_t)(((unsigned long long)key * 0x9E3779B97F4A7C15ULL) >> (64 - bits_));
        for (;;) {
            if (stamp_[h] != gen_) { found = false; return h; }
            if (keys_[h] == key) { found = true; return h; }
            h = (h + 1) & mask_;
        }
    }
    size_t value(size_t slot) const { return val_[slot]; }
    void set(size_t slot, long long key, size_t v) { keys_[slot] = key; val_[slot] = (uint32_t)v; stamp_[slot] = gen_; }
private:
    std::vector<long long> keys_;
    std::vector<uint32_t> val_, stamp_;
    uint32_t gen_ = 0;
    size_t mask_ = 0;
    int bits_ = 6;
};

// isotope_pattern of the neutral formula (m, probability) before normalisation;
// merge <= 0: nominal isotope peaks, else fine structure further apart than merge
std::vector<std::pair<double, double>> pattern(const Formula& g, double merge, double min_rel = 1e-4, bool pysum12 = false) {
    std::vector<double> dm = {0.0}, dp = {1.0};
    std::vector<double> nm, np_;
    KeyIndex idx;
    for (const auto& q : g) {
        if (q.second <= 0) continue;
        const std::vector<Iso>& iso = q.first->iso;
        for (long rep = 0; rep < q.second; rep++) {
            const size_t expect = dm.size() * iso.size();
            if (expect > 4000000000ULL) throw std::invalid_argument("The isotope pattern is too large to compute");
            nm.clear();
            np_.clear();
            nm.reserve(expect);
            np_.reserve(expect);
            idx.reset(expect);
            for (size_t a = 0; a < dm.size(); a++) {
                for (const Iso& i : iso) {
                    const double mm = dm[a] + i.mass;
                    const long long key = (long long)rne(mm * 1e4);   // 0.1 mDa fine structure
                    const double qq = dp[a] * i.p;
                    bool found;
                    const size_t slot = idx.find(key, found);
                    if (found) {
                        const size_t k = idx.value(slot);
                        const double om = nm[k], op = np_[k];
                        nm[k] = (om * op + mm * qq) / (op + qq);
                        np_[k] = op + qq;
                    } else {
                        idx.set(slot, key, nm.size());
                        nm.push_back(mm);
                        np_.push_back(qq);
                    }
                }
            }
            double top = np_[0];
            for (double v : np_) if (v > top) top = v;
            dm.clear();
            dp.clear();
            for (size_t a = 0; a < nm.size(); a++) {
                if (np_[a] >= top * min_rel * 0.01) { dm.push_back(nm[a]); dp.push_back(np_[a]); }
            }
        }
    }
    std::vector<std::pair<double, double>> dist(dm.size());
    for (size_t a = 0; a < dm.size(); a++) dist[a] = {dm[a], dp[a]};
    std::sort(dist.begin(), dist.end());
    const double mono_ref = pysum12 ? mono_mass_py(g, true) : mono_mass(g);
    // groups in order of first occurrence (dict insertion order)
    std::unordered_map<long long, size_t> gidx;
    std::vector<std::pair<double, double>> groups;   // (sum m p, sum p)
    for (const auto& d : dist) {
        const long long k = merge > 0 ? (long long)rne(d.first / merge)
                                      : (long long)rne((d.first - mono_ref) / 1.00235);
        auto it = gidx.find(k);
        if (it != gidx.end()) {
            groups[it->second].first += d.first * d.second;
            groups[it->second].second += d.second;
        } else {
            gidx[k] = groups.size();
            groups.push_back({d.first * d.second, d.second});
        }
    }
    std::vector<std::pair<double, double>> out;
    for (const auto& gr : groups) out.push_back({gr.first / gr.second, gr.second});
    std::sort(out.begin(), out.end());
    return out;
}

// rdb: ring and double bond equivalents
double rdb(const Formula& f) {
    const double tetra = count_of(f, "C") + count_of(f, "Si");
    const double tri = count_of(f, "N") + count_of(f, "P") + count_of(f, "B");
    long mono = 0;
    for (const char* e : {"H", "F", "Cl", "Br", "I", "Na", "K", "Li"}) mono += count_of(f, e);
    return tetra - mono / 2.0 + tri / 2.0 + 1.0;
}

struct Limit { const Element* e; long lo, hi; };

// "C0-80 H0-160 N0-10 O0-20 S0-3" (separators: spaces, commas, semicolons);
// C and H default to 0-80 and 0-160 when absent, other elements to 0
std::vector<Limit> parse_limits(const char* text) {
    std::vector<Limit> lim = {{element("C"), 0, 80}, {element("H"), 0, 160}};
    std::string s = text ? text : "";
    size_t pos = 0;
    while (pos < s.size()) {
        while (pos < s.size() && (s[pos] == ' ' || s[pos] == ',' || s[pos] == ';' || s[pos] == '\t' || s[pos] == '\n')) pos++;
        if (pos >= s.size()) break;
        std::string sym;
        while (pos < s.size() && isalpha((unsigned char)s[pos])) sym += s[pos++];
        const Element* e = element(sym);
        if (!e) throw std::invalid_argument("Unknown element in the limits: " + sym);
        char* end = nullptr;
        long lo = strtol(s.c_str() + pos, &end, 10);
        if (end == s.c_str() + pos) throw std::invalid_argument("Limits: number expected after " + sym);
        pos = end - s.c_str();
        long hi = lo;
        if (pos < s.size() && s[pos] == '-') {
            pos++;
            hi = strtol(s.c_str() + pos, &end, 10);
            if (end == s.c_str() + pos) throw std::invalid_argument("Limits: number expected after '-'");
            pos = end - s.c_str();
        }
        if (lo < 0 || hi < 0) throw std::invalid_argument("Limits: element counts start at 0 (" + sym + ")");
        if (lo > 100000 || hi > 100000) throw std::invalid_argument("Limits: at most 100000 atoms of an element (" + sym + ")");
        hi = std::max(lo, hi);
        bool done = false;
        for (Limit& l : lim) if (l.e == e) { l.lo = lo; l.hi = hi; done = true; }
        if (!done) lim.push_back({e, lo, hi});
    }
    return lim;
}

struct Candidate { std::string formula; double mass, ppm, score; };

// find_formulas for a neutral mass (Python adduct "M (neutral)")
std::vector<Candidate> find_formulas(double neutral, double ppm, const std::vector<Limit>& lim, long max_results,
                                     bool rules) {
    const double tol = neutral * ppm * 1e-6;
    const Element* C = element("C");
    const Element* H = element("H");
    Limit limC{C, 0, 80}, limH{H, 0, 160};
    std::vector<Limit> hetero;
    for (const Limit& l : lim) {
        if (l.e == C) limC = l;
        else if (l.e == H) limH = l;
        else if (l.hi > 0) hetero.push_back(l);
    }
    const double mC = mono_of(*C), mH = mono_of(*H);
    double combos = (double)(limC.hi - limC.lo + 1);
    for (const Limit& l : hetero) combos *= (double)(l.hi - l.lo + 1);
    if (combos > 1e9) throw std::invalid_argument("The element limits allow too many combinations (more than a billion); narrow them");
    std::vector<Candidate> out;
    std::vector<long> combo(hetero.size());
    for (size_t i = 0; i < hetero.size(); i++) combo[i] = hetero[i].lo;
    Formula f;
    for (;;) {
        double rest = 0.0;
        for (size_t i = 0; i < hetero.size(); i++) rest += mono_of(*hetero[i].e) * combo[i];
        if (rest <= neutral + tol) {
            const double cmx = (neutral - rest) / mC + 1;   // ((long) of a huge value is undefined)
            const long cmax = cmx < (double)limC.hi ? (long)cmx : limC.hi;
            for (long c = limC.lo; c <= cmax; c++) {
                const long h = (long)rne((neutral - rest - c * mC) / mH);
                if (h < limH.lo || h > limH.hi) continue;
                const double m = rest + c * mC + h * mH;
                if (std::fabs(m - neutral) > tol) continue;
                f.clear();
                f.push_back({C, c});
                f.push_back({H, h});
                for (size_t i = 0; i < hetero.size(); i++) if (combo[i]) f.push_back({hetero[i].e, combo[i]});
                if (rules) {
                    const double r = rdb(f);
                    if (r < -0.5 || std::fabs(r * 2 - rne(r * 2)) > 1e-6) continue;
                    if (c > 0) {   // plausibility (Kind and Fiehn, seven golden rules, loose)
                        // the lower bound counts the halogens with H (perfluorinated and perchlorinated compounds)
                        const long hx = h + count_of(f, "F") + count_of(f, "Cl") + count_of(f, "Br") + count_of(f, "I");
                        const double hc = (double)h / c, hxc = (double)hx / c;
                        if (!(0.1 <= hxc && hc <= 6.0) || (double)count_of(f, "N") / c > 4 || (double)count_of(f, "O") / c > 3)
                            continue;
                    } else if (h > 8 && !count_of(f, "N")) {
                        continue;
                    }
                    if (std::fabs(r - rne(r)) > 1e-6) continue;   // even electron neutral molecule
                }
                const double theo = mono_mass(f);
                const double err = (neutral - theo) / theo * 1e6;
                out.push_back({format_formula(f), theo, err, std::fabs(err) / std::max(ppm, 1e-9)});
            }
        }
        // next combination (itertools.product: the last element varies fastest)
        bool more = false;
        for (size_t k = hetero.size(); k-- > 0;) {
            if (++combo[k] <= hetero[k].hi) { more = true; break; }
            combo[k] = hetero[k].lo;
        }
        if (!more) break;
    }
    std::stable_sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) { return a.score < b.score; });
    if (max_results > 0 && (long)out.size() > max_results) out.resize(max_results);
    return out;
}

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_formula_mass(const char* formula, double* mono, double* average) {
    return ms::guarded([&] {
        if (!formula) throw std::invalid_argument("ms_formula_mass: null formula");
        const ms::Formula f = ms::parse(formula);
        if (mono) *mono = ms::mono_mass(f);
        if (average) *average = ms::average_mass(f);
    });
}

// Pattern of the formula with |z| charges: m/z = (M + |z| carrier) / |z| (z = 0:
// M + carrier). resolution > 0 merges isotopologues closer than M / resolution
// (the Python merge width), 0 groups them into the nominal isotope peaks.
MS_API int ms_isotope_pattern(const char* formula, int z, double carrier, double resolution,
                              double* mz, double* rel, long max_n, long* n) {
    return ms::guarded([&] {
        if (!formula || !mz || !rel || !n) throw std::invalid_argument("ms_isotope_pattern: null argument");
        if (max_n < 0) throw std::invalid_argument("ms_isotope_pattern: negative max_n");
        if (z < -100000 || z > 100000) throw std::invalid_argument("ms_isotope_pattern: charge out of range");
        if (!std::isfinite(carrier)) throw std::invalid_argument("ms_isotope_pattern: the carrier mass must be a finite number");
        const ms::Formula f = ms::parse(formula);
        long long atoms = 0;
        for (const auto& q : f) atoms += q.second > 0 ? q.second : 0;
        if (atoms > 100000) throw std::invalid_argument("The isotope pattern is computed for formulas of at most 100000 atoms");
        // (resolution above 1e9 or infinite: nominal peaks, as an infinite value always gave)
        const double merge = (resolution > 0 && resolution <= 1e9) ? ms::mono_mass(f) / resolution : 0.0;
        const std::vector<std::pair<double, double>> pat = ms::pattern(f, merge);
        double top = 0;
        for (const auto& p : pat) top = std::max(top, p.second);
        const int za = std::abs(z);
        long k = 0;
        for (const auto& p : pat) {
            const double r = p.second / top * 100.0;
            if (r < 1e-4 * 100.0) continue;
            if (k < max_n) {
                mz[k] = za ? (p.first + za * carrier) / za : p.first + carrier;
                rel[k] = r;
            }
            k++;
        }
        *n = std::min(k, max_n);   // truncated to max_n when the pattern is longer
    });
}

// Results as "formula\tmass\terror_ppm\n" (mass = monoisotopic mass of the
// candidate, error = (target - mass) / mass in ppm), best first; rdbe_check = 0
// skips the RDBE and plausibility rules.
MS_API const char* ms_find_formulas(double mass, double tol_ppm, const char* limits, int max_results, int rdbe_check) {
    static thread_local std::string out;
    out.clear();
    const int rc = ms::guarded([&] {
        if (!std::isfinite(mass) || !(mass > 0) || mass > 1e7) throw std::invalid_argument("ms_find_formulas: the mass must be a positive number (at most 10000000 Da)");
        if (!std::isfinite(tol_ppm) || !(tol_ppm >= 0)) throw std::invalid_argument("ms_find_formulas: the tolerance must be a positive number");
        const std::vector<ms::Limit> lim = ms::parse_limits(limits);
        const std::vector<ms::Candidate> res = ms::find_formulas(mass, tol_ppm, lim, max_results, rdbe_check != 0);
        char buf[128];
        for (const ms::Candidate& c : res) {
            snprintf(buf, sizeof buf, "\t%.8f\t%.6f\n", c.mass, c.ppm);
            out += c.formula;
            out += buf;
        }
    });
    return rc == MS_OK ? out.c_str() : nullptr;
}

}

// ---- 3.1 additions (agent P) --------------------------------------------------------
namespace ms {
namespace {

struct IonCand { std::vector<long> counts; double theo, err, rdb; };

// ms_formula.find_formulas, the enumeration of the candidates for the ion (before the isotope
// ranking and the scores): elements in the order of the Python limits dict (C and H among them),
// neutral = the neutral mass looked for, tol = its tolerance (Da), add / rem = the formulas the ion
// adds to / removes from M, electrons and zabs as in ms_formula.ADDUCTS
std::vector<IonCand> find_ion(double mz, double neutral, double tol, int zabs, double electrons,
                              const std::vector<const Element*>& els, const std::vector<long>& lo,
                              const std::vector<long>& hi, const std::string& add, const std::string& rem,
                              bool ion_itself, bool pysum12) {
    const Element* C = element("C");
    const Element* H = element("H");
    long cLo = 0, cHi = 80, hLo = 0, hHi = 160;
    std::vector<size_t> het;   // indices into els of the hetero elements (hi > 0), in order
    for (size_t i = 0; i < els.size(); i++) {
        if (els[i] == C) { cLo = lo[i]; cHi = hi[i]; }
        else if (els[i] == H) { hLo = lo[i]; hHi = hi[i]; }
        else if (hi[i] > 0) het.push_back(i);
    }
    double combos = 1.0;
    for (size_t k : het) combos *= (double)std::max(0L, hi[k] - lo[k] + 1);
    if (combos * (double)std::max(1L, cHi - cLo + 1) > 2e9)
        throw std::invalid_argument("The element limits allow too many combinations; narrow them");
    for (size_t k : het) if (hi[k] < lo[k]) return {};   // an empty range: itertools.product gives nothing
    Formula addf, remf;
    if (!add.empty()) addf = parse(add);
    if (!rem.empty()) remf = parse(rem);
    const double mC = mono_of(*C), mH = mono_of(*H);
    std::vector<IonCand> out;
    std::vector<long> combo(het.size());
    for (size_t i = 0; i < het.size(); i++) combo[i] = lo[het[i]];
    Formula f, ion;
    for (;;) {
        PySum rs(pysum12);   // sum(MONO[e] * n for e, n in zip(hetero, combo))
        for (size_t i = 0; i < het.size(); i++) rs.add(mono_of(*els[het[i]]) * (double)combo[i]);
        const double rest = het.empty() ? 0.0 : rs.value();
        if (!(rest > neutral + tol)) {
            const double v = (neutral - rest) / mC;   // int() truncates towards zero
            const double tv = std::trunc(v);
            const long cmax = (tv + 1.0 < (double)cHi) ? (long)tv + 1 : cHi;
            for (long c = cLo; c <= cmax; c++) {
                const double hv = rne((neutral - rest - (double)c * mC) / mH);
                if (!(hv >= (double)hLo && hv <= (double)hHi)) continue;
                const long h = (long)hv;
                const double m = rest + (double)c * mC + (double)h * mH;
                if (std::fabs(m - neutral) > tol) continue;
                f.clear();
                f.push_back({C, c});
                f.push_back({H, h});
                for (size_t i = 0; i < het.size(); i++) if (combo[i]) f.push_back({els[het[i]], combo[i]});
                const double r = rdb(f);
                if (r < -0.5 || std::fabs(r * 2 - rne(r * 2)) > 1e-6) continue;
                if (c > 0) {   // plausibility (Kind and Fiehn, seven golden rules, loose)
                    // the lower bound counts the halogens with H (perfluorinated and perchlorinated compounds)
                    const long hx = h + count_of(f, "F") + count_of(f, "Cl") + count_of(f, "Br") + count_of(f, "I");
                    const double hc = (double)h / (double)c, hxc = (double)hx / (double)c;
                    if (!(0.1 <= hxc && hc <= 6.0) || (double)count_of(f, "N") / (double)c > 4 ||
                        (double)count_of(f, "O") / (double)c > 3)
                        continue;
                } else if (h > 8 && !count_of(f, "N")) {
                    continue;
                }
                if (std::fabs(r - rne(r)) > 1e-6 && !ion_itself) continue;   // even electron neutral M
                ion = f;
                for (const auto& q : addf) {
                    bool done = false;
                    for (auto& w : ion) if (w.first == q.first) { w.second += q.second; done = true; break; }
                    if (!done) ion.push_back(q);
                }
                bool ok = true;
                for (const auto& q : remf) {
                    bool done = false;
                    for (auto& w : ion) if (w.first == q.first) { w.second -= q.second; ok = ok && w.second >= 0; done = true; break; }
                    if (!done) { ion.push_back({q.first, -q.second}); ok = false; }
                }
                if (!ok) continue;
                const double mono = mono_mass_py(ion, pysum12);
                const double theo = (mono - electrons * ELECTRON) / (double)zabs;
                IonCand cd;
                cd.counts.push_back(c);
                cd.counts.push_back(h);
                for (size_t i = 0; i < het.size(); i++) cd.counts.push_back(combo[i]);
                cd.theo = theo;
                cd.err = (mz - theo) / theo * 1e6;
                cd.rdb = r;
                out.push_back(std::move(cd));
                if (out.size() > 5000000) throw std::invalid_argument("More than five million candidate formulas; narrow the limits");
            }
        }
        bool more = false;   // itertools.product: the last element varies fastest
        for (size_t k = het.size(); k-- > 0;) {
            if (++combo[k] <= hi[het[k]]) { more = true; break; }
            combo[k] = lo[het[k]];
        }
        if (!more) break;
    }
    return out;
}

struct PatternOut { std::vector<double> mass, rel; };
struct FinderOut { std::vector<long> counts; std::vector<double> theo, err, rdb; };
thread_local PatternOut g_pattern;
thread_local FinderOut g_finder;

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_isotope_pattern2(const char* formula, double merge, double min_rel, int flags, const double** mass,
                               const double** rel, long* n) {
    return ms::guarded([&] {
        if (!formula || !mass || !rel || !n) throw std::invalid_argument("ms_isotope_pattern2: null argument");
        *mass = *rel = nullptr;
        *n = 0;
        if (!(min_rel >= 0) || !std::isfinite(min_rel)) throw std::invalid_argument("ms_isotope_pattern2: min_rel must be 0 or more");
        if (std::isnan(merge) || std::isinf(merge)) throw std::invalid_argument("ms_isotope_pattern2: merge must be a finite number");
        const ms::Formula f = ms::parse(formula);
        long long atoms = 0;
        for (const auto& q : f) atoms += q.second > 0 ? q.second : 0;
        if (atoms > 100000) throw std::invalid_argument("The isotope pattern is computed for formulas of at most 100000 atoms");
        const std::vector<std::pair<double, double>> pat = ms::pattern(f, merge > 0 ? merge : 0.0, min_rel, (flags & 1) != 0);
        double top = 0;
        bool first = true;
        for (const auto& p : pat) { if (first || p.second > top) top = p.second; first = false; }
        ms::g_pattern.mass.clear();
        ms::g_pattern.rel.clear();
        for (const auto& p : pat) {
            const double r = p.second / top * 100.0;   // arr[:, 1] / arr[:, 1].max() * 100.0
            if (!(r >= min_rel * 100.0)) continue;
            ms::g_pattern.mass.push_back(p.first);
            ms::g_pattern.rel.push_back(r);
        }
        *mass = ms::g_pattern.mass.data();
        *rel = ms::g_pattern.rel.data();
        *n = (long)ms::g_pattern.mass.size();
    });
}

MS_API int ms_find_formulas2(double mz, double neutral, double tol, int zabs, double electrons, const char* elements,
                             const long* lo, const long* hi, const char* add, const char* rem, int ion_itself, int flags,
                             const long** counts, const double** theo, const double** err, const double** rdbe,
                             long* ncand, int* ncols) {
    return ms::guarded([&] {
        if (!elements || !lo || !hi || !counts || !theo || !err || !rdbe || !ncand || !ncols)
            throw std::invalid_argument("ms_find_formulas2: null argument");
        *ncand = 0;
        *ncols = 0;
        if (!std::isfinite(mz) || !std::isfinite(neutral) || !std::isfinite(tol) || zabs < 1 || !std::isfinite(electrons))
            throw std::invalid_argument("ms_find_formulas2: bad numbers");
        std::vector<const ms::Element*> els;
        std::vector<long> vlo, vhi;
        std::string sym;
        const std::string text = elements;
        for (size_t i = 0; i <= text.size(); i++) {
            if (i < text.size() && text[i] != ' ' && text[i] != ',') { sym += text[i]; continue; }
            if (sym.empty()) continue;
            const ms::Element* e = ms::element(sym);
            if (!e) throw std::invalid_argument("Unknown element in the limits: " + sym);
            const size_t k = els.size();
            if (lo[k] < 0 || hi[k] < 0 || lo[k] > 100000 || hi[k] > 100000)
                throw std::invalid_argument("Limits: element counts from 0 to 100000 (" + sym + ")");
            els.push_back(e);
            vlo.push_back(lo[k]);
            vhi.push_back(hi[k]);
            sym.clear();
        }
        if (!ms::element("C") || !ms::element("H")) throw ms::Error("element table");
        const std::vector<ms::IonCand> res = ms::find_ion(mz, neutral, tol, zabs, electrons, els, vlo, vhi,
                                                         add ? add : "", rem ? rem : "", ion_itself != 0, (flags & 1) != 0);
        size_t nh = 0;
        for (size_t i = 0; i < els.size(); i++) if (els[i] != ms::element("C") && els[i] != ms::element("H") && vhi[i] > 0) nh++;
        ms::FinderOut& o = ms::g_finder;
        o.counts.clear(); o.theo.clear(); o.err.clear(); o.rdb.clear();
        for (const ms::IonCand& c : res) {
            o.counts.insert(o.counts.end(), c.counts.begin(), c.counts.end());
            o.theo.push_back(c.theo);
            o.err.push_back(c.err);
            o.rdb.push_back(c.rdb);
        }
        *counts = o.counts.data();
        *theo = o.theo.data();
        *err = o.err.data();
        *rdbe = o.rdb.data();
        *ncand = (long)res.size();
        *ncols = (int)(2 + nh);
    });
}

// ms_formula.isotope_match: agreement (0 to 100) of four peaks of a pattern with measured peaks; the
// measured peaks start at the monoisotopic peak and the pattern is aligned on its peak nearest to it
// (elements whose lightest isotope is not the most abundant one: peaks before the monoisotopic one)
MS_API double ms_isotope_match(const double* pat_mz_all, const double* pat_rel_all, long npat_all, const double* obs_mz,
                               const double* obs_rel, long nobs) {
    return ms::guarded_value(0.0, [&] {
        if (npat_all <= 0 || nobs <= 0 || !pat_mz_all || !pat_rel_all || !obs_mz || !obs_rel) return 0.0;
        long s = 0;
        {
            double best = std::fabs(pat_mz_all[0] - obs_mz[0]);
            for (long i = 1; i < npat_all && !std::isnan(best); i++) {   // np.argmin: first minimum (a NaN wins)
                const double d = std::fabs(pat_mz_all[i] - obs_mz[0]);
                if (std::isnan(d) || d < best) { best = d; s = i; }
            }
        }
        const double* pat_mz = pat_mz_all + s;
        const double* pat_rel = pat_rel_all + s;
        const long np4 = std::min(npat_all - s, 4L);
        double o[4] = {0, 0, 0, 0}, t[4] = {0, 0, 0, 0};
        for (long j = 0; j < np4; j++) {
            const double m = pat_mz[j];
            long k = 0;
            double best = std::fabs(obs_mz[0] - m);
            for (long i = 1; i < nobs && !std::isnan(best); i++) {   // np.argmin: first minimum (a NaN wins)
                const double d = std::fabs(obs_mz[i] - m);
                if (std::isnan(d) || d < best) { best = d; k = i; }
            }
            o[j] = best <= std::max(0.02, m * 10e-6) ? obs_rel[k] : 0.0;
        }
        if (!(o[0] > 0)) return 0.0;
        const double o0 = o[0], p0 = pat_rel[0];
        double num = 0.0, den = 0.0;   // np.sum of at most 4 values: in order
        for (long j = 0; j < np4; j++) {
            const double oo = o[j] / o0 * 100.0;
            t[j] = pat_rel[j] / p0 * 100.0;
            const double w = std::max(t[j], 5.0);
            num += std::fabs(oo - t[j]) / w * std::min(t[j], 100.0);
            den += std::min(t[j], 100.0);
        }
        const double diff = num / den;
        return std::max(0.0, 100.0 * (1.0 - diff));
    });
}

// ms_formula.measured_pattern: heights of the isotope peaks M, M + 1.00335 / z, ... in a spectrum
MS_API int ms_measured_pattern(const double* mz_p, const double* it_p, long n, long stride, double mono_mz, int z, int npk,
                               double ppm, double* out_mz, double* out_h) {
    return ms::guarded([&] {
        if (!out_mz || !out_h || n < 0 || npk < 0 || stride < 1 || (n > 0 && (!mz_p || !it_p)))
            throw std::invalid_argument("ms_measured_pattern: bad arguments");
        auto mz = [&](long i) { return mz_p[i * stride]; };
        auto it = [&](long i) { return it_p[i * stride]; };
        const int za = std::max(std::abs(z), 1);
        for (int k = 0; k < npk; k++) {
            const double m = mono_mz + (double)k * 1.00335 / (double)za;
            const double tol = std::max(m * ppm * 1e-6, 0.004);
            const double a = m - tol * 3, b = m + tol * 3;
            long best = -1;
            for (long i = 0; i < n; i++) {
                const double x = mz(i);
                if (!(x >= a && x <= b)) continue;
                if (best < 0) { best = i; continue; }
                if (std::isnan(it(best))) continue;                     // np.argmax: the first NaN wins
                if (std::isnan(it(i)) || it(i) > it(best)) best = i;
            }
            if (best < 0) { out_mz[k] = m; out_h[k] = 0.0; }
            else { out_mz[k] = mz(best); out_h[k] = it(best); }
        }
    });
}

}
