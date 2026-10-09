// Single precision logarithm and exponential for the UniDec engine of msengine (udengine.cpp).
//
// unidec.exe takes logf and expf from Intel's math library; the C runtimes of this library's two
// builds differ (glibc on Linux; on Windows the MinGW functions are about five times slower and
// made the charge smoothing of every iteration the slowest step). These evaluate the functions in
// double precision (a table and a short polynomial, relative error about 3e-16) and round once to
// float: the correctly rounded result except where the exact value lies within about 3e-16 of the
// midpoint between two floats, the same on every platform, with no call into the C runtime.
// (Written for msengine; not part of UniDec.)
#pragma once
#include <cstdint>
#include <cstring>
#include <limits>
#if defined(__SSE2__)
#include <emmintrin.h>
#endif

namespace ms {
namespace udengine {
namespace udmath {

// log: m in [0.75, 1.5) in 96 intervals of 1/128; 1/c of each interval centre rounded to 29 bits
// (so that m * inv is exact in a double), 1 for the two intervals next to 1; logc = -log(inv)
inline constexpr double kLogInv[96] = {
    1.3264248706400394, 1.3128205128014088, 1.299492385238409, 1.2864321619272232, 1.2736318409442902, 1.2610837444663048,
    1.2487804889678955, 1.2367149740457535, 1.2248803824186325, 1.213270142674446, 1.2018779329955578, 1.1906976737082005,
    1.1797235012054443, 1.1689497716724873, 1.1583710424602032, 1.1479820646345615, 1.137777779251337, 1.1277533024549484,
    1.1179039292037487, 1.1082251071929932, 1.0987124480307102, 1.0893617011606693, 1.0801687762141228, 1.0711297057569027,
    1.0622406639158726, 1.0534979440271854, 1.0448979586362839, 1.0364372469484806, 1.0281124487519264, 1.0199203193187714,
    1.011857707053423, 1.0, 1.0, 0.9884169884026051, 0.9808429125696421, 0.9733840301632881,
    0.9660377353429794, 0.9588014986366034, 0.951672863215208, 0.9446494467556477, 0.937728937715292, 0.9309090916067362,
    0.9241877254098654, 0.91756272315979, 0.9110320284962654, 0.9045936390757561, 0.8982456140220165, 0.8919860627502203,
    0.8858131486922503, 0.8797250855714083, 0.8737201374024153, 0.8677966110408306, 0.8619528617709875, 0.8561872914433479,
    0.8504983391612768, 0.8448844887316227, 0.8393442630767822, 0.833876222372055, 0.8284789640456438, 0.8231511246412992,
    0.8178913742303848, 0.8126984126865864, 0.8075709771364927, 0.8025078363716602, 0.7975077889859676, 0.7925696596503258,
    0.7876923084259033, 0.7828746177256107, 0.778115501627326, 0.773413896560669, 0.7687687687575817, 0.7641791049391031,
    0.759643916040659, 0.7551622427999973, 0.7507331371307373, 0.7463556844741106, 0.7420289851725101, 0.7377521619200706,
    0.7335243560373783, 0.729344729334116, 0.7252124641090631, 0.7211267612874508, 0.7170868348330259, 0.7130919229239225,
    0.7091412749141455, 0.7052341606467962, 0.7013698630034924, 0.6975476834923029, 0.6937669385224581, 0.6900269538164139,
    0.6863270774483681, 0.6826666668057442, 0.6790450923144817, 0.6754617411643267, 0.6719160098582506, 0.6684073098003864};
inline constexpr double kLogC[96] = {
    -0.28248725570564415, -0.27217788590126374, -0.2619737153195684, -0.251872620628185, -0.24187253653690205, -0.23197146593254026,
    -0.22216746627247688, -0.2124586497590019, -0.20284319222371316, -0.1933193114109496, -0.18388527770784754, -0.17453941575527115,
    -0.16528009000778035, -0.15610571464850975, -0.14701474446065693, -0.1380056746347063, -0.1290770435702628, -0.1202274256593836,
    -0.11145544009586365, -0.10275973302644636, -0.09413899244181301, -0.08559192944773668, -0.07711730319891213, -0.06871389128203519,
    -0.0603805110034594, -0.05211600269606895, -0.04391923341096654, -0.03578910783703337, -0.027724546996220793, -0.0197245059298552,
    -0.011787955300932868, 0.0, 0.0, 0.01165061723452719, 0.019342962137363046, 0.02697658796013655,
    0.03455238203052868, 0.04207121338954215, 0.0495339343219213, 0.056941376123652035, 0.06429435071994917, 0.07159365243758518,
    0.07884006194788262, 0.08603433827312573, 0.09317722482507947, 0.10026945371664793, 0.10731173580363997, 0.11430477124367884,
    0.12124924374200904, 0.1281458230775558, 0.1349951635261467, 0.1417979108561752, 0.1485546945341399, 0.15526612835815118,
    0.16193281992734324, 0.16855536069511262, 0.17513433119652658, 0.1816703020598968, 0.1881638328474645, 0.19461546862371829,
    0.20102574553672178, 0.2073951943606225, 0.2137243303654205, 0.2200136590765336, 0.2262636776172674, 0.23247487848115958,
    0.23864773691885244, 0.24478272643224283, 0.25088030614756623, 0.256940931828823, 0.26296504551543326, 0.26895308674159946,
    0.27490548702240053, 0.28082266169307885, 0.2867050337352769, 0.29255300356676833, 0.29836697300290665, 0.30414733473970096,
    0.30989447671878256, 0.31560877900085527, 0.32129061311584645, 0.3269403439917712, 0.3325583371618334, 0.3381449427208719,
    0.343700512900168, 0.3492253885556515, 0.35471990911748097, 0.3601844041934642, 0.3656191983313279, 0.37102461865174163,
    0.37640097560081054, 0.3817485812871215, 0.3870677437396998, 0.39235876098848965, 0.3976219316002889, 0.40285754605441165};
// exp: 2^(j / 32)
inline constexpr double kExp2[32] = {
    1.0, 1.0218971486541166, 1.0442737824274138, 1.0671404006768237, 1.0905077326652577, 1.1143867425958924,
    1.1387886347566916, 1.1637248587775775, 1.189207115002721, 1.215247359980469, 1.241857812073484, 1.2690509571917332,
    1.2968395546510096, 1.3252366431597413, 1.3542555469368927, 1.383909881963832, 1.4142135623730951, 1.4451808069770467,
    1.4768261459394993, 1.5091644275934228, 1.5422108254079407, 1.5759808451078865, 1.6104903319492543, 1.645755478153965,
    1.681792830507429, 1.718619298122478, 1.7562521603732995, 1.7947090750031072, 1.8340080864093424, 1.8741676341103,
    1.9152065613971474, 1.9571441241754002};
inline constexpr double kLn2 = 0.6931471805599453;
inline constexpr double kInvLn2_32 = 46.16624130844683;
inline constexpr double kLn2_32_hi = 0.021660849392446835;   // 38 significant bits: k * hi is exact for |k| < 2^15
inline constexpr double kLn2_32_lo = 5.145609244655338e-14;

inline float logf(float x) {
    uint32_t u;
    std::memcpy(&u, &x, 4);
    int e = 0;
    if (u - 0x00800000u >= 0x7f000000u) {   // not a positive normal number
        if (x != x) return x;
        if (x < 0) return std::numeric_limits<float>::quiet_NaN();
        if (x == 0) return -std::numeric_limits<float>::infinity();
        if (u == 0x7f800000u) return x;    // +inf
        x *= 8388608.0f;                   // subnormal: times 2^23 (exact)
        std::memcpy(&u, &x, 4);
        e = -23;
    }
    // m = 1.mantissa in [1, 2), halved (and e + 1) from 1.5 on: m in [0.75, 1.5) (no branch: the
    // test is a coin flip for the data and a mispredicted branch cost more than the rest)
    const uint32_t hi = (u & 0x7fffffu) >= 0x400000u ? 1u : 0u;
    e += (int)(u >> 23) - 127 + (int)hi;
    static constexpr double kHalf[2] = {1.0, 0.5};
    const double m = (1.0 + (double)(u & 0x7fffffu) * (1.0 / 8388608.0)) * kHalf[hi];
    const int i = (int)((m - 0.75) * 128.0);
    const double r = m * kLogInv[i] - 1.0;   // exact
    // log1p(r) = r - r^2/2 + r^3/3 - ... - r^8/8 for |r| <= 1/128 (Estrin's scheme: short dependency chains)
    const double r2 = r * r, r4 = r2 * r2;
    const double q = (-0.5 + r * (1.0 / 3)) + r2 * (-0.25 + r * 0.2) + r4 * ((-1.0 / 6 + r * (1.0 / 7)) + r2 * -0.125);
    return (float)(((double)e * kLn2 + kLogC[i]) + (r + r2 * q));
}

inline float expf(float x) {
    if (x != x) return x;
    if (x > 89.0f) return std::numeric_limits<float>::infinity();
    if (x < -104.0f) return 0.0f;   // below half the smallest float
    const double z = (double)x * kInvLn2_32;
    const double kd = (z + 6755399441055744.0) - 6755399441055744.0;   // round to nearest (|z| < 2^51)
    const int k = (int)kd;
    const double r = ((double)x - kd * kLn2_32_hi) - kd * kLn2_32_lo;   // |r| <= ln2 / 64
    // expm1(r) = r + r^2/2 + ... + r^6/720 (Estrin's scheme)
    const double r2 = r * r;
    const double p = r + r2 * ((0.5 + r * (1.0 / 6)) + r2 * ((1.0 / 24 + r * (1.0 / 120)) + r2 * (1.0 / 720)));
    const double t = kExp2[k & 31];
    const uint64_t sb = (uint64_t)((k >> 5) + 1023) << 52;   // 2^(k / 32, rounded down)
    double scale;
    std::memcpy(&scale, &sb, 8);
    return (float)((t + t * p) * scale);
}

// logf and expf of n values: the same results as the functions above, two at a time with SSE2 (every
// x86-64 processor has it) where all four values of a group are in the plain range, else one by one.
// The double operations are the same, in the same order, as in the functions above.
#if defined(__SSE2__)
inline void log_batch(const float* x, float* y, long n) {
    long i = 0;
    const __m128i one_m = _mm_set1_epi32(-1), lim = _mm_set1_epi32(0x7f000000), bias = _mm_set1_epi32(0x00800000);
    const __m128i mmask = _mm_set1_epi32(0x7fffff), half_m = _mm_set1_epi32(0x3fffff), e127 = _mm_set1_epi32(127);
    for (; i + 4 <= n; i += 4) {
        const __m128i u = _mm_loadu_si128((const __m128i*)(x + i));
        const __m128i v = _mm_sub_epi32(u, bias);
        const __m128i ok = _mm_and_si128(_mm_cmpgt_epi32(v, one_m), _mm_cmplt_epi32(v, lim));
        if (_mm_movemask_epi8(ok) != 0xffff) {
            for (int k = 0; k < 4; k++) y[i + k] = logf(x[i + k]);
            continue;
        }
        const __m128i mant = _mm_and_si128(u, mmask);
        const __m128i hi = _mm_cmpgt_epi32(mant, half_m);   // -1 where the mantissa is from 1.5 on
        const __m128i e = _mm_sub_epi32(_mm_sub_epi32(_mm_srli_epi32(u, 23), e127), hi);
        __m128 out[2];
        for (int h = 0; h < 2; h++) {
            const __m128i mh = h ? _mm_srli_si128(mant, 8) : mant;
            const __m128i hh = h ? _mm_srli_si128(hi, 8) : hi;
            const __m128i eh = h ? _mm_srli_si128(e, 8) : e;
            const __m128d mul = _mm_add_pd(_mm_set1_pd(1.0), _mm_mul_pd(_mm_set1_pd(0.5), _mm_cvtepi32_pd(hh)));   // 1 or 0.5
            const __m128d m = _mm_mul_pd(_mm_add_pd(_mm_set1_pd(1.0), _mm_mul_pd(_mm_cvtepi32_pd(mh), _mm_set1_pd(1.0 / 8388608.0))), mul);
            const __m128i idx = _mm_cvttpd_epi32(_mm_mul_pd(_mm_sub_pd(m, _mm_set1_pd(0.75)), _mm_set1_pd(128.0)));
            const int i0 = _mm_cvtsi128_si32(idx), i1 = _mm_cvtsi128_si32(_mm_srli_si128(idx, 4));
            const __m128d inv = _mm_set_pd(kLogInv[i1], kLogInv[i0]);
            const __m128d lc = _mm_set_pd(kLogC[i1], kLogC[i0]);
            const __m128d r = _mm_sub_pd(_mm_mul_pd(m, inv), _mm_set1_pd(1.0));
            const __m128d r2 = _mm_mul_pd(r, r), r4 = _mm_mul_pd(r2, r2);
            const __m128d a = _mm_add_pd(_mm_set1_pd(-0.5), _mm_mul_pd(r, _mm_set1_pd(1.0 / 3)));
            const __m128d b = _mm_mul_pd(r2, _mm_add_pd(_mm_set1_pd(-0.25), _mm_mul_pd(r, _mm_set1_pd(0.2))));
            const __m128d c = _mm_mul_pd(r4, _mm_add_pd(_mm_add_pd(_mm_set1_pd(-1.0 / 6), _mm_mul_pd(r, _mm_set1_pd(1.0 / 7))),
                                                       _mm_mul_pd(r2, _mm_set1_pd(-0.125))));
            const __m128d q = _mm_add_pd(_mm_add_pd(a, b), c);
            const __m128d res = _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_cvtepi32_pd(eh), _mm_set1_pd(kLn2)), lc),
                                           _mm_add_pd(r, _mm_mul_pd(r2, q)));
            out[h] = _mm_cvtpd_ps(res);
        }
        _mm_storeu_ps(y + i, _mm_movelh_ps(out[0], out[1]));
    }
    for (; i < n; i++) y[i] = logf(x[i]);
}

inline void exp_batch(const float* x, float* y, long n) {
    long i = 0;
    const __m128 lo = _mm_set1_ps(-104.0f), hi = _mm_set1_ps(89.0f);
    for (; i + 4 <= n; i += 4) {
        const __m128 xv = _mm_loadu_ps(x + i);
        // in [-104, 89] (false for NaN)
        const __m128 ok = _mm_and_ps(_mm_cmpge_ps(xv, lo), _mm_cmple_ps(xv, hi));
        if (_mm_movemask_ps(ok) != 15) {
            for (int k = 0; k < 4; k++) y[i + k] = expf(x[i + k]);
            continue;
        }
        __m128 out[2];
        for (int h = 0; h < 2; h++) {
            const __m128d xd = _mm_cvtps_pd(h ? _mm_movehl_ps(xv, xv) : xv);
            const __m128d z = _mm_mul_pd(xd, _mm_set1_pd(kInvLn2_32));
            const __m128d kd = _mm_sub_pd(_mm_add_pd(z, _mm_set1_pd(6755399441055744.0)), _mm_set1_pd(6755399441055744.0));
            const __m128i k = _mm_cvttpd_epi32(kd);
            const __m128d r = _mm_sub_pd(_mm_sub_pd(xd, _mm_mul_pd(kd, _mm_set1_pd(kLn2_32_hi))), _mm_mul_pd(kd, _mm_set1_pd(kLn2_32_lo)));
            const __m128d r2 = _mm_mul_pd(r, r);
            const __m128d inner = _mm_add_pd(_mm_add_pd(_mm_set1_pd(1.0 / 24), _mm_mul_pd(r, _mm_set1_pd(1.0 / 120))),
                                             _mm_mul_pd(r2, _mm_set1_pd(1.0 / 720)));
            const __m128d p = _mm_add_pd(r, _mm_mul_pd(r2, _mm_add_pd(_mm_add_pd(_mm_set1_pd(0.5), _mm_mul_pd(r, _mm_set1_pd(1.0 / 6))),
                                                                     _mm_mul_pd(r2, inner))));
            const int k0 = _mm_cvtsi128_si32(k), k1 = _mm_cvtsi128_si32(_mm_srli_si128(k, 4));
            const __m128d t = _mm_set_pd(kExp2[k1 & 31], kExp2[k0 & 31]);
            const __m128i ex = _mm_add_epi32(_mm_srai_epi32(k, 5), _mm_set1_epi32(1023));
            const __m128d scale = _mm_castsi128_pd(_mm_slli_epi64(_mm_unpacklo_epi32(ex, _mm_setzero_si128()), 52));
            out[h] = _mm_cvtpd_ps(_mm_mul_pd(_mm_add_pd(t, _mm_mul_pd(t, p)), scale));
        }
        _mm_storeu_ps(y + i, _mm_movelh_ps(out[0], out[1]));
    }
    for (; i < n; i++) y[i] = expf(x[i]);
}
#else
inline void log_batch(const float* x, float* y, long n) { for (long i = 0; i < n; i++) y[i] = logf(x[i]); }
inline void exp_batch(const float* x, float* y, long n) { for (long i = 0; i < n; i++) y[i] = expf(x[i]); }
#endif

}  // namespace udmath
}  // namespace udengine
}  // namespace ms
