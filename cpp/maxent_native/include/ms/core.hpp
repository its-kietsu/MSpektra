#pragma once
#include <cstdint>
#include <filesystem>
#include <optional>
#include <utility>
#include <vector>

namespace ms {
using Vector = std::vector<double>;
struct Point { double x, y; };
using Spectrum = std::vector<Point>;
using Range = std::pair<double,double>;
struct Peak { double mass, apex, height, area, relative, fraction; int isotope_count=0; };
constexpr double proton = 1.007276467;
double median(Vector v);
Spectrum restrict_spectrum(Spectrum s, std::optional<Range> range = {});
Spectrum engine_unique(const Spectrum& s);
Spectrum read_xy(const std::filesystem::path& path);
void write_xy(const std::filesystem::path& path, const Spectrum& s);
std::vector<Peak> pick_mass_peaks(const Spectrum& s, double window, double threshold);
Spectrum broaden_mass(const Spectrum& spectrum,double display_resolution,double model_resolution);
std::vector<Peak> group_isotopes(const Spectrum& spectrum,double threshold=.05,double spacing=1.00235);
Spectrum baseline(const Spectrum& s, double width);
double r_squared(const Vector& data, const Vector& fit);
Vector gaussian(const Vector& v, double sigma);
Vector cell_average(const Vector& x, const Vector& y, const Vector& centres, double width);
Vector grid(const std::vector<Range>& windows, double spacing, std::size_t limit);
Spectrum sum_profiles(const std::vector<Spectrum>& spectra);
double integrate(const Spectrum& trace, double baseline_start, double baseline_end);
Vector decode_pda(const std::vector<std::uint8_t>& raw, std::size_t records, std::size_t wavelengths);

// Exact discrete forward operator used by the current MaxEnt implementation.
// Optimizer migration is independent of this operator and its adjoint.
class ChargeModel {
    Vector coordinate_;
    std::vector<std::vector<std::size_t>> index_;
    std::vector<Vector> w0_, w1_;
    std::size_t masses_;
    double sigma_;
public:
    ChargeModel(Vector coordinate, double spacing, const Vector& masses,
                const Vector& charges, double carrier, double fwhm, bool logarithmic);
    Vector scatter(std::size_t charge, const Vector& mass) const;
    Vector forward(const Vector& mass, const Vector& envelope) const;
    Vector adjoint(const Vector& residual, const Vector& envelope) const;
    std::vector<Vector> columns(const Vector& mass) const;
};
}
