/*
 * msengine: the computational core of MSpektra in C++ (C interface).
 *
 * Phase 1 of the C++ conversion. Everything the LCMS / HRMS windows compute
 * lives here: reading data files (Bruker .d via the Bruker SDKs, mzML),
 * chromatograms, averaging, centroiding, peak picking and the maximum
 * entropy deconvolution. The interface is plain C so that the current
 * Python windows (ctypes) and the future C++ interface call the same code.
 *
 * Conventions
 *   - Functions that return a status code return 0 on success, a negative code
 *     on failure (functions that return a handle or a string return NULL on
 *     failure where noted). ms_last_error() gives the message of the call that
 *     failed last in the calling thread: every function that returns a status
 *     code (and ms_open, ms_find_formulas) clears it on entry, so after such a
 *     call it describes this call ("" on success). Functions that return a value
 *     (counts, strings, ms_apex, ms_resolution) leave it alone unless they fail.
 *   - No function lets an exception, an abort or a crash through for bad
 *     arguments: NULL pointers, negative lengths, NaN or infinite values and
 *     out of range indexes give MS_BAD_ARG (or 0, "" for functions without a
 *     code); damaged or truncated files give MS_ERROR with the reason; lack of
 *     memory gives MS_ERROR ("Not enough memory ..."). Arrays passed in must hold
 *     the number of elements given (the library cannot check that).
 *   - Arrays handed out by the library: see each section (data files: owned by
 *     the calling thread; deconvolution results: owned by the result object;
 *     spectrum tools: thread local scratch of the function).
 *   - Masses in Da, m/z in Th, times in minutes, intensities as recorded.
 *   - Progress callbacks: int cb(void* user, long i, long n), called in the
 *     calling thread. Returning 0 cancels (the function then returns
 *     MS_CANCELLED). A callback must not call functions on the same data file
 *     (MS_BAD_ARG); it may call ms_close.
 *   - Threads: different files, results and spectrum tools may be used from
 *     any threads at the same time; for one data file see below.
 */
#ifndef MSENGINE_H
#define MSENGINE_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#  if defined(MSENGINE_BUILD)
#    define MS_API __declspec(dllexport)
#  else
#    define MS_API __declspec(dllimport)
#  endif
#else
#  define MS_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MS_OK            0
#define MS_ERROR        -1
#define MS_CANCELLED    -2
#define MS_BAD_ARG      -3
#define MS_NOT_FOUND    -4
#define MS_UNSUPPORTED  -5

typedef int (*ms_progress_fn)(void* user, long i, long n);

MS_API const char* ms_version(void);
MS_API const char* ms_last_error(void);
MS_API int ms_set_threads(int n);     /* 0 = all cores but one (above 4); at most 4 x the cores (at least 8, at most 256) */
MS_API int ms_get_threads(void);
/* Memory (bytes) one deconvolution may use on this computer: half of the
 * physical memory, at least 1 GB. The deconvolutions refuse settings that
 * would need more (MS_BAD_ARG with an explanation). */
MS_API double ms_memory_budget(void);

/* ------------------------------------------------------------------ data */
/* A data file: Bruker .d (analysis.baf via baf2sql_c, analysis.tsf or
 * analysis.tdf via timsdata), Shimadzu .lcd or an mzML file.
 *
 * Threads: the functions of one file may be called from several threads at
 * once (the GUI averages in a worker thread while its main thread asks for
 * scans and chromatograms); the file serialises its reader, one scan at a time
 * for the averages of mzML and BAF files. Arrays and strings handed out for a
 * file belong to the calling thread: they stay valid until the same thread makes
 * the same call on the same file again (ms_event_scans: for the same event),
 * makes any call on a file after this file was closed, or ends.
 * Every other value than an open file (NULL, a closed file, a stale pointer)
 * makes the functions fail with MS_BAD_ARG (0 or "" for those without a code). */
typedef struct ms_file ms_file;

/* Opens a file (NULL on failure, the reason in ms_last_error). sdk_dir: folder
 * holding baf2sql_c.dll and timsdata.dll (or their .so on Linux), or their
 * baf2sql/ and timsdata/ sub folders; NULL = next to the library. */
MS_API ms_file* ms_open(const char* path, const char* sdk_dir, ms_progress_fn cb, void* user);
/* Closes a file at once: calls still running on it in other threads stop at
 * their next progress tick (MS_CANCELLED) and the last of them frees it.
 * NULL, a closed file or anything else that is not an open file: nothing. */
MS_API void     ms_close(ms_file* f);

MS_API const char* ms_file_kind(const ms_file* f);      /* "Bruker .d", "Bruker .d (TSF)", "Bruker .d (timsTOF)", "mzML" */
MS_API const char* ms_file_summary(const ms_file* f);   /* one line, as the window shows it */
MS_API int  ms_file_has_profile(const ms_file* f);
MS_API long ms_file_n_scans(const ms_file* f);          /* MS1 scans */
MS_API long ms_file_n_msms(const ms_file* f);
MS_API int  ms_file_n_events(const ms_file* f);         /* scan events (polarities) */
MS_API int  ms_file_event_polarity(const ms_file* f, int event);     /* +1 or -1 */
MS_API int  ms_file_event_mz_range(const ms_file* f, int event, double* lo, double* hi); /* MS_NOT_FOUND if unknown */
MS_API const char* ms_file_recalibration(const ms_file* f);          /* "" or e.g. "DataAnalysis, 15.09.2026" */

/* Scans of one event: times (min), TIC, base peak intensity; n = count. */
MS_API long ms_event_scans(const ms_file* f, int event, const double** rt, const double** tic, const double** bpc);

/* Mass chromatogram of one event: intensity between mz-|tol| and mz+|tol| per scan
 * (length as ms_event_scans); mz and tol must be finite. */
MS_API int ms_xic(ms_file* f, int event, double mz, double tol, const double** y, long* n);

/* Averaged spectrum of event between t0 and t1 (min, either order, +-inf allowed,
 * NaN not), minus the mean of the background ranges (bg: pairs t0,t1; n_bg
 * pairs, 0 = none). Returns the number of scans averaged in n_scans; the
 * spectrum as mz / intensity arrays of length n (0 when no scan is in the range). */
MS_API int ms_average(ms_file* f, int event, double t0, double t1, const double* bg, int n_bg,
                      const double** mz, const double** it, long* n, long* n_scans,
                      ms_progress_fn cb, void* user);

/* One scan (the nearest to t; t not NaN) as a spectrum; index may be NULL. */
MS_API int ms_scan(ms_file* f, int event, double t, const double** mz, const double** it, long* n, long* index);

/* ------------------------------------------------------- spectrum tools */
/* Centroids of a profile spectrum (local maxima above rel x the top, Gaussian
 * apex through the three top points). Output arrays owned by the library
 * (thread local scratch, valid until the next call of this function). */
MS_API int ms_centroid(const double* mz, const double* it, long n, double rel,
                       const double** cmz, const double** cit, long* nc);

/* Apex of the profile peak nearest to x (Gaussian through the points above
 * half height, as DataAnalysis labels it); returns x itself when x is not
 * at a maximum. */
MS_API double ms_apex(const double* mz, const double* it, long n, double x);

/* Resolving power m / FWHM measured on the tallest peaks. */
MS_API double ms_resolution(const double* mz, const double* it, long n);

/* --------------------------------------------------------- deconvolution */
typedef struct {
    int    z_lo, z_hi;            /* charge range */
    double mass_lo, mass_hi;      /* Da */
    double mass_step;             /* Da */
    double adduct_mass;           /* mass of the charge carrier, with its sign (H+ : 1.007276) */
    double mz_lo, mz_hi;          /* 0, 0 = the whole spectrum */
    double min_intensity;         /* 0 = none; peaks whose top is below are removed */
    int    min_intensity_pct;     /* 1: min_intensity is % of the base peak */
    double peak_width;            /* m/z FWHM, 0 = measured */
    double resolution;            /* resolving power for TOF data (log m/z grid), 0 = measured; <0 = linear m/z grid */
    int    resolved_isotopes;     /* 1: every isotope peak (needs a fine step); 0: envelopes (average masses) */
    int    rounds, min_rounds, iterations;   /* 10, 4, 50 */
    double peak_window;           /* Da, for the peak picking */
    double peak_threshold;        /* 0..1 of the tallest */
    int    baseline;              /* 1: subtract a rolling baseline first */
    double baseline_width;        /* 15 */
    /* since 2.8 (zero = automatic, so a zero filled struct of the old layout keeps working) */
    int    envelope_nodes;        /* charge envelope nodes across the mass range: 0 = automatic (log spaced to start
                                   * with, then at the ends of the mass range and at the species found, at most 12),
                                   * 1 = one envelope shared by all masses, n = n nodes log spaced */
    int    noise_model;           /* 0 or 1 = noise floor + counting noise, both measured on the data (sigma_i^2 =
                                   * s0^2 + y_i / g); 2 = one constant sigma (the 2.7 model); 3 = as 0 plus a term
                                   * proportional to the intensity for the peak shape mismatch (for comparisons) */
    double chi2_target;           /* chi squared per signal point aimed at by the alpha search, 0 = 1.0 */
} ms_maxent_params;

typedef struct ms_maxent_result ms_maxent_result;

MS_API int ms_maxent(const double* mz, const double* it, long n, const ms_maxent_params* p,
                     ms_maxent_result** out, ms_progress_fn cb, void* user);
MS_API void ms_maxent_free(ms_maxent_result* r);

/* Zero charge mass spectrum (mass, intensity), the fit on the m/z grid, the
 * charge envelope (one value per charge), peaks, and notes. */
MS_API long ms_maxent_mass(const ms_maxent_result* r, const double** mass, const double** it);
MS_API long ms_maxent_fit(const ms_maxent_result* r, const double** mz, const double** it);  /* NaN rows separate gaps */
MS_API long ms_maxent_zdist(const ms_maxent_result* r, const double** z, const double** it);
MS_API long ms_maxent_peaks(const ms_maxent_result* r, const double** mass, const double** height,
                            const double** area, const double** apex, const int** n_iso);
MS_API double ms_maxent_r2(const ms_maxent_result* r);
MS_API double ms_maxent_chi2(const ms_maxent_result* r);
MS_API int    ms_maxent_rounds(const ms_maxent_result* r);
MS_API const char* ms_maxent_notes(const ms_maxent_result* r);

/* Peak picking of any mass spectrum: local maxima within +-window above
 * threshold x the top; mass = centroid of the top half, apex = grid point. */
MS_API int ms_pick_peaks(const double* mass, const double* it, long n, double window, double threshold,
                         const double** pmass, const double** height, const double** area, long* np);

/* Isotope resolved mass spectrum -> one entry per species (average mass,
 * most abundant isotope, number of isotope peaks). */
MS_API int ms_group_isotopes(const double* mass, const double* it, long n, double threshold,
                             const double** avg, const double** apex, const double** height,
                             const double** area, const int** n_iso, long* np);

/* ---------------------------------------------------------------- phase 2 */
/* Shimadzu .lcd files are opened with ms_open as well (kind "Shimadzu .lcd");
 * events = scan events of the method (polarity switching), each with its
 * polarity and m/z range; spectra are centroid (unit resolution). */

/* PDA (UV/Vis) data of a Shimadzu .lcd: wavelengths (nm), times (min) and the
 * absorbance matrix (n_times x n_wl, row major, mAU). MS_NOT_FOUND if none. */
MS_API int ms_pda(ms_file* f, const double** wl, long* n_wl, const double** t, long* n_t, const double** a);
/* Sample information of a Shimadzu file as "key\tvalue\n" lines (sample name,
 * vial, volume, method, batch, tune file, operator, date). */
MS_API const char* ms_sample_info(ms_file* f);
/* Bin width (Da) of the Shimadzu averaged and single scan spectra (the
 * LCMS window's "Bin width"); 0 restores the default of 0.05 Da, otherwise
 * 0.000001 to 10 Da (MS_BAD_ARG outside). Other formats keep their native axis
 * and ignore it: MS_UNSUPPORTED. */
MS_API int ms_file_set_bin_width(ms_file* f, double binw);

/* ------------------------------------------------- LC peak integration */
typedef struct { double t0, t1, rt, height, area, base0, base1; } ms_lc_peak;
/* Automatic integration of a chromatogram: peaks above threshold_pct of the
 * tallest, at least min_width_s wide; baseline through the peak borders. */
MS_API int ms_integrate_auto(const double* t, const double* y, long n, double threshold_pct, double min_width_s,
                             const ms_lc_peak** peaks, long* np);
/* One peak between t0 and t1 (manual range, baseline between the borders). */
MS_API int ms_integrate_range(const double* t, const double* y, long n, double t0, double t1, ms_lc_peak* out);

/* ------------------------------------------------------ m/z calibration */
/* model: 0 linear, 1 quadratic, 2 cubic, 3 TOF (sqrt), 4 HPC quadratic + order k correction;
 * measured/reference arrays of n calibrants. Fills coefficients (up to 16) and the
 * errors (ppm) after calibration; rms and cross validated rms. */
MS_API int ms_calibration_fit(int model, int order, const double* measured, const double* reference, long n,
                              double* coef, int* ncoef, double* err_ppm, double* rms, double* cv_rms);
MS_API int ms_calibration_apply(int model, int order, const double* coef, int ncoef, double m_lo, double m_hi,
                                const double* mz, double* out, long n);
/* Finds the calibrant peaks (reference list) in a spectrum within tol_ppm, above min_rel of the tallest found. */
MS_API int ms_find_calibrants(const double* mz, const double* it, long n, int profile, const double* reference, long nref,
                              double tol_ppm, double min_rel, double* measured, double* relint, int* found);

/* -------------------------------------------------------- formulas */
/* Monoisotopic mass of a formula ("C6H12O6", charge carriers not included);
 * MS_BAD_ARG for an unknown element. */
MS_API int ms_formula_mass(const char* formula, double* mono, double* average);
/* Isotope pattern (m/z, relative intensity %) of a formula with charge z and
 * the carrier mass (per charge, signed); resolution for the merging of
 * isotopologues (0: nominal). n out of max_n. */
MS_API int ms_isotope_pattern(const char* formula, int z, double carrier, double resolution,
                              double* mz, double* rel, long max_n, long* n);
/* Formula finder: elements limits as "C0-80 H0-160 N0-10 O0-20 S0-3", target neutral
 * mass with tolerance (ppm); results as "formula\tmass\terror_ppm\n" lines, best first. */
MS_API const char* ms_find_formulas(double mass, double tol_ppm, const char* limits, int max_results, int rdbe_check);

/* ------------------------------------------------------------- UniDec */
/* Runs UniDec's engine on a spectrum (since 3.1 inside the library, see
 * ms_unidec_engine; engine_exe is not used and may be NULL; work_dir and name
 * are used only with write_files=1, which writes the engine input, conf.dat
 * and the engine outputs into work_dir/<name>_unidecfiles), takes mass, fit,
 * peaks and the mass x charge grid, and scales the
 * intensities as the Python code does. params: a "key=value\n" string with
 * the keys of ms_deconv.unidec_deconvolute (z_lo, z_hi, mass_lo, mass_hi,
 * mass_step, peak_width, peak_window, peak_threshold, adduct_mass, numit,
 * zzsig, psig, beta, psfun, msig, isotopemode, poolflag, smooth, subbuff,
 * subtype, peaknorm, binning (auto|none|linear|resolution), mzbins,
 * mz_lo, mz_hi, min_intensity, min_intensity_pct, baseline, baseline_width). */
typedef struct ms_unidec_result ms_unidec_result;
MS_API int ms_unidec(const double* mz, const double* it, long n, const char* engine_exe, const char* work_dir,
                     const char* name, const char* params, ms_unidec_result** out, ms_progress_fn cb, void* user);
MS_API void ms_unidec_free(ms_unidec_result* r);
MS_API long ms_unidec_mass(const ms_unidec_result* r, const double** mass, const double** it);
MS_API long ms_unidec_fit(const ms_unidec_result* r, const double** mz, const double** it);
MS_API long ms_unidec_zdist(const ms_unidec_result* r, const double** z, const double** it);
MS_API long ms_unidec_peaks(const ms_unidec_result* r, const double** mass, const double** height, const double** area,
                            const double** score, const double** apex, const int** n_iso);
MS_API double ms_unidec_r2(const ms_unidec_result* r);
MS_API double ms_unidec_uniscore(const ms_unidec_result* r);
MS_API const char* ms_unidec_notes(const ms_unidec_result* r);
MS_API const char* ms_unidec_folder(const ms_unidec_result* r);

/* IsoDec through isodeclib.dll (UniDec's C library): monoisotopic masses from
 * isotope distributions. params as for ms_unidec plus phaseres, matchtol,
 * minpeaks, css_thresh, knockdown_rounds, datathreshold, window, thresh. */
typedef struct ms_isodec_result ms_isodec_result;
MS_API int ms_isodec(const double* mz, const double* it, long n, const char* lib_dir, const char* params,
                     ms_isodec_result** out);
MS_API void ms_isodec_free(ms_isodec_result* r);
MS_API long ms_isodec_peaks(const ms_isodec_result* r, const double** mono, const double** height,
                            const double** mz_apex, const int** z, const char*** charges);
MS_API long ms_isodec_mass(const ms_isodec_result* r, const double** mass, const double** it);
MS_API const char* ms_isodec_notes(const ms_isodec_result* r);

/* ---- 2.8 additions (agent U: UniDec reliability) ---------------------------------------
 * Reliability of a UniDec result from UniDec's own scores (UniScore and the DScore of each peak,
 * Kostelic and Marty, Methods Mol. Biol. 2022, 2500, 159-180): returns the score (the UniScore,
 * 0 to 1); level: 0 good, 1 fair, 2 poor; text: one or two plain sentences on what is wrong and what
 * to try (owned by the result). The DScore of each peak is the score of ms_unidec_peaks. */
MS_API double ms_unidec_quality(const ms_unidec_result* r, int* level, const char** text);
/* ---- end of the 2.8 additions (agent U) ---- */

/* Since 2.8: flags of the scans of an event, in the order of ms_event_scans:
 * bit 0 = saturated detector readings were left out (Shimadzu .lcd). The array
 * belongs to the calling thread (valid until its next call of this function).
 * Returns the number of scans (0 for no such event). */
MS_API long ms_event_flags(const ms_file* f, int event, const unsigned char** flags);

/* ---- 3.1 additions (agent E: the UniDec engine in the library) -----------------------------
 * ms_unidec runs UniDec's engine inside this library: a C++ port of the one dimensional engine of
 * unidec.exe (UniDec 8.2.1, BSD license with a citation requirement), OpenMP on all logical processors
 * (params "threads" sets the count), progress through the callback (10 to 90 of 100 while it iterates),
 * cancellation within about 50 ms, every engine stop an error code with a plain message. Returns a
 * short description of the engine with its credit; its presence tells a caller that ms_unidec needs no
 * unidec.exe. */
MS_API const char* ms_unidec_engine(void);
/* ---- end of the 3.1 additions (agent E) ---- */

/* ---- 3.1 additions (agent P): the remaining computations of the windows ----------------------
 * Each function mirrors a Python function of MSpektra (named in its comment) and gives its
 * results; the Python code is the fallback when the library is off. Arrays handed out are thread
 * local scratch of the function (valid until the next call of the same function in the thread).
 * Functions that return a status give MS_NOT_FOUND where the Python function returns None. */

/* hrms_calib.Calibration: as ms_calibration_fit (same models and coefficient layout), but order < 0
 * chooses the HPC order by cross validation and order >= 0 (also 0) fixes it; cv_err (n values, may
 * be NULL) receives each calibrant predicted by the fit without it (ppm; NaN without a cross
 * validation, i.e. n <= the number of parameters or a fit failed); cv_rms NaN then too. */
MS_API int ms_calibration_fit2(int model, int order, const double* measured, const double* reference, long n,
                               double* coef, int* ncoef, double* err_ppm, double* rms, double* cv_err, double* cv_rms);
/* stride (ms_peak_position, ms_find_calibrants2, ms_measured_pattern): the values of mz and it are
 * stride doubles apart (2: the two columns of an n x 2 array, without a copy; 1: plain arrays). */
/* hrms_calib._peak_position: the tallest point between lo and hi (profile: the apex of a Gaussian
 * fit of its peak) and its height. */
MS_API int ms_peak_position(const double* mz, const double* it, long n, long stride, double lo, double hi, int profile,
                            double* pos, double* height);
/* hrms_calib.peak_near: of the peak maxima within target +- d at least min_rel of the tallest
 * there, the one nearest to target (mz sorted). */
MS_API int ms_peak_near(const double* mz, const double* it, long n, double target, double d, int profile,
                        double min_rel, double* pos, double* height);
/* hrms_calib.profile_apex: apex of the profile peak whose highest point is i (NaN out of range). */
MS_API double ms_profile_apex(const double* mz, const double* it, long n, long i);
/* hrms_calib.find_calibrants, the search of each reference within tol_ppm (the caller makes the
 * rows and gives a peak found by two references to the nearer one): measured m/z, height, height
 * in % of the tallest point and the found flag (0: absent, outside the spectrum or below min_rel). */
MS_API int ms_find_calibrants2(const double* mz, const double* it, long n, long stride, int profile,
                               const double* reference, long nref, double tol_ppm, double min_rel, double* measured,
                               double* height, double* relint, int* found);

/* flags of the formula functions: bit 0 = the caller's Python sums floats with compensation
 * (sum() of Python 3.12 and later), which ms_formula's mono_mass and finder use. */
#define MS_PY_SUM12 1
/* ms_formula.isotope_pattern of a formula whose elements come in the order given (it sets the order
 * of the sums): neutral masses and relative intensities (%, the largest 100) of the peaks of at least
 * min_rel x 100 %, merged into the nominal isotope peaks (merge <= 0) or merging isotopologues closer
 * than merge (Da). */
MS_API int ms_isotope_pattern2(const char* formula, double merge, double min_rel, int flags, const double** mass,
                               const double** rel, long* n);
/* ms_formula.find_formulas, the enumeration of the candidates (the isotope ranking and the scores are
 * the caller's): elements = symbols separated by spaces in the order of the limits (with C and H),
 * lo / hi their counts; neutral = mass of M looked for, tol its tolerance (Da), mz the ion m/z, zabs
 * and electrons as in ms_formula.ADDUCTS, add / rem the formulas the ion adds to and removes from M,
 * ion_itself: [M]+ / [M]- (half integer RDB allowed). Per candidate: ncols counts (C, H, then each
 * element with hi > 0 other than C and H, in order), the ion m/z, the error (ppm) and the RDB. */
MS_API int ms_find_formulas2(double mz, double neutral, double tol, int zabs, double electrons, const char* elements,
                             const long* lo, const long* hi, const char* add, const char* rem, int ion_itself, int flags,
                             const long** counts, const double** theo, const double** err, const double** rdbe,
                             long* ncand, int* ncols);
/* ms_formula.isotope_match: 0 to 100 (0 on bad arguments). */
MS_API double ms_isotope_match(const double* pat_mz, const double* pat_rel, long npat, const double* obs_mz,
                               const double* obs_rel, long nobs);
/* ms_formula.measured_pattern: m/z and height of the tallest point near M, M + 1.00335 / z, ... (npk values). */
MS_API int ms_measured_pattern(const double* mz, const double* it, long n, long stride, double mono_mz, int z, int npk,
                               double ppm, double* out_mz, double* out_h);

/* LC integration with the indexes of the peak borders (lcms_integrate peak dicts). */
typedef struct { double t0, t1, rt, height, area, base0, base1, apex_y; long i0, i1; } ms_lc_peak2;
/* lcms_integrate.smooth: Savitzky-Golay (points < 3 or n < points + 2: a copy). */
MS_API int ms_smooth(const double* y, long n, int points, double* out);
/* lcms_integrate.auto_integrate (peaks in the order made; finish() sorts them and adds the area %). */
MS_API int ms_lc_integrate(const double* t, const double* y, long n, double threshold_pct, double min_width_s,
                           int smooth_points, int use_range, double t_lo, double t_hi, double baseline_window_min,
                           const ms_lc_peak2** peaks, long* np);
/* lcms_integrate.manual_peak, peak_at and split_peak (out2: the two parts). */
MS_API int ms_lc_manual(const double* t, const double* y, long n, double t0, double t1, ms_lc_peak2* out);
MS_API int ms_lc_peak_at(const double* t, const double* y, long n, double x, double search_s, double baseline_window_min,
                         ms_lc_peak2* out);
MS_API int ms_lc_split(const double* t, const double* y, long n, long i0, long i1, double b0, double b1, double x,
                       ms_lc_peak2* out2);

/* lcms_pda.PDAData on its float32 absorbance matrix a (nt x nw, row major): chromatogram at wl0 averaged
 * over the band bw (nm), the max plot, and the spectrum at t0 (t1 NaN) or over t0 .. t1, minus a
 * background (bg_mode 0 none, 1 the spectrum nearest to b0, 2 over b0 .. b1). */
MS_API int ms_pda_chromatogram(const float* a, long nt, long nw, const double* wl, double wl0, double bw, float* out);
MS_API int ms_pda_max_plot(const float* a, long nt, long nw, float* out);
MS_API int ms_pda_spectrum(const float* a, long nt, long nw, const double* times, double interval_s,
                           double t0, double t1, int bg_mode, double b0, double b1, double* out);

/* ms_deconv helpers: isotope_peaks (label positions of a resolved mass spectrum), refine_peak_mass,
 * refine_species of one species (out[6]: checked 1/0, apex, apex2 (NaN: none), apex2_rel, source
 * 0 deconvolution / 1 m/z data, the isotope taken -1/0/1; nd = 0: no m/z data), group_isotopes with
 * the first and last isotope peak of each species, _restrict, peak_mask, subtract_baseline (width
 * 0: factor x the peak width) and the label maxima of unilcms.label_maxima (sep NaN: a fortieth of
 * the span; idx holds nlab values). */
MS_API int ms_isotope_peaks(const double* mass, const double* it, long n, double threshold, long max_n,
                            const double** apex, const double** height, long* np);
MS_API double ms_refine_peak_mass(const double* mass, const double* it, long n, double m0);
MS_API int ms_refine_species(const double* mass, const double* it, long n, const double* mz, const double* mz_it, long nd,
                             double carrier, int z_lo, int z_hi, double spacing, double apex, double* out);
MS_API int ms_group_isotopes2(const double* mass, const double* it, long n, double threshold, double spacing,
                              const double** avg, const double** apex, const double** height, const double** area,
                              const int** n_iso, const double** first, const double** last, long* np);
MS_API int ms_restrict(const double* mz, const double* it, long n, double lo, double hi, const double** omz,
                       const double** oit, long* nout);
MS_API int ms_peak_mask(const double* mz, const double* it, long n, double thr, unsigned char* mask, long* npeaks);
MS_API int ms_subtract_baseline(const double* mz, const double* it, long n, double width, double factor, double* out);
MS_API int ms_label_maxima(const double* x, const double* y, long n, int nlab, double sep, double rel, long* idx, long* nout);
/* Mass chromatograms of nw windows (mz[w] +- tol[w]) of one event in one pass over its scans: the
 * values of nw calls of ms_xic, window after window in out (n_out = the number of scans of the event,
 * as ms_event_scans). */
MS_API int ms_xic_multi(ms_file* f, int event, const double* mz, const double* tol, int nw, double* out, long n_out);
/* ---- end of the 3.1 additions (agent P) ---- */

/* ---- 3.31: the Compare view of LCMS Analysis (lcms_compare_core.py, src/compare.cpp) ----
 * y minus a blank run (bt rising, by) interpolated at t as np.interp; nothing subtracted outside the
 * blank's time range or where the blank is not finite (nb < 2: y unchanged). */
MS_API int ms_cmp_subtract_blank(const double* t, const double* y, long n, const double* bt, const double* by,
                                 long nb, double* out);
/* Time of the highest peak top of y (smoothed over smooth_points first; < 3: not) within tc +- win, refined
 * by a parabola (find_apex); *apex = NaN when there is no peak there. */
MS_API int ms_cmp_find_apex(const double* t, const double* y, long n, double tc, double win, int smooth_points,
                            double* apex);
/* y minus a baseline: kind 0 none, 1 offset (lowest point at zero), 2 drift (line between the medians of the
 * first and last 3 %), 3 rolling minimum over rolling_min minutes (ms_subtract_baseline). */
MS_API int ms_cmp_baseline(const double* t, const double* y, long n, int kind, double rolling_min, double* out);
typedef struct {
    int smooth;          /* smoothing points (< 3: off) */
    double shift;        /* minutes added to the times (alignment and the shift of the file) */
    double t0, t1;       /* time window (NaN: open) */
    int baseline;        /* kind of ms_cmp_baseline */
    double rolling_min;  /* window of the rolling minimum (min) */
    int scale;           /* 0 same scale, 1 each to its tallest peak (100), 2 each to its peak at ref_t (100) */
    double ref_t;        /* reference time (NaN: none) */
    double ref_win;      /* half width of the reference window (min, at least 0.02) */
} ms_cmp_opts;
/* One trace processed (process): smoothing, shift, time window, baseline, scale. t_out and y_out hold n
 * values, *n_out of them are used; y drawn = factor x signal; top = highest point (NaN ignored); ref_missing:
 * no signal at the reference time (scale 2). */
MS_API int ms_cmp_process(const double* t, const double* y, long n, const ms_cmp_opts* o, double* t_out,
                          double* y_out, long* n_out, double* factor, double* top, int* ref_missing);
/* Offsets (dx, dy) of n traces (stack): their tops and time spans; ok[i] = 0: left out of the height and
 * span. stacked = 0 (overlay): all zero. The first trace on top (reverse: at the bottom). */
MS_API int ms_cmp_stack(const double* tops, const double* spans, const int* ok, long n, double spacing,
                        int reverse, double skew, int stacked, double* dx, double* dy);
/* Peaks of a trace as drawn (y / factor, auto integration at thr_pct, 2 s minimum width), sorted by time,
 * and the numbers of its row in the table: index of the tallest peak (-1: none), its area % and the area %
 * of the peak at each of the ng guide times (NaN: none). The peaks stay valid until the next call. */
MS_API int ms_cmp_peaks(const double* t, const double* y, long n, double factor, double thr_pct,
                        const double* guides, long ng, const ms_lc_peak2** peaks, long* np, long* main,
                        double* main_pct, double* guide_pct);
/* Absorption maximum (nm, above 210) of the tallest peak of a PDA run: the peak of the chromatogram at wl0
 * (NaN: the max plot) smoothed over max(5, smooth) points, in the window t0 .. t1 of the shifted times (t0
 * NaN: after the first 8 % of the run), its spectrum minus the spectrum at the lowest point of the max plot
 * 0.15 to 1 min before it. MS_NOT_FOUND: no such peak or fewer than 3 wavelengths above 210 nm. */
MS_API int ms_cmp_lambda_max(const float* a, long nt, long nw, const double* wl, const double* times,
                             double interval_s, double wl0, double bw, int smooth, double shift, double t0,
                             double t1, double* lmax, double* t_apex, int* bg_used);
/* ---- end of the 3.31 additions ---- */

/* ---- 3.65: area of a region in a trace of the Compare view (lcms_compare_core._region_peak) ----
 * The positive area (signal x s, times in min) above the straight line between the signal at max(a, t[0]) and
 * min(b, t[n-1]), with the ends interpolated and the crossings of that line inserted. t_out, y_out and base_out
 * get the *n_out points (times, signal, baseline; cap at least 2 n); res: start, end, time of the highest point
 * above the line, its height, the area, the baseline at the start and at the end. MS_NOT_FOUND: no data (n < 2,
 * a value not finite, time not rising, or the region outside the trace). */
MS_API int ms_cmp_region(const double* t, const double* y, long n, double a, double b, double* t_out,
                         double* y_out, double* base_out, long cap, long* n_out, double* res);

/* ---- 3.33: suggested neutral mass of a compound from its ESI+ and ESI- spectra (unit resolution LC-MS;
 * ms_adducts.py, src/adducts.cpp; the algorithm is described at the top of src/adducts.cpp) ----
 * Adduct codes (ms_adduct_name; a lower code is the more common adduct): 0 [M+H]+, 1 [M-H]-, 2 [M+Na]+,
 * 3 [M+NH4]+, 4 [M+HCOO]-, 5 [M+K]+, 6 [M+Cl]-, 7 [M+CH3COO]-, 8 [M+2H]2+, 9 [M-2H]2-, 10 [2M+H]+,
 * 11 [2M+Na]+, 12 [2M-H]-, 13 assumed [M+H]+, 14 assumed [M-H]- (a single ion taken as [M+H]+ / [M-H]-). */
#define MS_ADDUCT_COUNT 15
typedef struct {
    int polarity;    /* +1 or -1 */
    int adduct;      /* adduct code */
    double mz;       /* m/z of the peak (apex refined) */
    double rel;      /* its height relative to the base peak of its spectrum */
    double mass;     /* the neutral mass this ion gives with this adduct */
} ms_mass_ion;
typedef struct {
    double mass;     /* suggested neutral mass (Da) */
    double score;    /* see src/adducts.cpp */
    int both;        /* 1: ions of both polarities support it */
    int assumed;     /* 1: a single ion taken as [M+H]+ or [M-H]- (no group of two ions anywhere) */
    long first;      /* index of its first ion in the ions array */
    long n;          /* number of its ions */
} ms_mass_cand;
/* The neutral masses (best first, at most top) suggested by the peaks of a positive (pmz, pit, np) and a
 * negative (nmz, nit, nn) spectrum; either may be empty (n = 0, pointers may then be NULL). m/z need not be
 * sorted; points with a NaN or infinite value are left out, negative intensities count as zero. tol: m/z
 * tolerance (0 < tol <= 5), min_rel: peaks at least this fraction of the base peak (0 to 1), max_ions: the
 * strongest peaks used per spectrum (1 to 1000), top: 1 to 1000. Results in thread local buffers valid until
 * the next call in the same thread; *ncands = 0 when there is nothing. */
MS_API int ms_neutral_masses(const double* pmz, const double* pit, long np, const double* nmz, const double* nit,
                             long nn, double tol, double min_rel, int max_ions, int top, const ms_mass_ion** ions,
                             long* nions, const ms_mass_cand** cands, long* ncands);
/* Name of an adduct code ("" for an unknown code). */
MS_API const char* ms_adduct_name(int code);
/* The same with everything the rules of src/adducts.cpp find (mobile phase additives, M+2 partners, base
 * peaks): the candidates of ms_neutral_masses in an extended struct; the ions of a candidate's M+2 partner
 * follow in the ions array (m2_first, m2_n). base_mz[0] / [1]: apex m/z of the strongest peak of the
 * positive / negative spectrum (NaN: none); explained[0] / [1]: 1 when the best candidate explains it (one
 * of its ions or of its M+2 partner's ions). ms_neutral_masses returns the same candidates. */
#define MS_ADDITIVE_COUNT 4
typedef struct {
    double mass, score;  /* as ms_mass_cand */
    int both, assumed;
    long first, n;
    int additive;        /* -1, or the mobile phase additive the mass matches (ms_additive_name) */
    double m2;           /* mass of its M+2 partner (Br, Cl), NaN: none */
    double m2_ratio;     /* height of the partner's strongest ion / height of its counterpart ion here */
    long m2_first, m2_n; /* the partner's ions in the ions array */
} ms_mass_cand2;
MS_API int ms_neutral_masses2(const double* pmz, const double* pit, long np, const double* nmz, const double* nit,
                              long nn, double tol, double min_rel, int max_ions, int top, const ms_mass_ion** ions,
                              long* nions, const ms_mass_cand2** cands, long* ncands, double* base_mz, int* explained);
/* Name of a mobile phase additive code: 0 TFA, 1 formic acid, 2 acetic acid, 3 DFA ("" for an unknown code). */
MS_API const char* ms_additive_name(int code);
/* ---- end of the 3.33 additions ---- */

/* ---- 3.35: mass shift finder on the masses of a deconvolution result (ms_shifts.py, src/shifts.cpp; the
 * algorithm is described at the top of src/shifts.cpp) ----
 * Which masses differ by one of the shifts (signed: a loss is negative), or by k x one of the tags plus up
 * to max_extra (0 to 2) shifts, within tol_da + tol_ppm x 1e-6 x the heavier mass; the differences from the
 * reference (ref: index; -1 automatic: the tallest, or with tags the lightest species of a tag series below
 * the tallest, whose tag counts are all <= 0, see step 11 of src/shifts.cpp) that nothing explains; the composition of every mass relative to
 * the reference; the degree of conjugation per tag. Masses below min_rel (0 to 1) x the tallest are left
 * out. Flags: MS_SHIFT_ISOTOPES (isotope resolved masses: an
 * explanation may differ by one isotope, 13C - 12C, in either direction), MS_SHIFT_COLLAPSE (the list holds
 * the isotope peaks of each species: a mass 1 or 2 isotopes above another one is its isotope peak). */
#define MS_SHIFT_ISOTOPES 1
#define MS_SHIFT_COLLAPSE 2
typedef struct {
    long from, to;      /* mass indexes (input order): the parent and the child (see src/shifts.cpp) */
    double delta;       /* mass[to] - mass[from] */
    double value;       /* the explanation (Da); NaN: unknown */
    double error;       /* delta - value; NaN: unknown */
    int kind;           /* 0 a shift, 1 tags (with up to max_extra shifts), 2 unknown (pairs with the reference
                           only), 3 isotope peak (MS_SHIFT_COLLAPSE) */
    int sign;           /* +1 as written, -1 reversed (the parent carries it) */
    int tag, k;         /* tag index (-1: none) and the number of tags */
    int s1, s2;         /* shift indexes (-1: none), s1 <= s2 */
    int iso;            /* isotope offset -1, 0, +1 (kind 3: the isotope peak's number) */
    int n_alt;          /* other explanations of another value (more than 1e-6 Da apart) within the tolerance */
    int alt_tag, alt_k, alt_s1, alt_s2, alt_sign, alt_iso;   /* the next best of them (alt_tag, alt_s1 -1 and
                           alt_sign 0 when n_alt is 0) */
    double alt_error;
    int ref;            /* 1: a pair with the reference */
    int link;           /* 1: drawn by default (the link of the child to its parent, see src/shifts.cpp) */
} ms_shift_pair;
typedef struct {
    int status;         /* 0 reference, 1 explained by its pair with the reference, 2 explained through
                           another mass, 3 unknown, 4 isotope peak of another mass (MS_SHIFT_COLLAPSE),
                           5 below min_rel (left out) */
    long parent;        /* the mass its composition comes through (status 4: the first mass of its chain), -1 none */
    long link;          /* index of its link in the pairs, -1 none */
    double value, error;   /* its composition (Da) and mass - reference - value (NaN: unknown) */
    double height, area;   /* as counted (MS_SHIFT_COLLAPSE: the largest height, the summed area of its chain) */
} ms_shift_species;
/* n: 0 to 2000 masses (positive, finite); heights and areas that are not finite or negative count as 0 (area
 * may be NULL); ns 0 to 64 shifts (not 0, |shift| <= 100000), nt 0 to 3 tags (0 < tag <= 1000000); kmax 1 to
 * 20. Results in thread local buffers valid until the next call in the same thread: npairs pairs, n species,
 * comp: n rows of nt + ns + 1 counts (each tag, each shift, the isotope offset), conj: nt rows of 2 kmax + 5
 * values (% by height for k = 0 .. kmax, % by area for k = 0 .. kmax, tags per molecule by height, by area,
 * the number of masses counted; NaN without data); ref_used: the reference (-1 for n = 0; with ref -1 the one
 * chosen). */
MS_API int ms_mass_shifts(const double* mass, const double* height, const double* area, long n,
                          const double* shifts, long ns, const double* tags, long nt, double tol_da, double tol_ppm,
                          long ref, int kmax, int max_extra, int flags, double min_rel, const ms_shift_pair** pairs,
                          long* npairs, const ms_shift_species** species, const int** comp, const double** conj,
                          long* ref_used);
/* ---- end of the 3.35 additions ---- */

#ifdef __cplusplus
}
#endif
#endif
