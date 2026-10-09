#pragma once
// ABI of the preserved msengine.dll, reconstructed from its supplied ctypes bridge.
// Never append fields to this structure; additional options belong to mx_run arguments.
using Progress = int (*)(void*,long,long);
struct Params {
    int z_lo,z_hi; double mass_lo,mass_hi,mass_step,adduct_mass,mz_lo,mz_hi,min_intensity;
    int min_intensity_pct; double peak_width,resolution;
    int resolved_isotopes,rounds,min_rounds,iterations; double peak_window,peak_threshold;
    int baseline; double baseline_width; int envelope_nodes,noise_model; double chi2_target;
};
static_assert(sizeof(Params)==152,"Installed engine parameter ABI changed");
