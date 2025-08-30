Overview of CRES.ipynb:
---------
This project processes cosmic ray air shower data from the LOFAR Radboud Air
Shower Array (LORA) and compares reconstructed energy spectra to the
Pierre Auger Observatory reference. The workflow spans multiple years and
seasons (May & November) to assess spectral stability and systematic offsets.

Workflow Summary: .dat files are used for input data here, .nur files are more suitable for NuRadio framework
-----------------
1. Event Data Processing (`process_lora_data`):
   - Read and parse LORA .dat files.
   - Extract zenith angle, charged particle number, and energy.
   - Apply quality cuts: zenith ≤ 35°, positive energy & particle counts.
   - Transform to log-space: log10(Ne) and log10(E).

2. Energy Reconstruction and Histogramming:
   - Fit log10(E) vs log10(Ne) using a linear model to reconstruct energy.
   - Compute uncertainties from fit parameters.
   - Bin reconstructed energies and calculate Poisson errors.
   - Produce histograms for May and November of selected years.

3. Auger Comparison (`auger_ratio.py`):
   - Load Auger reference flux and convert to GeV² units.
   - Apply LORA detector acceptance (m²·sr) from Dr. S. Thoudam’s paper:
     "Measurement of the cosmic-ray energy spectrum above 10^16 eV with the 
     LOFAR Radboud Air Shower Array."
   - Compute differential flux and E³·dI/dE for each year/month.
   - Plot LORA spectra vs Auger reference and calculate LORA/Auger ratio per bin.
   - Generate bin-averaged ratio summaries with statistical uncertainties.

4. Temporal Residual Analysis:
   - Compile mean ratios (⟨R⟩) and uncertainties (δR(stat)) from 2011–2021.
   - Compute residuals (⟨R⟩ − 1) to visualize systematic deviations.
   - Plot residuals vs. year for May (blue circles) and November (red squares)
     with error bars for statistical uncertainties.

Outputs:
--------
- Histograms of reconstructed energies per year/month.
- LORA vs Auger spectra with E³·dI/dE scaling.
- Bin-averaged LORA/Auger ratio summaries.
- Residual plot showing temporal stability of the flux ratio.
- Tabulated numerical summaries for further analysis.

Dependencies:
-------------
- numpy
- matplotlib
- scipy
- calendar
- NuRadioReco.utilities (for Auger flux reference)

Usage:
------
1. Run `process_lora_data` scripts to generate energy histograms for selected years.
2. Use `auger_ratio.py` to compute flux ratios and plot comparison to Auger data.
3. Run the residual plotting script to visualize stability over the decade.

Reference:
----------
- Thoudam, S., et al., "Measurement of the cosmic-ray energy spectrum above
  10^16 eV with the LOFAR Radboud Air Shower Array."

