Overview
--------
This project analyzes air-shower events detected by the LOFAR Radboud Air shower array (LORA).
The goal is to study shower core distributions, apply quality cuts, and compute
detector efficiency maps and yearly trends.

The pipeline is organized into several scripts, each focusing on one stage of the analysis.
.dat files are used for input data here, .nur files are more suitable for NuRadio framework.
-----------------------------------------------------------------------
1. Event Extraction and Core Distribution (script 1)
-----------------------------------------------------------------------
- Reads `.dat` files containing event parameters.
- Extracts core position (x, y), energy, Nch, and rM values.
- Applies quality cuts:
    • log10(Nch) > 6.25  
    • 10 m ≤ rM ≤ 200 m  
- Bins events into predefined energy ranges.
- Saves filtered event lists to a text file.
- Plots 2D histograms of shower core positions (all vs. filtered),
  overlaid with detector layout and a 150 m fiducial circle.

-----------------------------------------------------------------------
2. Efficiency Map Construction (script 2)
-----------------------------------------------------------------------
- Builds 2D efficiency maps by comparing detected shower densities to
  an idealized uniform distribution within the fiducial circle.
- Processes May and November data separately for each year (2011–2021).
- Produces efficiency maps per year and a combined map across all years.

-----------------------------------------------------------------------
3. Yearly Mean Efficiency with Error Bars (script 3)
-----------------------------------------------------------------------
- Computes mean detection efficiency inside the 150 m radius for each year.
- Separates results into May and November datasets.
- Estimates statistical (Poisson) uncertainties.
- Produces a time series plot showing:
    • May efficiency with error bars (blue)  
    • November efficiency with error bars (orange)  
    • Combined yearly average (black)

-----------------------------------------------------------------------
4. Efficiency Map Statistics Across Years (script 4)
-----------------------------------------------------------------------
- Stacks efficiency maps from all years to compute:
    • Mean efficiency map  
    • Standard deviation map (year-to-year variations)  
- Tracks yearly mean efficiencies (May, November, combined) with uncertainties.
- Plots long-term efficiency stability with error bars.

-----------------------------------------------------------------------
Key Parameters
-----------------------------------------------------------------------
- Fiducial radius: 150 m
- Cuts applied: log10(Nch) > 6.25, 10 m ≤ rM ≤ 200 m
- Efficiency definition: ratio of detected to expected showers
  (expected = uniform distribution normalized to total counts).

-----------------------------------------------------------------------
Outputs
-----------------------------------------------------------------------
- Filtered event list (text file with core position, energy, Nch, rM).
- Per-year and combined 2D histograms of core positions.
- Yearly efficiency maps (May, November, combined).
- Time series plots of mean detection efficiency with statistical errors.
- Mean and standard deviation efficiency maps across years.

-----------------------------------------------------------------------
Purpose
-----------------------------------------------------------------------
This workflow provides a systematic way to:
- Clean and filter LORA event data.
- Visualize shower core distributions.
- Quantify detector efficiency across time.
- Track stability and variations in detection efficiency from 2011–2021.


