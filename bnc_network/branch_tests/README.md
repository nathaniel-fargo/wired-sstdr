# Branch Tests: Reproducing Depth-of-Detection Analysis

This folder contains the pipeline that demonstrates SSTDR sensitivity vs. branch depth in a simple power grid testbed. The end result is the figure `Analysis/networks_analysis_12mhz_vs_networks_cal_analysis_12mhz_AreaNorm_intersection.png` which shows the exponential trend with an intersection beyond 18 branches deep.

## Data Layout

- `2025-05-30/LiveWire/LWS/` — Raw LiveWire `.lws` acquisitions, one per network ID
- `2025-05-30/LiveWire/CSV/` — Auto-generated `.csv` files parsed from `.lws`
- `2025-05-30/LiveWire/Plots/` — Waveform plots for quick visual inspection
- `2025-05-30/networks.csv` — Network catalogue; first column is the ID used to match LiveWire files
- `2025-05-30/Analysis/*.csv` — Pairwise comparison metrics and summary tables
- `scripts/` — MATLAB scripts to run the pipeline end-to-end

## One-Command Ingestion (LWS → CSV → Plots)

Use the unified utility in `livewire_scripts/`:

```
process_lws_folder('.../2025-05-30/LiveWire/Static', true, 4);
```

This will:
- Rename `Static/` → `LWS/` (idempotent)
- Convert all `.lws` → `.csv` with FFT interpolation (files saved to `../CSV/`)
- Generate waveform plots into `../Plots/`

The CSVs use a standard schema with waveform rows labeled as `DataType == "Waveform"` and include key fields like `MeasurementFrequency`, `UnitsPerSample`, and `ZeroIndex` for distance conversion in feet.

## Pairwise Network Comparison (metrics vs. depth)

Run comparisons to produce metrics for network pairs differing by up to one wire:

```
scripts/compare_networks( ...
  '2025-05-30/networks.csv', ...
  '2025-05-30/LiveWire/CSV', ...
  1, ...              % max difference level
  false, ...          % generate overlay plots per pair (optional, set true to save)
  "12 MHz" ...       % fixed frequency for analysis (recommended)
);
```

This writes `2025-05-30/Analysis/networks_analysis_12mhz.csv` (and `networks_cal_analysis_12mhz.csv` for calibration sets when run on calibration data). Columns include:
- `PeakMag`, `AreaNorm`, `AreaSquared` — difference metrics
- `NumWireDiff` — how many wires differ between the two networks
- `MinDiffDepth` — minimum depth at which the differing wires occur

## Trend vs Baseline and Intersection Plot

Create the “trend vs error level” figure and intersection point (branch-depth capability):

```
scripts/plot_data_with_baseline( ...
  '2025-05-30/Analysis/networks_analysis_12mhz.csv', ...
  '2025-05-30/Analysis/networks_cal_analysis_12mhz.csv', ...
  'AreaNorm' ...  % metric used in the published figure
);
```

This produces `2025-05-30/Analysis/networks_analysis_12mhz_vs_networks_cal_analysis_12mhz_AreaNorm_intersection.png` with:
- Blue points: metric vs. branch depth
- Red dashed curve: exponential fit, `y = a e^{bx}`
- Solid blue line: baseline (mean error level from calibration set)
- Black marker: intersection depth where trend crosses baseline (≥ 18 branches in our data)

## Helpful Scripts

- `scripts/plot_data_vs_depth.m` — Quick scatter + exponential fit for any metric vs. depth
- `scripts/plot_csv_difference.m` — Overlay two CSVs and show metrics
- `scripts/display_network_results.m` — Combine network diagram and waveform for a given ID
- `scripts/draw_network.m` — Draws the network tree with cable lengths (reads `cable_measurements/cable_logs.csv`)

## Notes

- The LiveWire processing is robust to older CSVs that included `VOP`; it is ignored for plotting.
- All LiveWire conversions are interpolated; the `_smooth` suffix is no longer used.

# Branch Tests

## Syntax Overview

Today we're starting to build larger branch networks and analyze the data that comes out of it. This file explains the notation used in files such as `./2025-05-21/branch_networks.csv`

Syntax is

`Network = { Wire1 { Wire2[Termination], Subnetwork } }`

For example this network

`{WIRE_A{WIRE_B[T],WIRE_C{WIRE_D[S],WIRE_E{WIRE_F[O]}}}}`

maps to the corresponding network

```ASCII
-> WIRE A -|-> WIRE_C -|-> WIRE_E -> WIRE_F -> [open]
	 |	         |
  WIRE_B 	   WIRE_D
	 |	         |
[terminated]  [short]
```

Where [X] is the termination type, which can be any of the following:

* O: Open
* S: Short
* T: Terminated
* R=_: Resistor value
* L=_: Inductor value
* C=_: Capacitance value

Series elements are not supported, and every unattached wire must specify it's ending type

A rib looks like

`{A00{A01{A02{A03{A04{B05[O],B04[O]},B03[O]},B02[O]},B01[O]},B00[O]}}`

or this,

`{A00{B00[O],A01{B01[O],A02{B02[O],A03{B03[O],A04{B04[O],B05[O]}}}}}}`

or this,

```
{
	A00 {
		B00[O],
		A01 {
			B01[O],
			A02 {
				B02[O],
				A03 {
					B03[O],
					A04 {
						B04[O],
						B05[O]
					}
				}
			}
		}
	}
}
```

They're all the same network

This format should be easily parseable by code for further analysis.

# Rules

* The top level network has to be contained within brackets `{}`
* Every wire must have a connection at its end
  * Other wires (1 or 2)
  * Termination type
* Each wire must have it's own .lws recording for proper simulation
