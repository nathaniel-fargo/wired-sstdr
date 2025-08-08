### LiveWire Scripts

These utilities convert LiveWire `.lws` acquisitions into analysis-ready `.csv` files and generate quick-look waveform plots. The pipeline now uses a single, interpolated path by default to avoid duplicate variants and file-name suffixes.

### Key commands

- Convert a folder of `.lws` files to `.csv` and produce plots:

```
process_lws_folder('/path/to/LiveWire/Static', true, 4);
```

This renames `Static/` → `LWS/` (idempotent), writes CSVs into `../CSV/`, and plots into `../Plots/`.

- Convert a single file with a custom interpolation factor:

```
lws_to_csv('T00_0.lws', 'T00_0.csv', 8);
```

### Unification changes

- `lws_to_csv.m` contains the unified, interpolated implementation and drops the historical `_smooth` suffix altogether.
- `lws_to_csv_interp.m` has been removed.
- `plot_smooth_livewire.m` is no longer needed; use `plot_livewire_csv.m`.

### CSV schema (relevant columns)

- Global: `SerialNumber`, `CableName`, `Units`, `SelectedFrequencyAtAcquisition`, `DistanceAtAcquisition`
- Waveform rows: `MeasurementFrequency`, `DataType == "Waveform"`, `DataIndex`, `Value`, `UnitsPerSample`, `ZeroIndex`

Distances are reported in feet inside plotting/analysis utilities using `UnitsPerSample` and `ZeroIndex`.


