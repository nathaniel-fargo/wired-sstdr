# SSTDR Dataset Visualization Tool

## Overview

The `visualize_sstdr_dataset.m` function provides comprehensive visualization and analysis capabilities for generated SSTDR datasets. It creates detailed plots showing network topology, signal analysis, correlation results, and statistical summaries.

## Features

### 📊 **Dataset Statistics**
- Network count and complexity metrics
- Segment distribution analysis
- Physical characteristics (length, impedance)
- Fault distribution and types
- Performance metrics

### 📈 **Signal Visualization**
- TX and RX signal plots
- Time-domain signal analysis
- Signal quality assessment

### 🔗 **Network Topology**
- Visual network diagrams
- Fault location and type visualization
- Segment and termination display
- Source and load representation

### 📡 **Cross-Correlation Analysis**
- Correlation magnitude plots
- Peak detection and marking
- Detailed correlation breakdown
- Configuration parameter display

### 💾 **Export Capabilities**
- Save plots in multiple formats (PNG, PDF, FIG)
- Batch processing for multiple networks
- Individual network analysis

## Usage Examples

### Basic Dataset Overview
```matlab
% Quick overview of entire dataset
[stats, figs] = visualize_sstdr_dataset('datasets/my_dataset_2025-07-17');
```

### Specific Network Analysis
```matlab
% Detailed analysis of network 3
[stats, figs] = visualize_sstdr_dataset('datasets/my_dataset_2025-07-17', ...
    'network_id', 3, ...
    'save_plots', true, ...
    'plot_format', 'png');
```

### Batch Visualization
```matlab
% Visualize all networks in small dataset
[stats, figs] = visualize_sstdr_dataset('datasets/my_dataset_2025-07-17', ...
    'show_all', true, ...
    'verbose', false);
```

### Direct File Analysis
```matlab
% Analyze specific network file
[stats, figs] = visualize_sstdr_dataset('datasets/my_dataset/network_001/data.mat');
```

## Function Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dataset_path` | string | *required* | Path to dataset folder or network .mat file |
| `network_id` | integer | 1 | Specific network ID to visualize |
| `show_all` | logical | false | Show all networks in dataset |
| `save_plots` | logical | false | Save plots to files |
| `plot_format` | string | 'png' | Format for saved plots ('png', 'pdf', 'fig') |
| `show_summary` | logical | true | Display dataset summary statistics |
| `verbose` | logical | true | Enable verbose output |

## Output Structure

### Statistics Structure
```matlab
stats = struct(
    'dataset_folder',     % Path to dataset
    'num_networks',       % Total number of networks
    'segments',           % Segment count statistics (min, max, mean, std)
    'length',             % Network length statistics  
    'faults',             % Fault count and distribution
    'fault_types',        % Series vs shunt fault counts
    'termination'         % Termination impedance statistics
);
```

### Figure Handles
Returns array of MATLAB figure handles for further manipulation or saving.

## Visualization Components

### Main Overview Figure (6 subplots)
1. **Network Topology** - Visual diagram with faults and components
2. **Network Statistics** - Text summary of key parameters
3. **SSTDR Signals** - TX/RX signal plots
4. **Cross-Correlation** - Correlation magnitude with peaks
5. **Fault Analysis** - Detailed fault location and type information

### Detailed Correlation Figure (4 subplots)
1. **Full Correlation** - Real and imaginary components
2. **Magnitude Spectrum** - Correlation magnitude with peak marking
3. **Peak Analysis** - Detected peaks visualization
4. **Configuration Info** - Correlation parameters and settings

## Network Diagram Symbols

| Symbol | Meaning |
|--------|---------|
| Green circle (SSTDR) | Signal source |
| Black line | Transmission line |
| Red zigzag | Series fault (resistance) |
| Red triangle | Shunt fault (conductance) |
| Black square | Termination load |
| White circles | Segment boundaries |

## Tips and Best Practices

### For Large Datasets (>10 networks)
- Set `verbose=false` to reduce output
- Use `show_summary=true` for overview without individual plots
- Process networks individually rather than `show_all=true`

### For Analysis and Documentation
- Use `save_plots=true` to preserve visualizations
- Choose appropriate `plot_format` for your needs:
  - PNG: Web/presentation use
  - PDF: Publication quality
  - FIG: MATLAB editing

### For Interactive Analysis
- Use returned figure handles for customization
- Modify plots after generation for specific requirements
- Combine with other MATLAB plotting functions

## Example Workflow

```matlab
% 1. Generate dataset
dataset_config = create_dataset_config('num_networks', 10);
[summary, failed] = generate_sstdr_dataset(dataset_config);

% 2. Quick overview
[stats, ~] = visualize_sstdr_dataset(summary.dataset_folder);

% 3. Detailed analysis of interesting networks
interesting_networks = [1, 5, 8];  % Networks with specific characteristics
for net_id = interesting_networks
    [~, figs] = visualize_sstdr_dataset(summary.dataset_folder, ...
        'network_id', net_id, ...
        'save_plots', true, ...
        'plot_format', 'pdf');
end

% 4. Generate summary report
fprintf('Dataset Analysis Summary:\n');
fprintf('Networks: %d\n', stats.num_networks);
fprintf('Average faults per network: %.1f\n', stats.faults.per_network_mean);
fprintf('Fault distribution: %d series, %d shunt\n', ...
    stats.fault_types.series, stats.fault_types.shunt);
```

## Integration with Dataset Generation

The visualization tool is designed to work seamlessly with the dataset generation system:

- Reads data saved by `save_sstdr_dataset.m`
- Compatible with all dataset configurations
- Handles various network sizes and fault patterns
- Works with different simulation parameters

## Troubleshooting

### Common Issues

**"No networks found"**: Check dataset folder structure
```
dataset_folder/
├── network_001/
│   └── sstdr_data_*.mat
├── network_002/
│   └── sstdr_data_*.mat
└── dataset_summary.mat
```

**"Network ID not found"**: Use `show_summary=true` to see available networks

**Empty plots**: Verify dataset contains required fields (network_config, signals, correlation_results)

### Performance Considerations

- Large datasets: Process incrementally
- Memory usage: Close figures when done (`close(fig_handles)`)
- Batch processing: Use `verbose=false` to reduce overhead

## Related Functions

- `generate_sstdr_dataset.m` - Dataset generation
- `build_and_simulate_network.m` - Individual network simulation  
- `save_sstdr_dataset.m` - Data persistence
- `create_dataset_config.m` - Configuration management 