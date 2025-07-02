# Quick Start Scripts for SSTDR Analysis

Three convenient scripts to get your SSTDR analysis running immediately!

## 🚀 **1. `simple_run.m` - Ultra Simple**
**Perfect for:** Quick testing with minimal setup

### How to use:
1. **Set your model name**: Change `MODEL_NAME = 'your_model_name'`
2. **Pick configuration**: Change `CONFIG_CHOICE = 1` (numbers 1-8)
3. **Hit Run!** 

### Configurations available:
```matlab
CONFIG_CHOICE = 1;  % Default (100 kHz carrier, 1 MHz sampling)
CONFIG_CHOICE = 2;  % Fast (250 kHz carrier, 2 MHz sampling)  
CONFIG_CHOICE = 3;  % High Resolution (100 kHz carrier, 5 MHz sampling)
CONFIG_CHOICE = 4;  % Unmodulated (no carrier, 1 MHz sampling)
CONFIG_CHOICE = 5;  % Custom Medium (500 kHz, 5 MHz)
CONFIG_CHOICE = 6;  % Custom High Frequency (1 MHz, 10 MHz)
CONFIG_CHOICE = 7;  % Custom Long Range (50 kHz, 2 MHz)
CONFIG_CHOICE = 8;  % Ultra High Resolution (unmodulated, 10 MHz)
```

---

## 🎛️ **2. `quick_sstdr_run.m` - Full Featured**
**Perfect for:** Comprehensive analysis with lots of options

### How to use:
1. **Set your model name**: Change `MODEL_NAME = 'your_model_name'`
2. **Uncomment ONE configuration section** (lines 18-34)
3. **Optionally uncomment additional features** (lines 100+)
4. **Hit Run!**

### Features:
- ✅ 9 different configuration options
- ✅ Automatic model configuration
- ✅ Complete error handling
- ✅ Results summary
- ✅ Optional parameter sweeps
- ✅ Performance comparisons
- ✅ Automatic result saving

### Example configurations:
```matlab
%% Uncomment ONE of these:
sstdr_config('default');                                    % Standard
% sstdr_config('sine_fast');                               % Fast
% sstdr_custom_config('carrier_freq', 500e3, 'fs', 5e6);  % Custom
% sstdr_custom_config();                                   % Interactive
```

### Optional features (uncomment to use):
```matlab
%% OPTION A: Parameter sweep different frequencies
%% OPTION B: Compare time vs frequency domain performance  
%% OPTION C: Save all results to file
```

---

## 📊 **3. `compare_configs.m` - Multi-Configuration Comparison**
**Perfect for:** Testing multiple setups and comparing results

### How to use:
1. **Set your model name**: Change `MODEL_NAME = 'your_model_name'`
2. **Hit Run!** (no other changes needed)
3. **View comparison table and plots**

### What it does:
- ✅ Tests 4 different configurations automatically
- ✅ Creates performance comparison table
- ✅ Generates side-by-side correlation plots
- ✅ Ranks configurations by speed
- ✅ Shows detected peaks for each config

### Output example:
```
=== COMPARISON SUMMARY ===
Configuration        Modulation      Carrier (kHz)   Sampling (kHz)  Time (ms)
--------------------------------------------------------------------------------
Default SSTDR        none            N/A             1000            2.1       
Fast Sine Wave       sine            250             2000            1.8       
High Resolution      sine            100             5000            4.2       
Unmodulated PN       none            N/A             1000            2.3
```

---

## 🎯 **Quick Setup for Any Script**

### Step 1: Set Your Model Name
```matlab
% In any script, change this line:
MODEL_NAME = 'my_sstdr_model';  % Your actual Simscape model name
```

### Step 2: If You Don't Have a Model Yet
```matlab
% Set to empty to use existing simulation data:
MODEL_NAME = '';  
```

### Step 3: Make Sure Your Model Has:
- **Signal Source**: `From Workspace` block with variable `pn_timeseries_interp`
- **Measurement**: `To Workspace` block saving to variable `simout`

---

## 💡 **Pro Tips**

### For Quick Testing:
```matlab
% Use simple_run.m and just change these two lines:
MODEL_NAME = 'my_model';
CONFIG_CHOICE = 3;  % Try different numbers 1-8
```

### For Exploring Different Frequencies:
```matlab
% Use quick_sstdr_run.m and try different custom configs:
sstdr_custom_config('carrier_freq', 750e3, 'fs', 8e6);
```

### For Research/Analysis:
```matlab
% Use compare_configs.m to see all options at once
% Then use quick_sstdr_run.m to dive deeper into the best option
```

### If You Get Errors:
1. **Model not found**: Check `MODEL_NAME` spelling
2. **No simulation data**: Run your Simscape model first, or set `MODEL_NAME = ''`
3. **Correlation fails**: Make sure your model saves data to variable `simout`

---

## 🚀 **Example Workflow**

```matlab
%% 1. Quick test to see if everything works
simple_run  % Uses CONFIG_CHOICE = 1 (default)

%% 2. Compare different options
compare_configs  % See all 4 configurations at once

%% 3. Deep dive into best option
% Edit quick_sstdr_run.m to use your favorite config
quick_sstdr_run  % Full analysis with all features
```

Just pick the script that matches your needs and hit run! 🎉 