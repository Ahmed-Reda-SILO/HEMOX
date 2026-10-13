<div align="center">

# 🫀 HEMOX

### Synchronous Dual-Site PPG–IMU Monitoring, Calibration, and Validation Platform

[![Android App](https://img.shields.io/badge/Android-v1.3.7.1-3DDC84?style=for-the-badge&logo=android&logoColor=white)](HEMOXV1.3.7.1.apk)
[![Explorer](https://img.shields.io/badge/Validation%20Explorer-v1.7.12-1261A0?style=for-the-badge)](index.html)
[![P01 Sample](https://img.shields.io/badge/Sample-P01%20%7C%204%20runs-8A2BE2?style=for-the-badge)](samples/P01/)
[![GitHub Pages](https://img.shields.io/badge/Open-HEMOX%20Explorer-222?style=for-the-badge&logo=github)](https://ahmed-reda-silo.github.io/HEMOX/)

**[📱 Download HEMOX Android App](HEMOXV1.3.7.1.apk)**  
**[🚀 Open the HEMOX Validation Explorer](https://ahmed-reda-silo.github.io/HEMOX/)**  
**[📦 Browse the P01 Sample Dataset](samples/P01/)**  
**[🐍 Browse the Python Analysis Scripts](scripts/)**

</div>

---

## 🔬 Overview

**HEMOX** is a wearable, wireless, dual-site photoplethysmography (PPG) and inertial measurement unit (IMU) platform developed for synchronized bilateral physiological monitoring. The system is designed to acquire red/infrared PPG and local motion information from two independent sensing sites, preserve a common temporal reference, and support motion- and signal-quality-aware physiological analysis and multi-site fusion.

The software workflow combines an Android acquisition/monitoring application, a browser-based validation explorer, and Python scripts for offline analysis and figure generation. Together, these resources support synchronized comparison against a commercial pulse-oximetry comparator, calibration development, independent validation, agreement analysis, and publication-oriented reporting.

> [!IMPORTANT]
> HEMOX is a research platform. Agreement with a commercial pulse-oximetry comparator does **not** constitute formal arterial-reference pulse-oximetry validation or clinical diagnostic certification.


## 📱 HEMOX Android application

Current APK: **HEMOX v1.3.7.1**

The Android application supports the HEMOX dual-node workflow, including:

- Bluetooth Low Energy (BLE) connection to two sensing nodes
- synchronized CH1/CH2 acquisition
- real-time red/infrared PPG visualization
- six-axis IMU monitoring
- synchronized session recording
- physiological and signal-quality visualization
- post-recording inspection and export

### Install

Download **[`HEMOXV1.3.7.1.apk`](HEMOXV1.3.7.1.apk)** and install it on a compatible Android device. Android may require permission to install applications from the browser or file manager used to open the APK.

> [!NOTE]
> The APK is provided for research and prototype evaluation. It is not distributed as a certified medical application.

## 📊 HEMOX Validation Explorer

The browser-based explorer supports synchronized analysis of HEMOX recordings against the commercial SP-20 pulse-oximetry comparator. All uploaded recordings are processed locally in the browser; no backend or package installation is required.

**[Open the hosted HEMOX Explorer](https://ahmed-reda-silo.github.io/HEMOX/)**

### Main analysis capabilities

- multi-run HEMOX and SP-20 data management
- synchronized one-to-one temporal matching
- separate calibration-development and independent-validation roles
- bias, MAE, RMSE, MAPE, accuracy, correlation, and limits of agreement
- Bland–Altman and regression analysis
- CH1, CH2, and fused physiological trends
- motion and signal-quality robustness analysis
- M0–M3 calibration-model comparison with fixed coefficients
- export of synchronized pairs, metrics, figures, reports, and calibration settings

## 📦 P01 sample dataset

The updated [`samples/P01/`](samples/P01/) dataset contains **four paired recording sessions**: two stationary runs (S01, S02) and two controlled-motion runs (NS01, NS02). Each HEMOX session includes all **nine export types** listed below, together with a matching SP-20 comparator CSV. The dataset contains **36 HEMOX export files and 4 SP-20 files**.

### Paired recordings

| Run | Condition | HEMOX exports | SP-20 comparator |
|---|---|---|---|
| **S01** | Stationary | [Session folder](samples/P01/HEMOX/HEMOX_DUAL_HEMOXP01S01_20260921_095135/) | [SP20P01S01.csv](samples/P01/SP20/SP20P01S01.csv) |
| **S02** | Stationary | [Session folder](samples/P01/HEMOX/HEMOX_DUAL_HEMOXP01S02_20260921_100338/) | [SP20P01S02.csv](samples/P01/SP20/SP20P01S02.csv) |
| **NS01** | Controlled motion | [Session folder](samples/P01/HEMOX/HEMOX_DUAL_HEMOXP01NS01_20260921_101149/) | [SP20P01NS01.csv](samples/P01/SP20/SP20P01NS01.csv) |
| **NS02** | Controlled motion | [Session folder](samples/P01/HEMOX/HEMOX_DUAL_HEMOXP01NS02_20260921_101755/) | [SP20P01NS02.csv](samples/P01/SP20/SP20P01NS02.csv) |

### Files included in each HEMOX session

| File | Description |
|---|---|
| `*_RAW.csv` | synchronized bilateral red/infrared PPG and IMU recording |
| `*_PARAMETER_TRENDS.csv` | CH1, CH2, and fused physiological, quality, and motion trends |
| `*_FUSED_TRENDS.csv` | fused physiological output, quality flags, and separately labelled calibrated values |
| `*_WINDOWS.csv` | window-level estimates, signal quality, motion burden, acceptance status, and SpO₂ eligibility |
| `*_MOTION_THRESHOLDS.csv` | time-varying CH1/CH2 acceleration and gyroscope thresholds |
| `*_SESSION_SUMMARY.csv` | compact session-level physiological and research-decision summary |
| `*_METADATA.txt` | recording, acquisition, and processing configuration |
| `*_REPORT.txt` | generated session report |
| `*_RECORDING.xlsx` | consolidated spreadsheet export |

The matching `SP20P01<run>.csv` files are stored in [`samples/P01/SP20/`](samples/P01/SP20/).

> [!NOTE]
> P01 is one participant from the calibration-development cohort. These recordings demonstrate the export format and single-participant analysis workflow. Participant-disjoint independent validation requires the additional development and held-out participant recordings.

## 🧭 Typical validation workflow

1. Add the required HEMOX recording runs.
2. Assign each run to **Calibration**, **Independent validation**, or **Excluded**.
3. Upload one SP-20 CSV/TXT file for each comparator recording.
4. Upload the corresponding HEMOX validated/fused trend files and optional supporting data.
5. Confirm synchronization settings and run the pooled analysis.
6. Review physiological agreement, calibration transfer, motion, signal quality, and bilateral-channel behavior.
7. Export synchronized data, performance metrics, figures, reports, and fixed calibration parameters.

### Quick start with P01

1. Clone or download the repository and locate `samples/P01/`.
2. Open the [HEMOX Explorer](https://ahmed-reda-silo.github.io/HEMOX/) and add a run, for example **P01/S01**.
3. Upload `SP20P01S01.csv` from the `SP20/` folder and `*_FUSED_TRENDS.csv` from the matching S01 HEMOX session folder in the table above.
4. Add the same session's `PARAMETER_TRENDS`, `RAW`, `WINDOWS`, `MOTION_THRESHOLDS`, `SESSION_SUMMARY`, `METADATA`, and `REPORT` files for supporting analyses. The `RECORDING.xlsx` workbook is available for spreadsheet inspection.
5. Check the common recording interval and HEMOX time offset before running the analysis. Repeat with S02, NS01, and NS02 to inspect all four sessions.

The `FUSED_TRENDS` files include native `fused_hr_bpm` and `fused_spo2_percent` fields, a `quality_valid` flag, and separately labelled calibrated fields. Use native estimates when reproducing post-fusion calibration so that the correction is applied once.

## 📥 Supported explorer inputs

| Source | Primary input | Optional supporting evidence |
|---|---|---|
| SP-20 | CSV, TXT, or TSV recording per run | — |
| HEMOX | `FUSED_TRENDS` or compatible validated-trends file | `PARAMETER_TRENDS`, `RAW`, `WINDOWS`, `MOTION_THRESHOLDS`, `SESSION_SUMMARY`, `METADATA`, and `REPORT` files |

The explorer automatically inspects compatible column names and uses synchronized HEMOX outputs for physiological comparison.

## 🧩 Analysis modules

| Module | Purpose |
|---|---|
| Data & synchronization | Verify uploads, run roles, common intervals, and temporal matching |
| Performance | Summarize per-run, pooled calibration, pooled validation, and overall metrics |
| Calibration | Estimate candidate models from calibration runs and freeze coefficients before validation |
| Filtered & fused | Inspect CH1, CH2, fused trends, and inter-channel disagreement |
| IMU motion | Relate acceleration/gyroscope evidence to signal reliability |
| Agreement | Generate regression and Bland–Altman analyses |
| SQI robustness | Examine error as a function of signal quality |
| Technical report | Produce structured findings, equations, limitations, and interpretation |
| Export | Generate reproducible tables, figures, synchronized pairs, and calibration settings |

## 🐍 Python analysis scripts

The [`scripts/`](scripts/) directory contains **13 Python analysis files** for offline HEMOX analysis and figure generation. They cover PPG preprocessing, adaptive motion thresholds, signal quality, fusion visualization, synchronization checks, calibration, and agreement with the SP-20 comparator.

### Script guide

| Script | Purpose |
|---|---|
| [`PPG_Preprocessing_Plots.py`](scripts/PPG_Preprocessing_Plots.py) | Visualize raw, artifact-suppressed, and filtered red/infrared PPG, including stationary and controlled-motion recordings. |
| [`AdaptiveMotionThresholds_ALL.py`](scripts/AdaptiveMotionThresholds_ALL.py) | Plot CH1/CH2 adaptive motion-threshold histories from exported thresholds or reconstruct them from raw IMU data. |
| [`FigS_5s_vs_Final_Thresholds.py`](scripts/FigS_5s_vs_Final_Thresholds.py) | Compare initial calibration motion limits reconstructed from raw IMU data with final thresholds recorded in session reports. |
| [`Reliability_Aware_Fusion.py`](scripts/Reliability_Aware_Fusion.py) | Visualize CH1, CH2, and fused HR/SpO₂ alongside motion and quality evidence, using Explorer-style window reconstruction or exported parameter trends. |
| [`CH1_CH2_CompleteHEMOX_Pairwise_Common.py`](scripts/CH1_CH2_CompleteHEMOX_Pairwise_Common.py) | Compare individual channels, simpler fusion methods, fixed/adaptive thresholds, and complete HEMOX using common paired observations and availability summaries. |
| [`HEMOX_R2_Synchronization_Validation.py`](scripts/HEMOX_R2_Synchronization_Validation.py) | Contains audit functions for bilateral timing, counter integrity, optical-data completeness, and comparator matching; the current main workflow is incomplete (see below). |
| [`bland_Altman_perf_befor_calibration.py`](scripts/bland_Altman_perf_befor_calibration.py) | Generate pooled Bland–Altman plots and performance metrics for native, uncalibrated SpO₂ and HR. |
| [`bland_Altman_perf_after_calibration.py`](scripts/bland_Altman_perf_after_calibration.py) | Fit the SpO₂ M1 offset on P01–P02 and report separate development/held-out agreement, retaining native HR M0. |
| [`Perf_Multi_calibration.py`](scripts/Perf_Multi_calibration.py) | Compare M0–M3 calibration models and calibration-cohort sizes using participant-disjoint combinations, with metric tables and figures. |
| [`RunLevel_Analysis_Figures.py`](scripts/RunLevel_Analysis_Figures.py) | Summarize held-out run-level SpO₂/HR errors with heatmaps and stationary-versus-motion comparisons. |
| [`HEMOX_vs_SP20_P03_P10_M0_vs_SelectedConfig.py`](scripts/HEMOX_vs_SP20_P03_P10_M0_vs_SelectedConfig.py) | Plot HEMOX and SP-20 time-domain traces under native or selected SpO₂ settings, retaining native HR. |
| [`HEMOX_Correlation_AgreementBand.py`](scripts/HEMOX_Correlation_AgreementBand.py) | Plot SpO₂/HR comparisons with configurable agreement bands and export associated agreement statistics. |
| [`Study_SQI_FullRange_Analysis.py`](scripts/Study_SQI_FullRange_Analysis.py) | Examine held-out error, agreement, and observation distribution across fused-SQI bins, including participant-bootstrap confidence intervals. |

**Synchronization-file status:** The current upload of `HEMOX_R2_Synchronization_Validation.py` ends inside its run-processing loop and has no script entry point. The main workflow needs completion before direct execution can produce the intended audit results.

### Requirements

Use **Python 3.10 or newer** with NumPy, pandas, and Matplotlib:

```bash
python -m pip install numpy pandas matplotlib
```

Run the examples below from the repository root.

### Input data and study configuration

The scripts read HEMOX session exports and, for comparator analyses, corresponding SP-20 recordings. Required files vary by script: `*_RAW.csv`, `*_WINDOWS.csv`, `*_PARAMETER_TRENDS.csv`, `*_FUSED_TRENDS.csv`, `*_MOTION_THRESHOLDS.csv`, and `*_REPORT.txt`.

Cohort analyses use participant identifiers **P01–P10** and run identifiers **S01/S02** (stationary) and **NS01/NS02** (controlled motion). Use extracted participant folders for the Bland–Altman and multi-model calibration scripts; several other scripts also support participant ZIP archives.

The selected-configuration analyses use **M1 offset-only calibration for SpO₂** and **M0 native output for HR**, with P01–P02 as the development cohort and P03–P10 held out. The calibration-size study evaluates additional participant-disjoint subsets. Review each script's cohort, synchronization, and calibration settings when preparing an analysis.

Full-cohort analyses require the corresponding participant/run recordings. The P01 sample folder alone cannot reproduce study-wide results.

### Example commands

Run the P01 examples from the repository root after installing the dependencies above. For other datasets, replace the input path and participant/run selection.

```bash
# PPG preprocessing figures using the P01/NS01 sample
python scripts/PPG_Preprocessing_Plots.py --input "samples" --participant P01 --run NS01 --out-dir "outputs/P01_preprocessing"

# CH1/CH2 and fused-trend visualization for the same sample run
python scripts/Reliability_Aware_Fusion.py --data-dir "samples" --participant P01 --run NS01 --out-dir "outputs/P01_fusion"

# Full-cohort agreement example: requires extracted P01–P10 study folders
python scripts/bland_Altman_perf_after_calibration.py --data-dir "path/to/full-study-data" --out-dir "outputs/agreement"

# View available options for a command-line script
python scripts/Study_SQI_FullRange_Analysis.py --help
```

`AdaptiveMotionThresholds_ALL.py` and `FigS_5s_vs_Final_Thresholds.py` use the **USER SETTINGS** block near the top of each file instead of command-line arguments. Set `DATA_ROOT` there to your extracted data directory.

Outputs include script-specific image files, CSV tables, and analysis summaries. Use `--out-dir` where available, or the output-directory setting in the script, to choose where results are saved.

## 📁 Repository contents

```text
HEMOX/
├── CHANGELOG.md              # Explorer version history
├── CITATION.cff              # Citation metadata
├── DEPLOYMENT.md             # GitHub Pages deployment notes
├── HEMOXV1.3.7.1.apk         # HEMOX Android application
├── README.md                 # Project documentation
├── index.html                # Self-contained HEMOX Validation Explorer
├── scripts/                 # Python analysis and figure-generation scripts
└── samples/
    └── P01/                 # Four paired sessions, 40 files
        ├── HEMOX/           # Four session folders, nine exports each
        └── SP20/            # Four matching comparator CSV files
```

## 🔐 Privacy and reproducibility

- Validation recordings loaded into the browser explorer remain in local browser memory.
- The static explorer does not upload research recordings to an external backend.
- Synchronized observations, performance metrics, and analysis settings can be exported for audit and reproducibility.
- The agreement convention used by the explorer is **HEMOX − SP-20**.
- Calibration-development and independent-validation data are handled separately to avoid fitting on held-out validation recordings.
- The `samples/P01/` directory provides four paired sessions with raw signals, fused and channel-level trends, quality/motion evidence, session reports, spreadsheet exports, and SP-20 comparator recordings.

## 📚 Citation

If you use the HEMOX software platform, Android application, validation explorer, analysis workflow, sample recordings, generated figures, or exported results in academic work, please cite the repository.

### IEEE-style reference

> A. R. Mohamed, A. H. Muqaibel, and M. A. Al-Absi, “HEMOX: Synchronous dual-site PPG–IMU monitoring and validation platform,” GitHub repository, 2026. [Online]. Available: https://github.com/Ahmed-Reda-SILO/HEMOX.

### BibTeX

```bibtex
@misc{Mohamed2026HEMOX,
  author       = {Ahmed Reda Mohamed and Ali Hussein Muqaibel and Muneer A. Al-Absi},
  title        = {{HEMOX}: Synchronous Dual-Site {PPG--IMU} Monitoring and Validation Platform},
  year         = {2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/Ahmed-Reda-SILO/HEMOX},
  note         = {Android application, validation explorer, sample recordings, and supporting software resources}
}
```

## 👨‍🔬 Authors

**Dr. Ahmed Reda Mohamed**  
Interdisciplinary Research Center for Communication Systems and Sensing (IRC-CSS)  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia

**Prof. Ali Hussein Muqaibel**  
Electrical Engineering Department  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia  
Interdisciplinary Research Center for Communication Systems and Sensing (IRC-CSS)  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia

**Prof. Muneer A. Al-Absi**  
Electrical Engineering Department  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia  
Interdisciplinary Research Center for Communication Systems and Sensing (IRC-CSS)  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia

### Contributions

- **Ahmed Reda Mohamed** — system design, wearable implementation, firmware, Android application, physiological-signal processing, validation framework, data analysis, and software development.
- **Ali Hussein Muqaibel** — research supervision, methodology, experimental validation, and technical review.
- **Muneer A. Al-Absi** — research supervision, system methodology, technical guidance, and manuscript/research review.

Project profile: [Ahmed-Reda-SILO](https://github.com/Ahmed-Reda-SILO)

---

<div align="center">

**HEMOX — synchronized multi-site physiological sensing for reproducible wearable-health research.**

</div>
