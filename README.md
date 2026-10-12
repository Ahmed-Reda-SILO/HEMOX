<div align="center">

# 🫀 HEMOX

### Synchronous Dual-Site PPG–IMU Monitoring, Calibration, and Validation Platform

[![Android App](https://img.shields.io/badge/Android-v1.3.7.1-3DDC84?style=for-the-badge&logo=android&logoColor=white)](HEMOXV1.3.7.1.apk)
[![Explorer](https://img.shields.io/badge/Validation%20Explorer-v1.7.12-1261A0?style=for-the-badge)](index.html)
[![Research Analyses](https://img.shields.io/badge/Python-Research%20Results-3776AB?style=for-the-badge&logo=python&logoColor=white)](research/README.md)
[![Sample Data](https://img.shields.io/badge/Sample-Dataset-8A2BE2?style=for-the-badge)](samples/HEMOX_DUAL_20121986L00_20260818_213825.zip)
[![GitHub Pages](https://img.shields.io/badge/Open-HEMOX%20Explorer-222?style=for-the-badge&logo=github)](https://ahmed-reda-silo.github.io/HEMOX/)

**[📱 Download HEMOX Android App](HEMOXV1.3.7.1.apk)**  
**[📦 Download Sample Dataset](samples/HEMOX_DUAL_20121986L00_20260818_213825.zip)**  
**[🚀 Open the HEMOX Validation Explorer](https://ahmed-reda-silo.github.io/HEMOX/)**

</div>

---

## 🔬 Overview

**HEMOX** is a wearable, wireless, dual-site photoplethysmography (PPG) and inertial measurement unit (IMU) platform developed for synchronized bilateral physiological monitoring. The system is designed to acquire red/infrared PPG and local motion information from two independent sensing sites, preserve a common temporal reference, and support motion- and signal-quality-aware physiological analysis and multi-site fusion.

The software workflow combines an Android acquisition/monitoring application with a browser-based validation explorer for synchronized comparison against a commercial pulse-oximetry comparator, calibration development, independent validation, agreement analysis, and publication-oriented reporting.

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

## 📦 Sample dataset

A complete example HEMOX dual-site recording session is provided to demonstrate the exported data structure and allow users to test the validation workflow without first acquiring their own recording.

**[Download `HEMOX_DUAL_20121986L00_20260818_213825.zip`](samples/HEMOX_DUAL_20121986L00_20260818_213825.zip)**

The archive contains one HEMOX session together with the corresponding SP-20 comparator recording:

| File | Description |
|---|---|
| `*_RAW.csv` | synchronized raw bilateral PPG and IMU recording |
| `*_PARAMETER_TRENDS.csv` | time-varying physiological and quality parameters |
| `*_FUSED_TRENDS.csv` | fused physiological output used for comparator analysis |
| `*_WINDOWS.csv` | window-level processing and quality information |
| `*_MOTION_THRESHOLDS.csv` | motion-related threshold information |
| `*_SESSION_SUMMARY.csv` | compact session-level summary |
| `*_METADATA.txt` | recording and configuration metadata |
| `*_REPORT.txt` | generated session report |
| `*_RECORDING.xlsx` | consolidated spreadsheet export |
| `sp20-L00.csv` | corresponding commercial pulse-oximetry comparator recording |

> [!NOTE]
> This archive is provided as an **example dataset** illustrating the HEMOX export format and analysis workflow. It should not be interpreted as the complete experimental dataset used for all calibration or validation results reported in associated research work.

## 🧭 Typical validation workflow

1. Add the required HEMOX recording runs.
2. Assign each run to **Calibration**, **Independent validation**, or **Excluded**.
3. Upload one SP-20 CSV/TXT file for each comparator recording.
4. Upload the corresponding HEMOX validated/fused trend files and optional supporting data.
5. Confirm synchronization settings and run the pooled analysis.
6. Review physiological agreement, calibration transfer, motion, signal quality, and bilateral-channel behavior.
7. Export synchronized data, performance metrics, figures, reports, and fixed calibration parameters.

For a quick test, extract the provided sample archive and load `*_FUSED_TRENDS.csv` together with `sp20-L00.csv` into the explorer. Supporting files from the same session can also be added for motion, SQI, raw-signal, and reporting analyses.

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

## 🧪 Research scripts, figures, and validation results

The 2026 research archive includes **12 analysis collections, 14 Python scripts, 117 CSV datasets/tables, and 203 figure assets** (103 JPEG, 57 PNG, and 43 SVG), with **343 original files** in total. Source scripts are retained with their corresponding result folders rather than separated from the evidence they generated.

**[Browse the research index](research/README.md)** · **[Reproduction guide](research/REPRODUCIBILITY.md)** · **[File provenance manifest](research/FILE_MANIFEST.csv)**

### Experimental protocol and selected models

| Property | Description |
|---|---|
| Participants | **N=10** (P01–P10) |
| Recordings | **40 × 5 min** (20 stationary; 20 controlled motion; **200 min**) |
| Development cohort | **P01–P02**, 8 runs and **589 valid calibration pairs** |
| Evaluation cohort | **P03–P10**, 32 runs, participant-disjoint from fitting |
| SpO₂ selection | **M1**, offset-only: `clip(native fused SpO2 − 0.9067538200 pp, 0, 100)` |
| HR selection | **M0**, native fused HR without additional calibration |
| Comparator | Creative SP-20; eligible one-to-one nearest observations within **±1 s** |
| Signed error | **HEMOX − SP-20**; SpO₂ differences expressed in **percentage points (pp)** |

> [!NOTE]
> P03–P10 are held out from the **M1 coefficient estimation**, but the evaluation observations also informed comparative characterization. They should not be represented as an untouched external test set for every analysis decision.

### Quantitative results from exported audit tables

| Analysis | Result | Supporting CSV |
|---|---|---|
| Held-out SpO₂, selected M1 | **2,202 pairs**; bias **−0.263 pp**; MAE **0.783 pp**; RMSE **0.968 pp** | [Matched calibration comparison](research/BA_after_calibration_outputs/Tables/heldout_validation_before_after.csv) |
| Same-pair SpO₂, native M0 | RMSE **1.132 pp**; applying M1 reduces matched RMSE by **14.5%** | [Matched calibration comparison](research/BA_after_calibration_outputs/Tables/heldout_validation_before_after.csv) |
| Held-out HR, native M0 | **2,255 pairs**; bias **+0.111 bpm**; MAE **1.869 bpm**; RMSE **2.495 bpm** | [Run-level pooled results](research/Run_Level_R2_Standalone_Output/pooled_heldout_descriptive_metrics.csv) |
| SpO₂ error within ±3 pp | **99.95%** of eligible held-out observations | [Agreement bands](research/Study_HEMOX_Correlation_outputs/heldout_agreement_band_metrics.csv) |
| HR error within ±5 bpm | **94.28%** of eligible held-out observations | [Agreement bands](research/Study_HEMOX_Correlation_outputs/heldout_agreement_band_metrics.csv) |
| Bilateral sync residual | **0.669 ms median**, **4.324 ms 95th percentile**, **5.189 ms maximum** over **480,000** records | [Synchronization summary](research/Study_HEMOX_Correlation_outputs/Synchronization_Validation_v2_Output/Synchronization_Validation_Summary.csv) |
| Bilateral optical completeness | **480,000 / 480,000** exported synchronized rows have CH1/CH2 red and IR readings; **4 missing counter increments** separately identified | [Synchronization summary](research/Study_HEMOX_Correlation_outputs/Synchronization_Validation_v2_Output/Synchronization_Validation_Summary.csv) |

All pooled paired-observation agreement statistics are **descriptive**, because repeated windows within participants are correlated. The default **zero-offset HEMOX–SP-20 timing assumption is provisional** where independent per-run synchronization verification is absent. These SP-20 comparator results are **not** arterial-reference accuracy certification.

### Featured figures

| Selected-calibration Bland–Altman | Bilateral synchronization QA |
|---|---|
| ![Bland–Altman SpO2 after selected calibration](research/BA_after_calibration_outputs/Figures/BA_after_SpO2.jpg) | ![HEMOX synchronization three-panel audit](research/Study_HEMOX_Correlation_outputs/Synchronization_Validation_v2_Output/Synchronization_Validation_3Panel.png) |
| [Underlying results](research/BA_after_calibration_outputs/) | [Underlying synchronization results](research/Study_HEMOX_Correlation_outputs/) |

| Run-level validation | Adaptive motion thresholds |
|---|---|
| ![HEMOX run-level RMSE heatmap](research/Run_Level_R2_Standalone_Output/Heldout_RunLevel_RMSE_Heatmap.jpg) | ![All-run adaptive motion thresholds](research/HEMOX_AdaptiveMotionThresholds_ALL/Adaptive_Motion_Thresholds_ALL_RUNS_COMBINED.png) |
| [Per-run tables](research/Run_Level_R2_Standalone_Output/) | [Per-run threshold data](research/HEMOX_AdaptiveMotionThresholds_ALL/) |

### Full analysis inventory

| Analysis collection | Main focus |
|---|---|
| [BA_without_calibration_outputs](research/BA_without_calibration_outputs/) | N=10 native M0 bias, RMSE, and Bland–Altman |
| [BA_after_calibration_outputs](research/BA_after_calibration_outputs/) | P01–P02 M1 development, P03–P10 selected-model agreement |
| [Study_perf_Multi_calibration_outputs](research/Study_perf_Multi_calibration_outputs/) | M1–M3 calibration model/size comparisons |
| [Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output](research/Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output/) | Bilateral channel comparative plots |
| [HEMOX_AdaptiveMotionThresholds_ALL](research/HEMOX_AdaptiveMotionThresholds_ALL/) | Motion-threshold runs and aggregate evidence |
| [HEMOX_5s_vs_Final_Thresholds](research/HEMOX_5s_vs_Final_Thresholds/) | Early-window versus final threshold study |
| [HEMOX_PPG_R2_Styled_outputs](research/HEMOX_PPG_R2_Styled_outputs/) | Stationary/motion optical preprocessing and AC/DC figures |
| [HEMOX_vs_SP20_AllParticipants_outputs](research/HEMOX_vs_SP20_AllParticipants_outputs/) | Selected-model SpO₂/HR time courses |
| [Study_Reliability_Aware_Fusion_outputs](research/Study_Reliability_Aware_Fusion_outputs/) | Dual-site reliability-aware fusion, IMU/PPG evidence |
| [Study_HEMOX_Correlation_outputs](research/Study_HEMOX_Correlation_outputs/) | Comparator agreement, timestamp and synchronization QA |
| [Run_Level_R2_Standalone_Output](research/Run_Level_R2_Standalone_Output/) | Held-out run-level RMSE, stationary vs motion |
| [Study_Absolute_Error_vs_Fused_SQI_outputs](research/Study_Absolute_Error_vs_Fused_SQI_outputs/) | Fused SQI, missingness, inclusion, accuracy by quality |

### Reproducing the analysis

```bash
python -m pip install -r research/requirements.txt
python research/BA_after_calibration_outputs/bland_Altman_perf_after_calibration.py \
  --data-dir "/path/to/HEMOX_and_SP20_sessions" \
  --out-dir "./reproduced/BA_after_calibration"
```

Original participant-level HEMOX and SP-20 recordings are **not included** in these derived research-output folders; complete regeneration requires access to those input exports, relevant manifests, and any verified clock-offset settings. See [reproduction details](research/REPRODUCIBILITY.md). No held/interpolated **display-only** trend values should be treated as valid paired physiological observations.

## 📁 Repository contents

```text
HEMOX/
├── CHANGELOG.md              # Explorer version history
├── CITATION.cff              # Citation metadata
├── DEPLOYMENT.md             # GitHub Pages deployment notes
├── HEMOXV1.3.7.1.apk         # HEMOX Android application
├── README.md                 # Project documentation + scientific results
├── index.html                # Self-contained HEMOX Validation Explorer
├── samples/
│   └── HEMOX_DUAL_20121986L00_20260818_213825.zip
│                             # Example HEMOX + SP-20 session
└── research/
    ├── README.md             # Index of all 12 analysis collections
    ├── REPRODUCIBILITY.md    # Data inputs and CLI reproductions
    ├── requirements.txt      # Python analysis requirements
    ├── FILE_MANIFEST.csv     # SHA-256 source artifact checksums
    └── [12 analysis directories with scripts, CSVs, figures]
```

## 🔐 Privacy and reproducibility

- Validation recordings loaded into the browser explorer remain in local browser memory.
- The static explorer does not upload research recordings to an external backend.
- Synchronized observations, performance metrics, and analysis settings can be exported for audit and reproducibility.
- The agreement convention used by the explorer is **HEMOX − SP-20**.
- Development and held-out-from-fitting recordings are handled separately. See the research notes for the scope of calibration-model selection and evaluation.
- The sample archive provides a reproducible example of the expected HEMOX session-export structure.

## 📚 Citation

If you use the HEMOX software platform, Android application, validation explorer, analysis workflow, sample data structure, generated figures, or exported results in academic work, please cite the repository.

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
  note         = {Android application, validation explorer, sample dataset, and supporting software resources}
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
