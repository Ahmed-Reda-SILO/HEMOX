<div align="center">

# 🫀 HEMOX

### Synchronous Dual-Site PPG–IMU Monitoring, Calibration, and Validation Platform

[![Android App](https://img.shields.io/badge/Android-v1.3.7.1-3DDC84?style=for-the-badge&logo=android&logoColor=white)](HEMOXV1.3.7.1.apk)
[![Explorer](https://img.shields.io/badge/Validation%20Explorer-v1.7.12-1261A0?style=for-the-badge)](index.html)
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

## 🎯 System concept

The HEMOX workflow follows an end-to-end synchronized multi-site architecture:

```text
CH1 PPG + IMU ─┐
               ├─> Synchronized wireless acquisition
CH2 PPG + IMU ─┘              │
                              ▼
                   Common bilateral timeline
                              │
              ┌───────────────┴───────────────┐
              ▼                               ▼
      Recording / export          Motion + signal quality
                                              │
                                              ▼
                                  Site-specific estimates
                                              │
                                              ▼
                                   Adaptive multi-site fusion
                                              │
                                              ▼
                                   Calibration correction
                                              │
                                              ▼
                                      HR / SpO₂ output
```

The architecture separates acquisition, storage, analysis, and calibration so that synchronized raw measurements are preserved independently of downstream processing.

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

## 📁 Repository contents

```text
HEMOX/
├── CHANGELOG.md              # Explorer version history
├── CITATION.cff              # Citation metadata
├── DEPLOYMENT.md             # GitHub Pages deployment notes
├── HEMOXV1.3.7.1.apk         # HEMOX Android application
├── README.md                 # Project documentation
├── index.html                # Self-contained HEMOX Validation Explorer
└── samples/
    └── HEMOX_DUAL_20121986L00_20260818_213825.zip
                              # Example HEMOX + SP-20 validation session
```

## 🔐 Privacy and reproducibility

- Validation recordings loaded into the browser explorer remain in local browser memory.
- The static explorer does not upload research recordings to an external backend.
- Synchronized observations, performance metrics, and analysis settings can be exported for audit and reproducibility.
- The agreement convention used by the explorer is **HEMOX − SP-20**.
- Calibration-development and independent-validation data are handled separately to avoid fitting on held-out validation recordings.
- The sample archive provides a reproducible example of the expected HEMOX session-export structure.

## 📚 Citation

If you use the HEMOX software platform, Android application, validation explorer, analysis workflow, sample data structure, generated figures, or exported results in academic work, please cite the repository.

### IEEE-style reference

> A. R. Mohamed, “HEMOX: Synchronous dual-site PPG–IMU monitoring and validation platform,” GitHub repository, 2026. [Online]. Available: https://github.com/Ahmed-Reda-SILO/HEMOX.

### BibTeX

```bibtex
@misc{Mohamed2026HEMOX,
  author       = {Ahmed Reda Mohamed},
  title        = {{HEMOX}: Synchronous Dual-Site {PPG--IMU} Monitoring and Validation Platform},
  year         = {2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/Ahmed-Reda-SILO/HEMOX},
  note         = {Android application, validation explorer, sample dataset, and supporting software resources}
}
```

## 👨‍🔬 Author

**Dr. Ahmed Reda Mohamed, Member, IEEE**  
Interdisciplinary Research Center for Communication Systems and Sensing (IRC-CSS)  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia

Project profile: [Ahmed-Reda-SILO](https://github.com/Ahmed-Reda-SILO)

---

<div align="center">

**HEMOX — synchronized multi-site physiological sensing for reproducible wearable-health research.**

</div>
