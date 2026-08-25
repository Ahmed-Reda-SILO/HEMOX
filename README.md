<div align="center">

# 🫀 HEMOX

### Dual-site PPG/IMU validation, calibration, and publication explorer

[![Version](https://img.shields.io/badge/version-1.7.12-1261A0?style=for-the-badge)](CHANGELOG.md)
[![Browser App](https://img.shields.io/badge/browser-app-18794E?style=for-the-badge&logo=googlechrome&logoColor=white)](index.html)
[![No Installation](https://img.shields.io/badge/installation-none-B64B2A?style=for-the-badge)](#-quick-start)
[![GitHub Pages](https://img.shields.io/badge/deployment-GitHub%20Pages-222?style=for-the-badge&logo=github)](https://ahmed-reda-silo.github.io/HEMOX/)

**[🚀 Open the HEMOX Explorer](https://ahmed-reda-silo.github.io/HEMOX/)**

</div>

---

## 🔬 Overview

HEMOX is a self-contained research explorer for comparing the dual-site HEMOX PPG/IMU platform with the commercial SP-20 pulse oximeter. It supports synchronized multi-run analysis, calibration-development and independent-validation separation, agreement assessment, motion and signal-quality evidence, and publication-ready exports.

All uploaded recordings are processed locally in the browser. The explorer does not require a backend, database, package installation, or build step.

> [!IMPORTANT]
> HEMOX–SP-20 results quantify agreement with a commercial comparator. They do not constitute formal arterial-reference pulse-oximetry validation.

## ✨ Main capabilities

- 📂 Multi-run SP-20 and HEMOX data management
- ⏱️ Strict one-to-one timestamp synchronization
- 🧪 Separate calibration and independent-validation roles
- 📊 Bias, MAE, RMSE, MAPE, ACC, correlation, and limits of agreement
- 📈 Comparator trends, per-run agreement, and Bland–Altman figures
- 🫁 SpO₂ and heart-rate analysis with independent availability metrics
- 🏃 Dual-IMU motion and signal-quality robustness evidence
- 🧮 M0–M3 calibration-model comparison and frozen coefficients
- 📱 Android-ready calibration parameter export
- 📝 Detailed technical reports and publication-ready SVG/PNG figures
- 🔒 Local browser processing for research-data privacy

## 🚀 Quick start

### Open the hosted explorer

Use [**Open the HEMOX Explorer**](https://ahmed-reda-silo.github.io/HEMOX/).

## 🧭 Typical workflow

1. Add the required recording runs.
2. Assign each run to **Calibration**, **Independent validation**, or **Excluded**.
3. Upload one SP-20 CSV/TXT file per run.
4. Upload the HEMOX `FUSED_TRENDS` file and any supporting export files.
5. Confirm synchronization settings and select **Analyze all runs**.
6. Review performance, calibration transfer, agreement, motion, SQI, and decision evidence.
7. Export synchronized pairs, metrics, figures, reports, and frozen calibration settings.

## 📥 Supported inputs

| Source | Required input | Optional supporting evidence |
|---|---|---|
| SP-20 | One CSV, TXT, or TSV recording per run | — |
| HEMOX | `FUSED_TRENDS` or a compatible validated-trends file | `PARAMETER_TRENDS`, `RAW`, `WINDOWS`, `MOTION_THRESHOLDS`, `SESSION_SUMMARY`, `METADATA`, and `REPORT` files |

The explorer automatically inspects column names and uses the validated HEMOX output for synchronized physiological comparison.

## 🧩 Analysis modules

| Module | Purpose |
|---|---|
| Data & synchronization | Verify uploads, roles, common intervals, and temporal matching |
| Performance | Summarize per-run, pooled calibration, pooled validation, and overall metrics |
| Calibration | Estimate models using calibration runs and freeze them before validation |
| Filtered & fused | Inspect CH1, CH2, fused trends, and inter-channel disagreement |
| IMU motion | Relate accelerometer/gyroscope evidence to recording quality |
| Agreement | Generate scatter, regression, and Bland–Altman results |
| SQI robustness | Examine error as a function of signal quality |
| Technical report | Produce a structured interpretation with equations, findings, and limitations |
| Export | Download reproducible data tables, figures, and calibration settings |

## 📁 Repository structure

```text
HEMOX/
├── .github/workflows/pages.yml   # Automatic GitHub Pages deployment
├── .nojekyll                     # Serve the static app without Jekyll processing
├── CHANGELOG.md                  # Version history
├── CITATION.cff                  # Citation metadata
├── DEPLOYMENT.md                 # Publishing and 404 troubleshooting
├── README.md                     # Project documentation
└── index.html                    # Complete HEMOX browser application
```

## 🔐 Privacy and reproducibility

- Recordings remain in browser memory during analysis.
- No research data are uploaded by this static application.
- The synchronized observations and metrics can be exported for audit and reproducibility.
- The difference convention is always **HEMOX − SP-20**.

## 📚 How to cite

If you use the HEMOX Explorer, its analysis workflow, generated figures, or exported results in an academic publication, please cite the software.

### IEEE reference format

> A. R. Mohamed, “HEMOX–SP-20 Validation Explorer,” version 1.7.12, Aug. 2026. [Online]. Available: https://github.com/Ahmed-Reda-SILO/HEMOX. Accessed: Aug. 25, 2026.

### IEEEtran LaTeX format

```latex
\bibitem{HEMOX}
A. R. Mohamed, ``HEMOX--SP-20 Validation Explorer,''
version 1.7.12, Aug. 2026. [Online]. Available:
\url{https://github.com/Ahmed-Reda-SILO/HEMOX}.
Accessed: Aug. 25, 2026.
```

### BibTeX format

```bibtex
@misc{Mohamed2026HEMOX,
  author       = {Ahmed Reda Mohamed},
  title        = {{HEMOX--SP-20 Validation Explorer}},
  year         = {2026},
  month        = aug,
  note         = {Version 1.7.12, accessed Aug. 25, 2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/Ahmed-Reda-SILO/HEMOX}
}
```

## 👨‍🔬 Author

**Dr. Ahmed Reda Mohamed, Member, IEEE**  
Postdoctoral Fellow, Interdisciplinary Research Center for Communication Systems and Sensing (IRC-CSS)  
King Fahd University of Petroleum and Minerals (KFUPM), Dhahran, Saudi Arabia

Project profile: [Ahmed-Reda-SILO](https://github.com/Ahmed-Reda-SILO)

---

<div align="center">

Developed for reproducible biomedical-signal validation and publication support. 🫀📊

</div>
