# Reproducing HEMOX analyses

## Requirements

Use Python 3.10 or later with NumPy, pandas, and Matplotlib (see [`requirements.txt`](requirements.txt)). Some plot scripts use `tkinter` for an optional interactive local-file chooser; headless runs can generally use CLI path flags. Arial/Arial Narrow may be used for manuscript-matched typography but are not bundled.

```bash
python -m pip install -r research/requirements.txt
```

## Distinguish retained outputs from source observations

This directory contains **14 source `.py` files and their 2026-generated figures/tables**. It does **not** contain the complete original participant-level `RAW`, `FUSED_TRENDS`, `PARAMETER_TRENDS`, `WINDOWS`, or SP-20 run exports for all P01–P10. The individual `.csv` summary and selected paired-observation files retained here do not reconstruct every underlying stream. Supply full original HEMOX/SP-20 exports to independently rerun analyses that read raw inputs.

Directory conventions expected by many scripts: `P01`, ..., `P10` folders or participant ZIP files containing labeled session exports for `S01`, `S02`, `NS01`, `NS02`, alongside their SP-20 counterparts. Many scripts accept `--data-dir`/`--input`/`--out-dir`; consult `--help` for each script. The following commands are *examples of invocation*, not claims that the source recordings are included.

## Representative invocations

From repository root, with `<DATA_ROOT>` pointing to the separately retained raw HEMOX and SP-20 exports:

```bash
# SpO2 M1 / HR M0 held-out Bland–Altman evaluation
python research/BA_after_calibration_outputs/bland_Altman_perf_after_calibration.py \
  --data-dir "<DATA_ROOT>" --out-dir "<OUTPUT_ROOT>/BA_after"

# Native, uncalibrated all-participant baseline
python research/BA_without_calibration_outputs/2-bland_Altman_perf_befor_calibration.py \
  --data-dir "<DATA_ROOT>" --out-dir "<OUTPUT_ROOT>/BA_before"

# Calibration size and model complexity analysis
python research/Study_perf_Multi_calibration_outputs/Study_perf_Multi_calibration_FIXED.py \
  --data-dir "<DATA_ROOT>" --out-dir "<OUTPUT_ROOT>/calibration_size"

# Held-out synchronized-record audit
python research/Study_HEMOX_Correlation_outputs/HEMOX_R2_Synchronization_Validation_FIXED.py \
  --data-root "<DATA_ROOT>" \
  --pairs-csv "research/Study_HEMOX_Correlation_outputs/native_and_selected_reference_pairs.csv" \
  --out-dir "<OUTPUT_ROOT>/synchronization"

# Run-level analysis, starting from exported aligned observations
python research/Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures.py \
  --pairs-csv "research/Study_HEMOX_Correlation_outputs/native_and_selected_reference_pairs.csv" \
  --out-dir "<OUTPUT_ROOT>/run_level"

# SQI acceptance / error analysis
python research/Study_Absolute_Error_vs_Fused_SQI_outputs/Study_SQI_FullRange_Zoom_Analysis.py \
  --data-dir "<DATA_ROOT>" --out-dir "<OUTPUT_ROOT>/sqi"
```

**Important:** Script defaults are not uniform. For complete reproduction, inspect each script's parameters and preserve any project-specific reference-offset CSV or selection manifest used at generation time. Two scripts (`HEMOX_AdaptiveMotionThresholds_ALL.py` and `HEMOX_FigS_5s_vs_Final_Thresholds.py`) use paths relative to their own file directory rather than a standard command-line `--data-dir`; they may need their input-path variables configured before rerunning. This package retains original script logic, without changing numerical methodology.

## Measurement model and safeguards

* Cohort: 10 participants, 40 sessions; 20 stationary and 20 controlled-motion recordings, approximately five minutes each.
* Calibration: frozen SpO₂ M1 offset estimated from P01–P02, using 589 valid pairs. Held-out P03–P10 results use this fixed coefficient and native HR M0.
* Comparator pairing: a ±1 s unique nearest-neighbor tolerance applied to eligible HEMOX/SP-20 windows. The absence of verified per-run clock offsets makes any default zero-offset assumption **provisional**.
* Different endpoint validity counts are normal: SpO₂ **2,202** and HR **2,255** held-out pairs. Do not reuse one denominator across both.
* Paired-point LoA use `bias ± 1.96 × sample SD` and are descriptive, not participant-independent 95% confidence intervals.
* CH1–CH2 sync metrics reflect exported `sync_delta_us` residuals; record-counter continuity is not a direct claim about individual BLE packet delivery. Optical-row completeness is distinct from counter continuity.
* The calibration-size sweep operates across many participant subsets. Its mean held-out performance by training-subset size **must not be presented as the exact P01–P02/P03–P10 result**.
* The browser explorer can display visually held values for continuity, but performance metrics must use original valid fused outputs, not display-only reconstructions.
* The SP-20 is not an invasive ground-truth clinical reference; no device certification is implied.

## Checksums and provenance

`FILE_MANIFEST.csv` lists SHA-256 hashes of **343 original extracted artifacts**; use it to verify the uploaded copies. Auxiliary `README.md`, `REPRODUCIBILITY.md`, and `requirements.txt` were added for repository usability. Separately uploaded ZIP archives have been verified identical to subfolders of the single consolidated archive and do not require duplicate storage.

## Privacy and repository visibility

Run-level CSVs and filenames refer to pseudonymous participants. A private GitHub repository is still an access-controlled sharing environment, so confirm institutional data-sharing approval before granting collaborators access or changing repository visibility to public. Avoid committing directly identifying participant records, credentials, and unrelated device dumps.
