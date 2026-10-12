# HEMOX · Analysis Scripts and Results

> **Repository research artifact index (2026).** The files originate from the user-supplied consolidated *Scripts and Results* archive. Eleven separately provided archive bundles were verified byte-for-byte identical to their counterparts in the consolidated archive and were not duplicated.

**Inventory:** 12 analysis collections; 14 Python scripts; 117 CSV tables/exports; 103 JPEG figures; 57 PNG figures; 43 SVG figures; five Markdown notes; three text notes; one TeX table. **343 original files** in total, plus documentation/manifest.

## Analysis modules

| Collection | Scope | Scripts and results |
|---|---|---|
| [`BA_without_calibration_outputs/`](BA_without_calibration_outputs/) | Native (M0) baseline agreement; all P01–P10; overall metrics, per-run tables, Bland–Altman | 11 files |
| [`BA_after_calibration_outputs/`](BA_after_calibration_outputs/) | Calibration-development P01–P02 / held-out P03–P10; selected SpO₂ M1 and HR M0, LoA | 11 files |
| [`Study_perf_Multi_calibration_outputs/`](Study_perf_Multi_calibration_outputs/) | Participant-subset calibration-size sweeps, M1–M3 model comparison | 25 files |
| [`Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output/`](Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output/) | Synchronized bilateral channel pairwise comparisons | 3 files |
| [`HEMOX_AdaptiveMotionThresholds_ALL/`](HEMOX_AdaptiveMotionThresholds_ALL/) | All-run and per-recording adaptive motion-threshold exports | 125 files |
| [`HEMOX_5s_vs_Final_Thresholds/`](HEMOX_5s_vs_Final_Thresholds/) | 5-second vs final motion-threshold illustration | 4 files |
| [`HEMOX_PPG_R2_Styled_outputs/`](HEMOX_PPG_R2_Styled_outputs/) | High-resolution PPG preprocessing plots, optical AC/DC traces | 61 files |
| [`HEMOX_vs_SP20_AllParticipants_outputs/`](HEMOX_vs_SP20_AllParticipants_outputs/) | P03–P10 trends for M0 versus selected model configuration | 2 files |
| [`Study_Reliability_Aware_Fusion_outputs/`](Study_Reliability_Aware_Fusion_outputs/) | P09/NS02 filtered channel plots, IMU evidence, fusion and display-hold safeguards | 41 files |
| [`Study_HEMOX_Correlation_outputs/`](Study_HEMOX_Correlation_outputs/) | Held-out correlation, agreement bands, strict synchronization QA | 15 files |
| [`Run_Level_R2_Standalone_Output/`](Run_Level_R2_Standalone_Output/) | Held-out run/participant metrics, heatmap, stationary versus controlled motion | 12 files |
| [`Study_Absolute_Error_vs_Fused_SQI_outputs/`](Study_Absolute_Error_vs_Fused_SQI_outputs/) | Fused SQI, native validity/exclusion accounting, agreement by SQI | 33 files |

## Python source files

Each script is preserved with its original filename and is listed here for direct review:

- [`BA_after_calibration_outputs/bland_Altman_perf_after_calibration.py`](BA_after_calibration_outputs/bland_Altman_perf_after_calibration.py)
- [`BA_without_calibration_outputs/2-bland_Altman_perf_befor_calibration.py`](BA_without_calibration_outputs/2-bland_Altman_perf_befor_calibration.py)
- [`Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output/HEMOX_R2_CH1_CH2_CompleteHEMOX_Pairwise_Common.py`](Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output/HEMOX_R2_CH1_CH2_CompleteHEMOX_Pairwise_Common.py)
- [`HEMOX_5s_vs_Final_Thresholds/HEMOX_FigS_5s_vs_Final_Thresholds.py`](HEMOX_5s_vs_Final_Thresholds/HEMOX_FigS_5s_vs_Final_Thresholds.py)
- [`HEMOX_AdaptiveMotionThresholds_ALL/HEMOX_AdaptiveMotionThresholds_ALL.py`](HEMOX_AdaptiveMotionThresholds_ALL/HEMOX_AdaptiveMotionThresholds_ALL.py)
- [`HEMOX_PPG_R2_Styled_outputs/HEMOX_PPG_Preprocessing_Plots_R2_Styled_2PanelMotion_Stationary.py`](HEMOX_PPG_R2_Styled_outputs/HEMOX_PPG_Preprocessing_Plots_R2_Styled_2PanelMotion_Stationary.py)
- [`HEMOX_vs_SP20_AllParticipants_outputs/HEMOX_vs_SP20_P03_P10_M0_vs_SelectedConfig.py`](HEMOX_vs_SP20_AllParticipants_outputs/HEMOX_vs_SP20_P03_P10_M0_vs_SelectedConfig.py)
- [`Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures.py`](Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures.py)
- [`Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures_Standalone.py`](Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures_Standalone.py)
- [`Study_Absolute_Error_vs_Fused_SQI_outputs/Study_SQI_FullRange_Zoom_Analysis.py`](Study_Absolute_Error_vs_Fused_SQI_outputs/Study_SQI_FullRange_Zoom_Analysis.py)
- [`Study_HEMOX_Correlation_outputs/HEMOX_R2_Synchronization_Validation_FIXED.py`](Study_HEMOX_Correlation_outputs/HEMOX_R2_Synchronization_Validation_FIXED.py)
- [`Study_HEMOX_Correlation_outputs/Study_HEMOX_Correlation_Heldout_AgreementBand_2Panel.py`](Study_HEMOX_Correlation_outputs/Study_HEMOX_Correlation_Heldout_AgreementBand_2Panel.py)
- [`Study_Reliability_Aware_Fusion_outputs/Study_Reliability_Aware_Fusion_R3_CleanLegend_MotionPeakReport.py`](Study_Reliability_Aware_Fusion_outputs/Study_Reliability_Aware_Fusion_R3_CleanLegend_MotionPeakReport.py)
- [`Study_perf_Multi_calibration_outputs/Study_perf_Multi_calibration_FIXED.py`](Study_perf_Multi_calibration_outputs/Study_perf_Multi_calibration_FIXED.py)

## Evidence and interpretation notes

* Held-out model-fit exclusion: P01–P02 develop the frozen SpO₂ offset; P03–P10 are evaluated after calibration, not used for fitting the M1 coefficient.
* Frozen SpO₂ M1 offset: **0.9067538200339558 pp**; selected HR model is uncalibrated M0.
* Key held-out evidence: [`BA_after_calibration_outputs/Tables/heldout_validation_before_after.csv`](BA_after_calibration_outputs/Tables/heldout_validation_before_after.csv), [`Study_HEMOX_Correlation_outputs/heldout_agreement_band_metrics.csv`](Study_HEMOX_Correlation_outputs/heldout_agreement_band_metrics.csv), and [`Run_Level_R2_Standalone_Output/pooled_heldout_descriptive_metrics.csv`](Run_Level_R2_Standalone_Output/pooled_heldout_descriptive_metrics.csv).
* Bilateral sync: [`Study_HEMOX_Correlation_outputs/Synchronization_Validation_v2_Output/Synchronization_Validation_Summary.csv`](Study_HEMOX_Correlation_outputs/Synchronization_Validation_v2_Output/Synchronization_Validation_Summary.csv).
* Pre/post calibration tables and figures differ in **cohort and sample selection**. Compare M0 versus M1 using **the same held-out pairs**, not by subtracting metrics from the all-participant native baseline.
* For reliability-aware fusion figures, `*_display_held_only` and any display-only interpolation are **not valid observations** for performance scoring.
* CSVs exported from analysis scripts may contain pseudo-participant/session-level records; review consent and privacy before making the repository public.
* Original participant source archives are **not present in this kit**. Some scripts require raw session files and SP-20 exports to regenerate tables. Reproducing the existing generated results requires those inputs and verified clock-offset settings.

## Reproduction documentation

Read [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for setup, script entry points, and methodological guardrails. SHA-256 file hashes for each supplied artifact are recorded in [`FILE_MANIFEST.csv`](FILE_MANIFEST.csv).
