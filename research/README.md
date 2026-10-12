# HEMOX Research Scripts and Results

The complete downloadable research upload kit contains **12 analysis collections, 14 Python scripts, 117 CSV tables, 103 JPG figures, 57 PNG figures, and 43 SVG figures**, representing 343 original files.

The repository documentation is staged in this review branch. **The full collection of large image/CSV/script files remains in the downloadable kit and has not yet been committed.** The files below become navigable when the research folder from the kit is uploaded.

| Collection | Analysis focus |
|---|---|
| BA_without_calibration_outputs | All-participant native Bland–Altman baseline |
| BA_after_calibration_outputs | P01–P02 calibration and held-out P03–P10 agreement |
| Study_perf_Multi_calibration_outputs | M1/M2/M3 participant subset and calibration-size comparisons |
| Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output | Bilateral common-time paired channels |
| HEMOX_AdaptiveMotionThresholds_ALL | Per-run and aggregate adaptive thresholds |
| HEMOX_5s_vs_Final_Thresholds | Five-second versus final thresholds |
| HEMOX_PPG_R2_Styled_outputs | Stationary and motion PPG preprocessing |
| HEMOX_vs_SP20_AllParticipants_outputs | Held-out model-specific trend overlays |
| Study_Reliability_Aware_Fusion_outputs | Reliability-aware bilateral fusion and IMU evidence |
| Study_HEMOX_Correlation_outputs | Correlation, agreement, synchronization |
| Run_Level_R2_Standalone_Output | Run-level validation and stationary/motion contrasts |
| Study_Absolute_Error_vs_Fused_SQI_outputs | SQI selection, error, and missingness |

[Reproduction instructions](REPRODUCIBILITY.md) · [Python dependencies](requirements.txt).

**Key original results:** matched held-out SpO2 M0 RMSE 1.132 pp versus M1 0.968 pp; HR M0 RMSE 2.495 bpm. The SpO2 coefficient fitted from 589 P01–P02 pairs was +0.906754 pp. Held-out P03–P10 observations are participant-disjoint from fitting, but are not an untouched external set for all analysis decisions.

Full-data reproduction requires original HEMOX participant recordings, the comparator SP-20 export, and verified time-offset settings where available.
