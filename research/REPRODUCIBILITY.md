# HEMOX analysis reproduction

Install dependencies (Python 3.10 or later):

```bash
python -m pip install -r research/requirements.txt
```

The supplied research-output folders preserve Python scripts, precomputed CSVs and figures. They **do not** include the complete set of participant `RAW`, `WINDOWS`, `FUSED_TRENDS`, `PARAMETER_TRENDS` or comparator SP-20 run-export source files. Original P01–P10 recordings and verified clock-offset settings are needed for full regeneration.

Example, from the repository root:

```bash
python research/BA_after_calibration_outputs/bland_Altman_perf_after_calibration.py --data-dir "/path/to/HEMOX_recordings" --out-dir "./reproduced/BA_after"
python research/BA_without_calibration_outputs/2-bland_Altman_perf_befor_calibration.py --data-dir "/path/to/HEMOX_recordings" --out-dir "./reproduced/BA_before"
python research/Study_perf_Multi_calibration_outputs/Study_perf_Multi_calibration_FIXED.py --data-dir "/path/to/HEMOX_recordings" --out-dir "./reproduced/multi_calibration"
python research/Run_Level_R2_Standalone_Output/HEMOX_R2_RunLevel_Analysis_Figures.py --pairs-csv "research/Study_HEMOX_Correlation_outputs/native_and_selected_reference_pairs.csv" --out-dir "./reproduced/run_level"
```

Run individual scripts with `--help` for their available options. Certain threshold scripts use fixed `DATA_ROOT` values, and optional verified-offset CSVs may be needed to replace provisional zero-offset comparator alignment.

**Selected evaluation scheme:** development participants P01–P02 (589 pairs) fit fixed SpO₂ M1 offset +0.9067538200 pp, applied as native-minus-offset; HR retains M0. Participants P03–P10 (32 recordings) were not used to fit that coefficient, but were also used for performance/model assessment and are not an untouched external test set for every analytical choice.

**Primary held-out measurements:** 2,202 eligible SpO₂ pairs (selected M1, RMSE 0.968 pp) and 2,255 HR pairs (M0, RMSE 2.495 bpm). Error sign is HEMOX minus SP-20. Paired-point Bland–Altman limits are descriptive since within-participant windows are correlated.

**Synchronization interpretation:** CH1–CH2 `sync_delta_us` and synchronized-record counter analysis must not be conflated with separate BLE packet-delivery audits. External SP-20 pairing accepts unique nearest eligible observations at ±1 s and uses a provisional zero offset where verified timestamps are unavailable.

**Statistical and ethical safeguards:** Never score visualized hold/interpolated estimates as native valid readings. SP-20 agreement is not invasive arterial-reference verification. Confirm participant-data sharing approval before changing repository visibility or sharing run-level exports.
