# R \<-\> SAS equivalents used in this project

Pharma teams work in both R and SAS. The logic here is language-neutral; this table shows how each element maps.

| This project (R) | SAS equivalent |
|------------------------------------|------------------------------------|
| `apply_spec()` - enforce names/order/labels/types from `spec/adam_spec.psv` | `ATTRIB` / `LENGTH` / `KEEP` driven from the spec; Pinnacle 21 conformance check |
| `qc_compare()` - production vs independent QC, tolerance 1e-6 | `PROC COMPARE ... CRITERION=1E-6 METHOD=ABSOLUTE` |
| `study_day()` - no day 0 | `ADY = ADT - TRTSDT + (ADT >= TRTSDT);` |
| `flag_last()` baseline flag | `PROC SORT; BY USUBJID PARAMCD ADT; ... IF LAST.PARAMCD` on pre-dose records |
| `nca_profile()` linear-up/log-down AUC | `PROC EXPAND` / macro; validated NCA tools (e.g. Phoenix WinNonlin) are the production standard |
| Pipeline run log + MD5 manifest | SAS log scan + `%PUT` audit macros; file hashes from the validated environment |
| `testthat` seeded-error tests | Negative test cases in the programming validation plan |
| ggplot2 swimmer / waterfall | `PROC SGPLOT` (`HBARPARM`, `SCATTER`, `VBARPARM`) |
| Shiny exploration app | SAS Viya / R Shiny front end on the same ADaM |
