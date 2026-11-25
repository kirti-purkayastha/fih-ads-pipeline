# Mock Statistical Analysis Plan (SAP) - ROC-FIH-SYN-001

**Status:** illustrative, written for a synthetic study. Section numbers below are the ones cited in `spec/adam_spec.psv` and in the code comments, so every derivation can be traced: *code -> spec -> SAP section*.

**Study:** open-label, multi-cohort, 3+3 / rolling-six dose escalation of DRUGX (IV, Q3W) in advanced solid tumours.
**Objectives:** safety and tolerability (primary), PK, pharmacodynamics (CXCL10), preliminary activity (best overall response), exploratory ctDNA (added by Protocol Amendment 1).

## 3. Analysis populations
| Ref | Population | Definition |
|---|---|---|
| 3.1 | Safety (SAFFL) | Received at least one dose of DRUGX. Age group: `>=65` / `<65`. |
| 3.2 | PK (PKFL) | At least one quantifiable (>= LLOQ) post-dose concentration. |
| 3.3 | DLT-evaluable (DLTEVLFL) | Safety population AND (a DLT in the 21-day window OR at least 21 days of follow-up from first dose). |

## 4. Treatment and study-day conventions
- 4.1 `TRTSDT`/`TRTSDTM` = earliest exposure record; `TRTEDT` = latest exposure date.
- 4.2 Study day = date - TRTSDT + 1 for dates on/after TRTSDT; date - TRTSDT before it (there is **no day 0**).
- 4.3 Ongoing subjects are followed to the data cut-off date (2026-03-31).

## 5. Pharmacokinetics
- 5.1 BLQ concentrations are set to 0 for listings and for AUC; excluded from geometric means. "Not done" samples stay missing.
- 5.2 Non-compartmental analysis, cycle 1 profile, actual elapsed time: linear-up/log-down trapezoidal AUC from 0 to the last quantifiable concentration; lambda_z by log-linear regression over the last 3 post-Tmax points, accepted if slope < 0 and adjusted R-squared >= 0.80; AUCinf = AUClast + Clast/lambda_z; CL = Dose / AUCinf. Parameters that fail the lambda_z criteria are reported as not estimable (NE), never imputed.
- 5.3 Dose proportionality: power model, log(parameter) = a + b * log(dose); slope b with 90% CI.

## 6. Safety
- 6.1 Treatment-emergent AE: onset on/after first dose and no later than last dose + 30 days.
- 6.2 Partial onset dates: `YYYY-MM` is imputed to the first of the month; if that month is the month of first dose, use first-dose date. Imputed records are flagged in `ASTDTF`.
- 6.3 DLT: treatment-emergent, grade >= 3, related to study drug, onset within study days 1-21.
- 6.4 AESI: cytokine release syndrome and infusion-related reaction.

## 7. Laboratory
- 7.1 Results converted to standard units before any derivation (haemoglobin g/dL x 10 -> g/L), with reference ranges converted identically. Original value and unit retained in `LBORRES`/`LBORRESU`.
- 7.2 Baseline = last non-missing value on or before first dose.

## 8. Exploratory biomarkers (semi-structured assay data)
- 8.1 Assay results arrive as nested JSON (one record per sample), flattened to one row per sample and analyte; ctDNA value = maximum variant allele fraction across called variants.
- 8.2 BLQ results are imputed as LLOQ/2 and flagged (`BLQFL`, `DTYPE`).
- 8.3 Change/ratio to baseline is **not** derived when baseline is missing or BLQ. Protocol Amendment 1 (2025-09-01) introduced ctDNA sampling, so subjects dosed earlier have no ctDNA baseline by design; they are retained for descriptive display but excluded from baseline-relative summaries.

## 9. Tumour response (simplified RECIST 1.1, target lesions only)
- 9.1 Visit response: PD if a new lesion is present or the sum of diameters rises >= 20% **and** >= 5 mm over the smallest prior value (nadir, baseline included); otherwise PR if change from baseline <= -30%; otherwise SD. CR cannot be assigned because non-target lesions are not collected.
- 9.2 Best overall response = best post-baseline visit response (PR > SD > PD); NE if no post-baseline assessment. Responses are unconfirmed. ORR = proportion PR; DCR = proportion PR or SD.

## 10. Quality control and traceability
- 10.1 Every analysis dataset conforms to the variable-level specification (`spec/adam_spec.psv`): names, order, labels, types.
- 10.2 ADSL, ADPP (NCA) and ADRS are independently double-programmed from raw data in base R and compared value-by-value (tolerance 1e-6).
- 10.3 62 automated integrity checks run on every execution (keys, referential integrity, raw-to-ADaM reconciliation, date logic, derivation consistency, controlled terminology, cross-dataset consistency). Any failure stops the pipeline with a non-zero exit status.
- 10.4 Each run writes a log, session information and an MD5 manifest of every input and output file.
