---
title: "Reanalyzing the Gautret et al. data on Covid-19 treatment. Asking for feedback."
output: html_document
---

Recently, a [paper](http://dx.doi.org/10.1016/j.ijantimicag.2020.105949) by Gautret et al. made waves with their strong support of using Hydroxychloroquine (HCQ) and Azithromycin (AZ) for Covid-19 treatment based on a small non-randomized study.  I wondered what would be a good way to analyze the data they reported? I made some attempts and am looking for feedback. 

Note that the paper was criticized on many methodological grounds. Most notably: patients with severe outcome (which were only in the HCQ group) were excluded, the treatment and control groups were recruited at different hospitals and obviously differed even before treatment and also some ethical and procedural considerations, see e.g. the [Pubpeer coments](https://pubpeer.com/publications/B4044A446F35DF81789F6F20F8E0EE) for more. 

In my analysis I am ignoring many of those shortcoming as most can't be solved by statistics. So wanted to note that this brings the analysis further from "practical" territory and more into "excercise in statistics" territory. 

## The dataset
The authors, in a rare and positive move (for a medical paper) provided reasonably detailed patient-level data. The data has a lot of structure:
- 42 patients,  16 served as controls, 20 got HCQ, 6 got HCQ + AZ,  6 were excluded because loss of followup, of the excluded 1 died and additional 3 were transfered to ICU (all on HCQ, not certain if also on AZ)
- For each patient we have 6-day time series of viral load (real-time PCR) measurements. Those could turn either negative or report Ct number between 0 and 35 - higher Ct means lower viral load, 35 is the detection limit.
  - For some patients some days are missing
  - For some patients in the control group the viral load is only reported as positive/negative
- The time series starts at treatment initiation (HCQ, HCQ + AZ) or inclusion in the study (Control)
 - The patient groups in the time between symptom onset and treatment initiation / inclusion in the study
 - We also know the age of the patients, and we know age is an important predictor of disease severity.
- For the excluded patients, we know their day-by-day history post-treatment, but we do not know their time from symptom onset.

This is one way to look at the data I've found helpful:

![image|500x500](upload://8asmobqs5sXcVrl40Vpxn3osBiK.png) 

Each line is a single patient, for the treatment groups, the line starts at first dose. Viral expression is computed as 35 - CT, so that 0 corresponds to the detection limit. The control group contains some patients where the PCR is reported only as positive/negative. For this figure, positive tests were imputed as 10 (CT = 25) plus tiny noise. The excluded cases (including 4 HCQ patients with severe outcomes) are excluded as in the paper.

## Previous analyses

The analysis in the paper is quite ad-hoc and IMHO suboptimal - it compares the proportions of positive at each day post-treatment with a Chi-square test. A [reanalysis using a frequentist survival model](https://github.com/andrewlover/HCQ_AZ_COVID_19/blob/master/Preprint_Update1_22_Mar_HCQ_AZ_COVID_19.pdf) (with "survival" as "stay infected") and [one with Bayes factors](https://osf.io/7ax9w/) was also made. But both reanalyses still group together the patients whose condition worsened (death/ICU) with those that stayed PCR-positive and ignore the signal provided by the actual viral load measurements beyond positive/negative - and as you see in the plot, for most patients that got healthy, the decrease in the viral load was apparent for a few days before.

## The challenges

I attempted to build a model that would take advantage of the full structure. This turned out to be more challenging than I expected. The biggest problem came from the fact that there are both continuous dynamics (viral load).

Also note that the viral load numbers are already on the log scale, with zero being the detection limit which is a somewhat arbitrary value.

## Attempt 1 - Discrete dynamics



As a bonus, the group should soon release a full report of their study with longer followup, so it should be possible to validate the model via predictions on out-of-sample data.