---
title: "Reanalyzing the Gautret et al. data on Covid-19 treatment. Asking for feedback."
output: html_document
---

Recently, a [paper](http://dx.doi.org/10.1016/j.ijantimicag.2020.105949) by Gautret et al. made waves with their strong support of using Hydroxychloroquine (HCQ) and Azithromycin (AZ) for Covid-19 treatment based on a small non-randomized study.  I wondered what would be a good way to analyze the data they reported? I made some attempts and am looking for feedback. 

Note that the paper was criticized on many methodological grounds. Most notably: patients with severe outcome (which were only in the HCQ group) were excluded, the treatment and control groups were recruited at different hospitals and obviously differed even before treatment and also some ethical and procedural considerations, see e.g. the nice review by [Darren L. Dahly et al.](https://doi.org/10.5281/zenodo.3724167) or [Pubpeer coments](https://pubpeer.com/publications/B4044A446F35DF81789F6F20F8E0EE) for more. 

In my analysis I am ignoring some of those shortcoming as they can't be solved by statistics. But I believe parts of them can be handled by a more custom-tailored analysis. But I wanted to note that the shortocomings I can't address brings the analysis further from "practical" territory and more into "excercise in statistics" territory. 

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

Each line is a single patient, for the treatment groups, the line starts at first dose. Viral expression is computed as 35 - CT, so that 0 corresponds to the detection limit as used in the paper. The control group contains some patients where the PCR is reported only as positive/negative. For this figure, positive tests were imputed as 10 (CT = 25) plus tiny noise. The excluded cases (including 4 HCQ patients with severe outcomes) are excluded as in the paper.

Also note that the viral load numbers are already on the log scale (because of the way rtPCR works), with zero being the detection limit which is a somewhat arbitrary value.


## Previous analyses

The analysis in the paper is quite ad-hoc and IMHO suboptimal - it compares the proportions of positive at each day post-treatment with a Chi-square test. A [reanalysis using a frequentist survival model](https://github.com/andrewlover/HCQ_AZ_COVID_19/blob/master/Preprint_Update1_22_Mar_HCQ_AZ_COVID_19.pdf) (with "survival" as "stay infected") and [one with Bayes factors for odds ratios](https://osf.io/7ax9w/) was also made. But both reanalyses still group together the patients whose condition worsened (death/ICU) with those that stayed PCR-positive and ignore the signal provided by the actual viral load measurements beyond positive/negative - and as you see in the plot, for most patients that got healthy, the decrease in the viral load was apparent for a few days before, so there is likely some signal to consider here.

Also both analyses ignore that groups differ in their timing from symptom onset to enrollment into the study, which is IMHO one of the deficiencies that can be somewhat accounted for statistically. 

## The challenges

I attempted to build a model that would take advantage of the full structure. This turned out to be more challenging than I expected. The biggest problem came from the fact that there are both continuous dynamics (viral load) and discrete dynamics - becoming healthy/severely ill.

Also, if we want to take into account time from symptom onset, we must handle the fact that for the originally excluded patients, we must treat this time as missing data.

## Attempt 1 - Discrete dynamics

In my first attempt (and the only model that currently completely works), I decided to model the underlying dynamic as discrete via a highly constrained Hidden markov model and treat the continous viral loads as emerging from this discrete model. The true states are ordered and consist of one _healthy_ state, an arbitrary number of _ill_ states and one _severe_ state (representing ICU or death as I don't think I can further distinguish those from the data). The number of _ill_ states is given as input.  

### Transition model

The _healthy_ and _severe_ states are terminal, i.e. once those states are entered, the model stays in this state for all future time points. So we define only transitions from _ill_ states. Since the states are ordered, we use the cumulative logistic model to define transitions. We have one global set of thresholds $c_1,...c_{k+1}$ where $ k $ is the number of _ill_ states and state-specific intercepts $s_1, ...,s_k}$. Of those, one is arbitrarily set to $0$ and the whole sequence is ordered. Given a time and patient dependent linear predictor $\mu_{p,t}$, the transition probabilities from and _ill_ state $i$ are 

$$
P(X_{p,t+1} | X_{p, t} = i) = CumulativeLogistic(X_{p,t+1} | \mu_{p,t} + s_i, c)
$$

Some (IMHO desirable) properties of this scheme:
- Doesn't have too much free parameters
- "Higher" _ill_ states are more likely to progress to _severe_ state and less likely to _healthy_ state, but the extent to which this happens is learned from data and can be negligible.

Initially, the model is assumed to be in any of the _ill_ states with equal probability.

### Observation model

First, the _severe_ state is considered to be directly observable, i.e. any observation of severe condition is treated as 100% evidence that the model is in _severe_ state and, conversely observation of any other condition is a 100% evidence that the state is different.

Each state is associated with a mean viral load $v_1, ..., v_k$, that is currently given as data (I simply split the observed range into $k$ equal parts). In addition to this, there are three parameters that drive the observation model: $\sigma, sens, spec$, corresponding to variations expected in viral load, sensitivity and specificity of the PCR tests respectively. 

When the model is in _healthy_ state, then a PCR negative is observed with $spec$ probability and with probability of $1 - spec$  a state $i$  is chosen uniformly from all the _ill_ states and the observed viral load is drawn from $N_0(v_i, \sigma)$ where $N_0$ is the normal distribution truncated at zero.

When the model is in one of the _ill_ states, then a PCR negative is observed with $1 - sens$ probability and with probability of $1 - sens$ the viral load is drawn from $N_0(v_i, \sigma)$.

For observations that only indicate positive PCR, the $N_0$ term is omitted from the likelihood.

I require $sens, spec \in [1/2,1]$ as otherwise the model can have an additional mode where the role if _ill_ and _healthy_ states is swapped.

### Missing data

For the patients with unknown time from symptom onset, I specify a maximum such time and then marginalize over the likelihood for each of this cases, treating their probabilities as a simplex parameter. The simplex is kept with the default flat prior. 

### Linear predictor

I have the usual dummy coding for treatments and also played with using time from symptom onset as a predictor. I have not so far incorporated any varying intercepts/effects as there are no obvious problems with the fits.


### Plots

Below are plots of fit to selected patients:

What do you think? Do you have any feedback?

The whole code and data can be found at https://github.com/martinmodrak/covid19-pharmacotherapy/, the model file is 
https://github.com/martinmodrak/covid19-pharmacotherapy/blob/master/hmm.stan, a simulator to create fake data is at https://github.com/martinmodrak/covid19-pharmacotherapy/blob/master/R/hmm_devel_tools.R, SBC and other checks are at https://github.com/martinmodrak/covid19-pharmacotherapy/blob/master/hmm_model_devel.Rmd and the code applying the model to the Gautret dataset is at https://github.com/martinmodrak/covid19-pharmacotherapy/blob/master/gautret.Rmd

## Attempt 2 - Continous dynamics

I also played with treating the disease progression as a simple dynamic system (basically a random walk + effect of the linear predictor), that is noisily observed via the viral load. I don't have this model working really well as it is IMHO implausible to model the severe/healthy outcomes as the viral load crossing a certain threshold, because we need the model to never return back... Also, it seems implausible to tie the probability of progression to severe state with the variability in viral load. So I am considering to use the dynamic system as a predictor in an separate model that would guide the transitions...  

I also tried to have a tri-state HMM and the linear model only loosely connected by using the same linear predictor (+a separate intercept) for both, but that could produce implausible simulated e.g. viral load increasing for everybody while everybody turned healthy.

Still working and thinking how to do this well, any feedback welcome, but implementing Kalman filter was fun.

## What's next?

I am working on additional posterior predictive checks and towards a multiverse analysis varying some aspects of the model (priors, number of ill states, whether the excluded cases got Azithromycin). I am also trying to use my (limited but nonzero) clinical contacts to get hands on additional time-course datasets with experimental treatments. I guess the model could be reasonably well extended to handle multiple clinics (via a varying intercept/effect) and other indirect measures of disease progression (e.g. radiological findings).

If possible, the progress of patients in the ICU should also be modelled, but that would need to be basically a separate model as it is implausible the disease dynamics are the same (especially, once the case progresses to severe, the relationship between viral load and outcomes becomes tenuous or completely non-existent).

As a bonus, the group should (in theory) soon release a full report of their study with longer followup, so it should be possible to validate the model via predictions on out-of-sample data.