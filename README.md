# Rule-Based Learning with Cardiorespiratory Physiology

This project investigates **rule-based learning** and its relationship with **cardiac cycle and respiratory rhythms** using behavioral performance, EEG, ECG, and respiration signals. The goal is to understand how physiological states influence learning performance and to develop **machine learning models** that predict learning dynamics.

This work is conducted at the **University of Jyväskylä**.

---

# Project Status

Current stage:

* Experiment design completed
* Rule-based learning paradigm implemented (Presentation)
* Easy, medium, and probabilistic difficulty blocks implemented
* Score-based motivation added
* Pilot measurement completed
* EEG, ECG, respiration recordings implemented
* MATLAB preprocessing pipeline in progress
* Machine learning pipeline planned

---

# Research Objectives

* Investigate rule-based learning performance
* Compare easy vs medium vs probabilistic rules
* Extract learning curves and learning trajectory
* Study cardiac phase effects (systole vs diastole)
* Study respiration phase effects (inspiration vs expiration)
* Analyze ERP components during learning
* Model learning performance using physiological signals
* Predict learning using machine learning

---

# Experiment Design

Participants perform a **rule-based learning task** with three difficulty levels:

### Difficulty Conditions

* Easy rule
* Medium rule
* Probabilistic (difficult) rule

### Measurements

Behavioral:

* Accuracy
* Reaction time
* Learning curve
* Trial-by-trial performance

Physiological:

* EEG
* ECG
* Respiration

Derived measures:

* Cardiac phase (systole / diastole)
* Respiration phase (inspiration / expiration)
* ERP components
* Learning rate
* Performance trajectory

---

# Repository Structure

```
rule-based-learning-physio/
│
├── experiment/              # Presentation experiment files
│   └── presentation/
│
├── matlab/                  # MATLAB analysis code
│   ├── preprocessing/
│   ├── analysis/
│   ├── stats/
│   └── plotting/
│
├── ml/                      # Machine learning models
│
├── scripts/                 # Run pipeline scripts
│
├── docs/                    # Documentation
│
├── results/                 # Figures and outputs
│
└── README.md
```

---

# Experiment Code

Location:

```
experiment/presentation/
```

Contains:

* Presentation (.sce / .pcl) files
* Stimuli
* Instructions
* Trigger definitions
* Timing parameters
* Rule difficulty implementation
* Score feedback system

---

# MATLAB Pipeline

## Preprocessing

Location:

```
matlab/preprocessing/
```

Includes:

* EEG preprocessing
* ECG R-peak detection
* Cardiac phase extraction
* Respiration phase extraction
* Trial segmentation
* Event alignment
* Trigger parsing

---

## Analysis

Location:

```
matlab/analysis/
```

Includes:

* Learning curves
* Rule difficulty comparison
* Correct vs incorrect trials
* ERP extraction
* Block-wise performance
* Trial-by-trial learning
* Behavioral modeling

---

## Physiology Analysis

Includes:

* Cardiac phase locking
* Respiration phase locking
* EEG locked to cardiac cycle
* EEG locked to respiration cycle
* Learning vs physiological phase
* Performance vs physiological phase

---

## Statistics

Location:

```
matlab/stats/
```

Includes:

* Mixed effects models
* Repeated measures analysis
* Cluster permutation tests
* Learning slope comparison

---

# Machine Learning (Planned)

Location:

```
ml/
```

Planned models:

* Predict learning speed
* Predict final performance
* Classify fast vs slow learners
* Predict rule difficulty
* Physiology-based prediction
* Multimodal prediction (EEG + ECG + respiration)

---

# Analysis Pipeline

1. Run experiment (Presentation)
2. Record EEG, ECG, respiration
3. Preprocess physiological signals
4. Extract cardiac phase
5. Extract respiration phase
6. Compute learning curves
7. ERP analysis
8. Statistical analysis
9. Machine learning prediction

---

# Data Organization

Raw data is not stored in this repository.

Expected structure:

```
data/
├── raw/
├── preprocessed/
└── derivatives/
```

---

# Signals Used

EEG
ECG
Respiration

Derived:

* Cardiac phase (systole / diastole)
* Respiration phase (inspiration / expiration)
* ERP components
* Learning curves
* Trial accuracy
* Reaction time

---

# Requirements

MATLAB toolboxes:

* FieldTrip
* Signal Processing Toolbox
* Statistics Toolbox

Optional:

* EEGLAB
* Machine Learning Toolbox

---

# Team
Praghajieeth Raajhen Santhana Goplan
Postdoctoral Researcher
University of Jyväskylä, Finland

Principal Investigators:
Prof. Tiina Parviainen
Prof. Raine Koskimaa

---

# Future Work

* Increase trial counts
* Additional pilot measurements
* Improve motivation using score feedback
* Machine learning prediction
* Physiology-specific learning analysis
* Paper preparation

---

# License

Research use only
