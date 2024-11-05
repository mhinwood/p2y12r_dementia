# P2Y12 Receptor Inhibitors and Post-Stroke Dementia Risk

## Overview
This repository contains the R code used to analyze the relationship between P2Y12 receptor inhibitor treatment and dementia risk in a nationwide Swedish cohort study of stroke survivors. The analysis uses target trial emulation methodology to investigate whether post-stroke treatment with P2Y12R inhibitors affects the risk of dementia or mild cognitive disorder (MCD).

## Summary
Animal studies suggest that antiplatelet drugs targeting the purinergic P2Y12 receptor (P2Y12R) may influence brain function through hypothesised effects on microglial activity. This study investigated whether post-stroke treatment with P2Y12R inhibitors is associated with dementia risk using nationwide Swedish registry data from 2006-2016, following individuals for up to 5 years after stroke diagnosis.

## Repository Structure
This repository contains four main R scripts:

1. **Data Management**
   * Creates wide and long datasets for analyses and prepares data structures for subsequent analytical steps

2. **Total Effect Analysis**
   * Implements Inverse Probability of Treatment Weighting (IPTW)
   * Analyses the total effect of treatment on dementia and death as separate outcomes

3. **Controlled Direct Effect Analysis**
   * Implements both IPTW and Inverse Probability of Censoring Weighting (IPCW)
   * Analyses the controlled direct effect of treatment on dementia accounting for death as a competing risk

4. **Utils**
   * Contains utility functions used across other scripts

## Methods
The study employed:
* Target trial emulation framework
* Inverse probability weighting techniques
* Competing risk analysis
* 5-year follow-up period

## Acknowledgments
* Code structure and analytical approach adapted from [Competing Risks Dementia Repository](https://github.com/palolili23/competing_risks_dementia)

## Contact
📧 madeleine.hinwood@newcastle.edu.au
🔗 https://orcid.org/0000-0002-2225-973X
