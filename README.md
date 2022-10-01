# Confounder adjustment using disease risk score analysis
Propensity score analysis is a common approach to addressing confounding in non-randomized studies. Its implementation, however, requires important assumptions (e.g., positivity). The disease risk score is an alternative balancing score that relax some of these assumptions. Like the propensity score, the disease risk score summarizes multiple confounders into a single score, on which conditioning (e.g., by subclassification or matching) allows the estimation of causal effects. However, conditioning by matching relies on arbitrary choices for pruning out data (e.g., matching ratio, algorithm and caliper width) and may be computationally demanding. Alternatively, weighting methods, common in propensity score analysis, are easy to implement and may entail fewer choices, yet none have been developed for the disease risk score. We present two weighting approaches: one derives directly from inverse probability weighting (IPW); the other named target distribution weighting (TDW) relates to importance sampling. We empirically show IPW and TDW display a performance comparable to matching techniques in terms of bias but outperform them in terms of efficiency (mean squared error) and computational speed (up to >1200 times faster in an illustrative study). We illustrate implementation of the methods in two case studies where we investigate placebo treatments for multiple sclerosis and the administration of Aspirin in stroke patients.

**Authors**: [Tri-Long Nguyen](https://orcid.org/0000-0002-6376-7212), [Thomas P.A. Debray](https://orcid.org/0000-0002-1790-2719), Bora Youn, Gabrielle Simoneau, Gary S. Collins
