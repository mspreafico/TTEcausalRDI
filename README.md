# TTEcausalRDI

Code for assessing the **causal effects** of chemotherapy **Received Dose Intensity (RDI)** on survival outcomes in osteosarcoma patients using a **Target Trial Emulation** approach.

- **Causal framework**:
  - Three RDI-based exposure strategies: 1) *standard*, 2) *reduced*, and 3) *highly-reduced* RDI (exposure variable *A*).
  - Outcome of interest is Event-Free Survival (EFS; variable *T*) defined as time from the end of therapy until the first event (local recurrence, evidence of new or progressive metastatic disease, second malignancy, death, or a combination of those events) or censoring at last contact.
  - Investigations are conducted between subgroups of patients characterised by *poor* or *good* Histological Responses (HRe; effect modifier *V*), i.e., the strongest known prognostic factor for survival in osteosarcoma.
  - Exposure-outcome confounders (*L*) are considered.
- **Research questions**: *Does reduced chemotherapy RDI lead to an improvement in EFS of patients with osteosarcoma? Does this effect vary among subjects characterized by different HRe?*
- **Methodology**: Inverse Probability of Treatment Weighting (IPTW) is first used to transform the original population into a pseudo-population which mimics the target randomized cohort. Then, a Marginal Structural Cox Model (Cox MSM) with effect modification is employed. Conditional Average Treatment Effects (CATEs) are ultimately measured as the difference between the Restricted Mean Survival Time of *reduced/highly-reduced* RDI strategy and the *standard* one. Confidence Intervals for CATEs are obtained using a novel IPTW-based bootstrap procedure.


### Reference

Spreafico M., Ieva F., Fiocco M. (2024). Causal effect of chemotherapy received dose intensity on survival outcome: a retrospective study in osteosarcoma. *BMC Medical Research Methodology*, 24, 296. https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-024-02416-x


### Data Availability

#### Original data
Data from the control arms of MRC BO03 and MRC BO06 (EORTC 80861 and 80931, respectively) clinical trials for osteosarcoma are analysed.
Data are not publicly available due to confidentiality and privacy restrictions.

#### Toy dataset
We provide a toy dataset in order to allow researchers who want to replicate the same analysis to properly get how the code has to be run and how results are displayed and should be read.


## Description

- Files:
  - **01 IPTW diagnostics.R**: Selection of the best IPTW model specifications (in terms of confounding covariates for the denominators) to determine the subject-specific standardized weights.
  - **02 Cox MSM.R**: Implementation of marginal structural Cox model with effect modification by IPTW (Cox MSM). Comparison with relative unweighted Cox model.
  - **03 CATEboot.R**: Computation of Conditional Average Treatment Effects (CATEs) with relative 95% bootstrap CIs.
- Sub-folder **./data/** contains the toy dataset to run the code with relative legend. [Original data are not publicy available due to privacy restrictions]
- Sub-folder **./functions/** contains utils functions
  - **tools_diagnostics.R**: Utils functions for IPTW diagnostics. Called in "01 IPTW diagnostics.R".
  - **tools_bootstrap.R**: Utils functions for estimating bootstrap CIs. Called in "03 CATEboot.R"
  - **tools_reshapes.R**:  Utils function for CATE result reshaping. Called in "03 CATEboot.R"
- Sub-folder **./results/** contains results generated by running the code. **Note that these represent the results for the fake toy dataset, NOT the real results.**



## Software
- R software.
- Packages: iptw, survival, splines, data.table, ggplot2, ggpubr


## Trial References
Details related to the primary analysis of the MRC BO03/EORTC 80861 Randomized Controlled Trial can be found in:

- Lewis I.J *et al.* (2000). Received Dose and Dose-Intensity of Chemotherapy and Outcome in Nonmetastatic Extremity Osteosarcoma. *Journal of Clinical Oncology*,18(24): 4028-4037. doi: 10.1200/JCO.2000.18.24.4028


Details related to the primary analysis of the MRC BO06/EORTC 80931 Randomized Controlled Trial can be found in:

- Lewis I.J. *et al.* (2007). Improvement in histologic response but not survival in osteosarcoma patients treated with intensified chemotherapy: a randomized phase III trial of the European Osteosarcoma Intergroup. *Journal of the National Cancer Institute*, 99(2):112-128. doi: 10.1093/jnci/djk015


(Last update: November 24th, 2023)
