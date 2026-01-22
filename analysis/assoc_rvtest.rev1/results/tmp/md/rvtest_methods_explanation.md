# Methods & Results: Rare Variant Association Testing

## 1. Overview
This document explains the **methodology** and **interpretation** of the rare variant association analysis for the CTEPH study.
The goal is to determine if rare genetic mutations in specific genes are associated with CTEPH status. We use the tool `rvtest` to perform this analysis using two complementary strategies: **Burden Test** and **SKAT-O**.

---

## 2. The Core Concept (Intuition before Math)

Before diving into formulas, it helps to understand the logical flow of the analysis:

1.  **Step 1: The "Baseline" Prediction (Null Model)**
    *   First, we ignore genetics. We try to predict who has CTEPH based *only* on their Sex and Ancestry (Principal Components).
    *   We calculate a "Residual" (Error) for each person:
        *   *Positive Residual*: The person has CTEPH, but our baseline model predicted they were healthy. **Something else (genetics?) must be causing it.**
2.  **Step 2: Check the Gene**
    *   We look at a specific gene. Does carrying rare mutations in this gene correlate with those "Positive Residuals"?
3.  **Step 3: Choose a Strategy**
    *   **Burden Test**: Assumes all mutations in the gene are **Bad**. It sums them up.
    *   **SKAT**: Assumes mutations are **Mixed** (some bad, some protective, some noise). It looks for *dispersion* (variance).
    *   **SKAT-O**: Tries **Both** and picks the winner.

---

## 3. Statistical Models

### 3.1 The Null Model (The Baseline)
We fit a Logistic Regression to remove the effect of confounders.

$$logit(\pi_i) = \alpha + \gamma Sex_i + \sum_{k=1}^{10} \delta_k PC_{ik}$$

*   **Result**: We get a **Residual** ($Y_i - \pi_i$) for every individual. This represents the "unexplained risk" that we hope the gene will explain.

### 3.2 The Three Testing Approaches (Detailed Mechanics)

Suppose we have a gene with **3 variants** in a patient "John".
*   Variant 1: Highly Damaging (Score +2), Very Rare ($w=1.0$)
*   Variant 2: Protective (Score -2), Less Rare ($w=0.5$)
*   Variant 3: Noise (Score 0), Very Rare ($w=1.0$)

#### A. Burden Test (CMC) $\rightarrow$ "The Cumulative Load"
*   **Analogy**: A **"Light Switch"**.
*   **Logic**: "Is the gene broken? Yes or No?"
*   **rvtest Weighting**: **None / Flat ($w=1$)**.
    *   `rvtest` assumes every variant contributes equally to the "broken" status.
*   **Detailed Calculation**:
    1.  **Check**: Does John have *any* of these variants? $\rightarrow$ Yes.
    2.  **Collapse**: $Score_{John} = 1$ (John is a carrier).
    3.  **Test**: Compare the rate of "Carriers" in Cases vs Controls.
    *   *Note*: The direction (+2 vs -2) and weight (1.0 vs 0.5) are **ignored**.
*   **Key Assumption**: **Homogeneity**. All variants act in the same direction (Risk).

#### B. SKAT $\rightarrow$ "The Variance Test"
*   **Analogy**: **"Scatter Detection"**.
*   **Logic**: "Are the effects in this gene dispersed?"
*   **rvtest Weighting**: **Beta Distribution ($Beta(1, 25)$)**.
    *   `rvtest` default: Weights are calculated as $w_j = Beta(MAF_j, 1, 25) \approx (1-MAF)^{24}$.
    *   Rare variants get $w \approx 1$, common variants get $w \approx 0$.
*   **Detailed Calculation**:
    1.  **Score**: Calculate signal for each variant $j$ by checking its correlation with the unexplained risk:
        $$U_j = \sum_{individuals} G_{ij} \times \text{Residual}_i$$
        *   **$G_{ij}$**: The **Genotype** (0, 1, or 2) from the VCF file.
        *   **Residual$_i$**: The **Unexplained Risk** from the Null Model (See Sec 3.1).
    2.  **Square**: Square the signal to remove direction ($U_j^2$).
    3.  **Weight**: Apply the Beta weights ($w_j^2 U_j^2$).
    4.  **Sum**: $Q_{SKAT} = (1.0 \times 2^2) + (0.5 \times -2^2) + \dots$
    *   *Note*: The +2 and -2 **accumulate** instead of canceling out.
*   **Key Assumption**: **Heterogeneity**. Variants can have different effects (+/-).

#### C. SKAT-O $\rightarrow$ "The Smart Adaptor"
*   **Analogy**: A **Hybrid Car**.
*   **Logic**: "Let's try mixing the two approaches."
*   **rvtest Weighting**: **Beta Distribution ($Beta(1, 25)$)**.
    *   Used for *both* the SKAT part and the Burden part inside SKAT-O.
*   **Detailed Calculation**:
    1.  Calculate $Q_{SKAT}$ ($\sum w_j^2 U_j^2$).
    2.  Calculate $Q_{Burden}$ (Weighted Sum: $(\sum w_j U_j)^2$).
    3.  **Mix**: Combine them using $\rho$ (0 to 1).
        $$Q_{\rho} = (1-\rho) Q_{SKAT} + \rho Q_{Burden}$$
    4.  **Optimize**: Find the $\rho$ that gives the strongest signal.
    *   If signals cancel out in Burden, SKAT wins ($\rho=0$).
    *   If signals add up effectively, Burden wins ($\rho=1$).


---

## 4. Methodological Details (The "How-To")

### 4.1 Comparison of Methods

| Feature | Burden Test (CMC) | SKAT | SKAT-O |
| :--- | :--- | :--- | :--- |
| **Statistical Model** | **Fixed Effect** | **Random Effect** | **Unified** |
| **Calculation** | **Correlation with Binary Score**<br>(Collapse to 0/1, then Test) | **Weighted Sum of Squares**<br>$Q = \sum w_j^2 U_j^2$<br>(Using Score $U_j$ from Sec 3.2) | Best of both worlds |
| **Constraint** | Assume all variants act in **Same Direction**. | Allows **Mixed Directions** (+/-). | **Adaptive**. |
| **Effect Size** | **YES** (Odds Ratio). | **NO**. | **NO**. |

### 4.2 Explicit Comparison of "Burden" Types

Why did we use **CMC** instead of **Zeggini** or **Madsen-Browning (MB)**?

| Feature | **CMC** (Our Choice) | **Zeggini** | **Madsen-Browning** | **SKAT-O Burden** |
| :--- | :--- | :--- | :--- | :--- |
| **Logic** | **Binary Switch** | **Counter** | **Weighted Sum** | **Weighted Sum of Scores** |
| **Rule** | "Do you have *any* mutation?"<br>(Yes/No) | "How *many* alleles do you have?"<br>(Unweighted Sum) | Sum, but rarer variants get **Inverse Variance Weight**. | Linear Sum of Scores $U_j$:<br>$(\sum w_j U_j)^2$ |
| **Best For** | **Dominant Diseases**.<br>(1 hit is enough to get sick) | **Additive Risk**.<br>(2 hits are 2x worse) | If you assume rarer variants are *much* stronger, based on control data. | If you want stable, theoretical weighting for rare variants. |
| **Weighting** | **Flat (Unweighted)**<br>$w_j = 1$<br>(A singleton counts as much as a 1% variant) | **Flat (Unweighted)**<br>$w_j = 1$<br>(Same as CMC) | **Inverse Variance**<br>$w_j \propto \frac{1}{\sqrt{MAF(1-MAF)}}$<br>(Data-Driven: Rare in *this* dataset = Higher Weight) | **Beta Distribution**<br>$w_j = Beta(MAF, 1, 25)$<br>$\approx (1-MAF)^{24}$<br>(Theoretical: MAF < 0.01 $\rightarrow$ Weight > 0.8) |

### 4.3 Technical Note: SKAT-O Components vs. Standalone Tests

It is common to ask: *"Is the Burden test inside SKAT-O the same as the Burden test I ran separately?"*

**1. The SKAT Component ($\rho=0$) $\rightarrow$ SAME**
*   **Standalone SKAT**: Uses the variance component score test with Beta[1,25] weights.
*   **SKAT-O ($\rho=0$)**: Uses the **exact same** calculation.
*   *Result*: If you run `rvtest --test skat`, the P-value will match the $\rho=0$ component inside SKAT-O.

**2. The Burden Component ($\rho=1$) $\rightarrow$ DIFFERENT**
*   **Standalone Burden (CMC)**:
    *   **Method**: **Binary Collapsing**. It collapses all variants to a simple 0/1 vector.
    *   **Weighting**: **None (Flat)**. All variants are treated equally.
    *   *Consequence*: A singleton mutation has the *same impact* as a mutation with MAF=1% or 5%.
*   **SKAT-O Burden ($\rho=1$)**:
    *   **Method**: **Weighted Sum of Scores**. It calculates the squared linear sum: $(\sum w_j U_j)^2$.
    *   **Weighting**: **Yes, Strongly MAF-Dependent**. It uses the Beta density function ($Beta(MAF, 1, 25)$).
        *   **Justification (Negative Selection)**: We assume that variants with lower MAF are more likely to be deleterious (damaging), because natural selection removes bad variants from the population.
        *   **Behavior (Beta[1,25] weights)**:
            *   **Singleton (MAF $\approx$ 0)**: Weight $\approx$ 1.0 (Maximum Impact).
            *   **Ultra-Rare (MAF = 0.001)**: Weight $\approx$ 0.98 (Still Very High).
            *   **Rare Threshold (MAF = 0.01)**: Weight $\approx$ 0.78 (Already dropped by ~22%).
    *   *Consequence*: SKAT-O's burden test is heavily biased towards **ultra-rare** variants. Even at the generic "rare" cutoff of 1% (MAF=0.01), the weight has already decreased significantly compared to a singleton. If the signal is driven by variants near 0.01, CMC (which treats them as 1.0) might perceive a stronger signal than SKAT-O.
*   *Result*: They often yield different P-values. The SKAT-O Burden is usually more sensitive to ultra-rare variants due to the weighting.

---

## 5. Result Interpretation Guide

### 5.1 Output Files and Columns

| Output File | Key Column | Meaning for Human Readers |
| :--- | :--- | :--- |
| **Burden (CMC)** | `Pvalue` | The probability that valid association is random. <br> **Small P** = Highly Significant (Carrying the gene burden increases risk). |
| | `NonRefSite` | The total number of mutations found in this gene across all people. |
| **SKAT-O** | `Pvalue` | The final significance *after* correcting for the complexity of the test. |
| | `rho` ($\rho$) | **The Diagnosis of the Gene**: <br> $\rho=1$: "This gene behaves like a classic Burden gene." <br> $\rho=0$: "This gene has complex/mixed effects." |

### 5.2 Where is the "Beta" (Effect Size)?

*   **Can ALL Burden Tests calculate a Beta?**
    *   **Yes, mathematically.** Since all Burden tests compress the gene into a single score variable ($X'$), you can always run a regression ($Y \sim \beta X'$) to get a Beta.
    *   **But interpretation varies wildly**:
        *   **CMC (Binary)**: $\beta$ = Log Odds Ratio of **Carriers vs Non-Carriers**. This is highly interpretable and clinically standard ("Risk increases 5-fold").
        *   **Zeggini (Count)**: $\beta$ = Log Odds Ratio **per Allele** ("Risk increases 2-fold for every extra mutation").
        *   **Weighted Burden (e.g. SKAT-O $\rho=1$, Madsen-Browning)**: $\beta$ = Log Odds Ratio **per unit of Weighted Score**. Since the score is an abstract sum of weights (e.g., 0.98 + 0.45 + ...), saying "Risk increases for every 1.0 score" is **hard for humans to interpret clinically**.
*   **SKAT / SKAT-O (Variance Tests)**:
    *   They do **NOT** have a single Beta.
    *   *Why?* Because they model the **variance** of effects, not the mean. They allow different variants to have different directions (+/-), so a single "Average Beta" would be meaningless (e.g., average of +2 and -2 is 0).

---

## 6. Methodological QA (Why is it like this?)

### 6.1 Strategy Rationale

#### 1. Why do we combine CMC Burden and SKAT-O?
We deliberately paired these two methods to balance **Discovery Power** with **Clinical Interpretability**:

*   **SKAT-O (The "Hunter")**: Its primary role is **Discovery/Sensitivity**. Since we do not know apriori if CTEPH genes behave as "Cumulative" (Burden) or "Dispersed" (Variance), SKAT-O scans all possibilities. It adapts to the data ($\rho$) to ensure we don't miss complex signals that simple tests might fail to see (Minimizing False Negatives).
*   **CMC (The "Measurer")**: Its primary role is **Quantification**. SKAT-O provides a P-value but no clear measure of risk magnitude. CMC forces the data into a binary model (Carrier vs Non-Carrier) specifically to calculate an **Odds Ratio (OR)**. Even if CMC is less significant than SKAT-O, the OR provides the critical clinical context (e.g., "Carriers are at 5x risk").

#### 2. Why did we NOT run a standalone SKAT?
*   **Redundancy**: SKAT is **already integrated** within SKAT-O.
    *   SKAT-O checks the SKAT model ($\rho=0$) as part of its optimization process.
    *   If a gene is best explained by the SKAT model (pure variance), SKAT-O will automatically select $\rho=0$ and yield a result equivalent to running SKAT alone.
*   **Statistical Efficiency**: Running SKAT as a third independent test would add no unique information but would effectively increase our Multiple Testing Burden (requiring a stricter P-value threshold). By using SKAT-O, we get the benefit of SKAT without the penalty of an extra test.

### 6.2 Why is the SKAT-O P-value different from Burden (even if $\rho=1$)?
You might see: *Burden P = 1.0e-5*, but *SKAT-O P = 2.0e-5*. Why?

*   **Reason 1: The "Search Penalty" (P-value $\uparrow$)**
    *   Burden Test asks **1 question**: "Is the Burden significant?"
    *   SKAT-O asks **11 questions**: "Is $\rho=0$ significant? ... Is $\rho=1$ significant?"
    *   Because SKAT-O asks more questions, statistics demands we **penalize** the final P-value to prevent cheating. This makes the P-value slightly larger (less significant).

*   **Reason 2: The "Weighting Boost" (P-value $\downarrow$)**
    *   CMC is **Unweighted** (1 mutation = 1 score).
    *   SKAT-O is **Weighted** (Beta Weights). Even at $\rho=1$, it gives extra points for *extremely rare* variants.
    *   If your gene is driven by ultra-rare variants, SKAT-O's weighting might make it **more significant** than CMC, despite the penalty.

### 6.3 Why does only Burden Test provide a Beta (Effect Size)?
It is often frustrating that SKAT/SKAT-O only gives a P-value, but no Odds Ratio. This is due to the mathematical design:

*   **Burden Test (CMC)**: 
    *   **Simplification**: It forces the gene into a **single variable** ($X'$: Carrier vs Non-Carrier).
    *   **Result**: Since there is only one variable, we can calculate **one Beta** ($\beta$). 
    *   **Meaning**: "On average, being a carrier increases risk by $\exp(\beta)$."

*   **SKAT / SKAT-O**:
    *   **Complexity**: They acknowledge that a gene contains many variants ($G_1, G_2, \dots, G_n$), and each variant has its **own** effect size ($\beta_1, \beta_2, \dots, \beta_n$).
    *   **Problem**: We cannot accurately estimate 50 different Betas for one gene.
    *   **Solution**: Instead of estimating the *Mean* effect (Beta), SKAT tests the **Variance** of the effects ($\tau$).
    *   **Result**: It answers "Is there significant genetic activity here?" but cannot give a single directional number because the variants might be acting in different directions (some increasing risk, some decreasing). A single Beta would be misleading (e.g., average of +2 and -2 is 0).

## 7. Mathematical Specification & Post-Analysis

### 7.1 The Regression Models (Formal Specification)

To ensure precision, we define the exact Generalized Linear Models (GLM) used for the **CMC Burden Test**, **SKAT**, and **SKAT-O**.

#### A. The Null Model (Baseline)
This model removes the effects of covariates to calculate the **Residuals** ($Y_i - \pi_i$). It is the foundation for all subsequent tests.
$$logit(\pi_i) = \ln\left(\frac{\pi_i}{1-\pi_i}\right) = \alpha + \gamma Sex_i + \sum_{k=1}^{10} \delta_k PC_{ik}^{(BBJ)}$$

#### B. Burden Test Models (Collapsing Strategies)
Burden tests compress the variants in a gene into a single "Meta-Genotype" predictor, denoted here as $M_i$.
$$logit(\pi_i) = \alpha + \mathbf{\beta_{burden} M_i} + \gamma Sex_i + \sum_{k=1}^{10} \delta_k PC_{ik}^{(BBJ)}$$
The definition of $M_i$ depends on the specific Burden method:
*   **CMC (Standard Burden)**: $M_i$ is **Binary**.
    $$M_i = I\left(\sum_{j=1}^{m} G_{ij} > 0\right)$$
    *(Interpretation: 1 if user carries ANY variant, 0 otherwise. $\beta_{burden}$ is the Log-OR of being a carrier.)*
*   **SKAT-O Burden ($\rho=1$)**: $M_i$ is a **Weighted Sum**.
    $$M_i = \sum_{j=1}^{m} w_j G_{ij}$$
    *(Interpretation: A linear sum where rare variants contribute more due to weight $w_j$.)*

#### C. SKAT (Variance Component Model)
Instead of collapsing, SKAT tests the joint distribution of all variant effects $\beta_j$ simultaneously.
$$logit(\pi_i) = \alpha + \sum_{j=1}^{m} \beta_j G_{ij} + \gamma Sex_i + \sum_{k=1}^{10} \delta_k PC_{ik}^{(BBJ)}$$
*   **Constraint**: We do not estimate each $\beta_j$ directly. Instead, we assume $\beta_j$ follows a distribution:
    $$\beta_j \sim N(0, \tau w_j^2)$$
*   **Hypothesis**: We test $H_0: \tau = 0$ (i.e., Is the variance of genetic effects significantly non-zero?).

#### D. SKAT-O (Unified)
SKAT-O does not have a separate regression equation. Instead, it constructs a test statistic $Q_{\rho}$ that linearly combines the squared scores from the **SKAT** approach and the **Weighted Burden** approach:
$$Q_{\rho} = (1-\rho) Q_{SKAT} + \rho Q_{Burden}$$

---

**Variable Definitions (Unified):**

| Symbol | Definition | Context & Notes |
| :--- | :--- | :--- |
| $\pi_i$ | **Probability of CTEPH** | $P(Y_i=1)$, estimated from the Null Model. |
| $\alpha, \gamma, \delta$ | **Covariate Coefficients** | Intercept, Sex effect, and PC effects. |
| $G_{ij}$ | **Raw Genotype** | Count of minor alleles (0, 1, or 2) for variant $j$ in person $i$. |
| $M_i$ | **Collapsed Score** | The single predictor variable used in Burden tests. <br> * **CMC**: Binary (0/1). <br> * **SKAT-O($\rho=1$)**: Weighted Sum ($\sum w G$). |
| $w_j$ | **Variant Weight** | * **CMC**: $w_j = 1$ (Flat). <br> * **SKAT/SKAT-O**: $w_j = Beta(MAF_j, 1, 25)$. |
| $\beta_{burden}$ | **Fixed Effect Size** | The scalar effect size we estimate in Burden tests. |
| $\beta_j$ | **Random Effect Size** | The random effect of variant $j$ in SKAT. |
| $\tau$ | **Variance Component** | Represents the magnitude of genetic variation in the gene. |

---

### 7.2 P-value Interpretation & Post-Processing

The raw P-values from `rvtests` (SKAT-O) already account for the internal multiple testing of different $\rho$ weights. However, they do **not** account for the fact that we tested thousands of genes across the genome. We apply the following two-step procedure to determine significance:

#### Step 1: Quality Control Filtering (Defining $N_{genes}$)
Before checking significance, we filter the gene list to ensure statistical power.
*   **Criterion**: A gene is only "tested" if it contains **at least 2 valid variants** in the dataset (after QC).
*   **Logic**: Genes with 0 or 1 variant have insufficient information for a variance-based test (SKAT) or association test, often yielding meaningless P-values (e.g., P=1.0 or unstable estimates). Including them artificially inflates the correction burden.
*   **Result**: This defines our final number of tests, $N_{genes}$ (e.g., reducing from 20,000 total genes to ~15,000 effective genes).

#### Step 2: Significance Thresholds

We employ two complementary approaches to identify candidate genes:

**A. Stringent Threshold (Bonferroni Correction)**
*   **Goal**: Strict control of Family-Wise Error Rate (probability of getting *any* false positive).
*   **Calculation**:
    $$\alpha_{Bonferroni} = \frac{0.05}{N_{genes}}$$
*   **Example**: If $N_{genes} \approx 15,000$, then $\alpha \approx 3.3 \times 10^{-6}$.
*   **Verdict**: Any gene with $P < \alpha_{Bonferroni}$ is considered **Significantly Associated**.

**B. Discovery Threshold (False Discovery Rate - FDR)**
*   **Goal**: Control the proportion of false positives among the top results (Benjamini-Hochberg method).
*   **Calculation**: We calculate the **q-value** (FDR adjusted P-value) for all $N_{genes}$.
*   **Verdict**: Genes with **FDR < 0.05** (or 0.10) are considered **Top Candidates** for validation. This approach is more sensitive than Bonferroni and helps ensure we don't miss true signals with moderate effect sizes.

