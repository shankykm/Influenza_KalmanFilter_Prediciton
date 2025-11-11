# Flu Data Simulation and Forecasting with Kalman Filter and SARIMA
This repository contains R code for simulating flu data and forecasting using **Kalman Filter (KF)** and **SARIMA** models. The code is designed for researchers, data scientists, and epidemiologists interested in time series forecasting for public health applications.

![GitHub](https://github.com/shankykm/Influenza_KalmanFilter_Prediciton.git)

---

## **Summary**
This project:
- Simulates synthetic weekly flu data (positive cases and total lab test volumes) for 4 years (2016–2019).
- Implements **Kalman Filter (KF)** models (ensemble and bivariate) to forecast flu positive cases (P) and lab test volumes (T).
- Uses **SARIMA** for comparison.
- Extracts **historical smoothed estimates** using a KF AR(1) model to serve as covariates for the Positive Case Model.
- Evaluates model performance using **Mean Absolute Deviation (MAD)** and **coverage probability**.
- Visualizes forecasts and performance metrics.

---

## **Breakdown of the R Code**

---

### **1. Setup and Data Simulation**
- **Libraries**: Loads required R packages (`tidyverse`, `lubridate`, `forecast`, `KFAS`, `rjags`, etc.).
- **Data Simulation**:
  - Simulates weekly flu data (`P_simulated` for positive cases, `T_simulated` for total volume) using a bivariate normal distribution.
  - Uses synthetic means and standard deviations, with a correlation of 0.93 between the two variables.
  - Ensures non-negativity of simulated values.

---

### **2. KF AR(1) Model to Extract Average Historical Smoothed Estimates**
- **Purpose**: Extracts smoothed estimates of historical total Lab test (T) to use as covariates for the Positive Case Model.
- **Process**:
  - Uses a **bivariate state-space model** (`SSModel`) to fit an AR(1) process to the simulated data.
  - Computes smoothed estimates for both total volume and detected cases using the Kalman Filter.
  - Aggregates these estimates by year, excluding the current year to avoid overfitting.
  - Visualizes the smoothed estimates for each year.
- **Output**: Adds smoothed estimates as new columns (`smoothed_total_vol.fit`, `smoothed_Detected.fit`) to the dataset.

---

### **3. Kalman Filter (KF) Models**
#### **A. Ensemble KF for Total Cases**
- **Model Definition** (`kf_ensemble_extractInformativePrior`):
  - Combines two KF models (one for total volume, one for detected cases) using a weight `w1`.
  - Extracts informative priors from historical data, excluding the current year to avoid overfitting.
- **Forecasting**:
  - Performs rolling 4-week-ahead forecasts for 2019 using informative priors.
  - Visualizes forecasts and performance metrics.

#### **B. Bivariate KF for Positive Cases**
- **Model Definition** (`bivar_kfar1gaussJags8_2`):
  - Models both total volume and detected cases, in difference scale and uses smoothed total volume (from step 2) as a covariate .
  - Extracts informative priors from historical data.
- **Forecasting**:
  - Forecasts positive cases and total volume, undoing the differencing to return to the original scale.
  - Evaluates performance using MAD and coverage probability.

---

### **4. SARIMA Model**
- **Model Fitting**:
  - Uses `auto.arima` to fit a seasonal ARIMA model to the data.
  - Generates rolling 4-week-ahead forecasts for 2019.
- **Evaluation**:
  - Computes MAD and coverage probability.
  - Visualizes forecasts and performance metrics.

---

### **5. Model Comparison**
- **Combines Results**:
  - Merges results from KF and SARIMA models.
  - Compares performance using MAD and coverage probability.
- **Visualization**:
  - Uses boxplots and faceted plots to compare absolute residuals and coverage probability across models.

---

## **⚠️ Important Cautions**

### **Step 3: Kalman Filter (KF) Models**
- **MCMC Sampling Warning**:
  - The **Kalman Filter models (Step 3)** use **Markov Chain Monte Carlo (MCMC) sampling** to generate posterior distributions for model parameters.
  - **This process is computationally intensive** and can take **several hours or more** to complete, depending on your system's resources.
  - **Recommendations**:
    - Test the code on a **small subset of data** (e.g., 1 year instead of 4) to estimate runtime.
    - Consider **parallelizing loops** or **reducing the number of MCMC iterations** for faster execution.
    - Use a **high-performance computing environment** if available.

---

## **How to Run the Code**
1. **Prerequisites**:
   - Install R (version 4.0 or later).
   - Install required R packages:
     ```r
     install.packages(c("tidyverse", "lubridate", "forecast", "KFAS", "rjags", "ggplot2", "dplyr", "knitr"))
     ```
   - Ensure `JAGS` is installed for Bayesian modeling.

2. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/yourrepo.git
   cd yourrepo
