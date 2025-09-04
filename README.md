# Leakage-in-Majorana-Qubits
This repository contains the MATLAB code for the numerical simulations in the paper "M. C. Goffage, A. Alase, M. C. Cassidy, S. N. Coppersmith, Leakage at Zero Temperature from Changes in Chemical Potential in Majorana Qubits."

Please note that this code uses the QPP_Library which is still work in progress, there are placeholders included for extended functionality that is yet to be added. All functionality used to generate results for our papers have been tested.

---

## Environment and Dependencies
Requires MATLAB 2024a or more recent.  

---

## How to Run the Code 
1. Download this entire repository as a ZIP. 
2. Unzip the contents into a local directory.  
3. There is **`.m`** run file corresponding to each figure in the main text and supplement. Open your chosen run file.  
4. In MATLAB, select the **Editor** tab and then select **Run**.  
5. View the generated figures in the "results" or "Supplement Plots/results" directory.

---

## Run Files 
There is **`.m`** run file corresponding to each figure in the main text and supplement (within the Supplement Plots directory). You may run the file corresponding to your figure of interest. Note that **`run_fig_S3.m`** requires **`run_fig_2_main.m`** to be run first. All other run files have no dependencies.  

Alternatively you may generate all figures using **`run_all_code.m`**, which has an expected run time of ~24 hours on a typical laptop.  

---

## Outputs
All figures (as .png and .fig files) and data (as .m files) generated for the main text figures are saved in the "results" directory, and the results for the Supplement figures are saved in the "Supplement Plots/results" directory. 

---
## arXiv Post
This work is currently posted on the arXiv at: https://arxiv.org/abs/2504.17485

---
## Troubleshooting

If you have any issues with this code, please contact Marcus Goffage at m.goffage@unsw.edu.au. 
