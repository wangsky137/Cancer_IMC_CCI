# CCI Score Calculator for Cell-Cell Interactions

## Overview

This repository provides the `CCI_score_Cal` function, an R tool designed to calculate a Cell-Cell Interaction (CCI) score matrix. The function evaluates the interaction strength between different cell types in a spatial context, which is vital for understanding cellular interactions within the tumor microenvironment.

The function uses k-nearest neighbors (KNN) and permutation testing to compute a z-score matrix based on observed and randomized cell type distributions. This approach is particularly valuable for identifying significant interactions between cancer and non-cancer cell types.

## Features

- **K-Nearest Neighbors (KNN):** Calculates distances between cells and identifies the nearest neighbors.
- **Permutation Testing:** Generates a null distribution by shuffling non-cancer cell types, allowing for robust z-score calculations.
- **Z-Score Matrix:** Quantifies interaction strength between each pair of cell types.
- **Parallel Processing:** Supports parallel processing for efficiency when working with large datasets.

## Installation

Ensure that R is installed on your system. Clone the repository and load the required data and functions as shown below:

```bash
git clone https://github.com/wangsky137/Cancer_IMC_CCI.git
```

## Usage

### Example Code

Here’s an example of how to use the `CCI_score_Cal` function with your data:

```r
# Clear environment
rm(list=ls())

# Load necessary library for parallel processing
library(parallel)

# Load the example data and CCI calculation script
load("./example_data.rda")
source("./CCI_Z_cal.R")

# Define cancer cell types
cancer_cells <- c(
  "panCKmed", "Vimentin+EMT", "CD56+NE", "PD-L1+GZMB+", "Helios+", 
  "CD15+", "pH2AX+DSB", "CKloGATA3+", "Basal", "CK8/18med", "TCF1+", 
  "MHCI&IIhi", "CA9+Hypoxia", "PD-L1+IDO+", "CKhiGATA3+", "AR+LAR", 
  "Apoptosis"
)

# Define immune cell types of interest
cell_type_of_interest <- c("CD8+T", "Treg", "CD8+GZMB+T", "CD8+PD1+T_Ex")

# Calculate z-score for the first image
z_score <- CCI_score_Cal(myList[[1]], cell_type_of_interest, cancer_cells)

# Apply CCI calculation across all image in the list in parallel
CCI<-mclapply(myList,CCI_score_Cal,CCI_cell=cell_type_of_interest,cancer_cell_types=cancer_cells)
df <- data.frame(matrix(unlist(CCI), nrow=length(CCI), byrow=T))
colnames(df)<- names(CCI[[1]])
rownames(df)<- names(CCI)
print(df)
```

### Parameters

- **`CCI_data`**: A data frame or matrix containing spatial coordinates of cells. The first two columns should be the x and y coordinates.
- **`CCI_cell`**: A vector of cell type labels that will be calcualted for CCI score to the cells in `CCI_data`.
- **`cancer_cell_types`**: A character vector specifying the names of cancer cell types in `CCI_cell`.
- **`Nknn`**: Integer specifying the number of nearest neighbors to consider (default is 10).
- **`nPermu`**: Integer specifying the number of permutations for the null distribution (default is 100).

## Output

The function returns a vector of z-scores, quantifying the interaction strength between each pair of cell types, with higher scores indicating stronger interactions.

## Example Output

Below is an example output from the `CCI_score_Cal` function, showing the z-score matrix that quantifies interaction strengths between cell types across different treatment conditions:


|                      |   CD8+T→CD8+T |   CD8+T→Treg |   CD8+T→CD8+GZMB+T |   CD8+T→CD8+PD1+T_Ex |   Treg→CD8+T |   Treg→Treg |   Treg→CD8+GZMB+T |   Treg→CD8+PD1+T_Ex |   CD8+GZMB+T→CD8+T |   CD8+GZMB+T→Treg |   CD8+GZMB+T→CD8+GZMB+T |   CD8+GZMB+T→CD8+PD1+T_Ex |   CD8+PD1+T_Ex→CD8+T |   CD8+PD1+T_Ex→Treg |   CD8+PD1+T_Ex→CD8+GZMB+T |   CD8+PD1+T_Ex→CD8+PD1+T_Ex |
|:---------------------|--------------:|-------------:|-------------------:|---------------------:|-------------:|------------:|------------------:|--------------------:|-------------------:|------------------:|------------------------:|--------------------------:|---------------------:|--------------------:|--------------------------:|----------------------------:|
| Baseline_C_1         |    -0.484659  |    0.0153345 |           0.638325 |             0.692064 |   -0.0555156 |  -0.0961158 |          0.125763 |          -0.829562  |          0.0941184 |        -0.245122  |                8.79459  |                  2.52505  |             0.771239 |           -0.178865 |                  2.89561  |                   -0.608101 |
| Baseline_C_2         |     1.97658   |    3.44946   |           0.32742  |             1.59357  |    2.19103   |   4.654     |         -0.49057  |           5.30062   |         -0.0861618 |        -0.23991   |               11.1959   |                  3.08653  |             1.02216  |            6.22291  |                  3.32259  |                    2.44522  |
| On-treatment_C_1     |     2.09853   |   -1.92603   |          -3.22885  |            -2.00298  |   -1.78545   |   1.56381   |         -2.56276  |           1.75427   |         -3.30442   |        -2.05823   |               15.111    |                  0.242271 |            -3.4004   |            0.744992 |                  0.370265 |                    7.95948  |
| On-treatment_C_2     |     2.4682    |   -0.760529  |           0.481257 |            -1.46504  |   -1.88556   |   4.48393   |         -1.3792   |          -0.0171329 |          0.292171  |        -1.78422   |               -0.151582 |                  1.45215  |            -1.92987  |           -0.266988 |                  2.20489  |                    5.75781  |
| Post-treatment_C&I_1 |     2.63283   |    3.49137   |           0.768165 |            -0.60093  |    3.40201   |  -0.875865  |          3.77011  |           5.11375   |          2.04723   |         3.93886   |                5.64968  |                  5.38102  |            -0.586812 |            5.28128  |                  6.52779  |                   14.382    |
| Post-treatment_C&I_2 |    -0.35463   |   -0.178652  |           1.33948  |             0.368687 |   -0.0417438 |   0.0495203 |          1.91757  |           3.92628   |          1.0428    |         1.19258   |                4.12903  |                  0.213503 |             0.925909 |            3.03776  |                  0.649114 |                   12.6199   |
| Post-treatment_C_1   |     0.0479624 |   -0.136458  |           0.421978 |             3.66491  |    0.341306  |   2.03716   |          0.428992 |          -0.69441   |         -0.489445  |        -0.370874  |                4.32161  |                  0.987199 |             0.668755 |           -0.658257 |                  0.435786 |                    1.44814  |
| Post-treatment_C_2   |    -0.728428  |    2.25241   |           0.5622   |             2.6457   |    1.56926   |   1.55658   |          1.57508  |           5.30864   |          1.93832   |         2.25809   |                4.65948  |                  0.193562 |             2.59155  |            5.37973  |                  1.25857  |                    1.99629  |
| Baseline_C&I_1       |     0.941425  |    0.805693  |          -0.518488 |             4.3738   |    0.732927  |   0.354487  |         -0.48129  |          -0.159274  |         -0.474246  |         0.0610132 |                2.18121  |                  2.30367  |             4.25094  |            0.245567 |                  1.67162  |                    1.26936  |
| Baseline_C&I_2       |    -0.251379  |    0.497494  |          -0.904204 |             0.673091 |    0.445368  |   0.508933  |         -0.827835 |          -1.06834   |         -0.905895  |        -0.797733  |                0.64121  |                 -2.69596  |             0.66932  |           -0.833162 |                 -1.80571  |                    1.28177  |
| On-treatment_C&I_1   |    -0.217773  |    0.994987  |           0.784718 |             2.99084  |    1.0345    |   6.5497    |          0.854086 |           1.60011   |          0.994987  |         0.22657   |                0.764966 |                  0.857542 |             1.75568  |            1.66304  |                  1.95274  |                    2.37914  |
| On-treatment_C&I_2   |     5.36496   |   -0.669847  |          -1.04883  |             3.0142   |    0.287884  |  -0.767019  |          0.100223 |           0.288679  |         -1.11991   |         0.194067  |                6.97218  |                 -0.504476 |             2.50183  |           -0.328637 |                 -0.471211 |                   -0.467361 |


Higher z-scores indicate stronger interactions, particularly between cancer cells and other cell types.

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests to enhance this project.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or further information, please contact [Xiang.Wang@bcm.edu].
