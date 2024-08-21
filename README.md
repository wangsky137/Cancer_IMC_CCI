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

# Identify unique cell types across all datasets
cel.typ = NULL
for(k in 1:length(myList)) {
  cat("
", k)
  data = myList[[k]]
  tmp = data$cel.typ
  cel.typ = unique(c(cel.typ, tmp))
}

# Define cancer cell types
cancer_cells <- c(
  "panCKmed", "Vimentin+EMT", "CD56+NE", "PD-L1+GZMB+", "Helios+", 
  "CD15+", "pH2AX+DSB", "CKloGATA3+", "Basal", "CK8/18med", "TCF1+", 
  "MHCI&IIhi", "CA9+Hypoxia", "PD-L1+IDO+", "CKhiGATA3+", "AR+LAR", 
  "Apoptosis"
)

# Define immune cell types
cc <- c("CD8+T", "Treg", "CD8+GZMB+T", "CD8+PD1+T_Ex")

# Calculate z-score for the first dataset
z_score <- CCI_score_Cal(myList[[1]], cc, cancer_cells)

# Apply CCI calculation across all datasets in parallel
xx <- mclapply(myList, CCI_score_Cal, CCI_cell = cc, cancer_cell_types = cancer_cells)
```

### Parameters

- ** `CCI_data`**: A data frame or matrix containing spatial coordinates of cells. The first two columns should be the x and y coordinates.
- **`CCI_cell`**: A vector of cell type labels that will be calcualted for CCI score to the cells in `CCI_data`.
- **`cancer_cell_types`**: A character vector specifying the names of cancer cell types in `CCI_cell`.
- **`Nknn`**: Integer specifying the number of nearest neighbors to consider (default is 10).
- **`nPermu`**: Integer specifying the number of permutations for the null distribution (default is 100).

## Output

The function returns a matrix of z-scores, quantifying the interaction strength between each pair of cell types, with higher scores indicating stronger interactions.

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests to enhance this project.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or further information, please contact [Xiang.Wang@bcm.edu].
