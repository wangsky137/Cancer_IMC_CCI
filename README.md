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

# Calculate z-score for the first dataset
z_score <- CCI_score_Cal(myList[[1]], cell_type_of_interest, cancer_cells)

# Apply CCI calculation across all datasets in parallel
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

```markdown
|    | CD8+T→CD8+T "CD8+T→Treg" "CD8+T→CD8+GZMB+T" "CD8+T→CD8+PD1+T_Ex" "Treg→CD8+T" "Treg→Treg" "Treg→CD8+GZMB+T" "Treg→CD8+PD1+T_Ex" "CD8+GZMB+T→CD8+T" "CD8+GZMB+T→Treg" "CD8+GZMB+T→CD8+GZMB+T" "CD8+GZMB+T→CD8+PD1+T_Ex" "CD8+PD1+T_Ex→CD8+T" "CD8+PD1+T_Ex→Treg" "CD8+PD1+T_Ex→CD8+GZMB+T" "CD8+PD1+T_Ex→CD8+PD1+T_Ex"   |
|---:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|  0 | Baseline_C_1 -0.484659417160015 0.0153345423219643 0.638325008092716 0.69206386126028 -0.055515570375648 -0.0961157652790334 0.125763428692838 -0.829562366176523 0.0941184026099537 -0.245122172442765 8.79458667449179 2.52504743130705 0.771238785529603 -0.178865364650215 2.89560780398735 -0.608101037747188      |
|  1 | Baseline_C_2 1.97658035689858 3.44946165657416 0.327420391786764 1.5935727331287 2.1910256385428 4.65399645026174 -0.490569984757674 5.30061768810627 -0.0861618459315523 -0.239910138138846 11.1959248506546 3.08652514945316 1.02216478461427 6.22290745101316 3.32259482579128 2.4452175845309                       |
|  2 | On-treatment_C_1 2.09852700549069 -1.92602613110541 -3.2288470976957 -2.00298066319461 -1.78545100029073 1.56381287337567 -2.56276303435132 1.7542729455055 -3.30442326507276 -2.0582313590896 15.1109722573049 0.24227055485477 -3.40040441020214 0.74499161156748 0.370265191733024 7.9594848351631                   |
|  3 | On-treatment_C_2 2.46820425966465 -0.760528897135702 0.481256654078399 -1.46504292431472 -1.88555956678701 4.48392940816328 -1.37919636757574 -0.017132855719005 0.292170863228588 -1.78422142303008 -0.151582184455513 1.45215105057029 -1.92987190587103 -0.266987988898902 2.20489354083254 5.75780834920729         |
|  4 | Post-treatment_C&I_1 2.63282932129703 3.49137398565268 0.768164728481173 -0.60092980065077 3.4020068727481 -0.875865130501177 3.77011103166458 5.11375476581328 2.047228705946 3.93886311721569 5.64967785988366 5.38101927921973 -0.586811576306191 5.28127526006957 6.52779301642847 14.3820268532547                 |
|  5 | Post-treatment_C&I_2 -0.354630282636209 -0.178651981964454 1.33948128976055 0.36868740284723 -0.0417437555128694 0.0495203044009818 1.91756856772699 3.92627553277636 1.04280085745965 1.19258195804425 4.12903499002992 0.213502727508516 0.925908912298289 3.03776107605555 0.649113960557953 12.6199164205555        |
|  6 | Post-treatment_C_1 0.0479624290958715 -0.136457647844203 0.421977515511518 3.66491307952701 0.341305697753096 2.03716435216709 0.428991956476223 -0.694409579542567 -0.489444826093362 -0.370873522661461 4.32161089124202 0.987198594690213 0.668754913307888 -0.658256870627521 0.435785720392048 1.4481382518224     |
|  7 | Post-treatment_C_2 -0.728427590386427 2.25240554233259 0.562199787960751 2.64570361795674 1.56926295020833 1.55657583643121 1.57507809215759 5.30864247161377 1.93831672910114 2.25809264068778 4.65947859108084 0.193561854260279 2.59155440235862 5.37973172020827 1.25857061780418 1.99628609638303                  |
|  8 | Baseline_C&I_1 0.941424540153286 0.805692535893954 -0.518488358612215 4.37379894228118 0.732926958351404 0.354487400481549 -0.48128957688279 -0.159274481283259 -0.47424607022541 0.061013229116527 2.18120543343747 2.30366635665496 4.25094339622641 0.245566789384248 1.67162301347356 1.26935795546972              |
|  9 | Baseline_C&I_2 -0.251379174478918 0.49749371855331 -0.904203666030459 0.673090519512515 0.445367967968345 0.508932529235945 -0.827834606294395 -1.06834392691747 -0.905894616004588 -0.797732779698271 0.641209517256826 -2.69595985088147 0.669319674606718 -0.833161618794051 -1.80570740394398 1.28177266957545      |
| 10 | On-treatment_C&I_1 -0.217773161111738 0.99498743710662 0.784718288728972 2.99083908537445 1.03449503503053 6.54970462733963 0.854086246337541 1.60010775499221 0.99498743710662 0.22656960463162 0.764966381151798 0.857542408975018 1.75567599233038 1.66304156639775 1.95273817311628 2.37914380582238                |
| 11 | On-treatment_C&I_2 5.36495731961771 -0.669846778232843 -1.04882995969421 3.01420143570274 0.287883637846879 -0.767019319321146 0.100222687667208 0.288679245283019 -1.11991120410831 0.194067201863847 6.97218302376739 -0.504475760392729 2.50183088585244 -0.328637086971624 -0.471211313528422 -0.467360998227699    |
```

Higher z-scores indicate stronger interactions, particularly between cancer cells and other cell types.

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests to enhance this project.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or further information, please contact [Xiang.Wang@bcm.edu].
