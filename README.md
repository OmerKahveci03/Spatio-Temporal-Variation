# Spatio-Temporal-Variation of Gene Expression across Different Pathologies Through the Lens of Aging
Below are the R scripts available, along with a step by step explanation of how to compile them.

## Prerequisites
- Download the raw tissue data from the gtex database, and leave it in a local directory named "data". File path: \data
- Download the supplementary file "supplementaryfile.xlsx" located in "data" directory as well. It is needed for all scripts that distinguish between AD genes and non-AD genes
- Download "9606.protein.links.v12.0.txt" and store it in the "data" directory

## Scripts
### auto_normalized.R
**This is the first script you should run**. Will create "normalized_{**tissue name**}.txt" and "normalized_{**tissue name**}.xlsx" for every tissue, in a directory named "output". File path: \output

### auto_threshhold.R
Creates a "plot_{**tissue name**}.png" file for every tissue, in a directory named "plots". File path: \plots
- Uses normalized files produced by ***auto_normalized.R***.
- The plots are # of switching genes x ***t*** value

### auto_threshhold_ad.R
Creates a "plot_{**tissue name**}.png" file for every tissue, in a directory named "AD" in the "plots" directory. File path: \plots\AD
- Needs normalized files produced by ***auto_normalized.R***
- Needs "supplementaryfile.xlsx" to be located in the "data" directory
- The plots are # of switching genes x ***t*** value, but only for the top 500 AD genes

### auto_bargraph.R
Uses normalized files produced by ***auto_normalized.R***. Creates bar graphs "ad_{**tissue name**}_bargraph_t0.5.png" and "all_{**tissue name**}_bargrapht0.5.png" in a directory named "graphs". File path: \graphs
- Needs "supplementaryfile.xlsx" to be located in the "data" directory
- t value can be configured on line 12 of the script to be any number between 0 and 1. The name of the files produced will change accordingly.

### auto_sets.R
Creates "{**tissue name**}_sets.txt" for every tissue in a directory "ad_sets" and "all_sets".
- Needs normalized files produced by ***auto_normalized.R***
- Needs "supplementaryfile.xlsx" to be located in the "data" directory
- Each line of the sets file represents an age group. Line 1 contains switching genes at ages 30-39. Line 2 contains switching genes at ages 40-49, and so on.

### Heatmap.R
Creates "heatmap_ad_age.png", "heatmap_ad_life.png", "heatmap_all_age.png", "heatmap_all_life.png" files in the "heatmaps" directory. File path is \heatmaps
- Needs sets files produced by ***auto_sets.R***

### connect.R
Creates "heatmap_connectivity.png" and "heatmap_connectivity.xlsx" in the "connectivity" directory. File path is \connectivity
- Needs "9606.protein.links.v12.0.txt" to be located in the "data" directory
- Needs sets files produced by ***auto_sets.R***

### Ad_Sets and All_Sets
- Contains a text file for every tissue, containing the switching genes
- ad_sets contains switching alzhiemers genes, all_sets contains all switching genes.
- Each line corresponds to an age group. Some lines may be empty, but there will only be up to five.
  - Line 1: Ages 30-39
  - Line 2: Ages 40-49
  - Line 3: Ages 50-59
  - Line 4: Ages 60-69
  - Line 5: Ages 70-79
