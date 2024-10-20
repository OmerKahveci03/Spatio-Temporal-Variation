# Spatio-Temporal-Variation of Gene Expression across Different Pathologies Through the Lens of Aging
Below are the R scripts available, along with a step by step explanation of how to compile them.

## Prerequisites
- Download the raw tissue data from the gtex database, and leave it in a local directory named "data". File path: \data
- Download the supplementary file located in "data" as well. It is needed for **auto_threshhold_ad.R**

## Scripts
### auto_normalized.R
**This is the first script you should run**. Will create "normalized_{**tissue name**}.txt" and "normalized_{**tissue name**}.xlsx" for every tissue, in a directory named "output". File path: \output

### auto_threshhold.R
Uses normalized files produced by ***auto_normalized.R***. Creates a "plot_{**tissue name**}.png" file for every tissue, in a directory named "plots". File path: \plots
- The plots are # of switching genes x ***t*** value

### auto_threshhold_ad.R
Uses normalized files produced by ***auto_normalized.R***. Creates a "plot_{**tissue name**}.png" file for every tissue, in a directory named "AD" in the "plots" directory. File path: \plots\AD.
- The plots are # of switching genes x ***t*** value, but only for the top 500 AD genes
- Needs "supplementaryfile.xcsl" to be located in the "data" directory
