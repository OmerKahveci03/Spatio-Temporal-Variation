# Spatio-Temporal-Variation of Gene Expression across Different Pathologies Through the Lens of Aging
Below are the R scripts available, along with a brief explanation of their functionality, and the order in which they must be run.

**1) auto_normalized.R**

This script processes the raw data into the format that every other script will be working with. 
- Raw data consists of text documents named: gene_tpm_{**tissue name**}.gctannotated.txt
- Raw data is stored in a local directory labeled "data"
- The script produces a normalized version for every tissue in the "output" directory
- The normalized files produced are named: normalized_{**tissue name**}.txt and normalized_{**tissue name**}.xlsx
