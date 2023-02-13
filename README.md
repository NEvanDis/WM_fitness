# Winter moth individual fitness and population growth
This folder contains all the scripts needed to reproduce the analysis of experimental and field winter moth data belonging to manuscript _Phenological mismatch affects individual fitness and population growth in the winter moth_.

**NB: The raw data, including both experimental data and field data, and the version of record analysis scripts can be found on Dryad [add DOI].**

&nbsp;

## Authors
Natalie E. van Dis, ORCID ID: 0000-0002-9934-6751

&nbsp;

## Analysis and visualization of Experimental data
### Script: ```2_scripts/1_CatFoodExp2021_analysis.R ```
R script to reproduce the analysis and visualization (incl. manuscript figures) of the 2021 winter moth caterpillar feeding experiment: (1) What are the fitness consequences of day to day timing (a)synchrony with budburst? and (2) Can food quality affect the timing of life stages?

See ```_src/env_CatFoodExp2021_analysis.txt``` for used R package versions.

&nbsp;

## Analysis and visualization of long-term Field data
### Script: ```2_scripts/2_prep_FieldData.R ```
R script to get trapping effort descriptives and to prep all the field data for analysis.

### Script: ```2_scripts/3_plot_popnum.R ```
R script to produce manuscript figure visualizing winter moth population dynamics and population phenological mismatch over time (1993-2021) at four locations in the Netherlands.

### Script: ```2_scripts/4_popdyn_analysis.R ```
R script to reproduce the analysis of winter moth population dynamics: how much variation in population growth can be explained by timing mismatch?

See ```_src/env_PopDyn_analysis.txt``` for used R package versions.
