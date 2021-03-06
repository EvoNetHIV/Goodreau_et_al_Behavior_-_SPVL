## Sexual network structure, HIV prevalence, and the evolution of set point viral load

This repository contains the code to reproduce the analysis of set point viral load (SPVL) evolution with varying network structures.  In this model, we investigate the changes to SPVL with four behavioral phenomena: coital frequency, relational duration, relational concurrency, and coital cessation.     

This model is written in the R programming language.  It requires [EpiModel](http://www.epimodel.org), the epidemic modeling software, to be installed first, as well as [EvoNetHIV](https://github.com/EvoNetHIV).      

Within R:  
`install.packages("EpiModel")`     
`if(!require(devtools)) { install.packages("devtools") }`     
`library(devtools)`         
`devtools::install_github("EvoNetHIV/Goodreau_et_al_Behavior_-_SPVL")`     

### Citation
> Goodreau SM, Stansfield SE, Murphy JT, Peebles KC, Gottlieb GS, Abernethy NF, Herbeck JT, Mittler JE. (submitted). Sexual network structure, HIV prevalence, and the evolution of set point viral load.     

### Abstract
HIV viral load (VL) predicts both transmission potential and rate of disease progression. For reasons that are still not fully understood, the set point viral load (SPVL) established after acute infection varies across individuals and populations. Previous studies have suggested that population mean SPVL (MSPVL) has evolved near an optimum that reflects a trade-off between transmissibility and host survival. Sexual networks affect rates of potential exposure during different within-host phases of infection marked by different transmission probabilities, and thus affect the number and timing of transmission events. Given this, we hypothesize that the sexual network can alter the fitness landscape for SPVL.     
      
To quantitatively test this hypothesis, we developed a dynamic, stochastic, network model of HIV transmission and evolution to assess the impact of key sexual network phenomena (coital frequency, relational concurrency, relational duration, and coital cessation in AIDS) on MSPVL evolution. Experiments were repeated in sensitivity runs that made different assumptions about transmissibility during acute infection, SPVL heritability, and the functional form of the relationship between VL and transmissibility.     
       
For our main transmission model, different scenarios yielded MSPVLs ranging from 4.4 to 4.75 log10 copies/mL, covering much of the observed empirical range. MSPVL evolved to be significantly higher in populations with high concurrency and shorter relational durations. For example, with mean relational duration of 750 days, MPSVL for 20% concurrency was 4.59 but for 0% concurrency was 4.41; fixing concurrency at 10%, MSPVL was 4.46 and 4.79 for mean relational durations of 1000 and 50 days, respectively. In regression analyses on these and other predictors, main effects were significant (p<0.05), as were interaction effects, indicating that effects are interdependent. We also noted a strong correlation between two key emergent properties measured at the end of the simulations—MSPVL and HIV prevalence—most clearly for phenomena that affect transmission networks early in infection. Interestingly, we observed significantly lower prevalence in runs in which SPVL heritability was zero, indicating the potential for viral evolution to exacerbate disease burden over time. Future efforts to understand empirical variation in MSPVL should consider local HIV burden and basic sexual behavioral and network structure.    
