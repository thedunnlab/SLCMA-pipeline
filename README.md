# IMPORTANT NOTE #
This SLCMA pipeline repository is no longer maintained, as we have transitioned to the slcma package, which is a streamlined and more accessible version of this initial set of scripts. Although this repository remains open for posterity's sake, we recommend you use the following slcma repository for future analyses: https://github.com/thedunnlab/slcma. 


# SLCMA-pipeline
Pipeline of the structured life course modeling approach analysis

## What is the structured life course modeling approach? 

The structured life course modeling approach (SLCMA) is an analytic method that formally compares a set of prespecified life course hypotheses to identify the one that best fits the data. To read more about the SLCMA: 

> Smith, A. D., Hardy, R., Heron, J., Joinson, C. J., Lawlor, D. A., Macdonald-Wallis, C., & Tilling, K. (2016). A structured approach to hypotheses involving continuous exposures over the life course. International Journal of Epidemiology, 45(4), 1271–1279. https://doi.org/10.1093/ije/dyw164

> Smith, A. D. A. C., Heron, J., Mishra, G., Gilthorpe, M. S., Ben-Shlomo, Y., & Tilling, K. (2015). Model Selection of the Effect of Binary Exposures over the Life Course. Epidemiology (Cambridge, Mass.), 26(5), 719–726. https://doi.org/10.1097/EDE.0000000000000348

> Mishra, G., Nitsch, D., Black, S., De Stavola, B., Kuh, D., & Hardy, R. (2009). A structured approach to modelling the effects of binary exposure variables over the life course. International Journal of Epidemiology, 38(2), 528–537. https://doi.org/10.1093/ije/dyn229

## Are there applied examples of the SLCMA? 

Yes, other researchers as well as our lab have used the SLCMA to study the effects of repeated exposures on various health outcomes, such as psychopathology, emotion recognition, or genome-wide DNA methylation. Here are some examples of our work:

> Lussier, A. A., Zhu Y., Smith B. J., Simpkin A. J., Smith A. D. A. C., Suderman M. J., et al. Updates to data versions and analytic methods influence the reproducibility of results from epigenome-wide association studies. bioRxiv. https://doi.org/10.1101/2021.04.23.441014

> Dunn, E. C., Soare, T. W., Zhu, Y., Simpkin, A. J., Suderman, M. J., Klengel, T., … Relton, C. L. (2019). Sensitive Periods for the Effect of Childhood Adversity on DNA Methylation: Results From a Prospective, Longitudinal Study. Biological Psychiatry, 85(10), 838–849. https://doi.org/10.1016/j.biopsych.2018.12.023

> Dunn, E. C., Crawford, K. M., Soare, T. W., Button, K. S., Raffeld, M. R., Smith, A. D. A. C., … Munafò, M. R. (2018). Exposure to childhood adversity and deficits in emotion recognition: Results from a large, population-based sample. Journal of Child Psychology and Psychiatry, 59(8), 845–854. https://doi.org/10.1111/jcpp.12881

> Dunn, E. C., Soare, T. W., Raffeld, M. R., Busso, D. S., Crawford, K. M., Davis, K. A., … Susser, E. S. (2018). What life course theoretical models best explain the relationship between exposure to childhood adversity and psychopathology symptoms: Recency, accumulation, or sensitive periods? Psychological Medicine, 1–11. https://doi.org/10.1017/S0033291718000181

## What are the scripts in this repository? 

- `SLCMA_completecase_analysis.R` 

A script to run the SLCMA / least angle regression (LARS) using complete data. This function does not accomodate imputed objects. The post-selection inference method, selective inference (`sI`), is used. 

- `SLCMA_toy_data.Rdata` 

A simulated toy data to use with the `SLCMA_completecase_analysis.R` script, if desired.

- `SLCMA_imputed_analysis.R`

A script to run the SLCMA / LARS using multiply imputed data. The post-selection inference method, selective inference (`sI`), is used. This script requires the `LARS-function-imputed_data.R` script below.

- `LARS-function-imputed_data.R`

A function to perform LARS that is sourced within the `SLCMA_imputed_analysis.R` script.

- `Imputation_example_code.R`

A script to use for performing multiple imputation using the `mice` R package.

- `imputation-function-v1.1.R`

A script that includes the do.impute() function for the imputation.

## Who should I contact if I have questions? 

Please email us at dunnlab@mgh.harvard.edu. 
