recode-transform
================

This repository contains R code to recode and transform variables in a R dataframe.

Transformations that are performed include:

- Blom rank based inverse normal transformation   
- Ordinal transformation using equal-quantile algorithm   
- Ordinal transformation using equal-interval algorithm   
  
The formula for the Blom transformation is:   

- y <- qnorm((r-c)/(N-2c+1))   
    - where qnorm is the standard normal quantile function   
    - r is the rank of your variable   
    - c = 3/8, the Blom constant   
    - N is the sample size   
    
Equal-quantile ordinal transformation divides the original variable into k categories based on k quantiles; the transformed variable will have a uniform distribution.

Equal-interval ordinal transformation divides the original variable into k categories such that the within category range, in original variable units, is equated across categories. The distribution of the transformed variable will approximate that of the original variable. Equal-interval will start with the number of specified categories (ncat) but will iteratively decrease this number until the minimum number of observed responses (nobs) is achieved in each category.

Equal-quantile is appropriate for continuous measures where the sample size divided by the
number of categories is greater than the minimum number of observations desired for each response category. Equal-quantile does not work well with ordinal scales that have highly skewed
distributions; equal-interval is recommended. It is important to check the distributions
of the recoded variables to assure that the recoding matches the purpose and goals for recoding.

Functions include:

recodeBlom(`df, varlist_orig, varlist_tr`)  

  - description: This function performs Blom transformation of selected variables   
  - input parameters   
    - df - label for the data frame that contains the variables to be recoded. Can be dataframe or dataframe name (in quotes).   
    - varlist_orig - list of labels for the original variables to be recoded   
    - varlist_tr - list of labels for the recoded variables   
  - output: data frame that consists of the input data frame with recoded variables   

recodeOrdinal(`df, varlist_orig, varlist_tr, type="interval", ncat=10, nobs=10`)    

  - description: This function performs equal-quantile or equal-interval recoding into ordinal response scale   
  - input parameters:   
    - df - label for the data frame that contains the variables to be recoded. Can be dataframe or dataframe name (in quotes).   
    - varlist_orig - list of labels for the original variables to be recoded   
    - varlist_tr - list of labels for the recoded variables   
    - type - "interval" for equal-interval recoding, "quantile" for equal-quantile   
    - ncat - maximum number of categories for recoded variables   
    - nobs - miminum number of observed responses in each category of the transformed variable   
  - output: data frame that consists of the input data frame with recoded variables   

recodeLookup(`df, varlist_orig, varlist_tr, type="continuous", lookup=NULL`)  

  - description: This funciton generates a reference table from recoded variables within a dataframe and uses this as a lookup table to recode vales that have not yet been recoded. The lookup values can from the "internal" input data frame or from an "external" table that has contains the original and recoded variables. "Internal" lookup is applicable where a subset of records are used for the original recode and this subset is then used as the reference to recode the other records. This function applies to continuous transformed variables (including but not limited to Blom) and to ordinal transformations.
  - input parameters:  
    -  df - label for the data frame that contains the variables to be recoded as well as recoded variables with recoded values for a subset of records. Can be dataframe or dataframe name (in quotes). 
    - varlist_orig - list of labels for the original variables to be recoded   
    - varlist_tr - list of labels for the recoded variables   
    - type - "continuous" for continuous (numeric) transformation, "ordinal" for ordinal   
    - lookup - NULL (no parameter input, the default) if the input table (df) will be used to lookup recode values; label for the data frame that contains the variables to be recoded if an external table is used. The external table has to be opened as a dataframe, and must contain all of the columns in varlist_orig and varlist_tr. Can be dataframe or dataframe name (in quotes).   
  - output: data frame that consists of the input dataframe with recoded values for all records for recoded variables  
  
  createLookupTable(`df, varlist_orig, varlist_tr, lookup=NULL`)  

  - description: This function generates a reference table from recoded variables within a dataframe that can be used as a lookup table to recode vales that have not yet been recoded. This is applicable for recoded variables that take on integer values.
  - input parameters:  
    -  df - label for the data frame that contains the variables to be recoded as well as recoded variables with recoded values for a subset of records. Can be dataframe or dataframe name (in quotes). 
    - varlist_orig - list of labels for the original variables to be recoded   
    - varlist_tr - list of labels for the recoded variables   
    - lookup - NULL (no parameter input, the default) if the input table (df) will be used to lookup recode values; label for the data frame that contains the variables to be recoded if an external table is used. The external table has to be opened as a dataframe, and must contain all of the columns in varlist_orig and varlist_tr. Can be dataframe or dataframe name (in quotes).   
  - output: data frame that can be used as a lookup table.


