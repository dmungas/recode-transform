recode-transform
================

This repository contains R code to recode and transform variables in a R dataframe.

Transformations that are performed include:
  Blom rank based inverse normal transformation
  Ordinal transformation using equal-quantile algorithm
  Ordinal transformation using equal-interval algorithm
  
The formula for the Blom transformation is: 
  y <- qnorm((r-c)/(N-2c+1))
    where qnorm is the standard normal quantile function
    r is the rank of your variable
    c = 3/8, the Blom constant
    N is the sample size
    
Equal-quantile ordinal transformation divides the original variable into k categories based on k quantiles; the transformed variable will have a uniform distribution.     

# ================================ blom transformation ============================================

# recodeBlom recodes original variables using a Blom rank based inverse normal transformation.
# recodeBlom parameters:
#   df - label for the data frame that contains the variables to be recoded (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables

# recodeOrdinal recodes original variables into ordinal variables using either equal-quantile 
#   or equal-interval recoding. Equal-quantile divides the original variable into k categories
#   based on k quantiles; the transformed variable will have a uniform distribution. Equal-interval
#   divides the original variable into k categories with the same range the original variable units.
#   Equal-interval will start with the number of specified categories (ncat) but will iterativel
#   decrease this number until the minimum number of observations (nobs) is achieved in each category
# recodeOrdinal parameters:
#   df - label for the data frame that contains the variables to be recoded (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables
#   type - "interval" for equal-interval recoding, "quantile" for equal quantile
#   ncat - maximum number of categories for recoded variables
#   nobs - miminum number of observed responses in each category of the transformed variable


# recodeLookup generates a reference table from recoded variables within a dataframe and uses
#   this as a lookup table to recode vales that have not yet been recoded. This is applicable
#   where a subset of records are used for the original recode and this subset is then used 
#   as the reference to recode the other records. It is applicable for continuous transformed variables
#   (including but not limited to Blom) and for ordinal transformations.
# recodeOrdinal parameters:
#   df - label for the data frame that contains the variables to be recoded (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables
#   type - "continuous" for continuous (numeric) transformation, "ordinal" for ordinal
