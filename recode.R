
# Blom formula
# y <- qnorm((r-c)/(N-2c+1))
# where qnorm is the standard normal quantile function
# r is the rank of your variable
# c = 3/8, the Blom constant
# N is the sample size

# ================================ blom transformation ============================================

# recodeBlom recodes original variables using a Blom rank based inverse normal transformation.
# recodeBlom parameters:
#   df - label for the data frame that contains the variables to be recoded (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables

recodeBlom <- function(df,varlist_orig,varlist_tr){
  tr <- eval(parse(text=df))
  for (j in 1:length(varlist_orig)){
      varlab <- varlist_tr[j]
      tr[,varlab] <- ifelse(is.na(tr[,varlist_orig[j]]),NA,qnorm((rank(tr[,varlist_orig[j]], ties="average")-0.375)/(sum(!is.na(tr[,varlist_orig[j]])) - 2*0.375 + 1))) 
  }
  return(tr)
}

# ------------------------------------- end blom ------------------------------------------------

# ------------------------------------ Example Code --------------------------------------------


# ================================ recode into ordinal variables ==============================

# ------------------------------- Combined Interval and Quantile -------------------------------

# recodeOrdinal recodes original variables into ordinal variables using either equal-quantile 
#   or equal-interval recoding. Equal-quantile divides the original variable into k categories
#   based on k quantiles; the transformed variable will have a uniform distribution. Equal-interval
#   divides the original variable into k categories with the same range the original variable units.
#   Equal-interval will start with the number of specified categories (ncat) but will iteratively
#   decrease this number until the minimum number of observations (nobs) is achieved in each category.
#   Equal-quantile is approrpiate for continuous measures where the sample size divided by the
#   number of categoriesis greater than the minimum number of observations for each response
#   category. Equal-quantile does not work well with ordinal scales that have highly skewed
#   distributions. It is important to check the distributions of the recoded variables to
#   assure that the recoding matches the purpose and goals for recoding.
# recodeOrdinal parameters:
#   df - label for the data frame that contains the variables to be recoded (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables
#   type - "interval" for equal-interval recoding, "quantile" for equal quantile
#   ncat - maximum number of categories for recoded variables
#   nobs - miminum number of observed responses in each category of the transformed variable

recodeOrdinal <- function(df,varlist_orig,varlist_tr,type="interval",ncat=10, nobs=10) {
  tro <- eval(parse(text=df))
  for (j in 1:length(varlist_orig)){
    cuts <- {}
    min <- min(tro[,varlist_orig[j]], na.rm=TRUE)
    max <- max(tro[,varlist_orig[j]], na.rm=TRUE)
    if (type == "interval"){
      go <- TRUE
      repeat{
        for (i in 1:(ncat+1)){
          cuts[i] <- min + ((i-1) * ((max - min)/ncat))
        }
        ordt <- cut(tro[,varlist_orig[j]], breaks=cuts, include.lowest=TRUE,
              labels=c(1:ncat))
        if (min(table(ordt)) >= nobs | ncat == 2) {go <- FALSE}
        if (!go) {break}
        ncat <- ncat-1
        cuts <- {}
      }
    } else {
        for (i in 1:(ncat+1)){
          if (i==1){
            cuts[i] <- min
          } else{
            cuts[i] <- quantile(tro[,varlist_orig[j]],(i-1) * (1.0/ncat), na.rm=TRUE)
          }  
        }
        cuts <- unique(cuts)
      }
    varlab <- varlist_tr[j]
    tro[,varlab] <- as.numeric(cut(tro[,varlist_orig[j]], breaks=cuts, 
            include.lowest=TRUE))
  }
  return(tro)
}

# -------------------------------- end ordinal recode ------------------------------------------- 


# =================================== Recode Lookup ================================================== 

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
#   lookup - NA if the input table (df) will be used to lookup recode values; label for the data frame 
#     that contains the variables to be recoded (in quotes) if an external table is used. The external
#     table has to be opened as a dataframe, and must contain all of the columns 
#     in varlist_orig and varlist_tr.

recodeLookup <- function(df,varlist_orig,varlist_tr,type="continuous",lookup=NA) {
  rcd <- eval(parse(text=df))
  if (is.na(lookup)){
    rcdlu <- rcd[,c(varlist_orig,varlist_tr)]
  } else{
    rcdlu <- eval(parse(text=lookup))
    rcdlu <- rcdlu[,c(varlist_orig,varlist_tr)]
  }
  for (j in 1:length(varlist_orig)){
    t5 <- unique(rcdlu[!is.na(rcdlu[,varlist_tr[j]]),c(varlist_orig[j],varlist_tr[j])])
    t5 <- t5[order(t5[varlist_orig[j]]),]
    luv <- as.data.frame(cbind(t5,rbind(t5[2:nrow(t5),],c(NA,NA))))
    colnames(luv) <- c("orig_min","tr_min","orig_max","tr_max")
    mino <- min(rcdlu[,varlist_orig[j]], na.rm=TRUE)
    mint <- min(rcdlu[,varlist_tr[j]], na.rm=TRUE)
    maxo <- max(rcdlu[,varlist_orig[j]], na.rm=TRUE)
    maxt <- max(rcdlu[,varlist_tr[j]], na.rm=TRUE)
    rcd[,varlist_tr[j]] <- ifelse(rcd[,varlist_orig[j]] <= mino,mint,rcd[,varlist_tr[j]])
    rcd[,varlist_tr[j]] <- ifelse(rcd[,varlist_orig[j]] >= maxo,maxt,rcd[,varlist_tr[j]])
    t3 <- rcd[,c(varlist_orig[j],varlist_tr[j])]
    sqlcd <- paste("SELECT * FROM t3 AS t3
        LEFT JOIN
        (SELECT * FROM luv)
        AS luv1 ON t3.",varlist_orig[j],
        " >= luv1.orig_min AND t3.", varlist_orig[j],
        " < luv1.orig_max",sep="")
    t4 <- sqldf(sqlcd)
    if (type == "continuous") {
      t4[,varlist_tr[j]] <- ifelse(!is.na(t4[,varlist_tr[j]]),t4[,varlist_tr[j]],( ( (t4[,varlist_orig[j]] - t4$orig_min) / (t4$orig_max - t4$orig_min) ) * (t4$tr_max - t4$tr_min) ) + t4$tr_min)
    } else {
      t4[,varlist_tr[j]] <- ifelse(!is.na(t4[,varlist_tr[j]]),t4[,varlist_tr[j]],t4$tr_min)
    }
    rcd[,varlist_tr[j]] <- t4[,varlist_tr[j]]
  }
 return(rcd) 
}
                                                      
# ------------------------------------ end recode lookup --------------------------------------------
