#' model_sig
#'
#' This function return TRUE if you model have all
#' coefficients significatives and FALSE otherwise.
#'
#' @param model R model object.
#' @param p_value p_value.
#' @param ignore_intercept_sig If TRUE function will ignore intercept significance.
#' @return TRUE or FALSE
#' @import stats 
#' @export
model_sig = function(model,
                     one_cat_sig = FALSE,
                     p_value=0.05,
                     ignore_intercept_sig = TRUE,
                     trim_chars = 2){

  df = data.frame(coef(summary(model)))
  names(df)[names(df) == 'Pr...z..'] = 'x'

  if(ignore_intercept_sig) df = df[-1,]
  if(one_cat_sig) df = aggregate(df$x, by=list(strtrim(x = rownames(df),width = nchar(rownames(df)) - trim_chars)), FUN=min)

  ifelse(sum(df$x>p_value),return(FALSE),return(TRUE))
}