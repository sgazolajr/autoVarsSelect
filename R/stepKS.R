#' stepKS
#'
#' This function select the best pool of variables of logistic regression
#' models based on KS metric and significative coefficients.
#'
#' @param data R data frame object.
#' @param y Name of you target column.
#' @param sample_col Name of you sample column.
#' @param train_sample Name of you sample that you trained the model.
#' @param include List with name of variables to be analyzed by the algorithm.
#' @param exclude If you do not want to pass the name of the Data Frame variables to be analyzed by the algorithm, you also have to pass the list of data frame variables that are not analyzed, for example the key variables, and everything else will be analyzed.
#' @param start List of variables that will already begin in the model.
#' @param sig_mode Parameter that defines whether the algorithm will only select models that have all significant coefficients, or if it can return models with non-significant coefficients.
#' @param direction "both", "forward" or "both_sophisticated".
#' @param sig_mode 'off', 'one_cat' or 'all'.
#' @param return_type "none","vars","model","formula","scored_data","ahead_vars" or "ahead_model".
#' @param vars_enable_both Number of vars in model to enble both method.
#' @param near_sample Defines a second sample in which the algorithm will always keep close to the main sample. The sample will need to be in the Data Frame.
#' @param pct_ks_dif "Defines how close the KS to the ks_comp need to be close to the ks_focus.
#' @param steps_ahead If in the current step when adding the best variable it does not get a KS higher than the KS of the previous step, this parameter indicates how many steps the algorithm will walk forward of the step with the best KS to try to find some variable that increases the KS.
#' @param progress_bar If TRUE function will show a progress_bar in every step.
#' @param show_time_elapsed If TRUE function will show the time_elapsed in every step.
#' @param ks_precision If TRUE function will score your dataset before compute the KS Score, for ervery model tested.
#' @param p_value p_value.
#' @param ignore_intercept_sig If TRUE function will ignore intercept significance.
#' @param warn Indicate if function will show warnings.
#' @param subsample_rows Proportion of your data that function will run.
#' @param seed Seed to subsample rows.
#' @param subsample_columns Proportion of your data that function will run.
#' @return A 'glm' model with the best pool of vars.
#' @import stats 
#' @export
stepKS <- function(data,
                   y,
                   sample_col,
                   val_sample,
                   train_sample = "DES",
                   exclude=NULL,
                   include=NULL,
                   start=NULL,
                   force_in_model=NULL,
                   return_type = "none",
                   sig_mode = "all",
                   direction = "both_sophisticated",
                   vars_enable_both = 5,
                   near_sample=NULL,
                   pct_ks_dif=0.05,
                   flag_bad=1,
                   steps_ahead = 5,
                   max_cat = 10,
                   ks_precision = FALSE,
                   progress_bar = TRUE,
                   show_time_elapsed = TRUE,
                   ignore_intercept_sig = TRUE,
                   p_value=0.05,
                   trim_chars = 2,
                   subsample_columns=1,
                   warn = NULL,
                   subsample_rows=1,
                   seed=NULL,
                   debug = FALSE,
                   debug2 = FALSE){

  start_time = Sys.time()

  epsilon = 0.000000001

  direction_list = c("both","forward","both_sophisticated")
  sig_mode_list = c("off","one_cat","all")
  return_type_list = c("none","vars","model","formula","scored_data","ahead_vars","ahead_model")

  if(!return_type  %in% return_type_list) stop('Direction must be "none","vars","model","formula","scored_data","ahead_vars" or "ahead_model".')
  if(!direction  %in% direction_list) stop("Direction must be 'forward', 'both' or 'both_sophisticated'.")
  if(!sig_mode %in% sig_mode_list) stop("sig_mode must be 'off', 'one_cat' or 'all'.")

  if(sig_mode == "off" & direction == "both_sophisticated") direction = "both"

  sig_analyze = ifelse(sig_mode != "off",TRUE,FALSE)
  one_cat_sig = ifelse(sig_mode == "one_cat",TRUE,FALSE)
  sophisticated_both = ifelse(direction == "both_sophisticated",TRUE,FALSE)


  if(debug2) progress_bar = FALSE
  if(!is.null(exclude) & !is.null(include)) stop("Just one param: 'include' or 'exclude', must be passed.")

  col_nameset <- colnames(data)

  if(!(y %in% col_nameset)) stop(paste(y," not found in data."))
  if(!(sample_col %in% col_nameset)) stop(paste(sample_col," not found in data."))
  if(!is.null(include)) if(sum(!(include %in% col_nameset))) stop("Some 'include' var not found in data.")
  if(!is.null(start)) if(sum(!(start %in% col_nameset))) stop("Some 'start' var not found in data.")
  if(!is.null(force_in_model)) if(sum(!(force_in_model %in% col_nameset))) stop("Some 'force_in_model' var not found in data.")

  data$sample_col = data[,sample_col]
  samples_names = unique(data$sample_col)

  if(!(val_sample %in% samples_names)) stop(paste(val_sample,"not found in",sample_col,"."))
  if(!is.null(near_sample)) if(!(near_sample %in% samples_names)) stop(paste(near_sample,"not found in",sample_col,"."))
  if(!(train_sample %in% samples_names)) stop(paste(train_sample,"not found in",sample_col,"."))
  
  if(!is.null(warn)){
    default_w <- getOption("warn")
    options(warn = warn)
  }

  if(is.null(include)){
    Xnames = col_nameset[!col_nameset %in% c(y,exclude,sample_col)]
  }else Xnames = include

  if (!is.null(near_sample)) {
    data = data[data$sample_col == train_sample | data$sample_col == val_sample | data$sample_col == near_sample,]
  }else data = data[data$sample_col == train_sample | data$sample_col == val_sample,]

  if(!is.null(seed)) set.seed(seed)
  if(subsample_rows<1-epsilon) data = data[sample(nrow(data),nrow(data)*subsample_rows),]



  for (var in Xnames) {

    if(debug) print(var)

    if(!(is.numeric(as.matrix(data[,var])))){
      if(length(unique(data[,var])) > max_cat | length(unique(data[,var])) == 1){
        Xnames = Xnames[!Xnames %in% var]
      }
    }
  }

  if(sum(is.na(data[,c(Xnames,y)]))) stop("NA found in data.")


  
  best_model = NULL
  best_fm = NULL
  vars_to_analize = NULL
  best_var = NULL
  best_ks = 0
  best_ks_iter = 0
  #ks_num = 0
  kscomp_num = 0
  count_steps_ahead = 0
  step = 0
  vars_in_model = unique(c(start,force_in_model))
  vars_remains = unique(Xnames)
  vars_remains = vars_remains[!vars_remains %in% vars_in_model]

  width_o = getOption("width")/3

  while (best_ks_iter + epsilon >best_ks | count_steps_ahead<steps_ahead) {

    best_model_iter = NULL
    best_ks_iter = 0
    best_fm_iter = NULL
    var_to_left_iter = NULL
    vars_remains_for = vars_remains

    if(subsample_columns<1-epsilon) vars_remains_for = vars_remains_for[sample(length(vars_remains_for),
                                                                               length(vars_remains_for)*subsample_columns)]

    step_var = 0
    if(progress_bar){
      len_vars_remains = length(vars_remains_for)
      pb <- txtProgressBar(min = 0, max = len_vars_remains, style = 3, width = width_o)
    }
    for (var in vars_remains_for) {
      var_to_left = NULL
      var_replace = NULL

      if(debug2) print(var)
      step_var = step_var+1
      if(progress_bar) setTxtProgressBar(pb, step_var)

      if(flag_bad == 1){
        fm <- paste("(1-",y,") ~ ",paste(vars_in_model,collapse = "+"),"+",var,sep="")
      }else fm <- paste(y," ~ ",paste(vars_in_model,collapse = "+"),"+",var,sep="")

      model <- glm(fm,data[data$sample_col == train_sample,],family="binomial")
      model$data=NULL

      if(sig_analyze & !model_sig(model = model,
                                  one_cat_sig = one_cat_sig,
                                  trim_chars = trim_chars,
                                  p_value = p_value,
                                  ignore_intercept_sig = ignore_intercept_sig)){

        if(sophisticated_both & length(vars_in_model) > vars_enable_both & !count_steps_ahead){

          ks_num_sb = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                               return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon

          var_replace = var

          vars_in_model_to_test = vars_in_model[!vars_in_model %in% force_in_model]

          ks_num = 0
          if(debug) if((ks_num_sb>best_ks) & (ks_num_sb>best_ks_iter)) print("Sophis. Back.")
          if((ks_num_sb>best_ks) & (ks_num_sb>best_ks_iter)) vars_to_analize = c(vars_to_analize,var_replace)
          if((ks_num_sb>best_ks) & (ks_num_sb>best_ks_iter)) for (var in vars_in_model_to_test){

            vars_to_test = vars_in_model[!vars_in_model %in% var]

            if(flag_bad == 1){
              fm_iter2 <- paste("(1-",y,") ~ ",paste(vars_to_test,collapse = "+"),"+",var_replace,sep="")
            }else fm_iter2 <- paste(y," ~ ",paste(vars_to_test,collapse = "+"),"+",var_replace,sep="")

            model_iter2 <- glm(fm_iter2,data[data$sample_col == train_sample,],family="binomial")
            model_iter2$data=NULL

            if(!model_sig(model = model,
                          one_cat_sig = one_cat_sig,
                          trim_chars = trim_chars,
                          p_value = p_value,
                          ignore_intercept_sig = ignore_intercept_sig)){
              ks_num_iter2 = -1
              kscomp_num_iter2 = -1

            }else if(!is.null(near_sample)){
              ks = ks_score(data = data,model = model_iter2,y = y,sample_col = sample_col,train_sample = train_sample,
                            return_data = F,ks_precision = ks_precision)
              ks_num_iter2 = as.numeric(ks[val_sample]) + epsilon
              kscomp_num_iter2 = as.numeric(ks[near_sample]) + epsilon
            }else ks_num_iter2 = ks_score(data = data,model = model_iter2,y = y,sample_col = sample_col,train_sample = train_sample,
                                          return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon

            if(is.null(near_sample)){
              if(ks_num_iter2 >  ks_num) {
                ks_num = ks_num_iter2
                model = model_iter2
                var_to_left = var
                fm = fm_iter2
              }
            }else if(ks_num_iter2 >  ks_num & ks_num_iter2 - kscomp_num_iter2 < pct_ks_dif) {
              ks_num = ks_num_iter2
              model = model_iter2
              var_to_left = var
              fm = fm_iter2
            }
          }
          if(debug) if((ks_num_sb>best_ks) & (ks_num_sb>best_ks_iter)) print(ks_num)

          var =var_replace

        }else{
          ks_num = -1
          kscomp_num = -1
        }

      }else if(!is.null(near_sample)){
        ks = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                      return_data = F,ks_precision = ks_precision)
        ks_num = as.numeric(ks[val_sample]) + epsilon
        kscomp_num = as.numeric(ks[near_sample]) + epsilon
      }else ks_num = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                              return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon
      if(debug) print(ks_num)

      if(is.null(near_sample)){
        if(ks_num >  best_ks_iter) {
          best_ks_iter = ks_num
          best_model_iter = model
          best_var = var
          best_fm_iter = fm
          var_to_left_iter = var_to_left
        }
      }else if(ks_num >  best_ks_iter & ks_num - kscomp_num < pct_ks_dif) {
        best_ks_iter = ks_num
        best_kscomp_iter = kscomp_num
        best_model_iter = model
        best_var = var
        best_fm_iter = fm
        var_to_left_iter = var_to_left
      }
      if(debug) print(best_ks_iter)
      if(debug) print(var_to_left_iter)

    }
    if(progress_bar) close(pb)

    vars_in_model_to_test = vars_in_model[!vars_in_model %in% force_in_model]

    if(debug) if((direction == "both" | direction == "both_sophisticated") & is.null(var_to_left_iter) &
                 length(vars_in_model) > vars_enable_both) print("Backward")
    if((direction == "both" | direction == "both_sophisticated") & is.null(var_to_left_iter) &
       length(vars_in_model) > vars_enable_both & !count_steps_ahead) for (var in vars_in_model_to_test){

         vars_to_test = vars_in_model[!vars_in_model %in% var]

         if(!is.null(best_var)){
           if(flag_bad == 1){
             fm <- paste("(1-",y,") ~ ",paste(vars_to_test,collapse = "+"),"+",best_var,sep="")
           }else fm <- paste(y," ~ ",paste(vars_to_test,collapse = "+"),"+",best_var,sep="")
         }else{
           if(flag_bad == 1){
             fm <- paste("(1-",y,") ~ ",paste(vars_to_test,collapse = "+"),sep="")
           }else fm <- paste(y," ~ ",paste(vars_to_test,collapse = "+"),sep="")
         }
         model <- glm(fm,data[data$sample_col == train_sample,],family="binomial")
         model$data=NULL

         if(sig_analyze & !model_sig(model = model,
                                     one_cat_sig = one_cat_sig,
                                     trim_chars = trim_chars,
                                     p_value = p_value,
                                     ignore_intercept_sig = ignore_intercept_sig)){
           ks_num = -1
           kscomp_num = -1
         }else if(!is.null(near_sample)){
           ks = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                         return_data = F,ks_precision = ks_precision)
           ks_num = as.numeric(ks[val_sample]) + epsilon
           kscomp_num = as.numeric(ks[near_sample]) + epsilon
         }else ks_num = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                                 return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon


         if(is.null(near_sample)){
           if(ks_num >  best_ks_iter) {
             best_ks_iter = ks_num
             best_model_iter = model
             var_to_left_iter = var
             best_fm_iter = fm
           }
         }else if(ks_num >  best_ks_iter & ks_num - kscomp_num < pct_ks_dif) {
           best_ks_iter = ks_num
           best_kscomp_iter = kscomp_num
           best_model_iter = model
           var_to_left_iter = var
           best_fm_iter = fm
         }
       }
    if(debug) if(!is.null(var_to_left_iter)) cat(paste("var_to_left_iter:",var_to_left_iter))
    if(!is.null(var_to_left_iter)) vars_to_analize = c(vars_to_analize,var_to_left_iter)

    if(is.null(best_model) & step > 1) stop("Not found significative model, change 'start' vars.")

    if(best_ks_iter == 0){
      model_best <<- best_model
      model_ahead <<- best_model_iter

      model_best_ks <<- best_ks
      model_ahead_ks <<- best_ks_iter

      best_ks_preciso = ks_score(data = data,model = best_model,y = y,sample_col = sample_col,train_sample = train_sample,
                                 return_ks_num = val_sample,return_data = F,ks_precision = TRUE)

      vars_to_analize = vars_to_analize[!vars_to_analize %in% best_model_vars]
      vars_to_analize = unique(vars_to_analize)

      best_vars_to_analize = NULL

      for (var in vars_to_analize){

        if(flag_bad == 1){
          fm <- paste("(1-",y,") ~ ",paste(best_model_vars,collapse = "+"),"+",var,sep="")
        }else fm <- paste(y," ~ ",paste(best_model_vars,collapse = "+"),"+",var,sep="")

        model <- glm(fm,data[data$sample_col == train_sample,],family="binomial")
        model$data=NULL

        ks_num = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                          return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon

        if(ks_num >  best_ks_preciso) best_vars_to_analize = c(best_vars_to_analize,var)
      }

      vars_to_analize <<- vars_to_analize
      best_vars_to_analize <<- best_vars_to_analize

      end_time = Sys.time()
      elapsed_time = as.character.Date(round(end_time - start_time,digits = 1))

      print(paste("Total time elapsed: ",elapsed_time,"| Best model found => KS",val_sample,":",
                  (round(best_ks_preciso,4)*100),"% | Formula:", best_fm))

      if(!is.null(warn)) options(warn = default_w)

      return_type_list = c("none","vars","model","formula","scored_data","ahead_vars","ahead_model")

      if(return_type == "vars"){
        return(best_model_vars)
      }else if(return_type == "model"){
        return(best_model)
      }else if(return_type == "ahead_model"){
        return(best_model_iter)
      }else if(return_type == "ahead_vars"){
        return(vars_in_model)
      }else if(return_type == "formula"){
        return(best_fm)
      }else if(return_type == "scored_data"){
        ks_score(data = data,model = best_model,y = y,sample_col = sample_col,train_sample = train_sample,
                 return_data = TRUE ,ks_precision = TRUE)
        return(scored_data)
      }else if(return_type == "none"){
        return(NULL)
      }else return(NULL)
    }


    vars_in_model = vars_in_model[!vars_in_model %in% var_to_left_iter]
    vars_in_model = c(vars_in_model,best_var)
    vars_in_model = unique(vars_in_model)

    vars_remains = vars_remains[!vars_remains %in% vars_in_model]
    vars_remains = c(vars_remains,var_to_left_iter)
    vars_remains =  unique(vars_remains)

    if(best_ks_iter >  best_ks){
      best_ks = best_ks_iter
      best_model = best_model_iter
      count_steps_ahead=0
      best_fm = best_fm_iter
      best_model_vars = vars_in_model

    }else count_steps_ahead = count_steps_ahead + 1

    if(show_time_elapsed){
      end_time = Sys.time()
      elapsed_time = as.character.Date(round(end_time - start_time,digits = 1))
      cat(paste("Time elapsed:",elapsed_time," "))
    }
    if(count_steps_ahead)  cat(paste("| Ahead:",count_steps_ahead," "))
    if(!is.null(near_sample)) cat(paste("| KS",near_sample,":",(round(best_kscomp_iter,4)*100),"% "))
    print(paste("KS",val_sample,":",(round(best_ks_iter,4)*100),"% | Formula:", best_fm_iter))
    # print(summary(best_model_iter))
    step = step+1

  }


  model_best <<- best_model
  model_ahead <<- best_model_iter

  model_best_ks <<- best_ks
  model_ahead_ks <<- best_ks_iter

  best_ks_preciso = ks_score(data = data,model = best_model,y = y,sample_col = sample_col,train_sample = train_sample,
                             return_ks_num = val_sample,return_data = F,ks_precision = TRUE)

  vars_to_analize = vars_to_analize[!vars_to_analize %in% best_model_vars]
  vars_to_analize = unique(vars_to_analize)

  best_vars_to_analize = NULL

  for (var in vars_to_analize){


    if(flag_bad == 1){
      fm <- paste("(1-",y,") ~ ",paste(best_model_vars,collapse = "+"),"+",var,sep="")
    }else fm <- paste(y," ~ ",paste(best_model_vars,collapse = "+"),"+",var,sep="")

    model <- glm(fm,data[data$sample_col == train_sample,],family="binomial")
    model$data=NULL

    ks_num = ks_score(data = data,model = model,y = y,sample_col = sample_col,train_sample = train_sample,
                      return_ks_num = val_sample,return_data = F,ks_precision = ks_precision) + epsilon

    if(ks_num >  best_ks_preciso) best_vars_to_analize = c(best_vars_to_analize,var)
  }

  vars_to_analize <<- vars_to_analize
  best_vars_to_analize <<- best_vars_to_analize

  end_time = Sys.time()
  elapsed_time = as.character.Date(round(end_time - start_time,digits = 1))

  print(paste("Total time elapsed: ",elapsed_time,"| Best model found => KS",val_sample,":",
              (round(best_ks_preciso,4)*100),"% | Formula:", best_fm))

  if(!is.null(warn)) options(warn = default_w)

  if(return_type == "vars"){
    return(best_model_vars)
  }else if(return_type == "model"){
    return(best_model)
  }else if(return_type == "ahead_model"){
    return(best_model_iter)
  }else if(return_type == "ahead_vars"){
    return(vars_in_model)
  }else if(return_type == "formula"){
    return(best_fm)
  }else if(return_type == "scored_data"){
    ks_score(data = data,model = best_model,y = y,sample_col = sample_col,train_sample = train_sample,
             return_data = TRUE ,ks_precision = TRUE)
    return(scored_data)
  }else if(return_type == "none"){
    return(NULL)
  }else return(NULL)
}
