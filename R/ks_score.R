#' KS score
#'
#' This function score (0-1000) a dataset with a binary
#' classification model, also compute the KS score of the model.
#'
#' @param data R data frame object.
#' @param model R model object.
#' @param score Name of you score column.
#' @param y Name of you target column.
#' @param sample_col Name of you sample column.
#' @param return_data If TRUE function will return the scored dataset.
#' @param return_conc_max If TRUE function will return the conc_max.
#' @param ks_precision If TRUE function will score your dataset before compute the KS Score.
#' @param score_var_name Name desired of computed score column.
#' @param return_ks_num Name of you sample that the function you return the KS Score numeric.
#' @param warn Indicate if function will show warnings.
#' @param subsample_data Proportion of your data that function will run.
#' @param seed Seed to subsample rows.#' @param train_sample Name of you sample that you trained the model.
#' @return scored dataset and the KS score.
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @import stats 
#' @export
ks_score <- function(data,
                     y,
                     sample_col,
                     model = NULL,
                     score = NULL,
                     return_data = TRUE,
                     return_conc_max = FALSE,
                     ks_precision = TRUE,
                     quebras = NULL,
                     score_var_name = "CON_stepKS_score",
                     return_ks_num = NULL,
                     warn = NULL,
                     subsample_data=1,
                     seed=NULL,
                     train_sample = "DES"){

  if((is.null(model) & is.null(score)) || (!is.null(model) & !is.null(score))) stop("Error: Only one param: 'model' or 'score', must be passed.")

  col_nameset <- colnames(data)

  if(!(y %in% col_nameset)) stop(paste("Error:",y," not found in data."))
  if(!(sample_col %in% col_nameset)) stop(paste("Error:",sample_col," not found in data."))

  data$sample_col = data[,sample_col]
  samples_names = unique(data$sample_col)

  if(!is.null(return_ks_num)) if(!(return_ks_num %in% samples_names)) stop(paste("Error:",
                                                                                 return_ks_num,"not found in",sample_col,"."))
  if(!(train_sample %in% samples_names)) stop(paste("Error:",train_sample,"not found in",sample_col,"."))

  if(!is.null(warn)){
    default_w <- getOption("warn")
    options(warn = warn)
  }

  if(!ks_precision) return_data = FALSE

  #if(!is.null(score)) score_var_name = score

  if(!is.null(seed)) set.seed(seed)
  if(subsample_data<0.99999) data = data[sample(nrow(data),nrow(data)*subsample_data),]

  # Flags para KS ----
  for (sample_name in samples_names) {
    data[[sample_name]] = ifelse(data$sample_col == sample_name, 1, 0)
  }

  if (!is.null(model)){
    if(!ks_precision){
      if(!is.null(return_ks_num)){
        data = data[data$sample_col == return_ks_num,]

        data$predict = predict(model, newdata = data)
        ks_return = ks.test(data$predict[data[,y]==0 & data[,return_ks_num]==1],
                            data$predict[data[,y]==1 & data[,return_ks_num]==1 ])

        if(!is.null(warn)) options(warn = default_w)
        return(as.numeric(ks_return$statistic))
      }
      data$score = predict(model, newdata = data)
    }else{
      # Processo de adicionar os Scores ----
      data$p0 = predict(model, newdata = data, type = "response")
      xbeta = -log( (1/(data$p0)) - 1 )
      aux = trunc(500 + xbeta * (100/log(2)) )

      # Faixas dos Scores ----
      if(is.null(quebras)) quebras = quantile(aux[data$sample_col == train_sample],
                                              probs = c(0.000,0.025,0.070,0.150,0.300,0.500,0.700,0.850,0.930,0.975,1.000))
      # Scores ----
      aux2 = ifelse( aux <= quebras[1]  , 0 ,
                     ifelse( aux <= quebras[2], (0   + ( ( aux - quebras[1] ) * ( 100 - 0    ) ) / ( quebras[2]  - quebras[1])),
                             ifelse( aux <= quebras[3], (101 + ( ( aux - quebras[2] ) * ( 200 - 101  ) ) / ( quebras[3]  - quebras[2])),
                                     ifelse( aux <= quebras[4], (201 + ( ( aux - quebras[3] ) * ( 300 - 201  ) ) / ( quebras[4]  - quebras[3])),
                                             ifelse( aux <= quebras[5], (301 + ( ( aux - quebras[4] ) * ( 400 - 301  ) ) / ( quebras[5]  - quebras[4])),
                                                     ifelse( aux <= quebras[6], (401 + ( ( aux - quebras[5] ) * ( 500 - 401  ) ) / ( quebras[6]  - quebras[5])),
                                                             ifelse( aux <= quebras[7], (501 + ( ( aux - quebras[6] ) * ( 600 - 501  ) ) / ( quebras[7]  - quebras[6])),
                                                                     ifelse( aux <= quebras[8], (601 + ( ( aux - quebras[7] ) * ( 700 - 601  ) ) / ( quebras[8]  - quebras[7])),
                                                                             ifelse( aux <= quebras[9], (701 + ( ( aux - quebras[8] ) * ( 800 - 701  ) ) / ( quebras[9]  - quebras[8])),
                                                                                     ifelse( aux <= quebras[10],(801 + ( ( aux - quebras[9] ) * ( 900 - 801  ) ) / ( quebras[10] - quebras[9])),
                                                                                             ifelse( aux <= quebras[11],(901 + ( ( aux - quebras[10]) * (1000 - 901  ) ) / ( quebras[11] - quebras[10])),
                                                                                                     1000 )))))))))))

      data$score =  trunc(aux2)
    }
  }else data$score = data[,score]

  if(!is.null(return_ks_num)){
    if(!is.null(warn)) options(warn = default_w)
    ks_return = ks.test(data$score[data[,y]==0 & data[,return_ks_num]==1],
                        data$score[data[,y]==1 & data[,return_ks_num]==1 ])
    return(as.numeric(ks_return$statistic))
  }

  ks = list()
  ks_num = list()

  for (sample_name in samples_names) {
    if(sum(is.na(data[data$sample_col == sample_name,][,y])) == 0){
      ks_tmp  = ks.test(data$score[data[,y]==0 & data[,sample_name]==1],
                        data$score[data[,y]==1 & data[,sample_name]==1 ])
      ks[[sample_name]] = ks_tmp
      ks_num[[sample_name]] = as.numeric(ks_tmp$statistic)
    }
  }

  ks$TOT = ks.test(data$score[data[,y]==0] , data$score[data[,y]==1])
  ks_num$TOT = as.numeric(ks$TOT$statistic)

  if(!is.null(warn)) options(warn = default_w)
  if(!return_data) return(ks_num)

  # Calcular Conc_max ----
  if(return_conc_max){
    library(dplyr)
    conc_max_tmp = group_by(data, score)
    conc_max_tmp = summarise(conc_max_tmp,Total = n(), Rep = n() /nrow(conc_max_tmp))
    conc_max_tmp = arrange(conc_max_tmp,desc(Rep))

    conc_max <<-  as.numeric(conc_max_tmp[1,3]*100)
  }

  colnames(data)[colnames(data) == 'score'] <- score_var_name

  ks_list <<- ks
  ks_list_num <<- ks_num
  if(is.null(score)){
    data_score <<- data
    quebras <<- quebras
  }
  return(ks_num)

}
