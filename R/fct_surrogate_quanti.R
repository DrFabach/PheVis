#' fct_surrogate_quanti
#' 
#' @description Compute the quantitative surrogate and then apply thresholds to get the qualitative surrogate.
#'
#' @param main_icd Character vector of the column names of the main ICD codes.
#' @param main_cui Character vector of the column names of the main CUIs.
#' @param df Dataframe containing all variables.
#' @param half_life Duration of accumulation. For a chronic disease you might chose Inf, for acute disease you might chose the duration of the disease.
#' @param date Column name of the time column. The time column should be numeric
#' @param patient_id Column name of the patient id column.
#' @param encounter_id Column name of the encounter id column.
#' @param omega Constant for the extrema population definition.
#' @param param param of a previous train_phevis() result.
#' @param best_encounter A bolean, if True Surrogate only on current Encounter
#' @param by_patient A Bolean, if True and best_encoutner T, select for each patient the best and worst encounter distribute among surrogate
#' @param    standardise
#' @param   other_icd
#' @param   other_cui 
#' @return A list
#' 
#' \itemize{
#'  \item table - Main result: \code{data.frame} with the rolling variables and the surrogates
#'  \item param - the parameters for the standardisation of ICD and CUI
#'  \item roll_all - a subset of table with the rolling variables only
#'  \item quantile_vec - the quantile defining the extrema populations
#' }
fct_surrogate_quanti <- function(main_icd,
                                 main_cui,
                                 df,
                                 half_life,
                                 date,
                                 patient_id,
                                 encounter_id,
                                 omega = 2,
                                 param = NULL,
                                 best_encounter=F,
                                 by_patient = F,
                                 standardise=F,
                                 other_icd=NULL,
                                 other_cui=NULL){
        ##### compute quantitative surrogate
        # icd

        tot_icd_all <- df[,colnames(df) %in% main_icd]
        if(length(main_icd) > 1) tot_icd_all <- rowSums(tot_icd_all)
        # nlp
        if(length(main_cui) > 1){
                tot_nlp_all <- rowSums(df[,colnames(df) %in% main_cui])
        } else {
                tot_nlp_all <- df[[main_cui]]
        }
        # normalized sum
        if(standardise) {
          if(is.null(other_icd)|is.null(other_cui)){
            stop("other_icd and other_cui should be existing")
          }else{
            nb_icd <- rowSums( df[,colnames(df) %in% other_icd])
            nb_cui <- rowSums( df[,colnames(df) %in% other_cui])
            tot_icd_all <- tot_icd_all/(nb_icd+1)
            tot_nlp_all <- tot_nlp_all/(nb_cui+1)
            }
        }
        sur_norm_both_all <- norm_var(tot_icd_all) + norm_var(tot_nlp_all)
        
        min_sur <- min(sur_norm_both_all)
        
        if(!is.null(param)){
                sur_norm_both_all <- (tot_icd_all - param$mean_icd)/param$sd_icd + (tot_nlp_all - param$mean_nlp)/param$sd_nlp
                min_sur <- param$min_all
        }
        
        list_stand <- list(mean_icd = mean(tot_icd_all),
                           sd_icd = sd(tot_icd_all),
                           mean_nlp = mean(tot_nlp_all),
                           sd_nlp = sd(tot_nlp_all),
                           min_all = min(sur_norm_both_all))
        
        sur <- sur_norm_both_all + abs(min_sur)
        ##### exponential correction
        cumul <- sur_exp_smooth(half_life = half_life,
                                sur = sur,
                                date = df[[date]],
                                patient_id = df[[patient_id]],
                                encounter_id = df[[encounter_id]])
        cumul <- cumul %>%
                rename(PATIENT_NUM = .data$patient_id,
                       ENCOUNTER_NUM = .data$encounter_id,
                       START_DATE = .data$date,
                       SUR_QUANTI = .data$Correct_value) %>%
                ungroup()
        quantile_vec <- build_qantsur(df = df, var.icd = main_icd, omega = omega)
        
        if(!best_encounter){
        ##### quantile thresholds
        
        ##### surrogate quali
          
          if(by_patient){
            cumul <- cumul %>% group_by(PATIENT_NUM) %>% mutate(surquantimax = max(SUR_QUANTI)) 
            sur_quali <- build_quali(x = cumul$surquantimax,
                                     p = quantile_vec["p.sur"],
                                     q = quantile_vec["q.sur"])
            
            cumul <- cumul %>% mutate(SUR_QUALI = sur_quali)
            
          
            
          }else{
        sur_quali <- build_quali(x = cumul$SUR_QUANTI,
                                 p = quantile_vec["p.sur"],
                                 q = quantile_vec["q.sur"])
        
        cumul <- cumul %>% mutate(SUR_QUALI = sur_quali)
          }
        }else{
          if(by_patient){
            
            
            sur_quali <- build_quali(x = cumul$SUR_QUANTI,
                                     p = quantile_vec["p.sur"],
                                     q = quantile_vec["q.sur"])
            
            cumul <- cumul %>% mutate(SUR_QUALI = sur_quali, 
                                      SUR = sur)
            
           cumul<- cumul %>% 
             group_by(PATIENT_NUM ) %>% arrange(sur) %>% 
             mutate(SUR_QUALI = c(rep(0, sum(SUR_QUALI==0)),rep(3, sum(SUR_QUALI==3)),rep(1, sum(SUR_QUALI==1)))) %>% 
             ungroup
            
            
            
          }else{

          sur_quali <- build_quali(x = jitter(sur),
                                   p = quantile_vec["p.sur"],
                                   q = quantile_vec["q.sur"])
          
          cumul <- cumul %>% mutate(SUR_QUALI = sur_quali)
          }
        }
        ##### rolling surrogate
        roll_all <- rolling_var(id = df[,patient_id],
                                var = sur,
                                start_date = df[,date],
                                id_encounter = df[,encounter_id])
        
        result <- roll_all %>%
                rename(PATIENT_NUM = .data$id,
                       ENCOUNTER_NUM = .data$id_encounter) %>%
                select(.data$PATIENT_NUM,
                       .data$ENCOUNTER_NUM,
                       .data$last_vis,
                       .data$last_5vis,
                       .data$cum,
                       .data$cum_month,
                       .data$cum_year) %>%
                right_join(cumul)
        
        return(list(table = result,
                    param = list_stand,
                    roll_all = roll_all,
                    quantile_vec = quantile_vec))
}