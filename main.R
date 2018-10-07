library(tercen)
library(dplyr)

calc_egfr <- function(crt, cyt, gender, age, race) {
  
  eGFR_0 <- {
    eGFR <- NULL 
    a <- ifelse (gender=="Female", -0.248, -0.207) 
    k <- ifelse (gender=="Female", 0.7, 0.9) 
    
    eGFR <-  135 * (min(crt/k, 1))^a * (max(crt/k, 1))^-0.601 * (min(cyt/0.8, 1))^-0.375 * (max(cyt/0.8, 1))^-0.711 * 0.995^age
    eGFR <- ifelse (gender=="Female",0.969 * eGFR, eGFR)
    eGFR <- ifelse (race=="Black",1.08 * eGFR, eGFR)
    eGFR
  }
  
  eGFR_1 <- {
    eGFR <- NULL 
    if (gender == "Female"){
      if   (crt <= 0.7) { eGFR <- 144 * (crt/ 0.7)^-0.329  * 0.995^age }
      if   (crt >  0.7) { eGFR <- 144 * (crt/ 0.7)^-1.209 * 0.995^age }
    }
    if (gender == "Male"){
      if   (crt <= 0.9) { eGFR <- 144 * (crt/ 0.9)^-0.411  * 0.995^age }
      if   (crt >  0.9) { eGFR <- 144 * (crt/ 0.9)^-1.209 * 0.995^age }
    }
    eGFR <- ifelse (race=="Black",1.159 * eGFR, eGFR)
    eGFR
  }
  
  eGFR_2 <- {
    eGFR <- NULL 
    if   (cyt <= 0.8) { eGFR <- 133 * (cyt/ 0.8)^-0.449 * 0.996^age }
    if   (cyt >  0.8) { eGFR <- 133 * (cyt/ 0.8)^-1.328 * 0.996^age }
    
    eGFR <- ifelse (gender=="Female",0.932 * eGFR, eGFR)
    eGFR
  }
  
  eGFR_3 <- {
    eGFR <- NULL  
    if (gender == "Female"){
      if   (crt <= 0.7 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.248 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt <= 0.7 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.248 * (cyt/ 0.8)^-0.711 * 0.995^age }
      if   (crt >  0.7 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt >  0.7 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.711 * 0.995^age }
    }
    if (gender == "Male"){
      if   (crt <= 0.9 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.207 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt <= 0.9 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.207 * (cyt/ 0.8)^-0.711 * 0.995^age }
      if   (crt >  0.9 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt >  0.9 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.711 * 0.995^age }
    }
    eGFR <- ifelse (race=="Black",1.159 * eGFR, eGFR)
    eGFR 
  }
  return(list(eGFR_0,eGFR_1, eGFR_2, eGFR_3))
}


do.eGFR = function(df) {
  
  ci <- df$.ci[1]
  ri <- df$.ri[1]
  col_annot <- ctx$cselect()
  age <- col_annot[[2]][ci+1]
  gender <- col_annot[[3]][ci+1]
  race = "Caucasian"
  
  crt_tbl <- df %>% group_by(.axisIndex) %>% 
    summarise(value = mean(.y)) %>% 
    filter(.axisIndex == 0)
  crt <- crt_tbl$value
  
  if (length(crt) == 0) crt = NaN
  
  cyt_tbl <- df %>% group_by(.axisIndex) %>% 
    summarise(value = mean(.y)) %>% 
    filter(.axisIndex == 1)
  cyt <- cyt_tbl$value
  
  if (length(cyt) == 0) cyt = NaN
  
  # convert units of crt
 if (as.logical(ctx$op.value('as uM/l'))) crt <- crt/ 88.4
  eGFR    <- calc_egfr(crt,cyt, gender, age, race)
  return(data.frame(
    .ri = df$.ri[1],
    .ci = df$.ci[1],
    eGFR_0 = eGFR[[1]],
    eGFR_1 = eGFR[[2]],
    eGFR_2 = eGFR[[3]],
    eGFR_3 = eGFR[[4]]
  ))
}

ctx = tercenCtx()

if (nrow(unique(ctx$select(c('.axisIndex')))) != 2)
  stop("Two layers are required, the first with the creatinine value and the second with the cystatin")

if (ncol(ctx$cselect()) < 3) stop ("Require age and gender as 2nd and 3rd factors on columns")
all_ages <- ctx$cselect()[[2]]
all_genders <- ctx$cselect()[[3]]

check_ages_within_range <- (all(unlist(lapply(all_ages, function(x) x %in% 1:100))))
check_gender_values <- (all(unlist(lapply(all_genders, function(x) x %in% c("Male", "Female")))))

if (!check_ages_in_range)  stop ("ages required as the second factor on columns")
if (!check_gender_values)  stop ("gender required as the third factor on columns")

df <- ctx %>%
  select(.ci, .ri, .y, .axisIndex) %>%
  group_by(.ci, .ri) %>%
  do(do.eGFR(.)) %>%
  ctx$addNamespace() %>%
  ctx$save()

