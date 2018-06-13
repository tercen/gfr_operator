library(tercen)
library(dplyr)

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
  crt <- crt/ 88.4
  
  a <- ifelse (gender=="Female",-0.248,-0.207) 
  k <- ifelse (gender=="Female",0.7,0.9) 
  
  
  eGFR <-
    135 * (min(crt/k, 1))^a *
    (max(crt/k, 1))^-0.601 *
    (min(cyt/0.8, 1))^-0.375 *
    (max(cyt/0.8, 1))^-0.711 *
    0.995^age
  
  eGFR <- ifelse (gender=="Female",0.969 * eGFR,eGFR)
  eGFR <- ifelse (race=="Black",1.08 * eGFR,eGFR)
  
  return(data.frame(
    .ri = df$.ri[1],
    .ci = df$.ci[1],
    eGFR = eGFR
  ))
}

ctx = tercenCtx()

if (nrow(unique(ctx$select(c('.axisIndex')))) != 2)
  stop("Two layers are required, the first with the creatinine value and the second with the cystatin")

df <- ctx %>%
  select(.ci, .ri, .y, .axisIndex) %>%
  group_by(.ci, .ri) %>%
  do(do.eGFR(.)) %>%
  ctx$addNamespace() %>%
  ctx$save()
