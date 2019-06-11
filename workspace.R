library(tercen)
library(tidyverse)
#http://127.0.0.1:5402/#ds/d30382066b71e6e7995cee981c001603/5-6
# http://127.0.0.1:5402/#ds/ee0cbc176ba159782635097a88009c5b/8-4
options("tercen.workflowId"= "d30382066b71e6e7995cee981c001603")
options("tercen.stepId"= "5-6")

ctx = tercenCtx()

core <- ctx %>% select() 
rows <- ctx %>% rselect()
cols <- ctx %>% cselect()

input <- list(
  alternative = "two.sided",
  mu = 0,
  paired = F,
  var.equal = F,
  conf.level= 0.95,
  abase=10,
  do.ttest = T,
NC.factor = "js0.compound",
  NC.annotation = "DMSO",
  value=".y",
  Ignore.negatives =F
)

ttest <- function(pop1, pop2, input){
  alternative <- input$alternative
  mu <- input$mu
  paired <- input$paired
  var.equal <- input$var.equal
  conf.level <- input$conf.level
  
  ds <- try(t.test(pop1, pop2, alternative, mu, paired, var.equal, conf.level))
  p  <-  if(!inherits(ds, "try-error")) {ds$p.value} else{1}
  return(p)
}

response1 <- function(df,pop1, input){
  
  abase <- input$abase
  pop2 <- df[[input$value]]
  df$binding <- round(mean(pop2), digits=0)
  df$n <- length(pop2)
  df$sd <- round(sd(pop2), digits=0)
  df$cv <- round(with(df, 100*sd/binding),digits=0)
  df$sem <- round(sd(pop2)/sqrt(length(pop2)-1), digits=0)
  dottest <- input$do.ttest 
  if(dottest){
    df$p <- ttest(pop1, pop2, input)
  }
  
  df$ratio <- round(mean(pop2)/mean(pop1), digits=2)
  df$MI <- round(log((mean(pop2)/mean(pop1)), abase), digits=2)
  
  
  return(df)
}

response2 <- function(coldf,pop1, input){
  abase <- input$ttest
  pop2 <- coldf[[input$value]]
  coldf$binding <- round(mean(pop2), digits=0)
  coldf$n <- length(pop2)
  coldf$sd <- round(sd(pop2), digits=0)
  coldf$cv <- round(with(coldf, 100*sd/binding),digits=0)
  coldf$sem <- round(sd(pop2)/sqrt(length(pop2)-1), digits=0)
  dottest <- input$do.ttest
  if(dottest=="Y"){
    coldf$p <- ttest(pop1, pop2, input)
  }
  
  coldf$ratio <- 1
  coldf$MI <- 0
  
  return(coldf)
}

responsefunction <- function(df, input){  
  NC.factor <- input$NC.factor
  NC.annot <-input$NC.annotation
  
  
  
  pop1 <- subset(df, df[[NC.factor]]==NC.annot)[[input$value]]
  
  # 
  
  if(input$Ignore.negatives){
    df <-df %>% plyr::ddply(NC.factor, function(x){response2(x, pop1, input)} )
  } else {
    df <-df %>% plyr::ddply(NC.factor, function(x){response1(x, pop1, input)} )
  }

  
  
  
  return(df)
}



core2 <- core

core2$.y <- with(core2, replace(.y, .y<1, 1))

q <- core2 %>%  plyr::ddply(~.ri, function(x){responsefunction(x, input)})



