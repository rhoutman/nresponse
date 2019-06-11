library(tercen)
library(tidyverse)

ctx <- tercenCtx() 

input <- list(
  alternative = "two.sided",
  mu = 0,
  paired = F,
  var.equal = F,
  conf.level = 0.95,
  abase = 10,
  do.ttest = T,
  NC.factor = ctx$colors[[1]],
  NC.annotation = "DMSO",
  value = ".y",
  Ignore.negatives =  as.logical(ctx$op.value('ignore negatives'))
)

ttest <- function(pop1, pop2, input) {
  alternative <- input$alternative
  mu <- input$mu
  paired <- input$paired
  var.equal <- input$var.equal
  conf.level <- input$conf.level

  ds <- try(t.test(pop1, pop2, alternative, mu, paired, var.equal, conf.level))
  p <- if (!inherits(ds, "try-error")) {
    ds$p.value
  } else {
    1
  }
  return(p)
}

dostats <- function(df, pop1, input) {
  abase <- input$abase
  pop2 <- df[[input$value]]
  df$binding <- round(mean(pop2), digits = 0)
  df$n <- length(pop2)
  df$sd <- round(sd(pop2), digits = 0)
  df$cv <- round(with(df, 100 * sd / binding), digits = 0)
  df$sem <- round(sd(pop2) / sqrt(length(pop2) - 1), digits = 0)
  dottest <- input$do.ttest
  if (dottest) {
    df$p <- ttest(pop1, pop2, input)
  }

  if (input$Ignore.negatives == F) {
    if (ctx %>% select() %>% pull(.y) %>% min() <= 0)
    stop("First, remove zeros and negative values from your measurement!")
    df$ratio <- round(mean(pop2) / mean(pop1), digits = 2)
    df$MI <- round(log((mean(pop2) / mean(pop1)), abase), digits = 2)
  } else {
    df$ratio <- 1
    df$MI <- 0
  }
  return(df)
}

responsefunction <- function(df, input) {
  NC.factor <- input$NC.factor
  NC.annot <- input$NC.annotation

  pop1 <- subset(df, df[[NC.factor]] == NC.annot)[[input$value]]
  df <- df %>% plyr::ddply(NC.factor, function(x) {
    dostats(x, pop1, input)
   })
  return(df)
}

ctx %>%
  select()  %>% 
  plyr::ddply(~.ri, function(x) {
    responsefunction(x, input)
  }) %>% 
  ctx$addNamespace() %>%
  ctx$save()
  

