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
  Ignore.negatives = as.logical(ctx$op.value('ignore negatives'))
)

# extract the data
ctxcore <- ctx %>%  
  select() 

# reformat column names for odd characters
colnames(ctxcore) <- colnames(ctxcore) %>% make.names()

# create a column for negative control annotation
if(ctx$labels %>% length() <1){
  stop("Drag the annotation factor for the negative control (NC) to 'Labels' first! \n Next, provide the annotion for the NC in the operator properties.")
}else{
  ctxcore$NCfactor = ctxcore[[ctx$labels[[1]] %>% make.names()]]
}

# create a column to split the data (if present)
if(ctx$colors %>% length() >0) {
  asplit <- ctx$colors %>% as.character() 
  ctxcore <- ctxcore %>% 
    unite("split" ,asplit, remove=F)
} else{ctxcore$split = ""} 

# filter for needed data
ctxcore <- ctxcore %>% 
  select(.ri, 
         .ci, 
         NCfactor,
         split,
         .y)

# if indicated in properties, convert negatives values to small positive
if(as.logical(ctx$op.value('convert negatives'))){
  ctxcore$.y <- with(ctxcore, replace(.y, .y<1, 1))
}


# ttest function  
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

# perform some basic calculations
dostats <- function(df, pop1, input) {
  abase <- input$abase
  pop2 <- df$.y
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
    if (ctxcore %>% pull(.y) %>% min() <= 0)
      stop("First, remove zeros and negative values from your measurement!")
    df$ratio <- round(mean(pop2) / mean(pop1), digits = 2)
    df$MI <- round(log((mean(pop2) / mean(pop1)), abase), digits = 2)
  } else {
    df$ratio <- 1
    df$MI <- 0
  }
  return(df)
}

# select the negative control and apply calculation
responsefunction <- function(df, input) {
  NC.annot <-  ctx$op.value('negative control') %>% as.character()
  if(!NC.annot %in% ctxcore$NCfactor){
    stop("This negative control is not found in the data. Please correct in the operator properties.")
    
  }
  pop1 <- subset(df, NCfactor == NC.annot)$.y
  
  df <- df %>% plyr::ddply(~NCfactor, function(x) {
    dostats(x, pop1, input)
  })
  return(df)
}


# perform actial calculations for the result
doper <- c(".ri", "split")

ctxcore <- ctxcore  %>% 
  plyr::ddply(~split+.ri, function(x){
    responsefunction(x, input) 
  }) %>%
  group_by(.ri, .ci) %>% 
  summarize(
    binding = mean(binding),
    sem = mean(sem),
    n = mean(n),
    sd = mean(sd),
    cv = mean(cv),
    ratio = mean(ratio),
    MI = mean(MI),
    p = mean(p)
  ) 


ctxcore$pcorr <- p.adjust(ctxcore$p, method='fdr', n=nrow(ctxcore))

ctxcore %>% 
  ctx$addNamespace() %>%
  ctx$save()
  

