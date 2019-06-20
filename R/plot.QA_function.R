plot.QA <- function(obj){
  # accuracy by prob
  xx <- seq(0.01, 0.99, 0.01)
  mean <- data.frame(xx, mean = apply(obj$prption_byProb, 2, mean), type = "QA")
  p <- ggplot(data=mean, aes(x=xx, y=mean/ncol(obj$val_stations),  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Relative probability interval") + 
    ylab("Proportion of validation stations included") + labs(title = "Assignment specificity")
  print(p)
  
  # accuracy by area
  mean <- data.frame(xx, mean = apply(obj$prption_byArea, 2, mean), type = "QA")
  p <- ggplot(data=mean, aes(x=xx, y=mean/ncol(obj$val_stations),  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Cumulative probability threshold") + 
    ylab("Proportion of validation stations included") + labs(title= "Assignment accuracy")
  print(p)
  
  # plot precision
  precision <- NULL
  for (i in 1:length(obj$precision)){
    precision <- append(precision, apply(obj$precision[[i]],1, median))
  }
  mean <- NULL
  for(i in 1:99){
    mean <- append(mean, mean(precision[0:(length(obj$precision)-1)*99+i]))
  }
  precision <- data.frame(xx,  mean = 1-mean, type = 'QA')
  p <- ggplot(data=precision, aes(x=xx, y=mean,  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Cumulative probability threshold") + 
    ylab("Proportion of area excluded") + labs(title = "Spatial precision")
  print(p)
  
  # plot pd
  pd <- data.frame(y=as.numeric(obj$pd_bird_val), type = "QA")
  p <- ggplot(pd, aes(x=type, y=y/obj$random_prob_density, fill=type)) + 
    geom_boxplot() +  guides(fill=FALSE) + xlab("") +
    ylab("Odds ratio (known origin:random)")
  print(p)
}

