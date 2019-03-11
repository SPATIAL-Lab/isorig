plot.QA <- function(obj){
  # accuracy by prob
  xx <- seq(0.01, 0.99, 0.01)
  mean <- data.frame(xx, mean = apply(obj$prption_byProb, 2, mean), type = "QA")
  p <- ggplot(data=mean, aes(x=xx, y=mean,  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Relative probability interval") + 
    ylab("Number of validation stations included") + labs(title = "Accuracy by top accumulated probability")
  print(p)
  
  # accuracy by area
  mean <- data.frame(xx, mean = apply(obj$prption_byArea, 2, mean), type = "QA")
  p <- ggplot(data=mean, aes(x=xx, y=mean,  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Relative probability interval") + 
    ylab("Number of validation stations included") + labs(title= "Accuracy by top percentage of area")
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
  precision <- data.frame(xx,  mean = mean, type = 'QA')
  p <- ggplot(data=precision, aes(x=xx, y=mean,  color = type)) + 
    geom_line(size=0.5) + 
    xlab("Relative probability interval") + 
    ylab("Proportion of the total surface") + labs(title = "Precision")
  print(p)
  
  # plot pd
  pd <- data.frame(y=as.numeric(obj$pd_bird_val), type = "QA")
  p <- ggplot(pd, aes(x=type, y=y, fill=type)) + geom_boxplot() +  guides(fill=FALSE)
  p + geom_hline(yintercept=obj$random_prob_density)
}

