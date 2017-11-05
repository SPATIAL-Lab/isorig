# bootstrap obs data
# resamaple (bootstrap with raplacement) calibration (known origin) data

# 10/20/2017: please ignore this function for now. It is used for bootstraping method in the future

resample_obs <- function(obs)
{
  t <- t(obs)
  t <- as.data.frame(t)
  t1 <- sample(t, replace = TRUE)
  t <- t(t1)
  resam <- cbind(as.numeric(t[,1]),as.numeric(t[,2]),as.numeric(t[,3]))

  return(resam)
}

