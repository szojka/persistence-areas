#-----------------------------------------------------------------------.
# Old Response Variable Calculation ####
# response is called 'sp' = average persistence over spatial scales (can be bw 0-4)

# set up dataframes
df <- dplyr::select(dat, names, run, scale, treatment, persist) %>% 
  # including names because otherwise n goes wonky
  distinct()
# scale has to be ascending
class(df$scale)
df <- arrange(df, run,treatment) 
df <- df[order(as.numeric(as.character(df$scale))), ]

n_sum <- df %>%
  dplyr::group_by(run, scale, treatment) %>%
  dplyr::distinct() %>%
  dplyr::add_count() %>%
  dplyr::select(run, scale, treatment, n)%>%
  dplyr::distinct() # create this so that frequencies are non repeating for each persistence   measure as the latter is on the plot level


nrows <- 15*2*nrun
mat <- matrix(data = NA, nrow = nrows, ncol = 4) 
colnames(mat) <- c("run", "scale", "treatment", "sp")

k <- 1 # counter for vector for matrix 'mat' rows

for (r in unique(df$run)) {
  temp <- subset(df, run == r)
  # because scale levels often change between runs
  
  n_temp <- subset(n_sum, run == r)
  # so that averages reflect the correct sum within each run
  
  for (i in unique(df$scale)) {
    temp2 <- subset(temp, scale == i)
    n_temp2 <- subset(n_temp, scale == i)
    
    for (t in levels(df$treatment)) {
      # calculates for A first then B (see arrange in df) as keeping an order is important for later
      
      sp <- # scale must be numeric for this to work
        sum(temp2$persist[temp2$scale <= i & temp2$treatment == t]) /
        sum(n_temp2$n[n_temp2$scale <= i & n_temp2$treatment == t])
      # dived by the sum of n for the given scale & treatment to get an appropriate mean at each group
      
      # fill the matrix with current values
      mat[k,1] <- r
      mat[k,2] <- i
      mat[k,3] <- t
      mat[k,4] <- sp
      
      k <- k + 1
    }
  }
} 
# may be experiencing an issue where scale does not ascent 1->15 in order

# fills to exactly the data that exists
mat.df <- as.data.frame(mat)
mat.df <- dplyr::filter(mat.df, !sp %in% "NaN")
str(mat.df)

mat.df$treatment <- as.factor(as.character(mat.df$treatment))
mat.df$scale <- as.numeric(as.character(mat.df$scale))
mat.df$sp <- as.numeric(as.character(mat.df$sp))

# combine mat.df & dat
dat <- left_join(dat, mat.df, by = c("run","scale","treatment"))

# rescale sp where 0 = 0% and 4=100%
dat$sp <- rescale(dat$sp, from = c(0,4), to = c(0,100)) 

rm(mat.df, n_sum, temp, n_temp, temp2, n_temp2)