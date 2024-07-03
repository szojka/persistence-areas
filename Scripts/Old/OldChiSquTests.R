# Old code from Paper1_stats objective 1 & 2 - chi squ has problems and doesn't properly answer my questins
# using df.scales and radar plot data

## Local scale

# overall
n_pops1 <- radat_plot_all[c(3,4),]
mat1 <- as.matrix(n_pops1)
chisq.test(mat1,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat1
# X-squared = 67.738, df = NA, p-value = 0.0004998

# plaere
n_pops1_p <-  radat_plot_plaere[c(3,4),]
mat1.p <- as.matrix(n_pops1_p)
chisq.test(mat1.p,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat1.p
# X-squared = 21.21, df = NA, p-value = 0.0004998

# miccal
n_pops1_m <- radat_plot_miccal[c(3,4),]
mat1.m <- as.matrix(n_pops1_m)
chisq.test(mat1.m,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat1.m
# X-squared = 10.241, df = NA, p-value = 0.01949

# vulmic
n_pops1_v <-  radat_plot_vulmic[c(3,4),]
mat1.v <- as.matrix(n_pops1_v)
chisq.test(mat1.v,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat1.v
# X-squared = 25.439, df = NA, p-value = 0.0004998

# brohor
n_pops1_b <-  radat_plot_brohor[c(3,4),]
mat1.b <- as.matrix(n_pops1_b)
chisq.test(mat1.b,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat1.b
# X-squared = 22.656, df = NA, p-value = 0.0004998

## Regional level
# overall
n_pops2 <- df.scales %>%
  dplyr::filter(scale%in%'regional') %>%
  dplyr::select(-scale,-species) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n) %>%
  distinct() 
n_pops2 <- pivot_wider(n_pops2, names_from = contingency, values_from = n_sum)
n_pops2$DL <- 0
n_pops2$SS_n <- 0
n_pops2 <- dplyr::select(n_pops2,-treatment)

mat2 <- as.matrix(n_pops2)
chisq.test(mat2,correct = F) # errror


# plaere
n_pops2_p <- df.scales %>%
  dplyr::filter(scale%in%'regional') %>%
  dplyr::filter(species%in%'plaere') %>%
  dplyr::select(-scale,-species) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n) %>%
  distinct()
n_pops2_p <- pivot_wider(n_pops2_p, names_from = contingency, values_from = n_sum)
n_pops2_p$DL <- 0
n_pops2_p$SS_n <- 0
n_pops2_p$ME <- 0
n_pops2_p <- dplyr::select(n_pops2_p,-treatment)

mat2.p <- as.matrix(n_pops2_p)
chisq.test(mat2.p,correct = T)


# miccal
n_pops2_m <- radat_plot_miccal[c(3,4),]

mat2.m <- as.matrix(n_pops1_m)
chisq.test(mat2.m,correct = F,simulate.p.value=TRUE)


# vulmic
n_pops2_v <-  radat_plot_vulmic[c(3,4),]

mat2.v <- as.matrix(n_pops2_v)
chisq.test(mat2.v,correct = F,simulate.p.value=TRUE)


# brohor
n_pops2_b <-  radat_plot_brohor[c(3,4),]

mat2.b <- as.matrix(n_pops2_b)
chisq.test(mat2.b,correct = F,simulate.p.value=TRUE)

## Global level

n_pops3 <- df.scales %>%
  dplyr::filter(scale%in%'global') %>%
  dplyr::select(-scale,-species) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n) %>%
  distinct() 
n_pops3 <- pivot_wider(n_pops3, names_from = contingency, values_from = n_sum)
n_pops3$DL <- 0
n_pops3$SS_n <- 0
n_pops3 <- dplyr::select(n_pops3,-treatment)

mat3 <- as.matrix(n_pops3)
chisq.test(mat3,correct = T) # error

# Obj 2 ####
# How does the relationship between fundamental and realized niches change with spatial scale?

# Test 2)  does the frequency of community statuses differ between spatial scales? (Fig 1.1)
# 4 factors (contingency outcomes) x 3 factors (scales) x 2 factors (treatment A or B)

# separate data to compare n of contingency groups between scales within A and within B

n_pops3_A <- dplyr::filter(df.scales, treatment%in%"A")
n_pops3_B <- dplyr::filter(df.scales, treatment%in%"B")

# chi squared test of independence
# Null hyp = The distribution of the outcome is independent of the groups

# first remove species

# treatment A
n_pops_AA <- dplyr::select(n_pops3_A, -species, -treatment)
n_pops_AA <- n_pops_AA %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n) %>%
  distinct()
n_pops_AA <- pivot_wider(n_pops_AA, names_from = "scale", values_from = "n_sum")
n_pops_AA[is.na(n_pops_AA)] <- 0

# observed data frame = n_pops_AA
n_pops_AA <- n_pops_AA %>%
  dplyr::ungroup() %>%
  dplyr::select(-contingency)
mat_A <- as.matrix(n_pops_AA) # contingency order is row1 =  SS_y, row2 = ME, row3 = DL, row4 = SS_n

chisq.test(mat_A,correct = F,simulate.p.value=TRUE)
# data:  mat_A
# X-squared = 128.31, df = NA, p-value = 0.0004998

# df = (r-1)*(c-1) = 2*3 = 6
1-pchisq(128.31, 6)  # 0?
# we use this command to find p-val because we want the porbability given by the upper tail of the chi-squared distribution

# treatment B
n_pops_BB <- dplyr::select(n_pops3_B, -species, -treatment)
n_pops_BB <- n_pops_BB %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n) %>%
  distinct()
n_pops_BB <- pivot_wider(n_pops_BB, names_from = "scale", values_from = "n_sum")
n_pops_BB[is.na(n_pops_BB)] <- 0

# observed data frame = n_pops_AA
n_pops_BB <- n_pops_BB %>%
  dplyr::ungroup() %>%
  dplyr::select(-contingency)
mat_B <- as.matrix(n_pops_BB) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_B,correct = F,simulate.p.value=TRUE)
# data:  mat_B
# X-squared = 130.91, df = NA, p-value = 0.0004998

1-pchisq(130.91, 6) #0?

#### Plaere ####
n_pops_pa <- dplyr::filter(df.scales, treatment%in%"A" & species%in%"plaere")
n_pops_pb <- dplyr::filter(df.scales, treatment%in%"B" & species%in%"plaere")

# plaere abiotic
n_pops_pa <- n_pops_pa %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_pa <- pivot_wider(n_pops_pa, names_from = "scale", values_from = "n_sum")
# n_pops_pa[4,1] <- "ME" #removed in an effort to fix zero marginal problem
n_pops_pa[is.na(n_pops_pa)] <- 0
# observed data frame = n_pops_AA
n_pops_pa <- n_pops_pa %>%
  dplyr::ungroup() %>%
  dplyr::select(-contingency)
mat_pa <- as.matrix(n_pops_pa) # scale order is row1 =  SS_y, row2 = DL, row3 = SS_n, row4 = ME

chisq.test(mat_pa,correct = T,simulate.p.value=T) #error because E-O / E where E = 0, and dividing by 0 causes problems
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_pa
# X-squared = 30.655, df = NA, p-value = 0.0004998

# plaere biotic
n_pops_pb <- n_pops_pb %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_pb <- pivot_wider(n_pops_pb, names_from = "contingency", values_from = "n_sum")
n_pops_pb[is.na(n_pops_pb)] <- 0
# observed data frame = n_pops_AA
n_pops_pb <- n_pops_pb %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_pb <- as.matrix(n_pops_pb) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_pb,correct = F,simulate.p.value=T) # error
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_pb
# X-squared = 326.67, df = NA, p-value = 0.0004998

#### Miccal ####
n_pops_ma <- dplyr::filter(df.scales, treatment%in%"A" & species%in%"miccal")
n_pops_mb <- dplyr::filter(df.scales, treatment%in%"B" & species%in%"miccal")

# miccal abiotic
n_pops_ma <- n_pops_ma %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_ma <- pivot_wider(n_pops_ma, names_from = "contingency", values_from = "n_sum")
n_pops_ma[is.na(n_pops_ma)] <- 0
# observed data frame = n_pops_AA
n_pops_ma <- n_pops_ma %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_ma <- as.matrix(n_pops_ma) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_ma, correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_ma
# X-squared = 52.951, df = NA, p-value = 0.0004998

# miccal biotic
n_pops_mb <- n_pops_mb %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_mb <- pivot_wider(n_pops_mb, names_from = "contingency", values_from = "n_sum")
n_pops_mb$SS_y <- 0
n_pops_mb[is.na(n_pops_mb)] <- 0
# observed data frame = n_pops_AA
n_pops_mb <- n_pops_mb %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_mb <- as.matrix(n_pops_mb) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_mb,correct = F,simulate.p.value=T) #Error

#### Vulmic ####
n_pops_va <- dplyr::filter(df.scales, treatment%in%"A" & species%in%"vulmic")
n_pops_vb <- dplyr::filter(df.scales, treatment%in%"B" & species%in%"vulmic")

#vulmic abiotic
n_pops_va <- n_pops_va %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_va <- pivot_wider(n_pops_va, names_from = "contingency", values_from = "n_sum")
n_pops_va[is.na(n_pops_va)] <- 0
# observed data frame = n_pops_AA
n_pops_va <- n_pops_va %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_va <- as.matrix(n_pops_va) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_va,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_va
# X-squared = 40.142, df = NA, p-value = 0.0004998

# vulmic biotic
n_pops_vb <- n_pops_vb %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_vb <- pivot_wider(n_pops_vb, names_from = "contingency", values_from = "n_sum")
n_pops_vb[is.na(n_pops_vb)] <- 0
# observed data frame = n_pops_AA
n_pops_vb <- n_pops_vb %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_vb <- as.matrix(n_pops_vb) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_vb,correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_vb
# X-squared = 63.709, df = NA, p-value = 0.0004998

#### Brohor ####
n_pops_ba <- dplyr::filter(df.scales, treatment%in%"A" & species%in%"brohor")
n_pops_bb <- dplyr::filter(df.scales, treatment%in%"B" & species%in%"brohor")

# brohor abiotic
n_pops_ba <- n_pops_ba %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_ba <- pivot_wider(n_pops_ba, names_from = "contingency", values_from = "n_sum")
n_pops_ba[is.na(n_pops_ba)] <- 0
# observed data frame = n_pops_AA
n_pops_ba <- n_pops_ba %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_ba <- as.matrix(n_pops_ba) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_ba, correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_ba
# X-squared = 59.31, df = NA, p-value = 0.0004998

# brohor biotic
n_pops_bb <- n_pops_bb %>%
  dplyr::group_by(scale, contingency) %>%
  dplyr::mutate(n_sum = sum(n)) %>%
  dplyr::select(-n, -treatment,-species) %>%
  distinct()
n_pops_bb <- pivot_wider(n_pops_bb, names_from = "contingency", values_from = "n_sum")
n_pops_bb[is.na(n_pops_bb)] <- 0
# observed data frame = n_pops_AA
n_pops_bb <- n_pops_bb %>%
  dplyr::ungroup() %>%
  dplyr::select(-scale)
mat_bb <- as.matrix(n_pops_bb) # contingency order is row1 =  ME, row2 = SS_Y, row3 = DL, row4 = SS_n

chisq.test(mat_bb, correct = F,simulate.p.value=TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  mat_bb
# X-squared = 31.592, df = NA, p-value = 0.0004998
