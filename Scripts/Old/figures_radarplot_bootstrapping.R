# Radarplot bootstrapping variability figures

load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_ga.RData")
load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_gb.RData")
load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_la.RData")
load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_lb.RData")
load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_sa.RData")
load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/avg_sb.RData")

# Local Abiotic

ggplot(data = avg_la, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) + # Make this work ####
theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")

# Local Biotic


ggplot(data = avg_lb, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")

# Grid Abiotic


ggplot(data = avg_ga, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")

# Grid Biotic


ggplot(data = avg_gb, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")

# Site Abiotic

ggplot(data = avg_sa, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")

# Site Biotic

ggplot(data = avg_sb, aes(x = contingency, y = true.prop)) + 
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_classic() +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
  ylim(0,1) +
  labs(x="Community type", y = "Proportion of total communities")


