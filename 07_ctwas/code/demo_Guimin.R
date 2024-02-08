library(tidyverse)
library(data.table)

tdf <- fread("../input/AD_Schwartz_GWASX.txt.gz")
tdf <- tdf %>% mutate(man_z = effect_size / standard_error)
tdf <- tdf %>% mutate(zscore_diff = abs(man_z - zscore)) %>% arrange(desc(man_z))
tdf %>% filter(zscore_diff > .01)
