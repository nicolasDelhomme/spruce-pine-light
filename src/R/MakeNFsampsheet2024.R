library(tidyverse)

# check sampsheet
pine <- read_csv("doc/pine.csv",show_col_types = F)
pineinf <- read_csv("doc/Pine-light-project-sample-info.csv",show_col_types = F)

spruce <- read_csv("doc/spruce.csv",show_col_types = F)
spruceinf <- read_csv("doc/Spruce-light-project-sample-info.csv",show_col_types = F)

# make sampsheet
pinefiles <- data.frame(path = list.files(path = "/mnt/picea/projects/spruce/facility/rgarcia-pine/raw",full.names = T)) %>%
  mutate(FileName = basename(path),
         fq = ifelse(grepl("_1.fastq.gz$",path),1,2)) %>%
  left_join(read_csv("doc/Pine-light-project-sample-info.csv",show_col_types = F), by = "FileName") %>%
  select(SampleName,path,fq)
sprucefiles <- data.frame(path = list.files(path = "/mnt/picea/projects/spruce/facility/rgarcia-spruce/raw",full.names = T)) %>%
  mutate(FileName = basename(path),
         fq = ifelse(grepl("_1.fastq.gz$",path),1,2)) %>%
  left_join(read_csv("doc/Spruce-light-project-sample-info.csv",show_col_types = F), by = "FileName") %>%
  select(SampleName,path,fq)

# separate fq1 and fq2
pinesamp <- pinefiles %>%
  group_split(fq) %>%
  bind_cols()
stopifnot(pinesamp$SampleName...1 == pinesamp$SampleName...4)
sprucesamp <- sprucefiles %>%
  group_split(fq) %>%
  bind_cols()
stopifnot(sprucesamp$SampleName...1 == sprucesamp$SampleName...4)

# save samsheet
pinesamp %>%
  select(SampleName...1,path...2,path...5) %>%
  setNames(c("sample","fastq_1","fastq_2")) %>%
  mutate(strandedness = "auto") %>%
  write_csv("doc/nf_samplesheet_pine.csv")
sprucesamp %>%
  select(SampleName...1,path...2,path...5) %>%
  setNames(c("sample","fastq_1","fastq_2")) %>%
  mutate(strandedness = "auto") %>%
  write_csv("doc/nf_samplesheet_spruce.csv")
