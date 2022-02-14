# Check for the necessary packages, and if they are missing, download it
library(tidyverse)

# 
# This script only identifies possible matches for a few mutational signatures:
# SBS1, SBS10a, SBS10b, SBS28, and ID1 & ID2
# It uses a nucleotide sequence as an input. Here we use  PTEN ORF from
# https://www.ncbi.nlm.nih.gov/nuccore/NM_000314.8

PTEN_ORF <- "ATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGA
CTTAGACTTGACCTATATTTATCCAAACATTATTGCTATGGGATTTCCTGCAGAAAGACTTGAAGGCGTA
TACAGGAACAATATTGATGATGTAGTAAGGTTTTTGGATTCAAAGCATAAAAACCATTACAAGATATACA
ATCTTTGTGCTGAAAGACATTATGACACCGCCAAATTTAATTGCAGAGTTGCACAATATCCTTTTGAAGA
CCATAACCCACCACAGCTAGAACTTATCAAACCCTTTTGTGAAGATCTTGACCAATGGCTAAGTGAAGAT
GACAATCATGTTGCAGCAATTCACTGTAAAGCTGGAAAGGGACGAACTGGTGTAATGATATGTGCATATT
TATTACATCGGGGCAAATTTTTAAAGGCACAAGAGGCCCTAGATTTCTATGGGGAAGTAAGGACCAGAGA
CAAAAAGGGAGTAACTATTCCCAGTCAGAGGCGCTATGTGTATTATTATAGCTACCTGTTAAAGAATCAT
CTGGATTATAGACCAGTGGCACTGTTGTTTCACAAGATGATGTTTGAAACTATTCCAATGTTCAGTGGCG
GAACTTGCAATCCTCAGTTTGTGGTCTGCCAGCTAAAGGTGAAGATATATTCCTCCAATTCAGGACCCAC
ACGACGGGAAGACAAGTTCATGTACTTTGAGTTCCCTCAGCCGTTACCTGTGTGTGGTGATATCAAAGTA
GAGTTCTTCCACAAACAGAACAAGATGCTAAAAAAGGACAAAATGTTTCACTTTTGGGTAAATACATTCT
TCATACCAGGACCAGAGGAAACCTCAGAAAAAGTAGAAAATGGAAGTCTATGTGATCAAGAAATCGATAG
CATTTGCAGTATAGAGCGTGCAGATAATGACAAGGAATATCTAGTACTTACTTTAACAAAAAATGATCTT
GACAAAGCAAATAAAGACAAAGCCAACCGATACTTTTCTCCAAATTTTAAGGTGAAGCTGTACTTCACAA
AAACAGTAGAGGAGCCGTCAAATCCAGAGGCTAGCAGTTCAACTTCTGTAACACCAGATGTTAGTGACAA
TGAACCTGATCATTATAGATATTCTGACACCACTGACTCTGATCCAGAGAATGAACCTTTTGATGAAGAT
CAGCATACACAAATTACAAAAGTCTGA"

PTEN_ORF <- str_replace_all(PTEN_ORF,"\\s","")

# To found mut. signatures you have to find locations with necessary 
# nucleodite context (e.g. "CG" for SBS1)
# https://cancer.sanger.ac.uk/signatures/sbs/sbs1/

PTEN_CG_coord <- str_locate_all(PTEN_ORF,"CG") %>% unlist()

# now, the coordinates for ID1 and ID2 signatures, combined 
# looking for polyA and polyT stretches
# https://cancer.sanger.ac.uk/signatures/id/id1/
# https://cancer.sanger.ac.uk/signatures/id/id2/

PTEN_IDT_coord <- str_locate_all(PTEN_ORF, c("A{4,}",  "T{4,}"))
# converting the list into a vector of possible coordinates 
PTEN_IDT_coord1 <- c()
for (i in 1:length(PTEN_IDT_coord)) {
  x = nrow(PTEN_IDT_coord[[i]])
  for (k in 1:x) {
    PTEN_IDT_coord1 <- c(PTEN_IDT_coord1, 
                         seq(PTEN_IDT_coord[[i]][k,][1], PTEN_IDT_coord[[i]][k,][2]))
  }
}

# Data set to be analyzed should be represented as .csv file, at least with the next columns:
# IDs, trinucleotideContext_SV, trinucleotideAlteration_SV, transcriptEffect_SV, Hugo_Symbol
# 
# age, MS_status (Microsatellite Stability), POLE_alt (Polymerase Epsilone alteration) and other variables
# can be used to extend and/or fine-tune the analysis

# The test data set is available as cancer_mutations_test_sample.csv
# It includes all the above data columns. Please note that this is a test set
# and represents randomized data, and thus is only suitable for training purposes.
# KRAS and TP53 in this test set are present for illustration purpose only

cancer_mutations <- read_csv("mut signatures PTEN CRC test intput.csv")

PTEN_variants <- filter(cancer_mutations, Hugo_Symbol == "PTEN")

PTEN_variants <- PTEN_variants %>% 
  mutate(CGfound=0, SBS1_subst = 0, SBS1 = 0, IDTfound = 0, IDT = 0) %>%
  as_tibble()

# Attach and detach as variant of working with column values

attach(PTEN_variants)

PTEN_variants$CGfound <-  str_extract(transcriptEffect_SV, "[:digit:]+") %in% PTEN_CG_coord + 0
PTEN_variants$SBS1_subst <- str_detect(transcriptEffect_SV, "C>T|G>A") + 0
PTEN_variants$SBS1 <- (PTEN_variants$CGfound * PTEN_variants$SBS1_subst +0) %>% replace_na(0) 

#  now, the ID1 and ID2 signatures, combined 
PTEN_variants$IDTfound <- str_extract(PTEN_variants$transcriptEffect_SV, "[:digit:]+") %in% PTEN_IDT_coord1 + 0
PTEN_variants$IDT <- PTEN_variants$IDTfound * (str_detect(PTEN_variants$transcriptEffect_SV, "insT$|insA$|delT$|delA$") + 0) %>% 
  replace_na(0)

detach(PTEN_variants)

# The other signatures studied here are described by the following nucleotide substitutions
# sbs10a - C>A | G>T (TCT | AGA)
# sbs10b - C>T | G>A (TCG, TCC | AGC, AGG)
# https://cancer.sanger.ac.uk/signatures/sbs/sbs10b/
# sbs28 TTT - A>C | T>G
# sbs28 ATT - A>C | T>G
# sbs28 CTT - A>C | T>G

PTEN_variants$SBS10a <- 0
PTEN_variants$SBS10a[which(str_detect(PTEN_variants$trinucleotideContext_SV, "[T][C][T]") & 
                             str_detect(PTEN_variants$trinucleotideAlteration_SV,"C>A"))] <- 1
PTEN_variants$SBS10b <- 0
PTEN_variants$SBS10b[which(str_detect(PTEN_variants$trinucleotideContext_SV, "[T][C][G]") & 
                             str_detect(PTEN_variants$trinucleotideAlteration_SV,"C>T"))] <- 1
PTEN_variants$SBS28 <- 0
PTEN_variants$SBS28[which(str_detect(PTEN_variants$trinucleotideContext_SV, "[T][T][T]") & 
                            str_detect(PTEN_variants$trinucleotideAlteration_SV,"T>G"))] <- 1
PTEN_variants$SBS28[which(str_detect(PTEN_variants$trinucleotideContext_SV, "[A][T][T]") & 
                            str_detect(PTEN_variants$trinucleotideAlteration_SV,"A>C"))] <- 1
PTEN_variants$SBS28[which(str_detect(PTEN_variants$trinucleotideContext_SV, "[C][T][T]") & 
                            str_detect(PTEN_variants$trinucleotideAlteration_SV,"C>T"))] <- 1

# Data aggregation to a convenient data format for further analysis.
mut_sig_count <- PTEN_variants %>% select(Sample_ID, SBS1, IDT, SBS10a, SBS10b, SBS28) %>% 
  group_by(Sample_ID) %>%
  summarise(SBS1_count = sum(SBS1), SBS1_presence = sign(SBS1_count),
            IDT_count = sum(IDT), IDT_presence = sign(IDT_count),
            SBS10a_count = sum(SBS10a), SBS10a_presence = sign(SBS10a_count), 
            SBS10b_count = sum(SBS10b), SBS10b_presence = sign(SBS10b_count),
            SBS28_count = sum(SBS28), SBS28_presence = sign(SBS28_count))

cancer_mutations_with_sig <- cancer_mutations %>%
  group_by(Sample_ID) %>%
  summarise(PTEN_count = sum(Hugo_Symbol == "PTEN", na.rm = T), PTEN_presence = sign(PTEN_count), 
            POLE_alt = sign(sum(POLE_alt)), MS_status = names(table(MS_status))[1], age = max(age)) %>%
  left_join(y = mut_sig_count, by = c("Sample_ID")) %>% 
  mutate(across(contains("count"), ~replace_na(.x, 0))) %>%
  mutate(across(contains("presence"), ~replace_na(.x, 0))) %>%
  mutate(PTEN_nonsig = PTEN_count-SBS1_count-IDT_count)

# To resolve the overlap between SBS1 and SBS10b signatures the POLE alteration status is required
cancer_mutations_with_sig$SBS1_presence[which(cancer_mutations_with_sig$SBS10b_presence == 1 & cancer_mutations_with_sig$SBS1_presence == 1
                                               & cancer_mutations_with_sig$POLE_alt == 1)] <- 0

# Rearrangement of data for plot the chart
cancer_mutations_with_sig <- cancer_mutations_with_sig %>%
  select(age, SBS1_count, IDT_count, SBS10a_count, SBS10b_count, SBS28_count, PTEN_nonsig, PTEN_count, MS_status) %>% 
  pivot_longer(-c(age,MS_status), names_to = "signature", values_to = "count") %>%
  mutate(signature = as.factor(signature))

# Inspecting and saving the output
cancer_mutations_with_sig %>%
  head()

write_csv(cancer_mutations_with_sig, "mut signatures PTEN CRC test output.csv")

# It can be read back for further analysis or plotting
samples_sig <- read_csv("mut signatures PTEN CRC test output.csv")

samples_sig %>%
  head()

# Chart plotting

# Here we examine the age trends of various 
# signatures . To do this, we first build a linear approximation (using method "lm"), with calculated 
# logistic regression coefficients, and then calculate statistical significance, in the box below.

ggplot(filter(cancer_mutations_with_sig, MS_status=="MT-L"), aes(x=age, y=count, color=signature))+
  geom_smooth(method = lm)
ggplot(filter(cancer_mutations_with_sig, MS_status=="MT-H"), aes(x=age, y=count, color=signature))+
  geom_smooth(method = lm)

# Chart statistics
# You can enter the appropriate output file name in the "file = " parameter

filter(cancer_mutations_with_sig, MS_status == "MT-L") %>%
  select(signature, count, age) %>%
  group_by(signature) %>%
  summarise(presence = sum(count), absence = sum(count == 0), 
            reg_coef = summary(lm(count ~ age))$coefficients[2],
            p_value = summary(lm(count ~ age))$coefficients[8]) %>%
  write_csv(file = "plot stats for mtl test.csv")

filter(cancer_mutations_with_sig, MS_status == "MT-H") %>%
  select(signature, count, age) %>%
  group_by(signature) %>%
  summarise(presence = sum(count), absence = sum(count == 0), 
            reg_coef = summary(lm(count ~ age))$coefficients[2],
            p_value = summary(lm(count ~ age))$coefficients[8]) %>%
  write_csv(file = "plot stats for mth test.csv")