###################### ACP proteomic Analysis ###


library(tidyverse)
library(seqinr)
library(ggpubr)


### Reading files
ACP_proteomics_MCOT_UniPep <- read_csv("path_to_csv_file")
ACP_proteomics_MCOT_TSC <- read_csv("path_to_csv_file")

### select unique peptides (unipep >=2) and total spectral count (TSC >=2) from adults and nymphs, removing decoys, joining files
ACP_adult_proteomics_MCOT_unipep2_TSC2 <- ACP_proteomics_MCOT_UniPep %>%
  select(Identified_Proteins:BBMV_Adult3)%>%
  filter(BBMV_Adult1 >=2 |  BBMV_Adult2 >=2 |  BBMV_Adult3 >=2)%>%
  filter(!str_detect(Identified_Proteins,"DECOY"))%>%
  select(Accession_Number)%>%
  left_join(ACP_proteomics_MCOT_TSC, by= "Accession_Number")%>%
  select(Accession_Number:BBMV_Adult3)%>%
  filter(BBMV_Adult1 >=2 |  BBMV_Adult2 >=2 |  BBMV_Adult3 >=2)%>%
  filter(!str_detect(Identified_Proteins,"DECOY"))
write_csv(ACP_adult_proteomics_MCOT_unipep2_TSC2 , "path_to_csv_file.csv")


ACP_nymph_proteomics_MCOT_unipep2_TSC2 <- ACP_proteomics_MCOT_UniPep %>%
  select(!BBMV_Adult1:BBMV_Adult3)%>%
  filter(BBMV_Nymph1 >=2 |  BBMV_Nymph2 >=2 |  BBMV_Nymph3 >=2)%>%
  filter(!str_detect(Identified_Proteins,"DECOY"))%>%
  select(Accession_Number)%>%
  left_join(ACP_proteomics_MCOT_TSC, by= "Accession_Number")%>%
  select(!BBMV_Adult1:BBMV_Adult3)%>%
  filter(BBMV_Nymph1 >=2 |  BBMV_Nymph2 >=2 |  BBMV_Nymph3 >=2)%>%
  filter(!str_detect(Identified_Proteins,"DECOY"))
write_csv(ACP_nymph_proteomics_MCOT_unipep2_TSC2 , "path_to_csv_file.csv")


### Removing duplicates accession number from the data
#adult
ACP_adult_proteomics_MCOT_unipep2_TSC2 <- read_csv("path_to_csv_file.csv")
ACP_adult_proteomics_unique_MCOT_unipep2_TSC2 <- ACP_adult_proteomics_MCOT_unipep2_TSC2 %>% distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_adult_proteomics_unique_MCOT_unipep2_TSC2,"path_to_csv_file.csv")
ACP_adult_proteomics_unique_MCOT_unipep2_TSC2_AccNum_only <- ACP_adult_proteomics_unique_MCOT_unipep2_TSC2 %>%
  select(Accession_Number)
write_csv(ACP_adult_proteomics_unique_MCOT_unipep2_TSC2_AccNum_only,"path_to_csv_file.csv")

#nymph
ACP_nymph_proteomics_MCOT_unipep2_TSC2 <- read_csv("path_to_csv_file.csv")
ACP_nymph_proteomics_unique_MCOT_unipep2_TSC2 <- ACP_nymph_proteomics_MCOT_unipep2_TSC2 %>% distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_nymph_proteomics_unique_MCOT_unipep2_TSC2,"path_to_csv_file.csv")
ACP_nymph_proteomics_unique_MCOT_unipep2_TSC2_AccNum_only <- ACP_nymph_proteomics_unique_MCOT_unipep2_TSC2 %>%
  select(Accession_Number)
write_csv(ACP_nymph_proteomics_unique_MCOT_unipep2_TSC2_AccNum_only,"path_to_csv_file.csv")



### join accession number from nymphs and adults and proceed with Cello2GO, DeepLoc and BUSCA analysis
ACP_adult_nymph_unique_MCOT_unipep2_TSC2_Acces_Number <- full_join(ACP_nymph_unipep2_TSC2_AccesNum,ACP_adult_unipep2_TSC2_AccesNum, by= "Accession_Number")
write_csv(ACP_adult_nymph_unique_MCOT_unipep2_TSC2_Acces_Number,"path_to_csv_file.csv")

### Extraction of aa sequences from  local database of the accession numbers identified by Mass-spectrometer analysis for nymphs and adults, and save in a fasta file using seqinr package
fastafile<- read.fasta(file = "path_to_fasta_file.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE) 
subsetlist<-read.csv("path_to_csv_file.csv",header = TRUE)
ACP_adult_nymph_unipep2_TSC2_Acces_Number_aaseq <- fastafile[c(which(names(fastafile) %in% subsetlist$Accession_Number))] 
write.fasta(ACP_adult_nymph_unipep2_TSC2_Acces_Number_aaseq, names = names(ACP_adult_nymph_unipep2_TSC2_Acces_Number_aaseq), nbchar = 80, file.out = "path_to_fasta_file.fasta")



### combine split files of BUSCA analysis
ACP_adult_nymph_BUSCA_altenative <- read_csv("path_to_csv_file.csv")

ACP_adult_nymph_BUSCA_final <- dir("path_to_csv_folder", 
                                   full.names = T) %>% map_df(read_csv) 
write_csv(ACP_adult_nymph_BUSCA_final,"path_to_csv_file.csv")
ACP_adult_nymph_BUSCA_final <- read_csv("path_to_csv_file.csv")

###select plasma membrane, GPI-anchored proteins and alternative localization prediction by BUSCA
ACP_adult_nymph_BUSCA_PM_GPI_alternative <- ACP_adult_nymph_BUSCA_final_alternative %>%
  filter(GOterms == "plasma membrane" | GOterms == "anchored component of plasma membrane" | Alternative_Localization == "plasma membrane")
write_csv(ACP_adult_nymph_BUSCA_PM_GPI_alternative,"path_to_csv_file.csv")
ACP_adult_nymph_BUSCA_PM_GPI_alternative_AccNum_Only <- ACP_adult_nymph_BUSCA_PM_GPI_alternative %>%
  select(Accession_Number)
write_csv(ACP_adult_nymph_BUSCA_PM_GPI_alternative_AccNum_Only,"path_to_csv_file.csv")


### combine split files of DeepLOC analysis
ACP_adult_nymph_deeploc_final <- dir("path_to_csv_folder", 
                                     full.names = T) %>% map_df(read_csv) 
write_csv(ACP_adult_nymph_deeploc_final,"path_to_csv_file.csv")

# select plasma membrane proteins from DeepLOC analysis 
ACP_adult_nymph_deeploc_final <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_deeploc_TM <- ACP_adult_nymph_deeploc_final %>%
  filter(Localization == "Cell membrane")
write_csv(ACP_adult_nymph_deeploc_TM,"path_to_csv_file.csv")

ACP_adult_nymph_deeploc_TM_AccNum_only <- ACP_adult_nymph_deeploc_TM %>%
  select(Accession_Number)
write_csv(ACP_adult_nymph_deeploc_TM_AccNum_only,"path_to_csv_file.csv")


### join cello2go, Deeploc and BUSCA accession numbers predicted to be plasma membrane
ACP_adult_nymph_cello2go <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_deeploc <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_BUSCA <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_cello2go_deeploc_BUSCA_AccNum_only <- full_join(ACP_adult_nymph_cello2go,ACP_adult_nymph_deeploc,ACP_adult_nymph_BUSCA, by = "Accession_Number") %>%
  distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_adult_nymph_cello2go_deeploc_BUSCA_AccNum_only,"path_to_csv_file.csv")

# Extract protein sequences of the accession numbers predicted to be localized to Plasma membrane using seqinr package 
fastafile<- read.fasta(file = "path_to_fasta_file.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE) 
subsetlist<-read.csv("path_to_csv_file.csv",header = TRUE)
ACP_adult_nymph_cello2go_deeploc_BUSCA_MCOT_aaseq <- fastafile[c(which(names(fastafile) %in% subsetlist$Accession_Number))] 
write.fasta(ACP_adult_nymph_cello2go_deeploc_BUSCA_MCOT_aaseq, names = names(ACP_adult_nymph_cello2go_deeploc_BUSCA_MCOT_aaseq), nbchar = 80, file.out = "path_to_fasta_file.fasta")


### combine all prediction files filtered for having transmembrane domain/s or Signal peptides 

ACP_adult_nymph_cello2go_deeploc_BUSCA_phobius <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_cello2go_deeploc_BUSCA_phobius_TM_SP <- ACP_adult_nymph_cello2go_deeploc_BUSCA_phobius %>%
  select(Accession_Number:SP) %>%
  filter(TM >= 1 | SP == "Y")
write_csv(ACP_adult_nymph_cello2go_deeploc_BUSCA_phobius_TM_SP  ,"path_to_csv_file.csv")
### combine all prediction files filtered for having transmembrane domain/s, GPI-anchor and Signal peptides from Phobius, PredGPI, GPISOM,Cello2GO, BUSCA and DeepLOC analyses
ACP_adult_nymph_phobius <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_BUSCA <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_cello2go <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_deeploc <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_PredGPI <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_GPISom <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_cello_deep_BUSCA_PreGPI_Phobius_GIPSom <- full_join(ACP_adult_nymph_PredGPI, ACP_adult_nymph_GPISom, by= "Accession_Number") %>%
  full_join(ACP_adult_nymph_phobius, by= "Accession_Number") %>%
  full_join(ACP_adult_nymph_deeploc, by= "Accession_Number") %>%
  full_join(ACP_adult_nymph_cello2go, by= "Accession_Number") %>%
  full_join(ACP_adult_nymph_BUSCA, by= "Accession_Number") %>%
  distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_adult_nymph_cello_deep_BUSCA_PreGPI_Phobius_GIPSom, "path_to_csv_file.csv")


### Statistical analysis
## finding protein amino acid length for calculation of distributed Normalized Spectral Abundance Factor (dNSAF)  

ACP_adul_nym_TSC_AccNum_stat <- ACP_proteomics_MCOT_TSC %>%
  filter(!str_detect(Identified_Proteins,"DECOY"))%>%
  select(Accession_Number) %>%
  mutate(Accession_Number= str_replace(Accession_Number,"\\s.*$", "")) %>% 
  distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_adul_nym_TSC_AccNum_stat, "path_to_csv_file.csv")

fastafile<- read.fasta(file = "path_to_fasta_file.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE) 
subsetlist<-read.csv("path_to_csv_file.csv",header = TRUE)
ACP_adul_nym_TSC_AccNum_stat_aaseq <- fastafile[c(which(names(fastafile) %in% subsetlist$Accession_Number))] 
write.fasta(ACP_adul_nym_TSC_AccNum_stat_aaseq, names = names(ACP_adul_nym_TSC_AccNum_stat_aaseq), nbchar = 80, file.out = "path_to_csv_file.csv")

ACP_adul_nym_TSC_AccNum_stat_aaseq_length <- summary(ACP_adul_nym_TSC_AccNum_stat_aaseq) 
write_csv(ACP_adul_nym_TSC_AccNum_stat_aaseq_length,"path_to_csv_file.csv")

ACP_proteomics_MCOT_TSC <- read_csv("path_to_csv_file.csv.csv") ## original data
ACP_adul_nym_TSC_AccNum_stat_aaseq_length_final <- left_join(ACP_adul_nym_TSC_AccNum_stat_aaseq_length, ACP_proteomics_MCOT_TSC, by= "Accession_Number") %>%
  distinct(Accession_Number, .keep_all = TRUE)
write_csv(ACP_adul_nym_TSC_AccNum_stat_aaseq_length_final,"path_to_csv_file.csv")

Unique_SpC <- read_csv("path_to_csv_file.csv.csv")
ACP_adul_nym_unique_SpC <- Unique_SpC %>%
  filter(!str_detect(Identified_Proteins,"DECOY"))%>%
  mutate(Accession_Number= str_replace(Accession_Number,"\\s.*$", "")) %>% # remove brackets 
  distinct(Accession_Number, .keep_all = TRUE) %>%
  select(Accession_Number, BBMV_Adult1_uniq:BBMV_Nymph3_uniq)
write_csv(ACP_adul_nym_unique_SpC, "path_to_csv_file.csv")

ACP_adul_nym_TSC_length_uniqueSpC_length_stat_final <- left_join(ACP_adul_nym_TSC_AccNum_stat_aaseq_length_final,ACP_adul_nym_unique_SpC, by = "Accession_Number")
write_csv(ACP_adul_nym_TSC_length_uniqueSpC_length_stat_final, "path_to_csv_file.csv")

## calculation of dNSAF from predicted plasma membrane proteins for nymphs (rep1, rep2, rep3) and adults (rep1, rep2 rep3). dNSAF should be calculated separate for each replicate

ACP_adul_nym_TSC_length_uSpC_length_stat_final <- read_csv("path_to_csv_file.csv")
ACP_adult_nymph_TM_Proteins <- read_csv("path_to_csv_file.csv")
ACP_adul_nym_TSC_length_uSpC_length_stat_final_PM <- left_join(ACP_adult_nymph_TM_Proteins,ACP_adul_nym_TSC_length_uSpC_length_stat_final, by="Accession_Number")

ACP_TSC_uSpC_length_adult_R1 <- ACP_adul_nym_TSC_length_uSpC_length_stat_final_PM %>%
  select(Accession_Number:BBMV_Adult1, BBMV_Adult1_uniq) %>%
  mutate(BBMV_Adult1 = ifelse(as.character(BBMV_Adult1) == "0", "0.29", as.character(BBMV_Adult1))) %>%
  mutate(BBMV_Adult1_uniq = ifelse(as.character(BBMV_Adult1_uniq) == "0", "0.29", as.character(BBMV_Adult1_uniq))) %>%
  rename(SpC = BBMV_Adult1, uSpC = BBMV_Adult1_uniq)
write_csv(ACP_TSC_uSpC_length_adult_R1, "path_to_csv_file.csv")

ACP_TSC_uSpC_length_adult_R1 <- read_csv("path_to_csv_file.csv")
ACP_TSC_uSpC_length_adult_R1_dNSAF <-ACP_TSC_uSpC_length_adult_R1 %>%
  mutate(sSpC = SpC - uSpC,
         dSAF = (uSpC + (uSpC/sum(uSpC))*sSpC)/Length,
         dNSAF =  dSAF/sum(dSAF),
         lndNSAF = log(dSAF))
write_csv(ACP_TSC_uSpC_length_adult_R1_dNSAF, "path_to_csv_file.csv.csv")

#NOTE! dNSAF should be calculated separate for every single replicate


## checking for normality of the lndNSAF values 

hist(ACP_TSC_uSpC_length_adult_nymph_R1_R3_dNSAF$lndNSAF)
qqnorm(ACP_TSC_uSpC_length_adult_nymph_R1_R3_dNSAF$lndNSAF)
ggqqplot(ACP_TSC_uSpC_length_adult_nymph_R1_R3_dNSAF$lndNSAF)







