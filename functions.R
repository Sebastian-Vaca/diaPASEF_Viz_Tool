library(shiny)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(UpSetR)
library(plotly)
library(shinyWidgets)
library(ggrepel)
library(heatmaply)
library(ggseqlogo)



####
###  Functions


load_spectronaut_data<-function(raw_input_file_path,metadata_filepath){
  metadata <- fread(metadata_filepath)
  
  raw_input_file <- fread(raw_input_file_path) %>%
    dplyr::select(-R.Condition, -R.Replicate) %>%
    left_join(metadata) %>%
    mutate(PEP.PeptideLength = nchar(PEP.StrippedSequence)) %>%
    mutate(PEP.IsProteotypic = as.logical(PEP.IsProteotypic))
  
  if(any(names(raw_input_file) %in% c("EG.IsImputed"))){
    
    raw_input_file<-raw_input_file %>%
      mutate(EG.IsImputed=as.logical(EG.IsImputed)) %>%
      filter(EG.IsImputed==F)
  }
  
  return(raw_input_file)
}

### load data from DIANN/TimsDIANN

format_DIANN_timsDIANN_toSpectronautFormat_longformat<-function(raw_input_file_path,metadata,software){
  if(software=="timsDIANN"){
    
    cols_selected=c(
      "File.Name","Protein.Group","Protein.Ids","Protein.Names","PG.Normalised",
      "Genes","Modified.Sequence","Precursor.Charge",
      "Q.Value","Protein.Q.Value","PG.Q.Value","GG.Q.Value",
      "Precursor.Quantity","Precursor.Normalised",
      "RT","Predicted.RT","CScore",
      "Precursor.Mz","Exp.1/K0","Precursor.FWHM","RT.Start","RT.Stop", "Proteotypic")
    
    diann<-fread(raw_input_file_path,
                 select=cols_selected)%>%
      rename(
        R.FileName=File.Name,
        PG.Genes=Genes,
        PG.ProteinAccessions=Protein.Ids,
        PG.ProteinGroups=Protein.Group,
        PG.ProteinNames=Protein.Names,
        PG.Qvalue=PG.Q.Value,
        PG.MS2Quantity=PG.Normalised,
        EG.IonMobility=`Exp.1/K0`,
        EG.ModifiedPeptide=Modified.Sequence,
        EG.Qvalue=Q.Value,
        EG.ApexRT=RT,
        EG.RTPredicted=Predicted.RT,
        EG.Cscore=CScore,
        FG.Charge=Precursor.Charge,
        FG.PrecMz=Precursor.Mz,
        FG.MS2Quantity=Precursor.Quantity,
        FG.MS2RawQuantity=Precursor.Normalised,
        EG.FWHM=Precursor.FWHM,
        PEP.IsProteotypic = Proteotypic)%>%
      left_join(metadata,by="R.FileName")%>%
      mutate(FG.MS2Quantity=as.numeric(FG.MS2Quantity))%>%
      mutate(EG.PeakWidth=(RT.Stop-RT.Start)*60)%>%
      mutate(PEP.IsProteotypic = as.logical(PEP.IsProteotypic))
  }
  
  if(software=="DIANN"){
    
    cols_selected=c(
      "File.Name","Protein.Group","Protein.Ids","Protein.Names","PG.Normalised",
      "Genes","Stripped.Sequence","Modified.Sequence","Precursor.Charge",
      "Q.Value","Protein.Q.Value","PG.Q.Value","GG.Q.Value",
      "Precursor.Quantity","Precursor.Normalised",
      "RT","Predicted.RT","CScore","IM","RT.Start","RT.Stop","Proteotypic")
    
    diann<-fread(raw_input_file_path, select=cols_selected) %>%
      rename(
        R.FileName = File.Name,
        PG.Genes = Genes,
        PG.ProteinAccessions = Protein.Ids,
        PG.ProteinGroups = Protein.Group,
        PG.ProteinNames = Protein.Names,
        PG.Qvalue = PG.Q.Value,
        PG.MS2Quantity = PG.Normalised,
        EG.IonMobility = IM,
        EG.ModifiedPeptide = Modified.Sequence,
        PEP.StrippedSequence = Stripped.Sequence,
        EG.Qvalue = Q.Value,
        EG.ApexRT = RT,
        EG.RTPredicted = Predicted.RT,
        EG.Cscore = CScore,
        FG.Charge = Precursor.Charge,
        #FG.PrecMz=Precursor.Mz,
        FG.MS2Quantity = Precursor.Quantity,
        FG.MS2RawQuantity = Precursor.Normalised,
        PEP.IsProteotypic = Proteotypic)%>%
      left_join(metadata,by="R.FileName")%>%
      mutate(FG.MS2Quantity = as.numeric(FG.MS2Quantity))%>%
      mutate(EG.PeakWidth = (RT.Stop - RT.Start)*60)%>%
      mutate(PEP.PeptideLength = nchar(PEP.StrippedSequence)) %>%
      mutate(PEP.IsProteotypic = as.logical(PEP.IsProteotypic))
  }
  
  return(diann)
}


load_diann_data_2<-function(raw_input_file_path,metadata_filepath,software){
  metadata<-fread(metadata_filepath)
  
  df<-format_DIANN_timsDIANN_toSpectronautFormat_longformat(raw_input_file_path,metadata,software)
  
  return(df)
}



load_data_2<-function(raw_input_file_path,metadata_filepath,software,needs_reassign){
  if(software%in%c("timsDIANN","DIANN")){
    df <- load_diann_data_2(raw_input_file_path = raw_input_file_path,
                            metadata_filepath = metadata_filepath,
                            software = software)
  }
  
  
  if(software=="Spectronaut"){
    df <- load_spectronaut_data(raw_input_file_path,metadata_filepath)
  }
  
  if(software=="null"){ df <- data.frame(Error="SelectSoftware") }
  if(is.null(software)){ df <- data.frame(Error="SelectSoftware") }
  
  return(df)
  
}

remove_data_based_on_metadata<-function(dt,remove_selected_runs=F){
  
  if(remove_selected_runs){
    dt2 <- dt %>%
      mutate(remove = ifelse(is.na(remove),0,remove)) %>%
      filter(remove != 1)
  }else{
    dt2 <- dt
  }
  
  return(dt2)
}

load_data_filter_data<-function(dt_unfiltered,
                                filter_EG_qValue = T,
                                EG_qValue_cutoff = 0.01,
                                filter_PG_qValue = T,
                                PG_qValue_cutoff = 0.01,
                                filter_PeptideLenght = F,
                                PeptideLenght_min_cutoff = 7,
                                PeptideLenght_max_cutoff = 35,
                                filter_isProteotypic = F
                                # filter_Protein_qValue=T,
                                # Protein_qValue_cutoff=0.01,
                                # filter_GG_qValue=T,
                                # GG_qValue_cutoff=0.01
){
  
  dt = dt_unfiltered %>%
    mutate(EG.Qvalue = as.numeric(EG.Qvalue)) %>%
    mutate(PG.Qvalue = as.numeric(PG.Qvalue)) %>%
    mutate(PEP.PeptideLength = as.numeric(PEP.PeptideLength))

  
  
  if(filter_EG_qValue){
    dt <- dt %>%
      filter(EG.Qvalue <= EG_qValue_cutoff)
  }
  
  if(filter_PG_qValue){
    dt <- dt %>%
      filter(PG.Qvalue <= PG_qValue_cutoff)
  }
  
  if(filter_PeptideLenght){
    dt <- dt %>%
      filter(PEP.PeptideLength >= PeptideLenght_min_cutoff) %>%
      filter(PEP.PeptideLength <= PeptideLenght_max_cutoff)
    
  }
  
  if(filter_isProteotypic){
    dt <- dt %>%
      filter(PEP.IsProteotypic)
  }
  
  # if(filter_Protein_qValue){
  # dt <- dt %>%
  # filter(Protein.Q.Value <= Protein_qValue_cutoff)
  # }
  # 
  # if(filter_GG_qValue){
  # dt <- dt %>%
  # filter(GG.Q.Value <= GG_qValue_cutoff)
  # }
  
  return(dt)
  
}

load_data_filter_text<-function(filter_EG_qValue = T,
                                EG_qValue_cutoff = 0.01,
                                filter_PG_qValue = T,
                                PG_qValue_cutoff = 0.01,
                                filter_PeptideLenght = F,
                                PeptideLenght_min_cutoff = 7,
                                PeptideLenght_max_cutoff = 35,
                                filter_isProteotypic = F,
                                remove_checkbox_input = F){
  
  x = "Data shown has:"
  
  if(filter_EG_qValue){
    x = c(x,paste0("- EG q-value <= ", EG_qValue_cutoff))
  }
  
  if(filter_PG_qValue){
    x = c(x, paste0( "- PG q-value <= ", PG_qValue_cutoff))
  }
  
  if(filter_PeptideLenght){
    x = c(x, paste0( "- ", PeptideLenght_min_cutoff, "<= Peptide lenght <= ", PeptideLenght_max_cutoff))
    
  }
  
  if(filter_isProteotypic){
    x = c(x, paste0( "- Only Proteotypic peptides"))
    
  }
  
  if(remove_checkbox_input){
    x = c(x, paste0( "- Selected MS runs removed"))
    
  }
  
  x = paste0(x, collapse = "<br>")
  return(x)
  
}


load_data_input_verification_afterLoad<-function(dt){
  
  report_names<- dt %>%
    names()
  report_names = data.frame(column_name = report_names, in_report =T)
  
  columns_needed <- c("R.FileName", "PG.ProteinGroups", "PG.ProteinAccessions",
                      "PG.ProteinNames", "PG.MS2Quantity", "PG.Genes",
                      "EG.ModifiedPeptide", "FG.Charge", "EG.Qvalue",
                      "Protein.Q.Value", "PG.Qvalue", "GG.Q.Value",
                      "FG.MS2Quantity", "FG.MS2RawQuantity", "EG.ApexRT",
                      "EG.RTPredicted", "EG.Cscore", "EG.IonMobility",
                      "RT.Start", "RT.Stop", "R.Condition",
                      "R.Replicate", "order", "remove", "Concentration",
                      "EG.PeakWidth","PEP.StrippedSequence","PEP.IsProteotypic","PEP.PeptideLength")
  
  columns_needed = data.frame(column_name = columns_needed, necessary_column = T) %>%
    # mutate(is_present = ifelse(column_name %in% report_names,1,0)) %>%
    full_join(report_names, by = "column_name")
  
  return(columns_needed)
}


## Metadata template
create_metadata_template <- function(input_filepath, software = NULL){
  
  if(software %in% c("timsDIANN", "DIANN")){
    df <- fread(input_filepath, select = c("File.Name")) %>% distinct()
    R.FileName = df$File.Name
    metadata_template = data.frame(R.FileName, R.Condition = "", R.Replicate="",	order="", remove ="", Concentration = "")
    
  }
  
  if(software == "Spectronaut"){
    df <- fread(input_filepath,select = c("R.FileName", "R.Condition","R.Replicate")) %>% distinct()
    
    metadata_template = data.frame(df, order="", remove ="", Concentration = "") %>%
      dplyr::select(R.FileName, R.Condition , 	R.Replicate,	order, remove, Concentration)
    
  }
  
  if(software == "null"){
    metadata_template = data.frame(R.FileName = "", R.Condition = "", 	R.Replicate="",	order="", remove ="", Concentration = "")
    }
  if(is.null(software)) {
    metadata_template = data.frame(R.FileName = "", R.Condition = "", 	R.Replicate="",	order="", remove ="", Concentration = "")
  }
  
  return(metadata_template)
  
}
# Upset Plots
upset_plot_from_spectronaut_precs <- function(dt){
  for_upset_precs <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge,EG.Qvalue) %>%
    filter(EG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(EG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge, observed) %>%
    mutate(prec_id = paste0(EG.ModifiedPeptide,"_", FG.Charge)) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0)
  rownames(for_upset_precs) <- for_upset_precs$prec_id
  
  upset_plot = upset(for_upset_precs, nsets = length(unique(dt$R.Condition)),
                     number.angles = 45, point.size = 3.5, line.size = 2,  
                     mainbar.y.label = "Precursors Intersections", sets.x.label = "# Precursors", 
                     text.scale = c(1.3, 1.3, 1), order.by = "freq") #550x330
  
  return(upset_plot)
}
upset_plot_from_spectronaut_precs_withConditionsFilter <- function(dt, conditions){
  
  if(!is.null(conditions)){
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  for_upset_precs <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge,EG.Qvalue) %>%
    filter(EG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(EG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge, observed) %>%
    mutate(prec_id = paste0(EG.ModifiedPeptide,"_", FG.Charge)) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0)
  rownames(for_upset_precs) <- for_upset_precs$prec_id
  
  upset_plot = upset(for_upset_precs, nsets = length(unique(dt$R.Condition)),
                     number.angles = 45, point.size = 3.5, line.size = 2,  
                     mainbar.y.label = "Precursors Intersections", sets.x.label = "# Precursors", 
                     text.scale = c(1.3, 1.3, 1), order.by = "freq") #550x330
  
  return(upset_plot)
}
upset_plot_from_spectronaut_precs_matrix <- function(dt, conditions){
  
  if(!is.null(conditions)){
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  for_upset_precs <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge,EG.Qvalue) %>%
    filter(EG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(EG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,EG.ModifiedPeptide,FG.Charge, observed) %>%
    mutate(prec_id = paste0(EG.ModifiedPeptide,"_", FG.Charge)) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0)
  rownames(for_upset_precs) <- for_upset_precs$prec_id
  
  
  return(for_upset_precs)
}

upset_plot_from_spectronaut_prots <- function(dt){
  for_upset_prots <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, PG.Qvalue) %>%
    filter(PG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(PG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, observed) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0) %>%
    distinct()
  rownames(for_upset_prots) <- for_upset_prots$PG.ProteinGroups
  
  upset_plot = upset(for_upset_prots, nsets = length(unique(dt$R.Condition)),
                     number.angles = 45, point.size = 3.5, line.size = 2,  
                     mainbar.y.label = "Proteins Intersections", sets.x.label = "# Proteins", 
                     text.scale = c(1.3, 1.3, 1), order.by = "freq") #550x330
  
  return(upset_plot)
}
upset_plot_from_spectronaut_prots_withConditionsFilter <- function(dt, conditions){
  
  if(!is.null(conditions)){
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  for_upset_prots <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, PG.Qvalue) %>%
    filter(PG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(PG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, observed) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0) %>%
    distinct()
  rownames(for_upset_prots) <- for_upset_prots$PG.ProteinGroups
  
  upset_plot = upset(for_upset_prots, nsets = length(unique(dt$R.Condition)),
                     number.angles = 45, point.size = 3.5, line.size = 2,  
                     mainbar.y.label = "Proteins Intersections", sets.x.label = "# Proteins", 
                     text.scale = c(1.3, 1.3, 1), order.by = "freq") #550x330
  
  return(upset_plot)
}
upset_plot_from_spectronaut_prots_matrix <- function(dt, conditions){
  
  if(!is.null(conditions)){
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  for_upset_prots <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, PG.Qvalue) %>%
    filter(PG.Qvalue<0.01) %>%
    mutate(observed = ifelse(is.na(PG.Qvalue), 0 , 1))%>%
    group_by(R.Condition) %>%
    dplyr::select(R.Condition,PG.ProteinGroups, observed) %>%
    top_n(n = 1, wt = observed) %>%
    distinct() %>%
    spread(key = R.Condition, observed) %>%
    data.frame() %>%
    replace(is.na(.), 0) %>%
    distinct()
  rownames(for_upset_prots) <- for_upset_prots$PG.ProteinGroups
  
  return(for_upset_prots)
}

## Metrics
create_counts_table <- function(dt){
  
  All_counts0 <- dt %>%
    dplyr::select(R.Condition) %>%
    distinct()%>%
    mutate(col1 = "Precursors", col2 = "Peptides", col3 =  "Proteins") %>%
    gather(key = cols, value = "type", which(str_detect(names(.),pattern = "col"))) %>%
    dplyr::select(-cols)
  ### Prots with >2peptides
  Prots_with_2_peptides = dt %>%
    dplyr::select(PG.ProteinGroups, EG.ModifiedPeptide) %>%
    distinct() %>%
    group_by(PG.ProteinGroups) %>%
    tally() %>%
    mutate(Prot_has_2_peptides = ifelse(n >= 2,1,0))
  
  
  
  counts <- dt %>%
    left_join(Prots_with_2_peptides) %>%
    filter(!is.na(FG.MS2Quantity)) %>%
    group_by(R.Condition) %>%
    summarise(n_precs = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[EG.Qvalue<0.01]),
              n_prots = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01 ]),
              n_pepts = n_distinct(EG.ModifiedPeptide[EG.Qvalue<0.01]),
              n_prots_with2Pepts = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01 & Prot_has_2_peptides == 1])) %>%
    gather(key = "type", value = "counts", which(str_detect(names(.), pattern = "n_")))%>%
    mutate(type = case_when(
      type == "n_precs" ~ "Precursors",
      type == "n_pepts" ~ "Peptides",
      type == "n_prots" ~ "Proteins",
      type == "n_prots_with2Pepts" ~ "Proteins_with2pepts"))
  
  counts_errorbars <- dt %>%
    left_join(Prots_with_2_peptides) %>%
    filter(!is.na(FG.MS2Quantity)) %>%
    group_by(R.FileName, R.Condition) %>%
    # summarise(n_precs = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)),
    #           n_prots = n_distinct(PG.ProteinGroups),
    #           n_pepts = n_distinct(EG.ModifiedPeptide)) %>%
    summarise(n_precs = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[EG.Qvalue<0.01]),
              n_prots = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01]),
              n_pepts = n_distinct(EG.ModifiedPeptide[EG.Qvalue<0.01]),
              n_prots_with2Pepts = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01 & Prot_has_2_peptides == 1])) %>%
    ungroup() %>%
    group_by(R.Condition) %>%
    summarise(mean_n_precs = mean(n_precs),sd_precs = sd(n_precs),
              mean_n_pepts = mean(n_pepts),sd_pepts = sd(n_pepts),
              mean_n_prots = mean(n_prots),sd_prots = sd(n_prots),
              mean_n_prots_with2Pepts = mean(n_prots_with2Pepts),sd_prots_with2Pepts = sd(n_prots_with2Pepts))%>%
    gather(key = "type", value = "value", which(str_detect(names(.), pattern = "mean|sd"))) %>%
    mutate(type2 = case_when(
      str_detect(type, "mean") ~ "mean",
      str_detect(type, "sd") ~ "sd")) %>%
    mutate(type = case_when(
      str_detect(type, "pepts") ~ "Peptides",
      str_detect(type, "precs") ~ "Precursors",
      str_detect(type, "prots$") ~ "Proteins",
      str_detect(type, "prots_with2Pepts") ~ "Proteins_with2pepts")) %>%
    spread(key = type2 , value = value)
  
  All_counts = All_counts0 %>%
    full_join(counts) %>%
    full_join(counts_errorbars) %>%
    replace(is.na(.), 0)
  
  return(All_counts)
  
}
create_counts_table_per_MSRun <- function(dt){
  
  All_counts0 <- dt %>%
    dplyr::select(R.Condition, R.Replicate) %>%
    mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
    distinct()%>%
    mutate(col1 = "Precursors", col2 = "Peptides", col3 =  "Proteins") %>%
    gather(key = cols, value = "type", which(str_detect(names(.),pattern = "col"))) %>%
    dplyr::select(-cols)
  
  ### Prots with >2peptides
  Prots_with_2_peptides = dt %>%
    dplyr::select(PG.ProteinGroups, EG.ModifiedPeptide) %>%
    distinct() %>%
    group_by(PG.ProteinGroups) %>%
    tally() %>%
    mutate(Prot_has_2_peptides = ifelse(n >= 2,1,0))
  
  
  
  counts <- dt %>%
    left_join(Prots_with_2_peptides) %>%
    filter(!is.na(FG.MS2Quantity)) %>%
    mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
    group_by(R.Condition,R.Replicate, run) %>%
    summarise(n_precs = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[EG.Qvalue<0.01]),
              n_prots = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01 ]),
              n_pepts = n_distinct(EG.ModifiedPeptide[EG.Qvalue<0.01]),
              n_prots_with2Pepts = n_distinct(PG.ProteinGroups[PG.Qvalue<0.01 & Prot_has_2_peptides == 1])) %>%
    gather(key = "type", value = "counts", which(str_detect(names(.), pattern = "n_")))%>%
    mutate(type = case_when(
      type == "n_precs" ~ "Precursors",
      type == "n_pepts" ~ "Peptides",
      type == "n_prots" ~ "Proteins",
      type == "n_prots_with2Pepts" ~ "Proteins_with2pepts"))
  
  
  All_counts = All_counts0 %>%
    full_join(counts) %>%
    replace(is.na(.), 0)
  
  return(All_counts)
  
}


plot_count_pepts_prots_precs_reordered <- function(counts, dt, reorder){
  
  if(reorder) {
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    
    counts$R.Condition <- factor(counts$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(counts>max(counts)/3, 1.1, -0.1))
  }
  
  counts <- counts %>%
    ungroup() %>%
    group_by(type) %>% 
    group_modify(~hjust_adjustment(.x))
  
  
  plot_counts = ggplot(counts, aes(x= R.Condition, y = counts, fill = R.Condition))+
    geom_bar(stat = "identity", color = "black")+
    geom_text(aes(label=counts), hjust=counts$hjust_value, angle=90, color= "black")+
    # scale_y_continuous (expand = c (0,0))+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~type, scales = "free_y")
  return(plot_counts)
  
}
plot_count_pepts_prots_precs_withErrorBars_reordered <- function(counts, dt, reorder){
  
  if(reorder) {
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    
    counts$R.Condition <- factor(counts$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(mean>max(mean)/3, 1.1, -0.1))
  }
  
  counts <- counts %>%
    ungroup() %>%
    group_by(type) %>% 
    group_modify(~hjust_adjustment(.x))
  
  plot_counts = ggplot(counts, aes(x= R.Condition, y = mean, fill = R.Condition))+
    geom_bar(stat = "identity", color= "black")+
    geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd),  width=0.25)+
    geom_text(aes(label=round(mean,0)), hjust=counts$hjust_value, angle=90, color= "black")+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~type, scales = "free_y")
  return(plot_counts)
  
}
plot_count_pepts_prots_precs_withErrorBars_reordered_filter <- function(counts, dt,plot_type = NULL, reorder){
  
  if(reorder) {
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    
    counts$R.Condition <- factor(counts$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(mean>max(mean)/3, 1.1, -0.1))
  }
  
  counts <- counts %>%
    ungroup() %>%
    group_by(type) %>% 
    group_modify(~hjust_adjustment(.x))
  
  if(!is.null(plot_type)){
    counts <- counts %>%
      filter(type %in% plot_type)
  }
  
  plot_counts = ggplot(counts, aes(x= R.Condition, y = mean, fill = R.Condition))+
    geom_bar(stat = "identity", color= "black")+
    geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd),  width=0.25)+
    geom_text(aes(label=round(mean,0)), hjust=counts$hjust_value, angle=90, color= "black")+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~type, scales = "free_y")
  return(plot_counts)
  
}

plot_count_per_MSRun_pepts_prots_precs_reordered <- function(counts_per_MSRun, dt, reorder){
  
  if(reorder) {
    order_dt <- dt %>%
      dplyr::select(R.Condition, R.Replicate, order) %>%
      distinct() %>%
      mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
      arrange(order)
    
    counts_per_MSRun$R.Condition <- factor(counts_per_MSRun$R.Condition, levels = unique(order_dt$R.Condition))
    counts_per_MSRun$run <- factor(counts_per_MSRun$run, levels = unique(order_dt$run))
    
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(counts>max(counts)/3, 1.1, -0.1))
  }
  
  counts_per_MSRun <- counts_per_MSRun %>%
    ungroup() %>%
    group_by(type) %>% 
    group_modify(~hjust_adjustment(.x))
  
  
  plot_counts = ggplot(counts_per_MSRun, aes(x= run, y = counts, fill = R.Condition))+
    geom_bar(stat = "identity", color = "black")+
    geom_text(aes(label=counts), hjust=counts_per_MSRun$hjust_value, angle=90, color= "black")+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~type, scales = "free_y")
  return(plot_counts)
  
}
plot_count_per_MSRun_pepts_prots_precs_reordered_filter <- function(counts_per_MSRun, dt, plot_type = NULL, reorder){
  
  if(reorder) {
    order_dt <- dt %>%
      dplyr::select(R.Condition, R.Replicate, order) %>%
      distinct() %>%
      mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
      arrange(order)
    
    counts_per_MSRun$R.Condition <- factor(counts_per_MSRun$R.Condition, levels = unique(order_dt$R.Condition))
    counts_per_MSRun$run <- factor(counts_per_MSRun$run, levels = unique(order_dt$run))
    
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(counts>max(counts)/3, 1.1, -0.1))
  }
  
  counts_per_MSRun <- counts_per_MSRun %>%
    ungroup() %>%
    group_by(type) %>% 
    group_modify(~hjust_adjustment(.x))
  
  
  if(!is.null(plot_type)){
    counts_per_MSRun <- counts_per_MSRun %>%
      filter(type %in% plot_type)
  }
  
  plot_counts = ggplot(counts_per_MSRun, aes(x= run, y = counts, fill = R.Condition))+
    geom_bar(stat = "identity", color = "black")+
    geom_text(aes(label=counts), hjust=counts_per_MSRun$hjust_value, angle=90, color= "black")+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~type, scales = "free_y")
  return(plot_counts)
  
}

plot_intensity_boxplot <- function(dt, reorder){
  
  if(reorder){
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    dt$R.Condition <- factor(dt$R.Condition, levels = unique(order_dt$R.Condition))
  }
  G = ggplot(dt, aes(x= R.Condition, y = log(FG.MS2Quantity,10), color = R.Condition))+
    geom_violin()+
    geom_boxplot(width = 0.25)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    labs(y = "Log10(Intensity)", x = "Condition")
  return(G)
  
}

plot_and_table_waterfall <- function(dt, condition, highlight_targets = F,
                                     list_of_targetProteinGroup = NULL,
                                     mean_or_median = c("mean", "median"),
                                     quant_value_to_use = c("PG.MS2Quantity", "FG.MS2Quantity")){
  
  if(mean_or_median == "mean") { mean_or_median_fct = mean}
  if(mean_or_median == "median"){ mean_or_median_fct = median}
  
  names(dt)[names(dt) == quant_value_to_use] <- 'Column_of_interest'
  
  
  dt2 <- dt %>%
    filter(!is.na(Column_of_interest)) %>%
    filter(Column_of_interest != 0) %>%
    ungroup() %>%
    filter(R.Condition == condition) %>%
    dplyr::select(PG.ProteinGroups, Column_of_interest) %>%
    distinct() %>%
    group_by(PG.ProteinGroups) %>%
    summarise(area = mean_or_median_fct(Column_of_interest, na.rm = T)) %>%
    mutate(rank = rank(-area, ties.method = "first")) %>%
    arrange(rank) %>%
    mutate(is_target = 
             case_when(
               highlight_targets == F ~ 0,
               highlight_targets == T & !(PG.ProteinGroups %in% list_of_targetProteinGroup) ~ 0,
               highlight_targets == T & PG.ProteinGroups %in% list_of_targetProteinGroup ~ 1))%>%
    mutate(label = ifelse(is_target ==0, "", PG.ProteinGroups))
  
  G = ggplot(dt2, aes(x= rank, y = log(area,10), label = label ))+
    geom_point(color = "grey50") +
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    labs(y = "Log10(Area)", x = "Rank")
  
  
  dt_fraction_plot <- dt2 %>%
    ungroup() %>%
    mutate(relative_fraction = 1-area/max(area,na.rm=T))
  
  dt_fraction_plot_counts = dt_fraction_plot %>%
    ungroup() %>%
    summarise(n_25 = n_distinct(PG.ProteinGroups[relative_fraction<0.25]),
              n_50 = n_distinct(PG.ProteinGroups[relative_fraction<0.5]),
              n_75 = n_distinct(PG.ProteinGroups[relative_fraction<0.75])) %>%
    gather() %>%
    mutate(relative_fraction = as.numeric(gsub(key, pattern = "n_",replacement = ""))/100,
           rank=max(dt_fraction_plot$rank)*0.1) %>%
    mutate(value = paste0(value, " Protein groups"))
  
  
  relative_fraction_plot = ggplot(dt_fraction_plot, aes(rank,relative_fraction ))+
    geom_point(color = "grey50")+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    ggtitle(paste0(condition, "   (", mean_or_median,"_",
                   gsub(quant_value_to_use,pattern = "\\.MS2Quantity",replacement = ""),")"))+
    geom_text(data =dt_fraction_plot_counts, aes(label = value),  hjust = "inward")+
    geom_segment(aes(x=-1,xend=max(dt_fraction_plot$rank)*0.09,y=0.25,yend=0.25), linetype= 2)+
    geom_segment(aes(x=-1,xend=max(dt_fraction_plot$rank)*0.09,y=0.5,yend=0.5), linetype= 2)+
    geom_segment(aes(x=-1,xend=max(dt_fraction_plot$rank)*0.09,y=0.75,yend=0.75), linetype= 2)
  
  
  
  if(highlight_targets){
    G = G +
      geom_text_repel(data =dt2[dt2$is_target !=0,], box.padding = 0.5, max.overlaps = Inf) +
      geom_point(data = dt2[dt2$is_target !=0,], color = "red")
    
    
    # dt_targets <-dt2 %>%
    #   filter(PG.ProteinGroups %in% list_of_targetProteinGroup) %>%
    #   dplyr::select(-is_target, -label)
  }
  
  plot_waterfall = cowplot::plot_grid(relative_fraction_plot, G, ncol=1, align = "v", rel_heights = c(1,4))
  return(plot_waterfall = plot_waterfall)
}

plot_mz_IM_map <- function(dt, condition, replicate, binwidth_mz = 50, binwidth_im = 0.01,
                           count_upper_limit = NULL,
                           mz_limits = c(300,1200),
                           im_limits = c(0.7,1.3)){
  
  map_data = dt %>%
    filter(R.Condition == condition,
           R.Replicate == replicate)
  
  density_2d <- ggplot(map_data, aes(x=FG.PrecMz, y=EG.IonMobility))+
    geom_bin2d(binwidth = c(binwidth_mz, binwidth_im)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    labs(y = "ion-mobility", x = "m/z")
  
  if (!is.null(count_upper_limit)) {
    
    density_2d <- density_2d +
      scale_fill_continuous(type = "viridis", limits = c(0,count_upper_limit), 
                            oob = scales::squish )
    
  }
  
  legend <- get_legend(density_2d)
  
  density_2d <- density_2d+
    theme(legend.position = "none",
          panel.grid = element_blank())
  
  hist_x <- axis_canvas(density_2d, axis = 'x') +
    geom_histogram(data = map_data, aes(x=FG.PrecMz), color= "black", binwidth = binwidth_mz)
  hist_y <- axis_canvas(density_2d, axis = 'y', coord_flip = TRUE) + 
    geom_histogram(data = map_data, aes(x=EG.IonMobility), color= "black",binwidth = binwidth_im)+
    coord_flip(xlim = im_limits)
  
  density_2d <- density_2d+
    coord_cartesian(xlim = mz_limits, ylim = im_limits)
  hist_x <- hist_x+
    coord_cartesian(xlim = mz_limits)
  
  
  
  G = density_2d %>% insert_xaxis_grob(hist_x, grid::unit(.2, "null"), position = "top") %>%
    insert_yaxis_grob(hist_y, grid::unit(.2, "null"), position = "right") %>%
    ggdraw()
  
  title <- ggdraw() + 
    draw_label(
      unique(paste0(map_data$R.Condition, "_", map_data$R.Replicate)),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  G2 = plot_grid(title,plot_grid(G, legend, rel_widths = c(1,0.2)), ncol = 1, rel_heights = c(0.1,1))
  
  
  return(G2)
}


data_for_one_metric <- function(dt,
                                metadata,
                                software,
                                metric_to_plot ){
  
  metadata2 = metadata %>%
    mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
    arrange(order)
    
    
  names(dt)[names(dt) == metric_to_plot] <- 'Column_of_interest'
  
  if(str_detect(metric_to_plot, pattern = "Quantity")){
    dt <- dt %>% 
      mutate(Column_of_interest =log(Column_of_interest,10))
  }
  
  if(str_detect(metric_to_plot, pattern = "FWHM|Width") & software == "Spectronaut"){
    dt <- dt %>% 
      mutate(Column_of_interest =Column_of_interest*60)
  }
  
  if(str_detect(metric_to_plot, pattern = "Qvalue")){
    dt <- dt %>% 
      mutate(Column_of_interest =-log(Column_of_interest,10))
  }
  
  dt_2 <- dt %>%
    dplyr::select(R.Condition, R.Replicate, EG.ModifiedPeptide, FG.Charge, Column_of_interest) %>%
    mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
    distinct()
  
  if(str_detect(metric_to_plot, pattern = "PG\\.")){
    dt_2 <- dt %>%
      dplyr::select(R.Condition, R.Replicate, PG.ProteinGroups, PG.ProteinAccessions,PG.ProteinNames, Column_of_interest) %>%
      mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
      distinct()
  }
  
  dt_2$run = factor(dt_2$run , levels = unique(metadata2$run) )
  
  
  return(dt_2)
  
}

create_list_matrices_for_each_metric <- function(list_matrices_for_each_metric,
                                                 dt,
                                                 metadata,
                                                 software,
                                                 metric_to_plot){
  
  dt2 = data_for_one_metric(dt, metadata, software, metric_to_plot)
  
  list_matrices_for_each_metric [[metric_to_plot]] = dt2
  
  return(list_matrices_for_each_metric)
}

Boxplot_OneMetric<- function(list_matrices_for_each_metric,
                               metric_to_plot,
                               conditions = NULL,
                               plot_type = c("freqpoly","histogram","boxplot"),
                               remove_outliers){
  
  dt_2 = list_matrices_for_each_metric[[metric_to_plot]]
  
  
  if(!is.null(conditions)){
    dt_2 <- dt_2 %>%
      filter(R.Condition %in% conditions)
  }
  
  if(remove_outliers){
    dt_2 = dt_2 %>%
      filter(Column_of_interest < quantile(dt_2$Column_of_interest, 0.95,na.rm = T))%>%
      filter(Column_of_interest > quantile(dt_2$Column_of_interest, 0.05,na.rm = T))
  }
  
  if(plot_type == "boxplot"){
    
    G = ggplot(dt_2, aes(x=run, y= Column_of_interest, fill= R.Condition))+
      geom_violin()+
      geom_boxplot(width=0.25)+
      theme_bw()+
      labs(y= metric_to_plot)+
      theme(axis.text.x = element_text(angle=90))+
      labs(y = metric_to_plot)
  }
  if(plot_type == "histogram"){
    
    G = ggplot(dt_2, aes(x=Column_of_interest, fill= R.Condition))+
      geom_histogram( color="black")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90))+
      labs(y= metric_to_plot)+
      facet_wrap(~run)+
      labs(x = metric_to_plot)
  }
  
  if(plot_type == "freqpoly"){
    
    G = ggplot(dt_2, aes(x=Column_of_interest, color= run))+
      geom_freqpoly()+
      theme_bw()+
      labs(x= metric_to_plot)
    
    
  }
  
  # G= ggplotly(G)
  
  return(G)
  
}

plot_charge_states <- function(dt,conditions = NULL){
  
  if(length(conditions)>=1) {
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  dt2 = dt %>%
    group_by(R.Condition) %>%
    select(R.Condition, EG.ModifiedPeptide, FG.Charge) %>% 
    distinct() %>%
    ungroup() %>%
    group_by(R.Condition, FG.Charge) %>% 
    summarise(num_peptides = n_distinct(EG.ModifiedPeptide)) %>%
    ungroup() %>%
    group_by(R.Condition) %>%
    mutate(Percent = paste0(round(num_peptides/sum(num_peptides)*100,1),"%"))
  
  
  dt2$FG.Charge <- factor(dt2$FG.Charge, unique(sort(dt2$FG.Charge)))
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(num_peptides>max(num_peptides)/3, 1.1, -0.1))
  }
  
  dt2 <- dt2 %>%
    ungroup() %>%
    group_by(R.Condition) %>% 
    group_modify(~hjust_adjustment(.x))
  
  
  ggplot(dt2, aes(x =FG.Charge, y = num_peptides, fill= FG.Charge))+
    geom_bar(stat = "identity", color= "black")+
    facet_wrap(~R.Condition)+
    theme_bw()+
    geom_text(aes(label=Percent), hjust=dt2$hjust_value, angle=90,size=3.5, color= "black")+
    scale_fill_viridis(discrete = T,option = "B",direction = -1)
  
}

plot_peptide_length <- function(dt,conditions = NULL,
                                type = c("histogram", "freqpoly")){
  
  if(length(conditions)>=1) {
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  dt1 = dt %>%
    select(R.Condition, EG.ModifiedPeptide) %>% 
    distinct() %>%
    ungroup() %>%
    mutate(Sequence = gsub("\\(([^()]|(?R))*\\)", "", EG.ModifiedPeptide, perl=TRUE)) %>%
    mutate(Sequence = gsub("\\[[^\\]]*\\]", "", Sequence, perl=TRUE))%>%
    mutate(Sequence = gsub("[^A-Za-z0-9]", ".", Sequence)) %>%
    mutate(Sequence = gsub("\\.", "", Sequence, perl=TRUE)) %>%
    mutate(PeptideLength = nchar(Sequence))
  
  dt2 <- dt1  %>%
    group_by(R.Condition, PeptideLength) %>% 
    summarise(num_peptides = n_distinct(Sequence)) %>%
    ungroup() %>%
    group_by(R.Condition) #%>%
  # mutate(Percent = paste0(round(num_peptides/sum(num_peptides)*100,1),"%"))
  
  
  dt2$PeptideLength <- factor(dt2$PeptideLength, unique(sort(dt2$PeptideLength)))
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(num_peptides>max(num_peptides)/3, 1.1, -0.1))
  }
  
  dt2 <- dt2 %>%
    ungroup() %>%
    group_by(R.Condition) %>% 
    group_modify(~hjust_adjustment(.x))
  
  if(type == "histogram"){
    G = ggplot(dt2, aes(x =PeptideLength, y = num_peptides, fill= PeptideLength))+
      geom_bar(stat = "identity", color= "black")+
      facet_wrap(~R.Condition)+
      theme_bw()+
      # geom_text(aes(label=Percent), hjust=dt2$hjust_value, angle=90,size=3.5, color= "black")+
      scale_fill_viridis(discrete = T)+
      theme(axis.text.x = element_text(angle=90),
            legend.position = "none")
  }
  
  
  if(type == "freqpoly"){
    G = ggplot(dt1, aes(x = PeptideLength, color= as.character(R.Condition)))+
      geom_freqpoly(binwidth=1)+
      theme_bw()+
      scale_fill_viridis(discrete = T)+
      theme(axis.text.x = element_text(angle=90))+
      guides(color=guide_legend(title="Condition"))
  }
  
  return(G)
  
  
  
  
}

Precursors_per_Protein <- function(dt){
  dt2 = dt %>%
    mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
    group_by(R.Condition, run, PG.ProteinGroups) %>%
    summarise(num_precursors = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)),
              num_peptides = n_distinct(paste0(EG.ModifiedPeptide)))
  
  return(dt2)
}

plot_number_of_peptides_hist <- function(Precursors_per_Protein_dt,
                                         item = c("Precursors", "Peptides"),
                                         conditions = NULL,
                                         remove_outliers = F){
  
  if(!is.null(conditions)){
    Precursors_per_Protein_dt = Precursors_per_Protein_dt %>%
      filter(R.Condition %in% conditions)
  }
  
  
  if(item == "Precursors"){
    G = ggplot( Precursors_per_Protein_dt, aes(num_precursors))+
      geom_histogram(binwidth = 1, color="#313331", fill="#008DEB")+
      facet_wrap(~run)+
      theme_bw()
    if(remove_outliers){
      G = G +
        coord_cartesian(xlim = c(0,
                                 quantile(Precursors_per_Protein_dt$num_precursors,0.95)))
    }
  }
  
  if(item == "Peptides"){
    G = ggplot( Precursors_per_Protein_dt, aes(num_peptides))+
      geom_histogram(binwidth = 1, color="#313331", fill="#008DEB")+
      facet_wrap(~run)+
      theme_bw()
    if(remove_outliers){
      G = G +
        coord_cartesian(xlim = c(0,
                                 quantile(Precursors_per_Protein_dt$num_peptides,0.95)))
    }
  }
  
  
  
  return(G)
  
}

plot_number_of_peptides_boxplot<- function(Precursors_per_Protein_dt,
                                           item = c("Precursors", "Peptides"),
                                           conditions = NULL,
                                           remove_outliers = F){
  
  if(!is.null(conditions)){
    Precursors_per_Protein_dt = Precursors_per_Protein_dt %>%
      filter(R.Condition %in% conditions)
  }
  
  
  if(item == "Precursors"){
    G = ggplot(Precursors_per_Protein_dt, aes(x = run , y=num_precursors, fill= R.Condition))+
      geom_violin()+
      geom_boxplot(color = "black", width = 0.25, alpha=0.2)+
      theme_bw()+
      scale_y_continuous(n.breaks = 20) +
      theme(axis.text.x = element_text(angle=90),
            legend.position = "none")
    if(remove_outliers){
      G = G +
        coord_cartesian(ylim = c(0,
                                 quantile(Precursors_per_Protein_dt$num_precursors,0.95)))
    }
  }
  
  if(item == "Peptides"){
    G = ggplot(Precursors_per_Protein_dt, aes(x = run , y=num_peptides, fill= R.Condition))+
      geom_violin()+
      geom_boxplot(color = "black", width = 0.25, alpha=0.2)+
      theme_bw()+
      scale_y_continuous(n.breaks = 20) +
      theme(axis.text.x = element_text(angle=90),
            legend.position = "none")
    if(remove_outliers){
      G = G +
        coord_cartesian(ylim = c(0,
                                 quantile(Precursors_per_Protein_dt$num_peptides,0.95)))
    }
  }
  
  
  
  return(G)
  
}

plot_number_of_peptides <- function(Precursors_per_Protein_dt,
                                    item = c("Precursors", "Peptides"),
                                    conditions = NULL,
                                    remove_outliers = F,
                                    plot_type = c("histogram", "boxplot")){
  if(plot_type == "histogram"){
    H = plot_number_of_peptides_hist(Precursors_per_Protein_dt,
                                     item , conditions, remove_outliers)
  }
  if(plot_type == "boxplot"){
    H = plot_number_of_peptides_boxplot(Precursors_per_Protein_dt,
                                        item , conditions, remove_outliers)
  }
  
  return(H)
  
}


### CVs
## CV precursor level
Calculate_CVs <- function(dt){
  CVs <- dt %>%
    ungroup %>%
    group_by(R.Condition, EG.ModifiedPeptide, FG.Charge) %>%
    summarise(mean_quant = mean(FG.MS2Quantity, na.rm = T),
              sd_quant = sd(FG.MS2Quantity, na.rm =T)) %>%
    mutate(cvs = sd_quant / mean_quant*100)%>%
    ungroup()
  
  return(CVs)
}
## CV Protein Level
Calculate_protein_CVs <- function(dt){
  CVs <- dt %>%
    ungroup %>%
    dplyr::select(R.Condition, R.Replicate, PG.ProteinGroups, PG.MS2Quantity) %>%
    distinct() %>%
    group_by(R.Condition, PG.ProteinGroups) %>%
    summarise(mean_quant = mean(PG.MS2Quantity, na.rm = T),
              sd_quant = sd(PG.MS2Quantity, na.rm =T)) %>%
    mutate(cvs = sd_quant / mean_quant*100)%>%
    ungroup()
  
  return(CVs)
}

plot_CVs <- function(dt){
  CVs <- dt %>%
    filter(!is.na(as.numeric(FG.MS2Quantity))) %>%
    ungroup %>%
    group_by(R.Condition, EG.ModifiedPeptide, FG.Charge) %>%
    summarise(mean_quant = mean(FG.MS2Quantity, na.rm = T),
              sd_quant = sd(FG.MS2Quantity, na.rm =T)) %>%
    mutate(cvs = sd_quant / mean_quant*100)
  
  plot_cvs1 = ggplot(CVs, aes(x= R.Condition, y = cvs, color = R.Condition))+
    geom_violin()+
    geom_boxplot(width = 0.25)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    labs(y = "CVs (%)", x = "Condition")
  
  plot_cvs2 = ggplot(CVs, aes(x= cvs, color = R.Condition))+
    geom_freqpoly(binwidth = 2)+
    theme_bw()+
    theme(legend.position = "none")+
    labs(x = "CVs (%)")
  plot_cvs3 = get_legend(ggplot(CVs, aes(x= R.Condition, y = cvs, fill = R.Condition))+
                           geom_violin()+
                           guides(fill=guide_legend(nrow=5, byrow=TRUE)))
  plot_cvs = cowplot::plot_grid(plot_grid(plot_cvs1, plot_cvs2),
                                plot_cvs3, nrow=2,rel_widths = c(2,1),rel_heights = c(2,1))
  
  return(plot_cvs)
  
}
plot_CVs_reordered_2<- function(CVs,
                                metadata,
                                reorder,
                                conditions = NULL,
                                plot_type = c("freqpoly","histogram","boxplot"),
                                remove_outliers){
  
  if(!is.null(conditions)){
    CVs <- CVs %>%
      filter(R.Condition %in% conditions)
  }
  
  if(reorder){
    order_dt <- metadata %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    CVs$R.Condition <- factor(CVs$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  if(remove_outliers){
    CVs = CVs %>%
      filter(cvs < quantile(CVs$cvs, 0.95,na.rm = T))
  }
  
  if(plot_type == "boxplot"){
    
    G = ggplot(CVs, aes(x= R.Condition, y = cvs, fill = R.Condition))+
      geom_violin()+
      geom_boxplot(width = 0.25)+
      geom_hline(yintercept = 20,linetype=20)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90))+
      labs(y = "CVs (%)", x = "Condition")
  }
  if(plot_type == "histogram"){
    
    G = ggplot(CVs, aes(x=cvs, fill= R.Condition))+
      geom_histogram( color="black")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90))+
      facet_wrap(~R.Condition)+
      labs(x = "CVs (%)")
  }
  
  if(plot_type == "freqpoly"){
    
    G = ggplot(CVs, aes(x= cvs, color = R.Condition))+
      geom_freqpoly(binwidth = 2)+
      theme_bw()+
      labs(x = "CVs (%)")
    
    
  }
  
  # G= ggplotly(G)
  
  return(G)
  
}
plot_CVs_reordered <- function(CVs, dt, reorder ){
  # CVs <- dt %>%
  #   ungroup %>%
  #   group_by(R.Condition, EG.ModifiedPeptide, FG.Charge) %>%
  #   summarise(mean_quant = mean(FG.MS2Quantity, na.rm = T),
  #             sd_quant = sd(FG.MS2Quantity, na.rm =T)) %>%
  #   mutate(cvs = sd_quant / mean_quant*100)
  
  if(reorder){
  order_dt <- dt %>%
    dplyr::select(R.Condition, order) %>%
    distinct() %>%
    arrange(order)
  CVs$R.Condition <- factor(CVs$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  CVs <- CVs %>%
    filter(!is.na(cvs))
    
  plot_cvs1 = ggplot(CVs, aes(x= R.Condition, y = cvs, color = R.Condition))+
    geom_violin()+
    geom_boxplot(width = 0.25)+
    geom_hline(yintercept = 20,linetype=20)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    labs(y = "CVs (%)", x = "Condition")
  
  plot_cvs2 = ggplot(CVs, aes(x= cvs, color = R.Condition))+
    geom_freqpoly(binwidth = 2)+
    theme_bw()+
    theme(legend.position = "none")+
    labs(x = "CVs (%)")
  plot_cvs_legend = get_legend(ggplot(CVs, aes(x= R.Condition, y = cvs, fill = R.Condition))+
                           geom_violin()+
                           guides(fill=guide_legend(nrow=10, byrow=TRUE)))
  plot_cvs = cowplot::plot_grid(plot_cvs1,
                                plot_grid(plot_cvs_legend, plot_cvs2),
                                nrow=2, rel_heights = c(2,1))
  
  return(plot_cvs)
  
}
plot_CVs_reordered_simple <- function(CVs, dt, reorder ){
  # CVs <- dt %>%
  #   ungroup %>%
  #   group_by(R.Condition, EG.ModifiedPeptide, FG.Charge) %>%
  #   summarise(mean_quant = mean(FG.MS2Quantity, na.rm = T),
  #             sd_quant = sd(FG.MS2Quantity, na.rm =T)) %>%
  #   mutate(cvs = sd_quant / mean_quant*100)
  
  if(reorder){
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    CVs$R.Condition <- factor(CVs$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  CVs <- CVs %>%
    filter(!is.na(cvs))
  
  plot_cvs1 = ggplot(CVs, aes(x= R.Condition, y = cvs, color = R.Condition))+
    geom_violin()+
    geom_boxplot(width = 0.25)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    labs(y = "CVs (%)", x = "Condition")+
    geom_hline(yintercept = 20, linetype=2)
  
  return(plot_cvs1)
  
}

plot_CVs_counts_Precursor <- function(CVs,dt = dt, reorder = T, conditions = NULL,
                                      type = c("Identified","Quantified","CV < 40","CV < 20","CV < 10")){
  
  if(length(conditions)>=1) {
    CVs <- CVs %>%
      filter(R.Condition %in% conditions)
  }
  
  CVs2 = CVs %>%
    # filter(!is.na(cvs)) %>%
    ungroup() %>%
    group_by(R.Condition) %>%
    summarise(Identified = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)),
              Quantified = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[!is.na(cvs)]),
              CV40 = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[cvs<40]),
              CV20 = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[cvs<20]),
              CV10 = n_distinct(paste0(EG.ModifiedPeptide, FG.Charge)[cvs<10])) %>%
    gather(key = "CV_cutoff", value = "value", which(str_detect(names(.),pattern = "CV|Quantified|Identified"))) %>%
    mutate(CV_cutoff = gsub(CV_cutoff, pattern = "CV",replacement = "CV < " ))
  
  CVs2$CV_cutoff <- factor(CVs2$CV_cutoff, levels = unique(CVs2$CV_cutoff))
  
  if(length(type)>=1) {
    CVs2 <- CVs2 %>%
      filter(CV_cutoff %in% type)
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(value>max(value)/3, 1.1, -0.1))
  }
  
  CVs2 <- CVs2 %>%
    ungroup() %>%
    group_by(R.Condition) %>% 
    group_modify(~hjust_adjustment(.x))
  
  CVs2_plot = ggplot(CVs2, aes(x= CV_cutoff, y = value, fill = CV_cutoff))+
    geom_bar(stat = "identity", color= "black")+
    geom_text(aes(label=round(value,0)#, y = ifelse(value>40000,value,value+40000)
    ), hjust=CVs2$hjust_value, angle=90,size=3.5, color= "black")+
    theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=90),
          axis.title.x = element_blank())+
    facet_grid(~R.Condition)+
    scale_fill_manual(values = c("#004E82","#0071BC","#E4EFF7", "#83929C", "#999999"))
  
  return(CVs2_plot)
  
}

plot_CVs_counts_Protein <- function(CVs_prot, dt = dt, reorder = T, conditions = NULL,
                                    type = c("Identified","Quantified","CV < 40","CV < 20","CV < 10")){
  
  if(length(conditions)>=1) {
    CVs_prot <- CVs_prot %>%
      filter(R.Condition %in% conditions)
  }
  
  CVs_count = CVs_prot %>%
    # filter(!is.na(cvs)) %>%
    ungroup() %>%
    group_by(R.Condition) %>%
    summarise(Identified = n_distinct(PG.ProteinGroups),
              Quantified = n_distinct(PG.ProteinGroups[!is.na(cvs)]),
              CV40 = n_distinct(PG.ProteinGroups[cvs<40]),
              CV20 = n_distinct(PG.ProteinGroups[cvs<20]),
              CV10 = n_distinct(PG.ProteinGroups[cvs<10])) %>%
    gather(key = "CV_cutoff", value = "value", which(str_detect(names(.),pattern = "CV|Quantified|Identified"))) %>%
    mutate(CV_cutoff = gsub(CV_cutoff, pattern = "CV",replacement = "CV < " ))
  
  CVs_count$CV_cutoff <- factor(CVs_count$CV_cutoff, levels = unique(CVs_count$CV_cutoff))
  
  if(length(type)>=1) {
    CVs_count <- CVs_count %>%
      filter(CV_cutoff %in% type)
  }
  
  hjust_adjustment = function(data){
    data %>%
      mutate(hjust_value = ifelse(value>max(value)/3, 1.1, -0.1))
  }
  
  CVs_count <- CVs_count %>%
    ungroup() %>%
    group_by(R.Condition) %>% 
    group_modify(~hjust_adjustment(.x))
  
  CVs_count_plot = ggplot(CVs_count, aes(x= CV_cutoff, y = value, fill = CV_cutoff))+
    geom_bar(stat = "identity", color= "black")+
    geom_text(aes(label=round(value,0)
    ), hjust=CVs_count$hjust_value, angle=90,size=3.5, color= "black")+
    theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=90),
          axis.title.x = element_blank())+
    facet_grid(~R.Condition)+
    scale_fill_manual(values = c("#004E82","#0071BC","#E4EFF7", "#83929C", "#999999"))
  
  return(CVs_count_plot)
  
}



plot_data_completeness <- function(dt){
  
  data.completeness = dt %>%
    filter(!is.na(FG.MS2Quantity))%>%
    group_by(PG.ProteinGroups) %>%
    summarise(n_runs = n_distinct(R.FileName)) %>%
    ungroup() %>%
    group_by(n_runs) %>%
    summarise(n_prots = n_distinct(PG.ProteinGroups)) %>%
    mutate(percent_prots = paste0(round(n_prots/length(unique(dt$PG.ProteinGroups))*100,1),"%"))
  
  a = ggplot(data.completeness, aes(x=n_runs, y= n_prots, label = percent_prots))+
    geom_bar(stat = "identity")+
    labs(x= "# of Proteins observed in X runs", y= "# Proteins")+
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle=90))+
    scale_x_continuous(breaks = seq(0,max(data.completeness$n_runs),1))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    geom_text(mapping = aes(x=n_runs, y= n_prots),nudge_y = max(data.completeness$n_prots)*0.10, angle= 90)
  
  
  total_number_of_runs =length(unique(dt$R.FileName))
  
  data.completeness = dt %>%
    filter(!is.na(FG.MS2Quantity))%>%
    group_by(PG.ProteinGroups) %>%
    summarise(n_runs = n_distinct(R.FileName[PG.Qvalue<0.01]),
              percent_runs = n_distinct(R.FileName)/total_number_of_runs*100) %>%
    ungroup() %>%
    arrange(-percent_runs) %>%
    mutate(rank = row_number())
  
  b = ggplot(data.completeness, aes(x=rank, y= percent_runs))+
    geom_line(stat = "identity")+
    labs(x= "# of Proteins observed", y= "% Ms runs")+
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle=90))+
    scale_x_continuous(breaks = scales::breaks_pretty(10))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       breaks = seq(0,100,10))
  
  plot_grid(a,b, ncol=1)
  
}


plot_TICS <- function(dt){
  dt2 <- dt%>%
    mutate(rep_id = paste0(R.Condition,"_", R.Replicate)) %>%
    dplyr::select(FG.MS2RawQuantity, EG.ApexRT, rep_id) %>%
    mutate(rt_bin = cut_interval(EG.ApexRT,length = 0.5)) %>%
    # mutate(rt_bin = cut_number(EG.ApexRT,n= 30)) %>%
    group_by(rt_bin,rep_id) %>%
    summarise(Int = sum(FG.MS2RawQuantity)) %>%
    mutate(RT = as.numeric(substr(rt_bin,start = 2, stop = str_locate(rt_bin, pattern = ",")-1)))
  
  a = ggplot(dt2, aes(x= RT, y = Int, color= rep_id, group = rep_id))+
    geom_line(size=1)+
    theme_bw()+
    # theme(legend.position = "none")+
    ggtitle("TIC of identified peptides")
  
  dt_sum = dt2 %>%
    ungroup() %>%
    group_by(rep_id) %>%
    summarise(Total_intensity = sum(Int))
  
  a2 = ggplot(dt_sum , aes(x = rep_id, y = log(Total_intensity), fill= rep_id))+
    geom_bar(stat = "identity")+
    # geom_text(aes(label=format(x = log(Total_intensity), digits = 3)),
    #           hjust=1.5, angle=90, color= "white")+
    theme_bw()+
    # theme(legend.position = "none")+
    ggtitle("Total intensity")
  
  G = subplot(ggplotly(a),ggplotly(a2), widths = c(0.7,0.3))
  return(G)
}


plot_SequenceMotif_by_PeptideLength <- function(dt, Peptide_Length, conditions = NULL){
  
  if(!is.null(conditions)){
    dt <- dt %>%
      filter(R.Condition %in% conditions)
  }
  
  dt2 <- dt %>%
    filter(PEP.PeptideLength == Peptide_Length)
  
  ggseqlogo::ggseqlogo(dt2$PEP.StrippedSequence, seq_type='aa')
  
  dt_list = dt2 %>%
    select(PEP.StrippedSequence, R.Condition) %>% 
    distinct() %>%
    group_by(R.Condition)%>% 
    group_map(~.x)
  
  names(dt_list) <- dt2 %>%
    select(PEP.StrippedSequence, R.Condition) %>% 
    distinct() %>%
    group_by(R.Condition)%>% 
    group_map(~.y) %>%
    unlist()
  
  dt_list <- lapply(dt_list, function(L){
    L$PEP.StrippedSequence})
  
  G = ggseqlogo(dt_list,seq_type='aa', font="roboto_bold")
  
  return(G)
  
}

#### QCs
#### RTs

plot_median_RT <- function(dt){

  dt2 <- dt %>%
    mutate(run = paste0(R.Condition,"_",R.Replicate))
  
  df_median = dt2 %>%
    group_by(run, R.Condition) %>%
    summarise(EG.ApexRT = median(EG.ApexRT, na.rm = T))
  
  ggplot(dt2, aes(x= run, y= EG.ApexRT, fill= R.Condition)) + 
    geom_violin()+
    geom_boxplot(width=0.25)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))+
    geom_line(data = df_median, 
              mapping = aes(x = run, y = EG.ApexRT, group=1))
  
} 

plot_median_RT_only_FullObservations <- function(dt){
  n_runs = length(unique(dt$R.FileName))
  
  Peptides_to_keep <- dt %>%
    group_by(EG.ModifiedPeptide, FG.Charge) %>%
    summarise(num_runs = n_distinct(R.FileName)) %>%
    filter(num_runs > n_runs*0.75) %>%
    mutate(keep = "keep")
  
  dt2 <- dt %>%
    left_join(Peptides_to_keep) %>%
    filter(keep == "keep") %>%
    mutate(run = paste0(R.Condition,"_",R.Replicate))
  
  df_median = dt2 %>%
    group_by(run, R.Condition) %>%
    summarise(EG.ApexRT = median(EG.ApexRT, na.rm = T))
  
  ggplot(dt2, aes(x= run, y= EG.ApexRT, fill= R.Condition)) + 
    geom_violin()+
    geom_boxplot(width=0.25)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))+
    geom_line(data = df_median, 
              mapping = aes(x = run, y = EG.ApexRT, group=1))
  
} 

plot_median_RT_2 <-function(dt, only_Full_observations = F){
  
  G = plot_median_RT(dt)
  
  if(only_Full_observations){
    G = plot_median_RT_only_FullObservations(dt)
  }
  return(G)
}
##### correlations plot

Longformat_to_matrixformat_precursors_OneMetric<- function(dt, metadata,
                                                           metric_to_plot =c("PG.MS2Quantity","FG.MS2Quantity","EG.ApexRT","EG.IonMobility")){
  
  names(dt)[names(dt) == metric_to_plot] <- 'Column_of_interest'
  
  if(str_detect(metric_to_plot, pattern = "Quantity")){
    dt <- dt %>% 
      mutate(Column_of_interest =log(Column_of_interest,10))
  }
  
  if(str_detect(metric_to_plot, pattern = "PG.MS2Quantity")){
    mat <- dt %>%
      dplyr::select(R.Condition, R.Replicate, PG.ProteinGroups, PG.ProteinAccessions,PG.ProteinNames, Column_of_interest) %>%
      group_by(R.Condition, R.Replicate, PG.ProteinGroups, PG.ProteinAccessions,PG.ProteinNames) %>%
      mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
      distinct() %>%
      ungroup()%>%
      select(-R.Condition, -R.Replicate) %>%
      tidyr::spread(key = run, value = Column_of_interest) %>%
      ungroup()
    
  } else {
    mat <- dt %>%
      dplyr::select(R.Condition, R.Replicate, EG.ModifiedPeptide, FG.Charge, Column_of_interest) %>%
      group_by(R.Condition, R.Replicate, EG.ModifiedPeptide, FG.Charge) %>%
      mutate(run = paste0(R.Condition,"_", R.Replicate)) %>%
      distinct() %>%
      ungroup()%>%
      select(-R.Condition, -R.Replicate) %>%
      tidyr::spread(key = run, value = Column_of_interest) %>%
      ungroup()
  }
  
  return(mat)
  
}



ggally_hexbin_SV <- function (data, mapping,  ...)  {
  p <- ggplot(data = data, mapping = mapping) +
    geom_hex(...)+
    theme_bw()+
    geom_abline(slope = 1, intercept = 0, color="black")+
    scale_fill_viridis(option="magma")
  p
}
upper_corr_plot <- function(data, mapping, method="p", use="pairwise", ...){
  
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  

  corr <- cor(x, y, method=method, use=use)
  
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  GGally::ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

# plot_correlations <- function(dt, conditions = NULL, metadata,
#                               metric_to_plot = c("FG.MS2Quantity","EG.ApexRT","EG.IonMobility"),
#                               fill_corrUpper = F){
#   
#   if(!is.null(conditions)){
#     dt <- dt %>%
#       filter(R.Condition %in% conditions)
#   }
#   
#   
#   mat_precs_RT = Longformat_to_matrixformat_precursors_OneMetric(dt, metadata,metric_to_plot)
#   mat_precs_RT = do.call(data.frame,lapply(mat_precs_RT, function(x) replace(x, is.infinite(x),NA)))
#   
#   mat_precs_RT = mat_precs_RT %>%
#     select(!which(names(mat_precs_RT) %in% c("R.Condition", "R.Replicate", "EG.ModifiedPeptide","FG.Charge")))
#   
#   if(fill_corrUpper){
#     GGally::ggpairs(mat_precs_RT, 
#                     upper= list(continuous = upper_corr_plot),
#                     ## use 'ggally_hexbin_SV' for continuous × continuous plots
#                     lower = list(continuous = "hexbin_SV",
#                                  ## use default plots for all other variable types
#                                  combo = "facethist", discrete = "facetbar", na =  "na"))+
#       # theme_bw()+
#       ggtitle(metric_to_plot)
#   } else {
#     GGally::ggpairs(mat_precs_RT, 
#                     ## use 'ggally_hexbin_SV' for continuous × continuous plots
#                     lower = list(continuous = "hexbin_SV",
#                                  ## use default plots for all other variable types
#                                  combo = "facethist", discrete = "facetbar", na =  "na"))+
#       theme_bw()+
#       ggtitle(metric_to_plot)
#   }
#   
#   
#   
# }


create_list_matrices_for_each_metric_2 <- function(list_matrices_for_each_metric_2,
                                                   dt,
                                                   metadata,
                                                   software,
                                                   metric_to_plot){
  
  dt2 = data_for_one_metric(dt, metadata, software, metric_to_plot)
  mat_precs = Longformat_to_matrixformat_precursors_OneMetric(dt2, metadata,metric_to_plot)
  mat_precs = do.call(data.frame,lapply(mat_precs, function(x) replace(x, is.infinite(x),NA)))
  
  list_matrices_for_each_metric_2 [[metric_to_plot]] = mat_precs
  
  return(list_matrices_for_each_metric_2)
}


plot_correlations_fromList <- function(list_matrices_for_each_metric_2,
                                       conditions = NULL,
                                       metadata,
                                       metric_to_plot = c("PG.MS2Quantity", "FG.MS2Quantity","EG.ApexRT","EG.IonMobility")){
  
  dt_2 = list_matrices_for_each_metric_2[[metric_to_plot]]
  
  
  if(!is.null(conditions)){
    
    selected_columns2 = which(str_detect(names(dt_2), pattern = paste0(conditions,collapse = "|")))
    selected_columns1 = which(str_detect(names(dt_2), pattern = paste0(c( "PG.ProteinGroups", "PG.ProteinAccessions","PG.ProteinNames","EG.ModifiedPeptide","FG.Charge"),collapse = "|")))
    selected_columns = c(selected_columns1, selected_columns2)
    
    dt_2 <- dt_2 %>%
      dplyr::select(all_of(selected_columns))
  }
  
  
  
  dt_2 = dt_2 %>%
    select(!which(names(dt_2) %in% c("PG.ProteinGroups", "PG.ProteinAccessions","PG.ProteinNames","R.Condition", "R.Replicate", "EG.ModifiedPeptide","FG.Charge")))
  
  GGally::ggpairs(dt_2, 
                  upper= list(continuous = upper_corr_plot),
                  ## use 'ggally_hexbin' for continuous × continuous plots
                  lower = list(continuous = "hexbin_SV",
                               ## use default plots for all other variable types
                               combo = "facethist", discrete = "facetbar", na =  "na"))+
    # theme_bw()+
    ggtitle(metric_to_plot)
  
}

#### PCAs

create_PCA_data <- function(list_matrices_for_each_metric_2, metadata){
  
  metadata2 = metadata %>% mutate(run = paste0(R.Condition,"_", R.Replicate))
  
  dt_1 = list_matrices_for_each_metric_2[["PG.MS2Quantity"]] 
  dt_2 = dt_1 %>%
    arrange(PG.ProteinGroups) %>%
    select(which(names(.) %in% metadata2$run))
  
  num_noNAs = apply(dt_2,MARGIN = 1,FUN = function(x) {sum(!is.na(x))})
  keep_rows = num_noNAs >= dim(dt_2)[2]*0.1
  
  dt_2 = dt_2[keep_rows,]
  
  imputation_value = quantile(dt_2,probs = 0.01,na.rm = T)
  
  dt_2 = do.call(data.frame,lapply(dt_2, function(x) replace(x, is.infinite(x),NA)))
  dt_2 = do.call(data.frame,lapply(dt_2, function(x) replace(x, is.na(x),imputation_value)))
  
  
  pca.data <- FactoMineR::PCA(t(dt_2), scale.unit = TRUE, graph = FALSE)
  
  return(pca.data)
}
plot_PCA <- function(pca.data, plot_type = c("scree", "pca"), metadata){
  
  metadata2 = metadata %>% mutate(run = paste0(R.Condition,"_", R.Replicate))
  
  # factoextra::fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))
  # 
  # 
  # factoextra::fviz_pca_ind(pca.data, col.ind = "cos2", 
  #              # gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
  #              repel = TRUE)
  # 
  # 
  # names_points= names(pca.data$ind$dist)
  # 
  # factoextra::fviz_pca_ind(pca.data, 
  #                          col.ind = names_points,
  #                          repel = TRUE,
  #                          addEllipses = TRUE)
  if(plot_type == "scree"){
    pca_Scree_plot = data.frame(pca.data$eig, Dimensions = rownames(pca.data$eig))
    pca_Scree_plot$Dimensions = factor(pca_Scree_plot$Dimensions, levels = unique(pca_Scree_plot$Dimensions))
    G = ggplot(pca_Scree_plot, aes(x = Dimensions, y = percentage.of.variance, group =1)) +
      geom_bar(stat = "identity", fill = "blue")+
      geom_line()+
      geom_point()+
      # geom_text(aes(label=paste0(round(percentage.of.variance,0),"%")), vjust=-0.25, hjust = -0.25,  color= "black")+
      theme_bw()+
      coord_cartesian(ylim = c(0, 100))
  }
  
  if(plot_type == "pca"){
    pca_for_plot <- data.frame(pca.data$ind$coord, run = rownames(pca.data$ind$coord)) %>%
      left_join(metadata2)
    
    # ggplot(pca_for_plot, aes(x = Dim.1, y = Dim.2 , color = run)) +
    #   geom_point(size = 3) +
    #   xlab("Principal Component 1") +
    #   ylab("Principal Component 2") +
    #   ggtitle("PCA plot")+
    #   theme_bw()+
    #   geom_hline(yintercept = 0)+
    #   geom_vline(xintercept = 0)
    
    
    G = ggplot(pca_for_plot, aes(x = Dim.1, y = Dim.2 ,color = R.Condition)) +
      geom_point(size = 3) +
      xlab("Principal Component 1") +
      ylab("Principal Component 2") +
      ggtitle("PCA plot")+
      theme_bw()+
      theme(panel.grid = element_blank())+
      geom_hline(yintercept = 0, color = "grey75")+
      geom_vline(xintercept = 0, color = "grey75")
  }
  
  ggplotly(G)
  
}



#PTMs
find_modifications <- function(dt, software){
  
  list_of_mods <- function(v){
    x = NULL
    for (i in 1:length(v)) {
      x[i] = str_match_all(v[i],
                           pattern = "\\$\\s*(.*?)\\s*\\#")[[1]][,2] %>%
        paste0(collapse = ",")
    }
    
    mods = data.frame(x = x) %>%
      separate_rows(x, sep = ",") %>%
      unique()
    
    return(mods)
  }
  print(software)
  if(software %in% c("timsDIANN", "DIANN")){
    Mods = dt %>%
      ungroup() %>%
      dplyr::select(EG.ModifiedPeptide) %>%
      distinct %>%
      mutate(EG.ModifiedPeptide = gsub(EG.ModifiedPeptide, pattern = "\\(", replacement = "$")) %>%
      mutate(EG.ModifiedPeptide = gsub(EG.ModifiedPeptide, pattern = "\\)", replacement = "#")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = "\\$"))
    }
  
  if(software == "Spectronaut"){
    Mods = dt %>%
      ungroup() %>%
      dplyr::select(EG.ModifiedPeptide) %>%
      distinct %>%
      mutate(EG.ModifiedPeptide = gsub(EG.ModifiedPeptide, pattern = "\\[", replacement = "$")) %>%
      mutate(EG.ModifiedPeptide = gsub(EG.ModifiedPeptide, pattern = "\\]", replacement = "#")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = "\\$"))%>%
      mutate(EG.ModifiedPeptide = gsub(EG.ModifiedPeptide, pattern = "_",replacement = ""))
  }
  
  Mods = list_of_mods(Mods$EG.ModifiedPeptide)
  Mods = Mods$x
  return(Mods)
}

PTM_counts_list_O_plots <- function(dt, reorder = F, software){
  mods = find_modifications(dt, software)
  list_plots = list()
  
  for (i in 1:length(mods)) {
    
    mod_string =str_replace_all(mods[i], "[^[:alnum:]]", "") 
    
    dt_mod = dt %>%
      mutate(EG.ModifiedPeptide = str_replace_all(EG.ModifiedPeptide, "[^[:alnum:]]", "")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = mod_string))
    
    counts_mod = create_counts_table(dt_mod)
    
    list_plots [[i]] =  plot_count_pepts_prots_precs_withErrorBars_reordered(counts_mod, dt_mod, reorder)+
      ggtitle(mods[i])
  }
  
  # plot_grid(plotlist=list_plots, ncol = 1)
  
  names(list_plots) <- mods
  list_plots
}

PTM_enrichment_yield_list_O_plots <- function(dt, counts,software,  reorder = F){
  mods = find_modifications(dt, software)
  list_plots = list()
  
  counts_all = counts %>%
    rename(counts_all = counts,
           mean_all = mean) %>%
    dplyr::select(-sd)
  
  for (i in 1:length(mods)) {
    mod_string =str_replace_all(mods[i], "[^[:alnum:]]", "") 
    
    dt_mod = dt %>%
      mutate(EG.ModifiedPeptide = str_replace_all(EG.ModifiedPeptide, "[^[:alnum:]]", "")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = mod_string))
    counts_mod = create_counts_table(dt_mod)
    
    counts_mod <- counts_mod %>%
      left_join(counts_all, by = c("R.Condition", "type")) %>%
      mutate(enrich_yield = counts/counts_all *100)
    
    if(reorder){
      order_dt <- dt %>%
        dplyr::select(R.Condition, order) %>%
        distinct() %>%
        arrange(order)
      counts_mod$R.Condition <- factor(counts_mod$R.Condition, levels = unique(order_dt$R.Condition))
    }
    
    
    hjust_adjustment = function(data){
      data %>%
        mutate(hjust_value = ifelse(enrich_yield>50, 1.5, -0.1))
    }
    
    counts_mod <- counts_mod %>%
      ungroup() %>%
      group_by(type) %>% 
      group_modify(~hjust_adjustment(.x))
    
    
    enrich_yield_plot = ggplot(counts_mod, aes(x= R.Condition, y = enrich_yield, fill = R.Condition))+
      geom_bar(stat = "identity", color= "black")+
      geom_text(aes(label=paste0(round(enrich_yield,1),"%")), hjust=counts_mod$hjust_value, angle=90, color= "black")+
      theme_bw()+
      theme(legend.position = "none", axis.text.x = element_text(angle=90))+
      coord_cartesian(ylim = c(0, 100))+
      facet_wrap(~type, scales = "free_y")+
      ggtitle(mods[i])
    
    list_plots[[i]] = enrich_yield_plot 
    
  }
  
  # plot_grid(plotlist=list_plots, ncol = 1)
  names(list_plots) <- mods
  list_plots
}
PTM_counts_list_O_plots_filter <- function(dt, reorder = F,plot_type = NULL, software){
  mods = find_modifications(dt, software)
  list_plots = list()
  
  for (i in 1:length(mods)) {
    
    mod_string =str_replace_all(mods[i], "[^[:alnum:]]", "") 
    
    dt_mod = dt %>%
      mutate(EG.ModifiedPeptide = str_replace_all(EG.ModifiedPeptide, "[^[:alnum:]]", "")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = mod_string))
    
    counts_mod = create_counts_table(dt_mod)
    
    list_plots [[i]] =  plot_count_pepts_prots_precs_withErrorBars_reordered_filter(counts_mod, dt_mod, plot_type, reorder)+
      ggtitle(mods[i])
  }
  
  # plot_grid(plotlist=list_plots, ncol = 1)
  
  names(list_plots) <- mods
  list_plots
}
PTM_enrichment_yield_list_O_plots_filter <- function(dt, counts,software, plot_type = NULL, reorder = F){
  mods = find_modifications(dt, software)
  list_plots = list()
  
  counts_all = counts %>%
    rename(counts_all = counts,
           mean_all = mean) %>%
    dplyr::select(-sd)
  
  for (i in 1:length(mods)) {
    mod_string =str_replace_all(mods[i], "[^[:alnum:]]", "") 
    
    dt_mod = dt %>%
      mutate(EG.ModifiedPeptide = str_replace_all(EG.ModifiedPeptide, "[^[:alnum:]]", "")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = mod_string))
    counts_mod = create_counts_table(dt_mod)
    
    counts_mod <- counts_mod %>%
      left_join(counts_all, by = c("R.Condition", "type")) %>%
      mutate(enrich_yield = counts/counts_all *100)
    
    if(reorder){
      order_dt <- dt %>%
        dplyr::select(R.Condition, order) %>%
        distinct() %>%
        arrange(order)
      counts_mod$R.Condition <- factor(counts_mod$R.Condition, levels = unique(order_dt$R.Condition))
    }
    
    
    hjust_adjustment = function(data){
      data %>%
        mutate(hjust_value = ifelse(enrich_yield>50, 1.5, -0.1))
    }
    
    counts_mod <- counts_mod %>%
      ungroup() %>%
      group_by(type) %>% 
      group_modify(~hjust_adjustment(.x))
    
    if(!is.null(plot_type)){
      counts_mod <- counts_mod %>%
        filter(type %in% plot_type)
    }
    
    
    enrich_yield_plot = ggplot(counts_mod, aes(x= R.Condition, y = enrich_yield, fill = R.Condition))+
      geom_bar(stat = "identity", color= "black")+
      geom_text(aes(label=paste0(round(enrich_yield,1),"%")), hjust=counts_mod$hjust_value, angle=90, color= "black")+
      theme_bw()+
      theme(legend.position = "none", axis.text.x = element_text(angle=90))+
      coord_cartesian(ylim = c(0, 100))+
      facet_wrap(~type, scales = "free_y")+
      ggtitle(mods[i])
    
    list_plots[[i]] = enrich_yield_plot 
    
  }
  
  # plot_grid(plotlist=list_plots, ncol = 1)
  names(list_plots) <- mods
  list_plots
}

PTM_CVs_list_O_plots<- function(CVs, dt, software, reorder = F){
  mods = find_modifications(CVs, software)
  list_plots = list()
  
  for (i in 1:length(mods)) {
    
    mod_string =str_replace_all(mods[i], "[^[:alnum:]]", "") 
    
    CVs_mod = CVs %>%
      filter(!is.na(cvs)) %>%
      mutate(EG.ModifiedPeptide = str_replace_all(EG.ModifiedPeptide, "[^[:alnum:]]", "")) %>%
      filter(str_detect(EG.ModifiedPeptide, pattern = mod_string))
    
    list_plots [[i]] =  plot_CVs_reordered_simple(CVs_mod, dt, reorder)+
      ggtitle(mods[i])
  }
  
  # plot_grid(plotlist=list_plots, ncol = 1)
  names(list_plots) <- mods
  list_plots
  
}

PTM_plot_counts <- function(list_plots_ptm, mod_string){
  list_plots_ptm[[mod_string]]
}

PTM_plot_enrichmentYield <- function(list_plots_ptm_enrich, mod_string){
  list_plots_ptm_enrich[[mod_string]]
}

PTM_plot_CVs<- function(list_plots_ptm_CVs, mod_string){
  list_plots_ptm_CVs[[mod_string]]
}




## Calibration curves
cal_curve_norm_to_single_condition <- function(dt, normalizing_condition,reorder){
  
  dt <- dt %>%
    filter(!is.na(FG.MS2Quantity))%>%
    mutate(FG.MS2Quantity = as.numeric(FG.MS2Quantity),
           Concentration = as.numeric(Concentration))
  
  
  condition_has_normCondition <- dt %>%
    dplyr::select(R.Condition,PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge ) %>%
    filter(R.Condition == normalizing_condition) %>%
    dplyr::select(-R.Condition) %>%
    mutate(keep= "keep") %>%
    distinct()
  
  
  mean_areas <- dt %>%
    filter(!is.na(FG.MS2Quantity))%>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge, R.Condition) %>%
    summarise(mean = mean(FG.MS2Quantity, na.rm =T),
              conc = mean(Concentration, na.rm =T),
              num_replicatesObserved = n_distinct(R.FileName)) %>%
    ungroup() %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge) %>%
    mutate(ratio = mean/mean[R.Condition == normalizing_condition],
           ratio_theo = conc/conc[R.Condition == normalizing_condition]) %>%
    ungroup()
  
  num_precursors <- dt %>%
    filter(!is.na(FG.MS2Quantity))%>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(R.Condition) %>%
    summarise(num_precursors = n_distinct(paste0(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge)))
  
  expected_ratios = mean_areas %>%
    ungroup() %>%
    dplyr::select(R.Condition, ratio_theo ) %>%
    distinct() %>%
    rename(ratio = ratio_theo)
  
  if(reorder){
    order_dt <- dt %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    mean_areas$R.Condition <- factor(mean_areas$R.Condition, levels = unique(order_dt$R.Condition))
    expected_ratios$R.Condition <- factor(expected_ratios$R.Condition, levels = unique(order_dt$R.Condition))
    num_precursors$R.Condition <- factor(num_precursors$R.Condition, levels = unique(order_dt$R.Condition))
  }
  
  a = ggplot(mean_areas , aes(x =R.Condition, y= ratio, color= R.Condition ))+
    geom_violin()+
    geom_boxplot(width = 0.25)+
    # coord_cartesian(ylim = c(0,6))+
    geom_point(data = expected_ratios,aes(x =R.Condition, y= ratio), color= "black" )+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),
          legend.position = "none")+
    scale_y_log10()
  
  b = ggplot(num_precursors, aes(x= R.Condition, y = num_precursors, fill = R.Condition))+
    geom_bar(stat = "identity", color= "black")+
    geom_text(aes(label=round(num_precursors,0)), hjust=-00.1, angle=90, color= "black")+
    theme_bw()+
    ggtitle(paste0("Ratios compared to condition ", normalizing_condition))+
    theme(legend.position = "none", axis.text.x = element_text(angle=90),
          axis.title.x = element_blank())+
    coord_cartesian(ylim = c(0, max(num_precursors$num_precursors)+max(num_precursors$num_precursors)/2))
  
  
  G = cowplot::plot_grid(b,a,ncol=1, align = "v" , rel_heights = c(1,2))
  
  return(G)
  
}


cal_curve_norm_to_single_condition_singlePept <- function(dt, normalizing_condition, reorder, PepSeq){
  
  dt2 <- dt %>%
    filter(EG.ModifiedPeptide == PepSeq) %>%
    filter(!is.na(FG.MS2Quantity))%>%
    mutate(FG.MS2Quantity = as.numeric(FG.MS2Quantity),
           Concentration = as.numeric(Concentration))
  
  
  condition_has_normCondition <- dt2 %>%
    dplyr::select(R.Condition,PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge ) %>%
    filter(R.Condition == normalizing_condition) %>%
    dplyr::select(-R.Condition) %>%
    mutate(keep= "keep") %>%
    distinct()
  
  
  mean_areas <- dt2 %>%
    filter(!is.na(FG.MS2Quantity))%>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge, R.Condition) %>%
    summarise(mean = mean(FG.MS2Quantity, na.rm =T),
              sd = sd(FG.MS2Quantity, na.rm =T),
              conc = mean(Concentration, na.rm =T),
              num_replicatesObserved = n_distinct(R.FileName)) %>%
    mutate(cv = sd/mean *100) %>%
    ungroup() %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge) %>%
    mutate(ratio = mean/mean[R.Condition == normalizing_condition],
           ratio_theo = conc/conc[R.Condition == normalizing_condition]) %>%
    ungroup()
  
  num_precursors <- dt2 %>%
    filter(!is.na(FG.MS2Quantity))%>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(R.Condition) %>%
    summarise(num_precursors = n_distinct(paste0(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge)))
  
  expected_ratios = mean_areas %>%
    ungroup() %>%
    dplyr::select(R.Condition, ratio_theo ) %>%
    distinct() %>%
    rename(ratio = ratio_theo)
  
  if(reorder){
    order_dt <- dt2 %>%
      dplyr::select(R.Condition, order) %>%
      distinct() %>%
      arrange(order)
    mean_areas$R.Condition <- factor(mean_areas$R.Condition, levels = unique(order_dt$R.Condition))
    expected_ratios$R.Condition <- factor(expected_ratios$R.Condition, levels = unique(order_dt$R.Condition))
    num_precursors$R.Condition <- factor(num_precursors$R.Condition, levels = unique(order_dt$R.Condition))
    dt2$R.Condition <- factor(dt2$R.Condition, levels = unique(order_dt$R.Condition))
    
  }
  
  # areas_plot = ggplot(mean_areas , aes(y =mean, x= conc, color = as.character(FG.Charge), group=FG.Charge ))+
  #   # geom_line()+
  #   geom_point(size= 2)+
  #   geom_smooth(method = "lm",se = F, size=0.7, linetype= 2)+
  #   theme_bw()+
  #   theme(axis.text.x = element_text(angle=90),
  #         legend.position = "none")+
  #   scale_x_log10()+
  #   scale_y_log10()+
  #   labs(x= "Concentration", y="Measured Area")+
  #   ggtitle(PepSeq)
  
  mean_areas <- mean_areas %>%
    rename(Concentration = conc,
           FG.MS2Quantity = mean)
  areas_plot2 = ggplot(dt2 , aes(y =FG.MS2Quantity, x= Concentration, color = as.character(FG.Charge), group=FG.Charge ))+
    # geom_line()+
    geom_point(size= 2)+
    geom_point(data = mean_areas,aes(y =FG.MS2Quantity, x= Concentration), color = "black", shape=3, size= 3)+
    geom_smooth(method = "lm",se = F, size=0.7, linetype= 2,color = "black")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),
          legend.position = "none")+
    scale_x_log10()+
    scale_y_log10()+
    labs(x= "Concentration", y="Measured Area")+
    ggtitle(PepSeq)
  
  
  
  
  
  charges <- mean_areas %>% dplyr::select(R.Condition,FG.Charge) %>% distinct()
  
  expected_ratios<- expected_ratios %>% left_join(charges)
  # ratio_plot = ggplot(mean_areas , aes(x =R.Condition, y= ratio,color = as.character(FG.Charge), group = FG.Charge))+
  #   geom_point()+
  #   geom_line()+
  #   geom_point(data = expected_ratios,aes(x =R.Condition, y= ratio), color= "black", shape= 3)+
  #   theme_bw()+
  #   ggtitle(paste0("Ratios compared to condition ", normalizing_condition))+
  #   theme(axis.text.x = element_text(angle=90),
  #         legend.position = "none",
  #         axis.title.x = element_blank())+
  #   scale_y_log10()
  ratio_plot2 = ggplot(mean_areas , aes(x= ratio_theo, y =ratio,  color = as.character(FG.Charge), group=FG.Charge ))+
    # geom_line()+
    geom_point(size= 2)+
    geom_smooth(method = "lm",se = F, size=0.7, linetype= 1, color="black")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),
          legend.position = "none")+
    labs(x= "Expected ratio", y="Measured ratio")+
    ggtitle(PepSeq)
  
  CV_plot = ggplot(mean_areas , aes(x =R.Condition, y= cv,color = as.character(FG.Charge), group = FG.Charge))+
    geom_point()+
    geom_line()+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),
          legend.position = "none",
          axis.title.x = element_blank())+
    geom_hline(yintercept = 20, linetype=2)+
    labs(x= "CV (%)")+
    ggtitle("CVs")
  
  num_replicatesObserved_plot = ggplot(mean_areas, aes(x= R.Condition, y = num_replicatesObserved,fill = as.character(FG.Charge)))+
    geom_bar(stat = "identity", color= "black")+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(angle=90))+
    facet_wrap(~FG.Charge)+
    ggtitle("Number of replicates observed")
  
  
  G = cowplot::plot_grid(areas_plot2, num_replicatesObserved_plot,ratio_plot2, CV_plot, ncol=2, align = "v" )
  
  return(G)
  
}

### GCT functions

Longformat_to_matrixformat <- function(dt, metadata){
  
  # metadata_columns = names(metadata)[-which(str_detect(names(metadata), "R.FileName"))]
  
  mat <- dt %>%
    dplyr::select(R.FileName, PG.ProteinGroups, EG.ModifiedPeptide, FG.Charge, FG.MS2Quantity) %>%
    # select(- which(names(.) %in% metadata_columns)) %>%
    # select(-PG.Qvalue, - EG.Qvalue) %>% 
    tidyr::spread(key = R.FileName, value = FG.MS2Quantity)
  
  return(mat)
  
}

Longformat_to_matrixformat_ProteinLevel<- function(dt, metadata){
  
  mat_prot <- dt %>%
    dplyr::select(R.FileName, PG.ProteinGroups, PG.ProteinAccessions, PG.ProteinNames,PG.Genes, PG.MS2Quantity ) %>%
    distinct() %>%
    tidyr::spread(key = R.FileName, value = PG.MS2Quantity)
  
  return(mat_prot)
  
}

create_gct_PG.Level <- function(dt, metadata){
  
  create_unique_sample_id <- function(metadata){
    
    metadata <- metadata %>%
      ungroup %>%
      distinct() %>%
      arrange(order, R.FileName) %>%
      mutate(index= row_number()) %>%
      mutate(sample_id = paste(R.Condition, R.Replicate, index,sep = "_"))
    
    return(metadata)
  } 
  create_unique_protein_id <- function(row_metadata){
    
    row_metadata <- row_metadata %>%
      ungroup %>%
      mutate(index= row_number()) %>%
      mutate(id = paste(PG.ProteinGroups, PG.Genes,index,sep = "_")) %>%
      dplyr::select(id, everything())
    
    return(row_metadata)
  } 
  
  metadata <- create_unique_sample_id(metadata)
  
  dt <- dt %>% left_join(metadata)
  
  mat_intermediate <- Longformat_to_matrixformat_ProteinLevel(dt)
  
  # mat
  mat <- mat_intermediate %>%
    select(which(names(.) %in% unique(dt$R.FileName))) %>%
    data.matrix()
  
  
  # column_metadata
  # col_metadata =  metadata%>%
  #   filter(R.FileName %in% unique(dt$R.FileName)) %>%
  #   select(id= sample_id, everything(),sample_id) %>%
  #   mutate(sample_id = id)
  # 
  # col_metadata = data.frame(R.FileName = colnames(mat)) %>%
  #   left_join(col_metadata)
  
  col_metadata =  metadata%>%
    filter(R.FileName %in% unique(dt$R.FileName)) %>%
    dplyr::select(id= sample_id, everything()) %>%
    mutate(sample_id = id)
  
  col_metadata = data.frame(R.FileName = colnames(mat)) %>%
    left_join(col_metadata) %>%
    dplyr::select(id, everything()) %>%
    mutate(sample_id = id)%>%
    rename(ii = id)
  
  #rename cid
  new_names = metadata %>% ungroup() %>%
    dplyr::select(R.FileName, sample_id)
  new_names = data.frame(R.FileName = colnames(mat)) %>%
    left_join(new_names)
  
  colnames(mat) = new_names$sample_id
  
  #row_metadata
  row_metadata <- mat_intermediate %>%
    dplyr::select(! which(names(.) %in% unique(dt$R.FileName))) %>%
    create_unique_protein_id()
  
  
  gct_dt = cmapR::GCT()
  gct_dt@mat = mat
  gct_dt@cid = colnames(mat)
  gct_dt@cdesc = col_metadata
  gct_dt@rid = row_metadata$id
  gct_dt@rdesc = row_metadata
  
  testthat::expect_equal(label = "GCT: same dimensions cdesc and cid",object = dim(gct_dt@cdesc)[1], expected = length(gct_dt@cid))
  testthat::expect_equal(label = "GCT: same dimensions rdesc and rid",object = dim(gct_dt@rdesc)[1], expected = length(gct_dt@rid))
  testthat::expect_equal(label = "GCT: uniques row ids",object = unique(length(gct_dt@rid)), expected = length(gct_dt@rid))
  testthat::expect_equal(label = "GCT: uniques col ids",object = unique(length(gct_dt@cid)), expected = length(gct_dt@cid))
  
  return(gct_dt)
  # cmapR::write_gct(gct_dt,
  #                  ofile= "C:/Users/Alvaro.Vaca_Jacome/OneDrive - Bruker Physik GmbH/Desktop/test_gct.gct",
  #                  precision = 3, appenddim = TRUE, ver = 3)
}


#### Targeted analysis
load_Targets <- function(Targets_file_path){
  
  Targets <- fread(Targets_file_path)
  
  return(Targets)
}

List_of_Observed_Proteins_n_Peptides <- function(dt){
  dt2 <- dt %>%
    dplyr::select(PG.ProteinGroups,PG.ProteinAccessions, PG.ProteinNames, PG.Genes, EG.ModifiedPeptide) %>%
    distinct()
  return(dt2)
}

plot_counts_Targets<- function(dt, Targets, reorder, pepts_or_prots = "pepts"){
  
  if(pepts_or_prots == "pepts") {
    list_of_targets = Targets$EG.ModifiedPeptide

    dt_mod = dt %>%
      filter(EG.ModifiedPeptide %in% list_of_targets)
    
    counts_mod = create_counts_table(dt_mod) %>%
      filter(!str_detect(type, "Prot"))
    
    g = plot_count_pepts_prots_precs_withErrorBars_reordered(counts_mod, dt_mod, reorder)
    
  }
  
  if(pepts_or_prots == "prots") {
    list_of_targets = Targets$PG.ProteinGroups

    dt_mod = dt %>%
      filter(PG.ProteinGroups %in% list_of_targets)
    
    counts_mod = create_counts_table(dt_mod)
    
    g= plot_count_pepts_prots_precs_withErrorBars_reordered(counts_mod, dt_mod, reorder)
  }
  
  return(g)
  
}

plot_HeatMap_Peptide_Targets<- function(dt, Targets){
  
  dt_mod = dt %>%
    filter(EG.ModifiedPeptide %in% Targets$EG.ModifiedPeptide) %>%
    mutate(FG.MS2RawQuantity = as.numeric(FG.MS2RawQuantity))%>%
    mutate(FG.MS2RawQuantity = log(FG.MS2RawQuantity,10))
  
  mat0 = Longformat_to_matrixformat(dt_mod)
  
  mat <- mat0 %>%
    select(-PG.ProteinGroups, -EG.ModifiedPeptide, -FG.Charge) %>%
    replace(is.na(.), 0) %>%
    data.matrix()
  
  mat0 <- mat0 %>%
    mutate(row_names= paste0(EG.ModifiedPeptide,"_",FG.Charge)) %>%
    select(row_names)
  row.names(mat) = mat0$row_names
  
  metadata = dt_mod %>%
    dplyr::select(R.FileName, R.Condition,R.Replicate) %>%
    distinct()
  
  metadata = data.frame(R.FileName = colnames(mat)) %>%
    left_join(metadata)
  
  
  H2 = heatmaply::heatmaply(
    mat,
    col_side_colors= metadata$R.Condition,
    seriate = "mean", 
    row_dend_left = TRUE,
    plot_method = "plotly"
  )
  
  return(H2)
}
plot_HeatMap_Peptide_Targets_ratios_by_condition<- function(dt, Targets, normalizing_condition){
  
  dt_mod = dt %>%
    filter(EG.ModifiedPeptide %in% Targets$EG.ModifiedPeptide) %>%
    filter(!is.na(FG.MS2Quantity))%>%
    mutate(FG.MS2RawQuantity = as.numeric(FG.MS2RawQuantity))
  
  condition_has_normCondition <- dt_mod %>%
    filter(!is.na(FG.MS2Quantity))%>%
    dplyr::select(R.Condition,PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge ) %>%
    filter(R.Condition == normalizing_condition) %>%
    dplyr::select(-R.Condition) %>%
    mutate(keep= "keep") %>%
    distinct()
  
  
  mean_areas <- dt_mod %>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge, R.Condition) %>%
    summarise(mean = mean(FG.MS2Quantity, na.rm =T)) %>%
    ungroup() %>%
    group_by(PG.ProteinGroups,EG.ModifiedPeptide, FG.Charge) %>%
    mutate(ratio = mean/mean[R.Condition == normalizing_condition]) %>%
    ungroup()%>%
    mutate(ratio = log(ratio,2))# %>%
  # filter(R.Condition != normalizing_condition)
  
  
  
  mat0 <- mean_areas %>%
    dplyr::select(PG.ProteinGroups, EG.ModifiedPeptide, FG.Charge,R.Condition, ratio) %>%
    tidyr::spread(key = R.Condition, value = ratio)
  
  mat <- mat0 %>%
    select(-PG.ProteinGroups, -EG.ModifiedPeptide, -FG.Charge) %>%
    # replace(is.na(.), 0) %>%
    data.matrix()
  
  row_names <- mat0 %>%
    mutate(row_names= paste0(EG.ModifiedPeptide,"_",FG.Charge)) %>%
    select(row_names)
  row.names(mat) = row_names$row_names
  
  
  gradient_col <- ggplot2::scale_fill_gradient2(
    low = "blue", high = "red",mid = "white", 
    midpoint = 0, limits = c(-4, 4)
  )
  
  heatmaply::heatmaply(
    mat,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
      low = "blue", 
      high = "red", 
      midpoint = 0, 
      limits = c(-max(abs(mat)),max(abs(mat)))
      ),
    seriate = "mean"
  )
  
  
  # mat_cor  = cor(mat, use = "pairwise.complete.obs", method = "pearson")
  # mat_cor[is.na(mat_cor)] <- 0
  # 
  # heatmaply::heatmaply_cor(
  #   mat_cor,
  #   xlab = "Condition",
  #   ylab = "Condition"
  # )
  
}
plotly_Peptide_Intensity_accross_runs<- function(dt, Targets, reorder){
  
  dt_mod = dt %>%
    filter(EG.ModifiedPeptide %in% Targets$EG.ModifiedPeptide) %>%
    filter(!is.na(FG.MS2Quantity))%>%
    mutate(FG.MS2RawQuantity = as.numeric(FG.MS2RawQuantity)) %>%
    mutate(ms_run = paste0(R.Condition,"_",R.Replicate))
  
  if(reorder) {
    order_dt <- dt %>%
      mutate(ms_run = paste0(R.Condition,"_",R.Replicate))%>%
      dplyr::select(ms_run, R.Condition,R.Replicate,order) %>%
      distinct() %>%
      arrange(order)
    
    dt_mod$ms_run <- factor(dt_mod$ms_run, levels = unique(order_dt$ms_run))
  }
  
  g = ggplot(dt_mod, aes(x = ms_run, y=  log(FG.MS2Quantity,10),  color= paste0(EG.ModifiedPeptide, FG.Charge),group= paste0(EG.ModifiedPeptide, FG.Charge)))+
    geom_point()+
    geom_line()+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))+
    labs(x="MS run", y= "log10(Quantity)")
  
  ggplotly(g)
  
}

plot_HeatMap_Protein_Targets <- function(dt, Targets, clustering=F){
  
  
  dt_mod = dt %>%
    filter(PG.ProteinGroups %in% Targets$PG.ProteinGroups) %>%
    dplyr::select(R.FileName, R.Condition,R.Replicate, PG.ProteinGroups, PG.ProteinAccessions, PG.ProteinNames,PG.Genes, PG.MS2Quantity ) %>%
    distinct()%>%
    mutate(PG.MS2Quantity = as.numeric(PG.MS2Quantity))%>%
    mutate(PG.MS2Quantity = log(PG.MS2Quantity,10))
  
  mat0 = Longformat_to_matrixformat_ProteinLevel(dt_mod)
  
  mat <- mat0 %>%
    select(!which(names(.) %in% c("PG.ProteinGroups", "PG.ProteinAccessions", "PG.ProteinNames", "PG.Genes" ))) %>%
    # replace(is.na(.), 0) %>%
    data.matrix()
  
  mat0 <- mat0 %>%
    mutate(row_names= PG.ProteinGroups) %>%
    select(row_names) %>%
    group_by(row_names) %>%
    mutate(id = row_number()) %>%
    mutate(row_names  = paste0(row_names,"_g",id))
  row.names(mat) = mat0$row_names
  
  colnames_mat = dt %>%
    dplyr::select(R.FileName, R.Condition, R.Replicate) %>%
    distinct() %>%
    mutate(colnames_mat_values = paste0(R.Condition,"_", R.Replicate))
  
  colnames_mat <- data.frame(R.FileName = colnames(mat)) %>%
    left_join(colnames_mat)
  colnames(mat) <-colnames_mat$colnames_mat_values
  
  
  if(clustering){
    H2 = heatmaply::heatmaply(
      mat,
      col_side_colors= colnames(mat),
      seriate = "mean", 
      row_dend_left = TRUE,
      plot_method = "plotly")
  } else {
    H2 = heatmaply::heatmaply(
      mat,
      col_side_colors= colnames(mat),
      seriate = "mean", 
      row_dend_left = TRUE,
      plot_method = "plotly",
      dendrogram = "none")
  }
  
  
  # ComplexHeatmap::Heatmap(
  #   mat )
  # 
  return(H2)
}
plot_HeatMap_Protein_Targets_ratios_by_condition<- function(dt, Targets, normalizing_condition){
  
  dt_mod = dt %>%
    filter(PG.ProteinGroups %in% Targets$PG.ProteinGroups) %>%
    dplyr::select(R.FileName, R.Condition,R.Replicate, PG.ProteinGroups, PG.ProteinAccessions, PG.ProteinNames,PG.Genes, PG.MS2Quantity ) %>%
    distinct()%>%
    mutate(PG.MS2Quantity = as.numeric(PG.MS2Quantity))
  
  condition_has_normCondition <- dt_mod %>%
    filter(!is.na(PG.MS2Quantity))%>%
    dplyr::select(R.Condition,PG.ProteinGroups ) %>%
    filter(R.Condition == normalizing_condition) %>%
    dplyr::select(-R.Condition) %>%
    mutate(keep= "keep") %>%
    distinct()
  
  
  mean_areas <- dt_mod %>%
    left_join(condition_has_normCondition) %>% 
    filter(keep == "keep") %>%
    group_by(PG.ProteinGroups, R.Condition) %>%
    summarise(mean = mean(PG.MS2Quantity, na.rm =T)) %>%
    ungroup() %>%
    group_by(PG.ProteinGroups) %>%
    mutate(ratio = mean/mean[R.Condition == normalizing_condition]) %>%
    ungroup()%>%
    mutate(ratio = log(ratio,2))
  
  
  
  mat0 <- mean_areas %>%
    dplyr::select(PG.ProteinGroups, R.Condition, ratio) %>%
    tidyr::spread(key = R.Condition, value = ratio)
  
  mat <- mat0 %>%
    select(-PG.ProteinGroups) %>%
    # replace(is.na(.), 0) %>%
    data.matrix()
  
  row_names <- mat0 %>%
    mutate(row_names= PG.ProteinGroups) %>%
    select(row_names)
  row.names(mat) = row_names$row_names
  
  
  # gradient_col <- ggplot2::scale_fill_gradient2(
  #   low = "blue", high = "red",mid = "white", 
  #   midpoint = 0, limits = c(-4, 4)
  # )
  mat<- data.matrix(mat)
  
  heatmaply::heatmaply(
    mat,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
      low = "blue", 
      high = "red", 
      midpoint = 0, 
      limits = c(-max(abs(mat)),max(abs(mat)))
    )
  )

}

plotly_Protein_Intensity_accross_runs<- function(dt, Targets, reorder){
  
  dt_mod = dt %>%
    filter(PG.ProteinGroups %in% Targets$PG.ProteinGroups) %>%
    dplyr::select(R.FileName, R.Condition,R.Replicate, PG.ProteinGroups, PG.ProteinAccessions, PG.ProteinNames,PG.Genes, PG.MS2Quantity ) %>%
    distinct()%>%
    mutate(PG.MS2Quantity = as.numeric(PG.MS2Quantity))%>%
    mutate(ms_run = paste0(R.Condition,"_",R.Replicate))
  
  if(reorder) {
    order_dt <- dt %>%
      mutate(ms_run = paste0(R.Condition,"_",R.Replicate))%>%
      dplyr::select(ms_run, R.Condition,R.Replicate,order) %>%
      distinct() %>%
      arrange(order)
    
    dt_mod$ms_run <- factor(dt_mod$ms_run, levels = unique(order_dt$ms_run))
  }
  
  g = ggplot(dt_mod, aes(x = ms_run, y=  log(PG.MS2Quantity,10),  color= PG.ProteinGroups,group= PG.ProteinGroups))+
    geom_point()+
    geom_line()+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))+
    labs(x="MS run", y= "log10(Quantity)")
  
  ggplotly(g)
  
}


#### PRM_Tools
format_diaData_to_PRMmethod <- function(dt,
                                        Targets,
                                        Isolation_Width = 3,
                                        RT_Range_seconds=300,
                                        IM_range = 0.07,
                                        experiment_name = "Targeted_peptides"){
  
  
  dt2 <- dt %>%
    filter(EG.ModifiedPeptide %in% Targets$EG.ModifiedPeptide) %>%
    ungroup() %>%
    group_by(EG.ModifiedPeptide,FG.Charge ) %>%
    select(EG.ModifiedPeptide,FG.Charge,R.Condition, R.Replicate, EG.RTPredicted, EG.IonMobility, FG.PrecMz)%>%
    distinct() %>%
    summarise(EG.RTPredicted = mean(EG.RTPredicted, na.rm =T)*60,
              EG.IonMobility= mean(EG.IonMobility, na.rm =T),
              FG.PrecMz= unique(FG.PrecMz, na.rm =T))
  
  
  dt_formatted <- dt2 %>%
    ungroup %>%
    mutate(Isolation_Width_col = Isolation_Width,
           RT_Range_seconds_col = RT_Range_seconds,
           IM_range_col = IM_range,
           "Start IM [1/K0]" = EG.IonMobility - IM_range_col/2,
           "End IM [1/K0]" = EG.IonMobility + IM_range_col/2,
           "CE [eV]" = NA,
           "External ID"= paste0(EG.ModifiedPeptide,"_",FG.Charge),
           "Description" = experiment_name) %>%
    select("Mass [m/z]"  = FG.PrecMz,
           "Charge" = FG.Charge,
           "Isolation Width [m/z]" = Isolation_Width_col,
           "RT [s]" = EG.RTPredicted,
           "RT Range [s]" = RT_Range_seconds_col,
           "Start IM [1/K0]" ,
           "End IM [1/K0]" ,
           "CE [eV]",
           "External ID",
           "Description")%>%
    mutate(`RT [s]` = ifelse(`RT [s]` - `RT Range [s]`/2 < 0,  RT_Range_seconds  ,`RT [s]`))#%>%
  # mutate(`RT [s]` = ifelse(`RT [s]` - `RT Range [s]`/2 < 0, `RT [s]` = RT_Range_seconds  ,`RT [s]`))%>%
  
  return(dt_formatted)
}


