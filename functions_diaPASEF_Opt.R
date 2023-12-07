library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(UpSetR)
library(plotly)
library(digest)
library(tidygenomics)

redo_DIA_cycles <-function(dt){
  
  
  dt_OG = dt %>%
    ungroup() %>%
    arrange(desc(end_im))%>%
    mutate(is_recycled = 0) %>%
    mutate(id = row_number()) %>%
    mutate(new_cycle = 0)
  
  dt_list = list()
  
  cycle_n = 1
  dt_intern = dt_OG
  ids_used = NULL
  
  while (dim(dt_intern)[1]> 0) {
    dt_working = dt_intern
    
    while(dim(dt_working)[1]>0) {
      x = dt_working[1,]
      x$new_cycle = cycle_n
      x$is_recycled = 1
      dt_list[[x$id]] = x
      
      # dt_working <- dt_working %>% filter(`End IM [1/K0]` < x$`Start IM [1/K0]`-0.006)
      dt_working <- dt_working %>% filter(end_im < x$start_im )
      
      ids_used = append(ids_used,values = x$id)
    }
    
    cycle_n = cycle_n+1
    dt_intern = dt_OG %>%
      filter(!id %in% ids_used)
    
  }
  
  dt2 = bind_rows(dt_list) %>%
    arrange(new_cycle, start_im ,end_im  ) %>%
    mutate(cycle = new_cycle) %>%
    select(-is_recycled, - id, - new_cycle)
  
  
  return(dt2)
}
redo_DIA_cycles_distribute_windows_to_have_same_number <-function(dt){
  
  max_num_cycles = length(unique(dt$cycle))
  number_of_windows = dim(dt)[1]
  dt2 = dt %>%
    arrange(desc(end_im))
  
  dt2$cycle = rep_len(1:max_num_cycles, number_of_windows)
  
  
  return(dt2)
}
redo_DIA_cycles_distribute_windows_to_have_same_number_plusONEcycle <-function(dt){
  
  max_num_cycles = length(unique(dt$cycle))+1
  number_of_windows = dim(dt)[1]
  dt2 = dt %>%
    arrange(desc(end_im))
  
  dt2$cycle = rep_len(1:max_num_cycles, number_of_windows)
  
  
  return(dt2)
}

extend_windows_to_upper_lower_limit <- function(dt, im_upper_limit, im_lower_limit){
  dt2 <- dt %>%
    ungroup() %>%
    group_by(cycle) %>%
    mutate(end_im = ifelse(end_im == max(end_im), im_upper_limit, end_im))%>%
    mutate(start_im = ifelse(start_im == min(start_im), im_lower_limit, start_im))
  return(dt2)
    
}

extend_windows_to_upper_lower_limit_2 <- function(dt, im_upper_limit, im_lower_limit, im_value_to_add =0.05){
  dt2 <- dt %>%
    ungroup() %>%
    group_by(cycle) %>%
    mutate(end_im = ifelse(end_im == max(end_im), end_im + im_value_to_add, end_im))%>%
    mutate(start_im = ifelse(start_im == min(start_im), start_im - im_value_to_add, start_im)) %>%
    mutate(end_im = ifelse(end_im > im_upper_limit, im_upper_limit,end_im ),
           start_im = ifelse(start_im < im_lower_limit, im_lower_limit, start_im))
  return(dt2)
  
}

extend_windows_to_meet_halfway <- function(dt){
  dt2 <- dt %>%
    ungroup() %>%
    group_by(cycle)%>%
    arrange(desc(end_im)) %>%
    arrange(cycle) %>%
    mutate(end_im_next = lead(end_im),
           start_in_previous = lag(start_im)) %>%
    mutate(star_im_2 = ifelse(!is.na(end_im_next), (start_im + end_im_next)/2, start_im),
           end_im_2 = ifelse(!is.na(start_in_previous), (end_im + start_in_previous)/2, end_im)) %>%
    mutate(start_im = star_im_2,
           end_im = end_im_2)
  
  return(dt2)
  
}

cut_IM_windows_in_two = function(dt){
  x <- dt$center
  y <- dt$center_im 
  
  mod <- lm(y ~ x)
  cf <- coef(mod)
  
  Intercept <- cf[1]
  Slope <- cf[2]
  
  dt2 <- dt %>%
    mutate(new_lim = Slope*center + Intercept)
  dt2_1 <- dt2 %>%
    mutate(start_im = ifelse(start_im < new_lim,new_lim,start_im)) %>%
    mutate(im_window = "up")
  dt2_2 <- dt2 %>%
    mutate(end_im = ifelse(end_im > new_lim, new_lim,end_im))%>%
    mutate(im_window = "down")
  
  dt3= rbind(dt2_1,dt2_2) %>%
    distinct()%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) 
  
  # median_widths = dt3 %>%
  #   group_by(im_window) %>%
  #   summarise(median_im_width = median(width_im,na.rm = T))
  # 
  # dt4 = dt3 %>%
  #   left_join(median_widths) %>%
  #   mutate(end_im = ifelse(im_window == "up", start_im + median_im_width,end_im))%>%
  #   mutate(start_im = ifelse(im_window == "down",end_im - median_im_width,start_im ))
  # 
  
  
}

cut_IM_windows_in_two_with_slope_lower_bound = function(dt){
  dt_filt = dt%>% filter(end_mz<800)
  
  x <- dt_filt$center
  y <- dt_filt$start_im  
  
  mod <- lm(y ~ x)
  cf <- coef(mod)
  
  Intercept <- cf[1]
  Slope <- cf[2]
  
  dt2 <- dt %>%
    mutate(new_lim = Slope*center + Intercept)
  dt2_1 <- dt2 %>%
    mutate(start_im = ifelse(start_im < new_lim,new_lim,start_im)) %>%
    mutate(im_window = "up")
  dt2_2 <- dt2 %>%
    mutate(end_im = ifelse(end_im > new_lim, new_lim,end_im))%>%
    mutate(im_window = "down")
  
  dt3= rbind(dt2_1,dt2_2) %>%
    distinct()%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) 
  
  # median_widths = dt3 %>%
  #   group_by(im_window) %>%
  #   summarise(median_im_width = median(width_im,na.rm = T))
  # 
  # dt4 = dt3 %>%
  #   left_join(median_widths) %>%
  #   mutate(end_im = ifelse(im_window == "up", start_im + median_im_width,end_im))%>%
  #   mutate(start_im = ifelse(im_window == "down",end_im - median_im_width,start_im ))
  # 
  
  
}


cut_IM_windows_in_two_above_mzvalue = function(dt){
  x <- dt$center
  y <- dt$center_im 
  
  mod <- lm(y ~ x)
  cf <- coef(mod)
  
  Intercept <- cf[1]
  Slope <- cf[2]
  
  dt2 <- dt %>%
    mutate(new_lim = Slope*center + Intercept)
  dt2_1 <- dt2 %>%
    mutate(start_im = ifelse(start_im < new_lim,new_lim,start_im)) %>%
    mutate(im_window = "up")
  dt2_2 <- dt2 %>%
    mutate(end_im = ifelse(end_im > new_lim, new_lim,end_im))%>%
    mutate(im_window = "down")
  
  dt3= rbind(dt2_1,dt2_2) %>%
    distinct()%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) 
  
  dt_finale = rbind(dt2 %>%
    filter(end_mz<800)%>%
      mutate(center = (end_mz+start_mz)/2,
             center_im = (end_im+start_im)/2,
             width = end_mz - start_mz,
             width_im = end_im - start_im)%>%
      mutate(im_window = "normal"),
    dt3 %>%
      filter(end_mz>800))
  
  # median_widths = dt3 %>%
  #   group_by(im_window) %>%
  #   summarise(median_im_width = median(width_im,na.rm = T))
  # 
  # dt4 = dt3 %>%
  #   left_join(median_widths) %>%
  #   mutate(end_im = ifelse(im_window == "up", start_im + median_im_width,end_im))%>%
  #   mutate(start_im = ifelse(im_window == "down",end_im - median_im_width,start_im ))
  # 
  
  
}

cut_IM_windows_in_three = function(dt,Slope =0.0006, Intercept1=0.60,Intercept2=0.5){
  # x <- dt$center
  # y <- dt$center_im 
  # 
  # mod <- lm(y ~ x)
  # cf <- coef(mod)
  # 
  # Intercept <- cf[1]
  # Slope <- cf[2]
  
  p = ggplot() + 
    scale_x_continuous(name="mz") +
    scale_y_continuous(name="im",n.breaks = 20) +
    # geom_point(data = HLA_data, aes(x= mz, y=im, alpha=0.2))+
    ggpointdensity::geom_pointdensity(data = HLA_data, aes(x= mz, y=im))+
    viridis::scale_color_viridis()
  p +
    coord_cartesian(xlim = c(0,1200), ylim = c(0,1.3))+
    geom_abline(slope = 0.0006, intercept = 0.60)+
    geom_abline(slope = 0.0006, intercept = 0.45)
  
  
  dt2 <- dt %>%
    mutate(new_lim_up = Slope*center + Intercept1,
           new_lim_down = Slope*center + Intercept2)
  dt2_up <- dt2 %>%
    mutate(im_window = "down") %>%
    mutate(start_im = start_im,
           end_im = new_lim_down)
  dt2_middle <- dt2 %>%
    mutate(im_window = "middle") %>%
    mutate(start_im = new_lim_down,
           end_im = ifelse(end_im< new_lim_up,end_im, new_lim_up))
  dt2_down <- dt2 %>%
    mutate(im_window = "up") %>%
    mutate(start_im = new_lim_up,
           end_im = end_im)
  
  
  dt3= rbind(dt2_up,dt2_middle,dt2_down) %>%
    filter(!end_im<= start_im) %>%
    distinct()%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im)
  
  # median_widths = dt3 %>%
  #   group_by(im_window) %>%
  #   summarise(median_im_width = median(width_im,na.rm = T))
  # 
  # dt4 = dt3 %>%
  #   left_join(median_widths) %>%
  #   mutate(end_im = ifelse(im_window == "up", start_im + median_im_width,end_im))%>%
  #   mutate(start_im = ifelse(im_window == "down",end_im - median_im_width,start_im ))
  
  
  
}

count_datapoints_per_window = function(data, dia_method){
  
  data2 = data %>%
    ungroup() %>%
    group_by(bin_mz) %>%
    filter(precursor_mz < max(dia_method$end_mz),
           precursor_mz > min(dia_method$start_mz)) %>%
    ungroup() %>%
    select(im, bin_mz )
  
  
  list_dia_method = list()
  
  
  
  for(i in 1:length(unique(intersect(dia_method$bin_mz, data2$bin_mz)))){
    
    bin = unique(data2$bin_mz)[i]
    
    dia_method2 = dia_method %>%
      filter(bin_mz == bin)
    
    list_dia_method_intermediate = list()
    
    for(j in 1: dim(dia_method2)[1]){
      
      dia_method3 = dia_method2[j,]
      
      n_precursor = data2 %>%
        filter(bin_mz == bin) %>%
        filter(im> dia_method3$start_im,
               im< dia_method3$end_im) %>%
        ungroup %>% 
        tally() %>%
        data.frame()
      dia_method3$n_precursors = n_precursor$n
      
      list_dia_method_intermediate[[j]] = dia_method3
      
    }
    
    list_dia_method[[i]] = bind_rows(list_dia_method_intermediate)
  }
  
  dia_method_2 = bind_rows(list_dia_method)
  
}
count_datapoints_per_MZ_and_RT_window = function(data, dia_method){
  
  data2 = data %>%
    ungroup() %>%
    group_by(bin_mz, bin_RT) %>%
    filter(precursor_mz < max(dia_method$end_mz),
           precursor_mz > min(dia_method$start_mz)) %>%
    ungroup() %>%
    select(im, bin_mz, bin_RT)
  
  list_dia_method_RT = list()
  
  for(r in length(unique(intersect(dia_method$bin_RT, data2$bin_RT)))){
    list_dia_method = list()
    
    for(i in 1:length(unique(intersect(dia_method$bin_mz, data2$bin_mz)))){
      
      bin = unique(data2$bin_mz)[i]
      
      dia_method2 = dia_method %>%
        filter(bin_mz == bin)
      
      list_dia_method_intermediate = list()
      
      for(j in 1: dim(dia_method2)[1]){
        
        dia_method3 = dia_method2[j,]
        
        n_precursor = data2 %>%
          filter(bin_mz == bin) %>%
          filter(im> dia_method3$start_im,
                 im< dia_method3$end_im) %>%
          ungroup %>% 
          tally() %>%
          data.frame()
        dia_method3$n_precursors = n_precursor$n
        
        list_dia_method_intermediate[[j]] = dia_method3
        
      }
      
      list_dia_method[[i]] = bind_rows(list_dia_method_intermediate)
    }
    
    list_dia_method_RT[[r]] = bind_rows(list_dia_method)
  }
  
  dia_method_RT = bind_rows(list_dia_method_RT)
  
}

redo_DIA_cycles_withCounts <-function(dt){
  
  
  dt_OG = dt %>%
    ungroup() %>%
    arrange(desc(n_precursors))%>%
    mutate(is_recycled = 0) %>%
    mutate(id = row_number()) %>%
    mutate(new_cycle = 0)
  
  dt_list = list()
  
  cycle_n = 1
  dt_intern = dt_OG
  ids_used = NULL
  
  while (dim(dt_intern)[1]> 0) {
    dt_working = dt_intern
    
    while(dim(dt_working)[1]>0) {
      x = dt_working[1,]
      x$new_cycle = cycle_n
      x$is_recycled = 1
      dt_list[[x$id]] = x
      
      # dt_working <- dt_working %>% filter(end_im < x$start_im-0.006)
      dt_working_1 <- dt_working %>%
        filter(end_im < x$start_im)
      dt_working_2 <- dt_working %>%
        filter(start_im> x$end_im)
      
      dt_working = rbind(dt_working_1, dt_working_2) %>%
        arrange(desc(n_precursors))
      
      
      ids_used = append(ids_used,values = x$id)
    }
    
    cycle_n = cycle_n+1
    dt_intern = dt_OG %>%
      filter(!id %in% ids_used)
    
  }
  
  dt2 = bind_rows(dt_list) %>%
    arrange(new_cycle, start_im,end_im  ) %>%
    mutate(cycle = new_cycle) %>%
    select(-is_recycled, - id, - new_cycle)
  
  
  return(dt2)
}
redo_DIA_cycles_withCounts_per_RT_window <-function(dt){
  dia_method_cycles_Redone = list()
  
  for (r in 1:length(unique(dt$bin_RT))) {
    print(r)
    selected_bin_RT = unique(dt$bin_RT) [r]
    
    dt_OG = dt %>%
      filter(bin_RT == selected_bin_RT) %>%
      ungroup() %>%
      arrange(desc(n_precursors))%>%
      mutate(is_recycled = 0) %>%
      mutate(id = row_number()) %>%
      mutate(new_cycle = 0)
    
    dt_list = list()
    
    cycle_n = 1
    dt_intern = dt_OG
    ids_used = NULL
    
    
    
    while (dim(dt_intern)[1]> 0) {
      dt_working = dt_intern
      
      while(dim(dt_working)[1]>0) {
        x = dt_working[1,]
        x$new_cycle = cycle_n
        x$is_recycled = 1
        dt_list[[x$id]] = x
        
        # dt_working <- dt_working %>% filter(end_im < x$start_im-0.006)
        dt_working_1 <- dt_working %>%
          filter(end_im < x$start_im)
        dt_working_2 <- dt_working %>%
          filter(start_im> x$end_im)
        
        dt_working = rbind(dt_working_1, dt_working_2) %>%
          arrange(desc(n_precursors))
        
        
        ids_used = append(ids_used,values = x$id)
      }
      
      cycle_n = cycle_n+1
      dt_intern = dt_OG %>%
        filter(!id %in% ids_used)
      
    }
    
    dia_method_cycles_Redone[[r]] = bind_rows(dt_list) %>%
      arrange(new_cycle, start_im,end_im  ) %>%
      mutate(cycle = new_cycle) %>%
      select(-is_recycled, - id, - new_cycle)
  }
  
  dia_method_cycles_Redone = bind_rows(dia_method_cycles_Redone)
  
  return(dia_method_cycles_Redone)
}

Add_mz_overlaps <- function (dt, overlap_mz){
  dt2 <- dt %>%
    mutate(end_mz = end_mz + overlap_mz,
           start_mz = start_mz - overlap_mz)
  return(dt2)
}

cut_SlicePASEF_mz_windows_in_two = function(dt){
  
  dt2_1 <- dt %>%
    mutate(end_mz = center_mz )
  dt2_2 <- dt %>%
    mutate(start_mz = center_mz )
  
  dt3= rbind(dt2_1,dt2_2) %>%
    distinct()%>%
    mutate(center_mz = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) 
  
  
  
  
}

calc_number_of_windows <- function(data_mz_im, mz_width = 25){
  
  bin_mz_borders = seq(mz_range_lower, mz_range_upper, mz_width)
  bin_mz_borders_2 = c(bin_mz_borders[-length(bin_mz_borders)],mz_range_upper)
  
  data_mz_im <- data_mz_im %>%
    filter(precursor_mz > mz_range_lower, precursor_mz<mz_range_upper) %>%
    mutate(bin_mz = cut(precursor_mz, breaks = bin_mz_borders_2)) 
  
  dia_method = data_mz_im%>%
    group_by(bin_mz)%>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    summarise(im_limit_low = ifelse(max<= 850, quantile(im, mz_window_lower_percentile_limit), quantile(im, 0.15)),
              im_limit_high = ifelse(max<= 850, quantile(im, mz_window_upper_percentile_limit), quantile(im, 0.95)),
              num_precursors = n()) %>%
    distinct() %>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    rename(start_im = im_limit_low,
           end_im = im_limit_high,
           start_mz = min,
           end_mz = max) %>%
    mutate(`Cycle Id` = 1)%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) %>%
    mutate(end_im = center_im+width_im/2+im_offset_lower,
           start_im = center_im - width_im/2-im_offset_upper) %>%
    count_datapoints_per_window(data = data_mz_im, dia_method = .) %>%
    filter(n_precursors>20) %>%
    Add_mz_overlaps(overlap_mz = 0.5) %>%
    redo_DIA_cycles %>%
    redo_DIA_cycles_distribute_windows_to_have_same_number %>%
    extend_windows_to_upper_lower_limit_2(im_upper_limit = 1.5,
                                          im_lower_limit = 0.6,
                                          im_value_to_add = 0.05) %>%
    extend_windows_to_meet_halfway
  
  number_of_cycles = length(unique(dia_method$cycle))
  return(number_of_cycles)
}

format_data_input <- function(dt, im_column_name,
                              precursor_mz_column_name,
                              retention_time_column_name = NULL,
                              sequence_column_name = NULL,
                              charge_column_name = NULL,
                              score_wt_column_name = NULL){
  
  names(dt)[names(dt) == im_column_name] <- 'im'
  names(dt)[names(dt) == precursor_mz_column_name] <- 'precursor_mz'
  if(!is.null(retention_time_column_name)) {names(dt)[names(dt) == retention_time_column_name] <- 'retention_time'}
  if(!is.null(sequence_column_name)) {names(dt)[names(dt) == sequence_column_name] <- 'sequence'}
  if(!is.null(charge_column_name)) {names(dt)[names(dt) == charge_column_name] <- 'charge'}
  if(!is.null(score_wt_column_name)) {names(dt)[names(dt) == score_wt_column_name] <- 'score_wt'}
  
  dt2 <- dt %>%
    dplyr::select(which(names(.) %in% c("im","precursor_mz", "retention_time", "sequence","charge","score_wt"))) %>%
    distinct()
    
  return(dt2)
}

format_dia_method_input <- function(dia_method){
  dia_method2 <- dia_method %>%
    filter(`#MS Type` != "MS1") %>%
    rename(start_mz = `Start Mass [m/z]`,
           end_mz = `End Mass [m/z]`,
           start_im = `Start IM [1/K0]`,
           end_im = `End IM [1/K0]`) %>%
    mutate(start_mz = as.numeric(start_mz),
           end_mz = as.numeric(end_mz),
           start_im = as.numeric(start_im),
           end_im = as.numeric(end_im),
           `Cycle Id` = as.numeric(`Cycle Id`)) %>%
    mutate(width = end_mz - start_mz) %>%
    mutate(width_im = end_im - start_im)
    
  return(dia_method2)
  }


format_DIA_PASEF_method <- function(dt){
  dt2 <- dt %>%
    select(cycle,
           start_im,
           end_im,
           start_mz,
           end_mz) %>%
    distinct() %>%
    rename("Cycle Id" = cycle,
           "Start IM [1/K0]" = start_im,
           "End IM [1/K0]" = end_im,
           "Start Mass [m/z]" = start_mz,
           "End Mass [m/z]" = end_mz) %>%
    mutate(`#MS Type` = "PASEF",
           `CE [eV]` = "-",
           `Start IM [1/K0]` = round(`Start IM [1/K0]`,3),
           `End IM [1/K0]` = round(`End IM [1/K0]`,3)) %>%
    select("#MS Type","Cycle Id","Start IM [1/K0]","End IM [1/K0]","Start Mass [m/z]", "End Mass [m/z]", "CE [eV]") 
  
  dt_MS1 = tibble(`#MS Type` = "MS1",`Cycle Id`=0,`Start IM [1/K0]` = "-",
                      `End IM [1/K0]`= "-",`Start Mass [m/z]`= "-", `End Mass [m/z]`= "-", `CE [eV]`= "-")
  
  dia_method = rbind(dt_MS1, dt2)
  
  return(dia_method)
  }

format_SlicePASEF_method <- function(dt, n_repetitions){
  dt2 <- dt %>%
    select(
           start_im,
           end_im,
           start_mz,
           end_mz) %>%
    distinct() %>%
    rename("Start IM [1/K0]" = start_im,
           "End IM [1/K0]" = end_im,
           "Start Mass [m/z]" = start_mz,
           "End Mass [m/z]" = end_mz) %>%
    mutate(`#MS Type` = "PASEF",
           `CE [eV]` = "-",
           `Start IM [1/K0]` = round(`Start IM [1/K0]`,3),
           `End IM [1/K0]` = round(`End IM [1/K0]`,3),
           `Start Mass [m/z]` = round(`Start Mass [m/z]`,3),
           `End Mass [m/z]` = round(`End Mass [m/z]`,3)) 
  
  dt_list = list()
 for(i in 1:n_repetitions){
   dt2$`Cycle Id` = i
   dt_list[[i]] = dt2
 } 
  
  dt2 = bind_rows(dt_list) %>%
    select("#MS Type","Cycle Id","Start IM [1/K0]","End IM [1/K0]","Start Mass [m/z]", "End Mass [m/z]", "CE [eV]") 
  
  dt_MS1 = tibble(`#MS Type` = "MS1",`Cycle Id`=0,`Start IM [1/K0]` = "-",
                  `End IM [1/K0]`= "-",`Start Mass [m/z]`= "-", `End Mass [m/z]`= "-", `CE [eV]`= "-")
  
  dia_method = rbind(dt_MS1, dt2)
  
  return(dia_method)
}



do_optimization_graph_data <- function(optimization_params, optimization_values_to_try){
  optimization_params
  number_of_cycles_against_DaWindow = list()
  for(i in optimization_values_to_try){
    
    optimization_params$mz_window_width = i
    dt_opti = dt <- do.call(create_DIA_PASEF_method, optimization_params)

    
    num_cycles = length(unique(dt_opti$cycle))
    number_of_cycles_against_DaWindow[[i]] = data.frame(mz_window = i, num_cycles)
    print(i)
  }
  
  number_of_cycles_against_DaWindow_results = bind_rows(number_of_cycles_against_DaWindow)
  
  # ggplot(number_of_cycles_against_DaWindow_results, aes(x= mz_window, y= num_cycles, group =1, label = num_cycles))+
  #   geom_line()+
  #   geom_point()+
  #   geom_text(nudge_x = 2,nudge_y = 1)+
  #   theme_bw()+
  #   scale_x_continuous(breaks = number_of_cycles_against_DaWindow_results$mz_window)+
  #   scale_y_continuous(n.breaks = 10)
  # 
  # estimated_cycle_times = number_of_cycles_against_DaWindow_results %>%
  #   mutate(Acc_50 = num_cycles*(50+5.9),
  #          Acc_75 = num_cycles*(75+5.9),
  #          Acc_100 = num_cycles*(100+5.9),
  #          Acc_125 = num_cycles*(125+5.9),
  #          Acc_166 = num_cycles*(166+5.9)) %>%
  #   gather(key = "Acc_time", value = "cycle_time_ms", which(str_detect(names(.),pattern = "Acc_"))) %>%
  #   mutate(Acc_time = as.numeric(gsub(Acc_time, pattern = "Acc_",replacement = "")))
  # 
  # ggplot(estimated_cycle_times, aes(x= mz_window, y= cycle_time_ms, group = Acc_time, color = as.character(Acc_time)))+
  #   geom_line()+
  #   geom_point()+
  #   theme_bw()+
  #   scale_x_continuous(breaks = estimated_cycle_times$mz_window)+
  #   scale_y_continuous(n.breaks = 10)
  
  return(number_of_cycles_against_DaWindow_results)
  
}
do_optimization_graph_data_VW <- function(optimization_params, optimization_values_to_try){
  optimization_params
  number_of_cycles_against_DaWindow = list()
  for(i in optimization_values_to_try){
    
    optimization_params$number_of_MZcuts = i
    dt_opti = dt <- do.call(create_DIA_PASEF_method_with_VW, optimization_params)
    
    
    num_cycles = length(unique(dt_opti$cycle))
    number_of_cycles_against_DaWindow[[i]] = data.frame(mz_window = i, num_cycles)
    print(i)
  }
  
  number_of_cycles_against_DaWindow_results = bind_rows(number_of_cycles_against_DaWindow)
  
  # ggplot(number_of_cycles_against_DaWindow_results, aes(x= mz_window, y= num_cycles, group =1, label = num_cycles))+
  #   geom_line()+
  #   geom_point()+
  #   geom_text(nudge_x = 2,nudge_y = 1)+
  #   theme_bw()+
  #   scale_x_continuous(breaks = number_of_cycles_against_DaWindow_results$mz_window)+
  #   scale_y_continuous(n.breaks = 10)
  # 
  # estimated_cycle_times = number_of_cycles_against_DaWindow_results %>%
  #   mutate(Acc_50 = num_cycles*(50+5.9),
  #          Acc_75 = num_cycles*(75+5.9),
  #          Acc_100 = num_cycles*(100+5.9),
  #          Acc_125 = num_cycles*(125+5.9),
  #          Acc_166 = num_cycles*(166+5.9)) %>%
  #   gather(key = "Acc_time", value = "cycle_time_ms", which(str_detect(names(.),pattern = "Acc_"))) %>%
  #   mutate(Acc_time = as.numeric(gsub(Acc_time, pattern = "Acc_",replacement = "")))
  # 
  # ggplot(estimated_cycle_times, aes(x= mz_window, y= cycle_time_ms, group = Acc_time, color = as.character(Acc_time)))+
  #   geom_line()+
  #   geom_point()+
  #   theme_bw()+
  #   scale_x_continuous(breaks = estimated_cycle_times$mz_window)+
  #   scale_y_continuous(n.breaks = 10)
  
  return(number_of_cycles_against_DaWindow_results)
  
}

plot_optimization_results<-function(number_of_cycles_against_DaWindow_results, ylim = 3000){
  a = ggplot(number_of_cycles_against_DaWindow_results, aes(x= mz_window, y= num_cycles, group =1, label = num_cycles))+
    geom_line()+
    geom_point()+
    geom_text(nudge_x = 2,nudge_y = 1)+
    theme_bw()+
    scale_x_continuous(breaks = number_of_cycles_against_DaWindow_results$mz_window)+
    scale_y_continuous(n.breaks = 10)

  estimated_cycle_times = number_of_cycles_against_DaWindow_results %>%
    mutate(Acc_30 = num_cycles*(30+5.9),
           Acc_50 = num_cycles*(50+5.9),
           Acc_75 = num_cycles*(75+5.9),
           Acc_100 = num_cycles*(100+5.9),
           Acc_125 = num_cycles*(125+5.9),
           Acc_166 = num_cycles*(166+5.9)) %>%
    gather(key = "Acc_time", value = "cycle_time_ms", which(str_detect(names(.),pattern = "Acc_"))) %>%
    mutate(Acc_time = as.numeric(gsub(Acc_time, pattern = "Acc_",replacement = "")))

  estimated_cycle_times$Acc_time <- factor(estimated_cycle_times$Acc_time, levels =  c(30,50,75,100,125,166))
  
  b= ggplot(estimated_cycle_times, aes(x= mz_window, y= cycle_time_ms, group = Acc_time, color = Acc_time))+
    geom_line()+
    geom_point()+
    theme_bw()+   
    theme(legend.position = "top",
          panel.grid.minor = element_blank())+
    labs(color='Acc. time (ms)') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))+
    scale_x_continuous(breaks = estimated_cycle_times$mz_window)+
    scale_y_continuous(n.breaks = 10)+
    coord_cartesian(ylim = c(0,ylim))
  
  plot_grid(a,b, nrow = 1)
  
}

## DIA_PASEF

create_DIA_PASEF_method <- function(data_mz_im,
                                    im_range_lower = 0.7,
                                    im_range_upper = 1.3,
                                    mz_range_lower = 300,
                                    mz_range_upper = 1200,
                                    mz_window_width = 15,
                                    mz_window_lower_percentile_limit = 0.1,
                                    mz_window_upper_percentile_limit = 0.95,
                                    im_offset_lower = 0.01,
                                    im_offset_upper = 0.01,
                                    num_precursors = 10,
                                    Add_mz_overlaps_width = 1){
  bin_mz_borders = seq(mz_range_lower, mz_range_upper, mz_window_width)
  bin_mz_borders_2 = c(bin_mz_borders[-length(bin_mz_borders)],mz_range_upper)
  
  data_mz_im <- data_mz_im %>%
    filter(precursor_mz > mz_range_lower, precursor_mz<mz_range_upper) %>%
    # mutate(bin_mz = cut(precursor_mz, breaks = seq(200,1600, mz_window_width))) %>%
    mutate(bin_mz = cut(precursor_mz, breaks = bin_mz_borders_2, dig.lab = 5)) #%>%
  # mutate(bin_mz = cut_number(precursor_mz, n = 22)) 
  
  
  dt = data_mz_im%>%
    group_by(bin_mz)%>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    summarise(im_limit_low = ifelse(max<= 850, quantile(im, mz_window_lower_percentile_limit), quantile(im, 0.15)),
              im_limit_high = ifelse(max<= 850, quantile(im, mz_window_upper_percentile_limit), quantile(im, 0.95)),
              num_precursors = n()) %>%
    distinct() %>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    rename(start_im = im_limit_low,
           end_im = im_limit_high,
           start_mz = min,
           end_mz = max) %>%
    mutate(`Cycle Id` = 1)%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) %>%
    mutate(end_im = center_im+width_im/2+im_offset_lower,
           start_im = center_im - width_im/2-im_offset_upper) %>%
    # filter(num_precursors>20) %>%
    # cut_IM_windows_in_two %>%
    # cut_IM_windows_in_three %>%
    # cut_IM_windows_in_two_with_slope_lower_bound %>%
    count_datapoints_per_window(data = data_mz_im, dia_method = .) %>%
    # filter(n_precursors>=num_precursors) %>%
    # redo_DIA_cycles_withCounts %>%
    Add_mz_overlaps(overlap_mz = 0.5) %>%
    redo_DIA_cycles %>%
    redo_DIA_cycles_distribute_windows_to_have_same_number %>%
    # redo_DIA_cycles_distribute_windows_to_have_same_number_plusONEcycle %>%
    # extend_windows_to_upper_lower_limit(im_upper_limit = 1.5,
    #                                     im_lower_limit = 0.6) %>%
    extend_windows_to_upper_lower_limit_2(im_upper_limit = im_range_upper,
                                          im_lower_limit = im_range_lower,
                                          im_value_to_add = 1) %>%
    extend_windows_to_meet_halfway
}

create_DIA_PASEF_method_with_VW <- function(data_mz_im,
                                            im_range_lower = 0.75,
                                            im_range_upper = 1.3,
                                            mz_range_lower = 380,
                                            mz_range_upper = 980,
                                            # mz_window_width = 15,
                                            number_of_MZcuts = 20,
                                            mz_window_lower_percentile_limit = 0.1,
                                            mz_window_upper_percentile_limit = 0.95,
                                            im_offset_lower = 0.01,
                                            im_offset_upper = 0.01,
                                            num_precursors = 0,
                                            Add_mz_overlaps_width = 1){
  
  
  data_mz_im <- data_mz_im %>%
    filter(precursor_mz > mz_range_lower, precursor_mz<mz_range_upper) %>%
    # mutate(bin_mz = cut(precursor_mz, breaks = seq(200,1600, mz_window_width))) %>%
    # mutate(bin_mz = cut(precursor_mz, breaks = bin_mz_borders_2)) #%>%
    mutate(bin_mz = cut_number(precursor_mz, n = number_of_MZcuts))
  
  
  dt = data_mz_im%>%
    group_by(bin_mz)%>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    summarise(im_limit_low = ifelse(max<= 850, quantile(im, mz_window_lower_percentile_limit), quantile(im, 0.15)),
              im_limit_high = ifelse(max<= 850, quantile(im, mz_window_upper_percentile_limit), quantile(im, 0.95)),
              num_precursors = n()) %>%
    distinct() %>%
    mutate(x_tmp = str_sub(bin_mz, 2, -2)) %>% 
    separate(x_tmp, c("min", "max"), sep = ",") %>% 
    mutate_at(c("min", "max"), as.double) %>%
    rename(start_im = im_limit_low,
           end_im = im_limit_high,
           start_mz = min,
           end_mz = max) %>%
    mutate(`Cycle Id` = 1)%>%
    mutate(center = (end_mz+start_mz)/2,
           center_im = (end_im+start_im)/2,
           width = end_mz - start_mz,
           width_im = end_im - start_im) %>%
    mutate(end_im = center_im+width_im/2+im_offset_lower,
           start_im = center_im - width_im/2-im_offset_upper) %>%
    # filter(num_precursors>20) %>%
    # cut_IM_windows_in_two %>%
    # cut_IM_windows_in_three %>%
    # cut_IM_windows_in_two_with_slope_lower_bound %>%
    count_datapoints_per_window(data = data_mz_im, dia_method = .) %>%
    filter(n_precursors>20) %>%
    # redo_DIA_cycles_withCounts %>%
    Add_mz_overlaps(overlap_mz = 0.5) %>%
    redo_DIA_cycles %>%
    redo_DIA_cycles_distribute_windows_to_have_same_number %>%
    # redo_DIA_cycles_distribute_windows_to_have_same_number_plusONEcycle %>%
    # extend_windows_to_upper_lower_limit(im_upper_limit = 1.5,
    #                                     im_lower_limit = 0.6) %>%
    extend_windows_to_upper_lower_limit_2(im_upper_limit = im_range_upper,
                                          im_lower_limit = im_range_lower,
                                          im_value_to_add = 0.15) %>%
    extend_windows_to_meet_halfway
}


#### read spectral library

read_column_names <- function(filepath){
  dt<- fread(file = filepath, nrows = 1)
  
  dt = data.frame(columnnames = names(dt))
  dt
}



read_and_format_Spectral_library <- function(spectralLib_filePath,
                                             im_column_name = NULL,
                                             precursor_mz_column_name = NULL,
                                             retention_time_column_name = NULL,
                                             sequence_column_name = NULL,
                                             charge_column_name = NULL,
                                             score_wt_column_name = NULL){
  
  Spectral_Lib <- fread(spectralLib_filePath)
  
  data_mz_im <- Spectral_Lib %>%
    format_data_input(im_column_name = im_column_name ,
                      precursor_mz_column_name = precursor_mz_column_name,
                      retention_time_column_name = retention_time_column_name,
                      sequence_column_name = sequence_column_name,
                      charge_column_name = charge_column_name,
                      score_wt_column_name = score_wt_column_name) %>%
    group_by(sequence, charge) 
  
    if(!is.null(score_wt_column_name)) {data_mz_im <- data_mz_im %>% top_n(1, score_wt)}
  
  data_mz_im
    
}


Optimize_dia_PASEF_method_from_SpectralLib <- function(data_mz_im, DIA_params){
  
  # DIA_params = list(data_mz_im = data_mz_im,
  #                   im_range_lower = 0.70,
  #                   im_range_upper = 1.3,
  #                   mz_range_lower = 300,
  #                   mz_range_upper = 1200,
  #                   mz_window_width = 35,
  #                   mz_window_lower_percentile_limit = 0.1,
  #                   mz_window_upper_percentile_limit = 0.95,
  #                   im_offset_lower = 0.01,
  #                   im_offset_upper = 0.01,
  #                   Add_mz_overlaps_width = 1)
  
  
  
  
  dt <- do.call(create_DIA_PASEF_method, DIA_params)
  
  
  dt
}
plot_dia_PASEF_method_over_SpectralLib <- function(data_mz_im, diaPASEF_method, DIA_params){
  # DIA_params = list(data_mz_im = data_mz_im,
  #                   im_range_lower = 0.70,
  #                   im_range_upper = 1.3,
  #                   mz_range_lower = 300,
  #                   mz_range_upper = 1200,
  #                   mz_window_width = 35,
  #                   mz_window_lower_percentile_limit = 0.1,
  #                   mz_window_upper_percentile_limit = 0.95,
  #                   im_offset_lower = 0.01,
  #                   im_offset_upper = 0.01,
  #                   Add_mz_overlaps_width = 1)
  
  
  d <- ggplot() +
    scale_x_continuous(name="precursor_mz",n.breaks = 10) +
    scale_y_continuous(name="im",n.breaks = 20) +
    geom_density_2d_filled(data = data_mz_im, aes(precursor_mz, im, alpha = after_stat(level)), contour_var = "ndensity") + 
    # facet_wrap(vars(cut), nrow = 1)+
    scale_discrete_manual("alpha", values = c(0, rep(1, 50)))+
    theme_bw()
  
  G= d +
    geom_rect(data=diaPASEF_method, mapping=aes(xmin=start_mz, xmax=end_mz,
                                                ymin=start_im, ymax=end_im),fill="transparent", color="black") +
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle(paste0("#cycles = ", max(diaPASEF_method$cycle)))+
    theme(legend.position = "none")+
    coord_cartesian(xlim = c(DIA_params$mz_range_lower,DIA_params$mz_range_upper),
                    ylim = c(DIA_params$im_range_lower,DIA_params$im_range_upper))
  
  G
}

