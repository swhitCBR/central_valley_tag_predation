##### Implement expert-opinion rule-based filter (version 1_1) with Phase 1 data

### Author: Rebecca Buchanan, UW


  ## load in data
  {
    # User will need to update file paths
    
    # DEDF with metrics
    x1 <- readr::read_rds(file.path("data","sample_data_sets","Phase_1","v1.3","Phase_1_sds_v_1.3_METRICS.rds"))
  }
  
  
  # extract data
  {
    Nsmolt15.df <- x1[["N_smolts_2015"]]
    Nsmolt16.df <- x1[["N_smolts_2016"]]
    Ssmolt15.df <- x1[["S_smolts_2015"]]
    Ssmolt16.df <- x1[["S_smolts_2016"]]
    Npred15.df <- x1[["N_predators_2015"]]
    Npred16.df <- x1[["N_predators_2016"]]
    Spred15.df <- x1[["S_predators_2015"]]
    Spred16.df <- x1[["S_predators_2016"]]    
  }

  
  ### Run rule-based filter
  source("functions/RULE_BASED_FILTER_FXNS.R")
    # Define critical thresholds
    {
      # South Delta smolt
      {
        rb.th.ls.s<-list(time_since_release = 15*24,  # hours (max)
                              time_since_last_at_station = 12*24, # hours (max)
                              migration_rate = c(0.2, 1.5),    # km/hr (min, max)
                              BLPS = c(0.4, 5),                # (min, max) 
                              resid_time_near = 5,             # hours (max)
                              resid_time_mid = 39,             # hours (max)
                              cumul_resid_time = 15,           # hours (max)
                              
                              cumul_visit_count = 3,           # max
                              consec_visit_count = 3,          # max
                              
                              cumul_time_hrz = 117,            # hours (max)
                              cumul_down_up_switch = 3,        # max
                              cumul_against_flow = 1           # max
        )    
      }
      
      # North Delta smolt
      {
        rb.th.ls.n<-list(time_since_release = 45*24,      # hours (max)
                              time_since_last_at_station = 12*24, # hours (max)
                              migration_rate = c(0.1, 0.5),    # km/hr (min, max)
                              BLPS = c(0.2, 5),                # (min, max) 
                              resid_time_near = 10,            # hours (max)
                              resid_time_mid = 54,             # hours (max)
                              cumul_resid_time = 30,           # hours (max)
                              
                              cumul_visit_count = 3,           # max
                              consec_visit_count = 2,          # max
                              
                              cumul_time_hrz = 152,            # hours (max)
                              cumul_down_up_switch = 3,        # max
                              cumul_against_flow = 1           # max
        )    
      }
    }
    
    # Run for the smolt data 
    {
      # run
      {
        # with elapsed time measured
        {
          ss15_start<-Sys.time()
          ss15<-rb.filter(dat=Ssmolt15.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          ss15_end<-Sys.time(); ss15_elapsed<-ss15_end - ss15_start
          
          ss16_start<-Sys.time()
          ss16<-rb.filter(dat=Ssmolt16.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          ss16_end<-Sys.time(); ss16_elapsed<-ss16_end - ss16_start
          
          ns15_start<-Sys.time()
          ns15<-rb.filter(dat=Nsmolt15.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          ns15_end<-Sys.time(); ns15_elapsed<-ns15_end - ns15_start
          
          ns16_start<-Sys.time()
          ns16<-rb.filter(dat=Nsmolt16.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          ns16_end<-Sys.time(); ns16_elapsed<-ns16_end - ns16_start          
        }
        
        # without elapsed time (not implemented)
        {
          #ss15<-rb.filter(dat=Ssmolt15.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          #ss16<-rb.filter(dat=Ssmolt16.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          #ns15<-rb.filter(dat=Nsmolt15.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          #ns16<-rb.filter(dat=Nsmolt16.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
        }
        
      }
    }
    
    # Run on the predator data sets
    {
      # run
      {
        # with elapsed time
        {
          sp15_start<-Sys.time()
          sp15<-rb.filter(dat=Spred15.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          sp15_end<-Sys.time(); sp15_elapsed<-sp15_end - sp15_start
          
          sp16_start<-Sys.time()
          sp16<-rb.filter(dat=Spred16.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          sp16_end<-Sys.time(); sp16_elapsed<-sp16_end - sp16_start
          
          np15_start<-Sys.time()
          np15<-rb.filter(dat=Npred15.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          np15_end<-Sys.time(); np15_elapsed<-np15_end - np15_start
          
          np16_start<-Sys.time()
          np16<-rb.filter(dat=Npred16.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          np16_end<-Sys.time(); np16_elapsed<-np16_end - np16_start          
        }
        
        
        # without elapsed time (not implemented)
        {
          #sp15<-rb.filter(dat=Spred15.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          #sp16<-rb.filter(dat=Spred16.df,filter_criteria = rb.th.ls.s,pred_score_crit = 2, detailed_output = TRUE)
          #np15<-rb.filter(dat=Npred15.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
          #np16<-rb.filter(dat=Npred16.df,filter_criteria = rb.th.ls.n,pred_score_crit = 2, detailed_output = TRUE)
        }
        
      }
    }
    
    
    ### Summarize results
    {
      ## Remove study ID from gen_loc values for the South Delta smolt data sets
      {
        ss15$Data$gen_loc<-gsub("_SouthDelta_2015","",ss15$Data$gen_loc)
        ss15$details$gen_loc<-gsub("_SouthDelta_2015","",ss15$details$gen_loc)
        ss15$flag$gen_loc<-gsub("_SouthDelta_2015","",ss15$flag$gen_loc)
        ss15$first_pred$gen_loc<-gsub("_SouthDelta_2015","",ss15$first_pred$gen_loc)
        
        ss16$Data$gen_loc<-gsub("_SouthDelta_2016","",ss16$Data$gen_loc)
        ss16$details$gen_loc<-gsub("_SouthDelta_2016","",ss16$details$gen_loc)
        ss16$flag$gen_loc<-gsub("_SouthDelta_2016","",ss16$flag$gen_loc)
        ss16$first_pred$gen_loc<-gsub("_SouthDelta_2016","",ss16$first_pred$gen_loc)
      }
      
      ## Plot distribution of gen_loc of first_pred
      {
        library(ggplot2)
        
        # South Delta smolts
        {
          ss15_fp.df<-as.data.frame(table(ss15$first_pred$gen_loc)); names(ss15_fp.df)<-c("gen_loc","Count")
          ss15_fp.df$gen_loc<-as.character(ss15_fp.df$gen_loc)
          ss15_fp.df$year<-2015
          
          ss16_fp.df<-as.data.frame(table(ss16$first_pred$gen_loc)); names(ss16_fp.df)<-c("gen_loc","Count")
          ss16_fp.df$gen_loc<-as.character(ss16_fp.df$gen_loc)
          ss16_fp.df$year<-2016
          
          ss_fp.df<-rbind(ss15_fp.df,ss16_fp.df)
          ss_fp.df$gen_loc<-as.factor(ss_fp.df$gen_loc)
          ss_fp.df$year<-as.factor(ss_fp.df$year)
          
          p_ss <- ggplot(data=ss_fp.df, aes(x=gen_loc, y=Count, fill=year)) + geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            ggtitle("Site of first predator assignment: South Delta smolts") + 
            labs(x = "General Location", y = "Count") +
            scale_fill_brewer(palette = "Set1")
          p_ss    
        }
        
        # North Delta smolts
        {
          ns15_fp.df<-as.data.frame(table(ns15$first_pred$gen_loc)); names(ns15_fp.df)<-c("gen_loc","Count")
          ns15_fp.df$gen_loc<-as.character(ns15_fp.df$gen_loc)
          ns15_fp.df$year<-2015
          
          ns16_fp.df<-as.data.frame(table(ns16$first_pred$gen_loc)); names(ns16_fp.df)<-c("gen_loc","Count")
          ns16_fp.df$gen_loc<-as.character(ns16_fp.df$gen_loc)
          ns16_fp.df$year<-2016
          
          ns_fp.df<-rbind(ns15_fp.df,ns16_fp.df)
          ns_fp.df$gen_loc<-as.factor(ns_fp.df$gen_loc)
          ns_fp.df$year<-as.factor(ns_fp.df$year)
          
          p_ns <- ggplot(data=ns_fp.df, aes(x=gen_loc, y=Count, fill=year)) + geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            ggtitle("Site of first predator assignment: North Delta smolts") + 
            labs(x = "General Location", y = "Count") +
            scale_fill_brewer(palette = "Set1")
          p_ns    
        }
        
        # South Delta predators
        {
          sp15_fp.df<-as.data.frame(table(sp15$first_pred$gen_loc)); names(sp15_fp.df)<-c("gen_loc","Count")
          sp15_fp.df$gen_loc<-as.character(sp15_fp.df$gen_loc)
          sp15_fp.df$year<-2015
          
          sp16_fp.df<-as.data.frame(table(sp16$first_pred$gen_loc)); names(sp16_fp.df)<-c("gen_loc","Count")
          sp16_fp.df$gen_loc<-as.character(sp16_fp.df$gen_loc)
          sp16_fp.df$year<-2016
          
          sp_fp.df<-rbind(sp15_fp.df,sp16_fp.df)
          sp_fp.df$gen_loc<-as.factor(sp_fp.df$gen_loc)
          sp_fp.df$year<-as.factor(sp_fp.df$year)
          
          p_sp <- ggplot(data=sp_fp.df, aes(x=gen_loc, y=Count, fill=year)) + geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            ggtitle("Site of first predator assignment: South Delta predators") + 
            labs(x = "General Location", y = "Count") +
            scale_fill_brewer(palette = "Set1")
          p_sp    
        }
        
        # North Delta predators
        {
          np15_fp.df<-as.data.frame(table(np15$first_pred$gen_loc)); names(np15_fp.df)<-c("gen_loc","Count")
          np15_fp.df$gen_loc<-as.character(np15_fp.df$gen_loc)
          np15_fp.df$year<-2015
          
          np16_fp.df<-as.data.frame(table(np16$first_pred$gen_loc)); names(np16_fp.df)<-c("gen_loc","Count")
          np16_fp.df$gen_loc<-as.character(np16_fp.df$gen_loc)
          np16_fp.df$year<-2016
          
          np_fp.df<-rbind(np15_fp.df,np16_fp.df)
          np_fp.df$gen_loc<-as.factor(np_fp.df$gen_loc)
          np_fp.df$year<-as.factor(np_fp.df$year)
          
          p_np <- ggplot(data=np_fp.df, aes(x=gen_loc, y=Count, fill=year)) + geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            ggtitle("Site of first predator assignment: North Delta predators") + 
            labs(x = "General Location", y = "Count") +
            scale_fill_brewer(palette = "Set1")
          p_np    
        }
        
        ## plot smolts and predators on same plot
        {
          # South Delta
          {
            ss_fp.df$focal_group<-"Smolts"
            sp_fp.df$focal_group<-"Predators"
            s_fp.df<-rbind(ss_fp.df,sp_fp.df)
            s_fp.df$gen_loc<-gsub("_"," ",s_fp.df$gen_loc)
            
            # equate gen_loc between data sets (use best judgement here)
            s_fp.df$gen_loc[s_fp.df$gen_loc %in% c("Old River South","Old River at Middle River Head")]<-"Old River at Middle River Head" # 
            
            # group sites (use best judgement here)
            s_fp.df$gen_loc[s_fp.df$gen_loc %in% c("Chipps","Mallard")]<-"Chipps" 
            s_fp.df$gen_loc[s_fp.df$gen_loc %in% c("Radial Gate Downstream","Radial Gate Upstream")]<-"Radial Gates" 
            
            s_fp.df$gen_loc<-as.factor(s_fp.df$gen_loc)
            s_fp.df$focal_group<-as.factor(s_fp.df$focal_group)
            s_fp.df$group_year<-as.factor(paste(s_fp.df$focal_group,s_fp.df$year))
            
            # convert to proportion
            s_fp_prop.df<- (s_fp.df %>%
                                   group_by(year,focal_group) %>%
                                   mutate(Proportion = Count/sum(Count)))
            
            # plot count distribution
            {
              p_s <- ggplot(data=s_fp.df, aes(x=gen_loc, y=Count, fill=group_year)) + geom_bar(stat="identity") +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                ggtitle("Site of first predator assignment: South Delta") + 
                labs(x = "General Location", y = "Count") +
                scale_fill_brewer(palette = "Set1")
              p_s       
            }
            
            # plot proportion distribution
            {
              pp_s <- ggplot(data=s_fp_prop.df, aes(x=gen_loc, y=Proportion, fill=group_year)) + geom_bar(stat="identity") +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                ggtitle("Site of first predator assignment: South Delta") + 
                labs(x = "General Location", y = "Proportion") +
                scale_fill_brewer(palette = "Set1")
              pp_s   
            }
            
          }
          
          # North Delta
          {
            ns_fp.df$focal_group<-"Smolts"
            np_fp.df$focal_group<-"Predators"
            n_fp.df<-rbind(ns_fp.df,np_fp.df)
            n_fp.df$gen_loc<-gsub("_"," ",n_fp.df$gen_loc)
            
            # equate gen_loc between data sets (use best judgement here)
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("CHIPP","Chipps")]<-"Chipps" 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("Decker Is","DeckerIsland")]<-"Decker Is"
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("AbvTisdale","SR AbvTisdale")]<-"SR Above Tisdale" # 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("BlwChinaBend","SR BLWCHIBEND")]<-"SR Below ChinaBend" # 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("Blw Sutter","SR BlwSutter")]<-"SR Below Sutter" # close
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("Blw FRConf","SR BLWFEATHER")]<-"SR Below Feather" # 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("I80-50 Br","SR I-80/50Br")]<-"SR I-80/50 Br" # 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("KnightsLandingBr","SR KinghtsLanding")]<-"SR Knights Landing" # 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("RioVistBr","SR RioVista")]<-"SR Rio Vista" # 
            
            # group sites (use best judgement here)
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("Chipps","Mallard")]<-"Chipps" 
            n_fp.df$gen_loc[n_fp.df$gen_loc %in% c("Radial Gate Downstream","Radial Gate Upstream")]<-"Radial Gates" 
            
            n_fp.df$gen_loc<-as.factor(n_fp.df$gen_loc)
            n_fp.df$focal_group<-as.factor(n_fp.df$focal_group)
            n_fp.df$group_year<-as.factor(paste(n_fp.df$focal_group,n_fp.df$year))
            
            # convert to proportion
            n_fp_prop.df<- (n_fp.df %>%
                                   group_by(year,focal_group) %>%
                                   mutate(Proportion = Count/sum(Count)))
            
            # plot count distribution
            {
              p_n <- ggplot(data=n_fp.df, aes(x=gen_loc, y=Count, fill=group_year)) + geom_bar(stat="identity") +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                ggtitle("Site of first predator assignment: North Delta") + 
                labs(x = "General Location", y = "Count") +
                scale_fill_brewer(palette = "Set1")
              p_n       
            }
            
            # plot proportion distribution
            {
              pp_n <- ggplot(data=n_fp_prop.df, aes(x=gen_loc, y=Proportion, fill=group_year)) + geom_bar(stat="identity") +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                ggtitle("Site of first predator assignment: North Delta") + 
                labs(x = "General Location", y = "Proportion") +
                scale_fill_brewer(palette = "Set1")
              pp_n 
            }
            
          }  
          
        }
        
      }
      

    }
    
    ### Run code to compile results for report
    {
      source("scripts/phase_1/P1_summary_fxns.R")
      
      # report how metrics were used (they were used the same way for all data sets)
      {
        rb_metrics_tab <- get_default_metrics_tab()
        rb_metrics_tab$Level="Event"
        rb_metrics_tab$Used[!(rb_metrics_tab$Metric %in% names(rb.th.ls.s))]<-FALSE
        
        # check
        check_metrics_tab(rb_metrics_tab)
        # --> PASS
        
      }
      
      # define "metrics_flag" character vector for use in "extra_cols" (used the same metrics for all data sets)
      metrics_flag<-base::intersect(paste(names(rb.th.ls.s),"flag",sep="_"),names(ss15$details))
      
      # South Delta 2015 smolts: ss15$Data
      {
        # define column indicating whether tag was ever classified as predated
        ss15$details$pred_tag<-c(unlist(sapply(unique(ss15$details$tagID),
                                                    get_pred_tag,
                                                    dat=ss15$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="S_smolts_2015",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=ss15$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
      # South Delta 2016 smolts: ss16$details
      {
        # define column indicating whether tag was ever classified as predated
        ss16$details$pred_tag<-c(unlist(sapply(unique(ss16$details$tagID),
                                                    get_pred_tag,
                                                    dat=ss16$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="S_smolts_2016",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=ss16$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
      # North Delta 2015 smolts: ns15$details
      {
        # define column indicating whether tag was ever classified as predated
        ns15$details$pred_tag<-c(unlist(sapply(unique(ns15$details$tagID),
                                                    get_pred_tag,
                                                    dat=ns15$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="N_smolts_2015",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=ns15$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
        
        
      }
      
      # North Delta 2016 smolts: ns16$details
      {
        # define column indicating whether tag was ever classified as predated
        ns16$details$pred_tag<-c(unlist(sapply(unique(ns16$details$tagID),
                                                    get_pred_tag,
                                                    dat=ns16$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="N_smolts_2016",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=ns16$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
      
      # South Delta 2015 predators: sp15$details
      {
        # define column indicating whether tag was ever classified as predated
        sp15$details$pred_tag<-c(unlist(sapply(unique(sp15$details$tagID),
                                                    get_pred_tag,
                                                    dat=sp15$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="S_predators_2015",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=sp15$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
      # South Delta 2016 predators: sp16$details
      {
        # define column indicating whether tag was ever classified as predated
        sp16$details$pred_tag<-c(unlist(sapply(unique(sp16$details$tagID),
                                                    get_pred_tag,
                                                    dat=sp16$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="S_predators_2016",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=sp16$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
      # North Delta 2015 predators: ns15$details
      {
        # define column indicating whether tag was ever classified as predated
        np15$details$pred_tag<-c(unlist(sapply(unique(np15$details$tagID),
                                                    get_pred_tag,
                                                    dat=np15$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="N_predators_2015",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=np15$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
        
        
      }
      
      # North Delta 2016 predators: np16$details
      {
        # define column indicating whether tag was ever classified as predated
        np16$details$pred_tag<-c(unlist(sapply(unique(np16$details$tagID),
                                                    get_pred_tag,
                                                    dat=np16$details)))
        
        # create and write RDS file for output
        summarize_PF_res(
          data_ver="1.3",               
          sdsID="N_predators_2016",        
          Analyst="RAB_1",                 
          PF_method="Expert-Opinion Rule-Based Event-Level Filter",   
          metric_descript_tab = rb_metrics_tab, 
          DEDF_PF=np16$details,  
          PRED_TAG_col = "pred_tag",
          output_pth = "data/tmp_data/P1_res_tmp_out",
          PRED_EVENT_col = "predator",     
          extra_cols = c("pred_score","first_pred",metrics_flag))
      }
      
    }
    
    ### Summarize effort
    {

      el_vec<-c(el_vec_s,el_vec_p)
      # Time differences in secs
      # [1] 14.990534 21.674810 10.412777 10.778809  1.231374  1.602282  1.361755  1.172178
      
      n_vec<-c(n_vec_s,n_vec_p)
      
      summary(as.numeric(el_vec/n_vec)) # secs per row
      #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
      # 0.0005437 0.0006013 0.0010179 0.0010896 0.0015834 0.0017568 
      
      sum(el_vec) # Time difference of 63.22452 secs
      
      # South Delta 2015
      indices<-c(1,5)
      el_vec[indices] # 14.990534  1.231374
      sum(el_vec[indices]) # 16.22191
      sum(el_vec[indices])/sum(n_vec[indices]) # 0.001502307
      
      # South Delta 2016
      indices<-c(2,6)
      el_vec[indices] # 21.674810  1.602282
      sum(el_vec[indices]) # 23.27709
      sum(el_vec[indices])/sum(n_vec[indices]) # 0.00141176
      
      # North Delta 2015
      indices<-c(3,7)
      el_vec[indices] # 10.412777  1.361755
      sum(el_vec[indices]) # 11.77453
      sum(el_vec[indices])/sum(n_vec[indices]) # 0.001351685      
      
      # North Delta 2016
      indices<-c(4,8)
      el_vec[indices] # 10.778809  1.172178
      sum(el_vec[indices]) # 11.95099
      sum(el_vec[indices])/sum(n_vec[indices]) # 0.0009337438      
      
      # all data sets combined
      sum(el_vec) # 63.22452 secs
      sum(el_vec)/sum(n_vec) # 0.001295691 secs per row
      
    }
    
    
    ### develop some summary stats for results
    {
      
      # South Delta 2015 smolts
      {
        # number of tags detected
        print(nrow(ss15$first_pred))
        
        # summarize first_pred
        ss15_ddp<-summarize_firstpred(dlist=ss15, min_det=5)
        View(ss15_ddp[,c("gen_loc","print")])
      }
      
      # South Delta 2015 predators
      {
        # number of tags detected
        print(nrow(sp15$first_pred))
        
        # summarize first_pred
        sp15_ddp<-summarize_firstpred(dlist=sp15, min_det=5)
        View(sp15_ddp[,c("gen_loc","print")])
      }  
      
      # South Delta 2016 smolts
      {
        # number of tags detected
        print(nrow(ss16$first_pred))
        
        # summarize first_pred
        ss16_ddp<-summarize_firstpred(dlist=ss16, min_det=5)
        View(ss16_ddp[,c("gen_loc","print")])
      }
      
      # South Delta 2016 predators
      {
        # number of tags detected
        print(nrow(sp16$first_pred))
        
        # summarize first_pred
        sp16_ddp<-summarize_firstpred(dlist=sp16, min_det=5)
        View(sp16_ddp[,c("gen_loc","print")])
      }  
      
      
      
      # North Delta 2015 smolts
      {
        # number of tags detected
        print(nrow(ns15$first_pred))
        
        # summarize first_pred
        ns15_ddp<-summarize_firstpred(dlist=ns15, min_det=5)
        View(ns15_ddp[,c("gen_loc","print")])
      }
      
      # North Delta 2015 predators
      {
        # number of tags detected
        print(nrow(np15$first_pred))
        
        # summarize first_pred
        np15_ddp<-summarize_firstpred(dlist=np15, min_det=5)
        View(np15_ddp[,c("gen_loc","print")])
      }  
      
      # North Delta 2016 smolts
      {
        # number of tags detected
        print(nrow(ns16$first_pred))
        
        # summarize first_pred
        ns16_ddp<-summarize_firstpred(dlist=ns16, min_det=5)
        View(ns16_ddp[,c("gen_loc","print")])
      }
      
      # North Delta 2016 predators
      {
        # number of tags detected
        print(nrow(np16$first_pred))
        
        # summarize first_pred
        np16_ddp<-summarize_firstpred(dlist=np16, min_det=5)
        View(np16_ddp[,c("gen_loc","print")])
      }  
      
    }
    
    








