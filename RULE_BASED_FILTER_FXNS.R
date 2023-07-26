###### Functions for Rule-Based Filter

### Author: Rebecca Buchanan, UW


###### Contents

# rb.filter:     Rule-based predator filter for data set with multiple tags
# rb.filter.tag: Rule-based predator filter for single tag
# test_crit:     Test whether observed metric in x exceeds threshold value defined in limit_vec of type limit_type
# get_pred_tag:  Identify whether tag was ever classified as predated from event-level classification
# summarize_firstpred:   Summarize the detection counts and number of tags diagnosed as first predated at each site

####################################################################################################################

rb.filter<-function(dat,filter_criteria,criteria_min=NULL,pred_score_crit=1,detailed_output=FALSE,extra_cols=c("studyID","year","species"))
{
  ## Rule-based predator filter for data set with multiple tags
  
  # Necessary packages:
  library(dplyr)  
  
  ## Arguments:
  
  # dat = data.frame or data.table of detection event data on general location spatial scale; see Details
  # filter_criteria = named list of threshold criteria for metrics in data
  # criteria_min = character vector naming the metrics for which the threshold in filter_criteria is a minimum for smolt-like behavior; see Details
  # pred_score_crit = minimum predator score required for predator classification; default value = 1 (i.e., violation of any criterion results in predator classification)
  # detailed_output = logical; if TRUE, identifies metrics that violated the threshold criteria and where they were violated; see Value
  # extra_cols = additional columns from dat to be included in output
  
  ## Value:
  
  # if detailed_output = FALSE:
  #    data.frame identical to dat but with additional fields recording the predator score and predator designation for each detection event per tag; 
  #
  #    additional fields = pred_score, smolt, predator
  #    where 
  #           pred_score = numeric count of number of metrics that violated criteria for detection event,
  #           smolt = !predator
  #           predator = logical indicator of whether tag is classified as predated by the start of the detection event
  #
  # if detailed_output = TRUE:
  #    a list with components:
  # 
  #    'Data' = the output for detailed_output=FALSE
  #    'details' = data frame with identical structure to 'Data' with a second field for each metric ('[metric]_flag') indicating whether it exceeded of relevant f_crit criterion (TRUE/FALSE); 
  #              omits metrics not tested in filter_crit;
  #              includes field first_pred = logical; indicates whether tag was first classified as predated at current detection event
  #    'flag' = data frame in 'details' limited to detection events in which at least one criterion in filter_crit is exceeded (i.e., pred_score>0);
  #           tags with no filter_crit violation are included with detection fields assigned NA;
  #           limited to fields: tagID, gen_loc, [extra_cols], [metric_flag fields limited to the metrics included in f_crit], pred_score, smolt, predator, first_pred
  #    'first_pred' = data frame in 'flag' limited to first detection event in which number of criteria exceeded is >= pred_score_crit (i.e., pred_score>=pred_score_crit);
  #                 tags that have no detection event with pred_score>=pred_score_crit are included with detection fields assigned NA 
  
  ## Function type: exposed to user
  
  ## Details:
  
  # Metrics defined for each general location detection event are compared to user-specified threshold criteria.
  
  # If an observed metric exceeds the criterion of that metric, the predator score (pred_score) is incremented by 1 for that detection event.
  
  # When the value of pred_score reaches pred_score_crit, the tag is classified as in a predator and the 'predator' indicator variable is set to TRUE for the 
  # current detection event and all future detection events.
  
  # The filter_criteria argument is a named list that defines threshold criteria for predator designation for any of the metrics included as fields in dat
  
  # The threshold values in f_crit are either numeric scalars or else vectors of length 2. Scalar values are assumed by default to be 
  # maximum values consistent with smolt-like behavior. Any scalar thresholds that instead represent minimum values must be named in crit_min. 
  # 2-element vectors are assumed to be (min, max).
  
  # Required fields in dat:
  #    tagID = alphanumeric code representing tag ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    ... metrics named in the filter_criteria argument
  
  
  
  # check filter_criteria versus available metrics in dat
  missing_metrics<-names(filter_criteria)[!(names(filter_criteria) %in% names(dat))]
  if(length(missing_metrics)>0) stop(paste("dat is missing metrics used in filter_criteria ",paste(missing_metrics,sep=", "),sep=": "))
  
  
  # define vector of tag ID
  tags<-unique(dat$tagID)
  
  
  # apply rule-based filter to each tag
  x<-sapply(tags, rb.filter.tag, Data=dat, f_crit=filter_criteria, crit_min=criteria_min, pred_score_min=pred_score_crit, detailed_output=detailed_output, extra_cols=extra_cols)  
  
  
  # compile into single data frame(s)
  {
    if(!detailed_output) 
    {
      Data<-as.data.frame(do.call(rbind,x))
      row.names(Data)<-1:nrow(Data)
      out<-Data
      
    } else {
      
      # identify fields to define for all tags
      base_fields<-c("tagID","det_event_ID","gen_loc",extra_cols)
      
      
      ## Data
      {
        Data<-as.data.frame(do.call(rbind,x["Data",]))
        row.names(Data)<-1:nrow(Data)
        
        # create sub-dataframe that is only the first row for each tag 
        first.rows.df<-do.call(rbind,lapply(split(Data,Data$tagID),function(x) x[1,]))
        
        # identify the tags that were never flagged
        tag.noflag<-tags[!(tags %in% unique(Data$tagID[Data$pred_score>0]))]
        
        # identify the tags that were never diagnosed as predators
        tag.nopred<-tags[!(tags %in% unique(Data$tagID[Data$predator]))]        
      }
      
      
      ## details
      details.df<-do.call(rbind,x["details.df",])
      row.names(details.df)<-1:nrow(details.df)
      
      
      ## flag
      flag.df<-do.call(rbind,x["flag.df",])
      if(ncol(flag.df)==1)
      {
        # no tags were flagged for violating any criteria
        flag.df<-first.rows.df[,intersect(base_fields,names(details.df))]
      } else 
      {
        flag.df[is.na(flag.df$tagID),intersect(base_fields,names(flag.df))]<-first.rows.df[first.rows.df$tagID %in% tag.noflag,intersect(base_fields,names(flag.df))]
        
        # limit the "metrics" fields to the TRUE/FALSE values indicating being flagged by predator filter
        flag.df<-flag.df[,intersect(c(base_fields,"pred_score","smolt","predator","first_pred",names(flag.df)[grep("flag",names(flag.df))]),names(flag.df))]
        names(flag.df)[grep("flag",names(flag.df))]<-gsub("_flag","",names(flag.df)[grep("flag",names(flag.df))])
      }  
      
      row.names(flag.df)<-1:nrow(flag.df)
      
      ## first_pred
      first_pred.df<-do.call(rbind,x["first_pred.df",])
      if(ncol(first_pred.df)==1)
      {
        # no tags were diagnosed as predators
        first_pred.df<-first.rows.df[,intersect(base_fields,names(details.df))]
      } else {
        
        first_pred.df[is.na(first_pred.df$tagID),intersect(base_fields,names(first_pred.df))]<-first.rows.df[first.rows.df$tagID %in% tag.nopred,intersect(base_fields,names(first_pred.df))]
        
        # limit the "metrics" fields to the TRUE/FALSE values indicating being flagged by predator filter
        first_pred.df<-first_pred.df[,intersect(c(base_fields,"pred_score","smolt","predator","first_pred",names(first_pred.df)[grep("flag",names(first_pred.df))]),names(first_pred.df))]
        names(first_pred.df)[grep("flag",names(first_pred.df))]<-gsub("_flag","",names(first_pred.df)[grep("flag",names(first_pred.df))])
        
      }
      row.names(first_pred.df)<-1:nrow(first_pred.df) 
      
      ## define output list
      out<-list(Data=Data,
                details=details.df,
                flag=flag.df,
                first_pred=first_pred.df)
    }
  }
  
  # return output
  return(out)
}


rb.filter.tag<-function(tag,Data,f_crit,crit_min=NULL,pred_score_min=1,detailed_output=FALSE,extra_cols=c("studyID","year","species"))
{
  ## Rule-based predator filter for single tag
  
  # Necessary packages: NA
  
  ## Arguments:
  # tag = tag ID to implement; alphanumeric code 
  # Data = data.frame or data.table of detection event data on general location spatial scale; may include multiple tags
  # f_crit = named list of threshold criteria for metrics in data
  # crit_min = character vector naming the metrics for which the threshold in f_crit is a minimum for smolt-like behavior; see Details
  # pred_score_min = minimum predator score required for predator classification
  # detailed_output = logical; if TRUE, identifies metrics that violated the threshold criteria and where they were violated; see Details
  # extra_cols = additional columns from dat to be included in output
  
  
  ## Value:
  # data.frame recording the predator score and predator designation for each detection event per tag; 
  # see rb.filter for more information
  
  #  where pred_score = numeric count of number of metrics that violated criteria for detection event,
  #  smolt = !predator
  #  and predator = logical indicator of whether tag is classified as predated by the start of the detection event
  
  ## Function type: internal
  #  used in rb.filter
  
  ## Details: 
  # method used in function rb.filter
  # see rb.filter for more information  
  
  # The threshold values in f_crit are either numeric scalars or else vectors of length 2. Scalar values are assumed by default to be 
  # maximum values consistent with smolt-like behavior. Any scalar thresholds that instead represent minimum values must be named in crit_min. 
  # 2-element vectors are assumed to be (min, max).
  
  # Names of f_crit elements must correspond to the metrics
  
  # detailed output includes
  # (1) 'Data': input data frame with additional fields: pred_score, smolt, predator
  # (2) 'details': data frame with identical structure to Data with a second field for each metric ('[metric]_flag') indicating whether it exceeded of relevant f_crit criterion (TRUE/FALSE); 
  #     omits metrics not tested in f_crit;
  #     includes field first_pred = logical; indicates whether tag was first classified as predated at current detection event
  # (3) 'flag': data frame in (2) limited to detection events in which at least one criterion in f_crit is exceeded (i.e., pred_score>0);
  #     value = NA if none exceeded; 
  #     limited to fields: tagID, gen_loc, [extra_cols], [metric_flag fields limited to the metrics included in f_crit], pred_score, smolt, predator, first_pred
  # (4) 'first_pred': data frame in (3) limited to first detection event in which number of criteria exceeded is >= pred_score_min (i.e., pred_score>=pred_score_min);
  #     value = NA if pred_score < pred_score_min for all detection events
  
  
  # Restrict data to tag
  Data<-Data[Data$tagID==tag,,drop=F]  
  
  # Set up data structure
  nhits<-nrow(Data) 
  pred_score<-rep(0,length=nhits)  # count of metrics that violate criteria for each detection event
  smolt<-rep(TRUE,length=nhits)    # tag is assumed to be un-predated at the start of the detection history
  
  # identify metrics and thresholds to be used
  metrics<-names(f_crit)
  threshold_crit<-unlist(f_crit) # vector of numeric threshold criteria; if observed metric > threshold_crit value (or < for metrics in crit_min), pred_score will be incremented by 1 for detection event
  
  #threshold_type<-rep("max",length(f_crit)); threshold_type[metrics %in% crit_min]<-"min"
  threshold_type<-rep("max",sum(c(unlist(lapply(f_crit,length)))))
  threshold_type[names(threshold_crit) %in% crit_min]<-"min"
  n_threshold<-c(unlist(lapply(f_crit,length)))
  threshold_type[which(n_threshold==2)+c(1:sum(n_threshold==2)-1)]<-"min"
  
  metrics<-rep(names(f_crit),times=n_threshold)
  
  # change metrics fields to numeric
  Data<-Data %>% mutate_at(metrics,as.numeric)
  
  # identify metrics and detection events that exceed threshold criteria
  exceed_crit<-apply(Data[,metrics],1,test_crit,limit_vec=threshold_crit,limit_type=threshold_type,simplify=FALSE)
  if(length(exceed_crit)==0) exceed_crit<-vector(mode='list',length=nhits)
  
  # calculate predator score for each detection event
  pred_score<-unlist(lapply(exceed_crit,length))
  
  # identify detection event where predator score exceeds maximum allowable for smolt
  first_pred<-ifelse(length(which(pred_score>=pred_score_min))>0,which(pred_score>=pred_score_min),NA)
  
  # update smolt classification and define predator classification: is the tag predated at the start of each detection event
  if(!is.na(first_pred)) smolt[first_pred:nhits]<-FALSE
  pred<-!smolt                    
  
  # store predator filter output in Data  
  Data$pred_score<-pred_score
  Data$smolt<-smolt
  Data$predator<-pred
  
  out<-list(Data=Data)
  
  if(detailed_output)
  {
    # Identify fields to include in all output
    base_fields<-c("tagID","det_event_ID","gen_loc",extra_cols)
    base_fields<-intersect(base_fields,names(Data))
    
    # Identify detection event where tag was first classified as predated
    Data$first_pred<-FALSE; if(!is.na(first_pred)) Data$first_pred[first_pred]<-TRUE
    
    # Identify metrics that were flagged as violating criteria (includes events with pred_score < minimum predator score)
    details.df<-Data; #details.df<-details.df[,c(base_fields,unique(metrics),"pred_score","smolt","predator","first_pred")]
    details.df[,unique(metrics)]<-FALSE
    for(i in 1:nhits) details.df[i,unlist(lapply(exceed_crit,names)[i])]<-TRUE
    names(details.df)[names(details.df) %in% unique(metrics)]<-paste(unique(metrics),"flag",sep="_")
    
    #metrics.obs.df<-Data[,unique(metrics)]
    details.df<-cbind(details.df,Data[,unique(metrics)])
    #details.df<-details.df[,c(base_fields,unique(metrics),"pred_score","smolt","predator","first_pred",paste(unique(metrics),"flag",sep="_"))]
    details.df<-details.df[,c(names(Data),paste(unique(metrics),"flag",sep="_"))]
    
    # Limit details.df to only detection events with flagged metrics (includes events with pred_score < minimum predator score)
    if(max(pred_score)==0) flag.df<-NA else flag.df<-details.df[pred_score>0,]
    #flag.df<-flag.df[,c(base_fields,unique(metrics),"pred_score","smolt","predator","first_pred",paste(unique(metrics),"flag",sep="_"))]
    
    # Limit details.df to the first detection event with pred_score >= minimum predator score
    if(is.na(first_pred)) first_pred.df<-NA else first_pred.df<-details.df[details.df$first_pred,,drop=F]
    #first_pred.df<-first_pred.df[,c(base_fields,unique(metrics),"pred_score","smolt","predator","first_pred",paste(unique(metrics),"flag",sep="_"))]
    
    
    out<-c(out,list(details.df=details.df,flag.df=flag.df,first_pred.df=first_pred.df))
  }
  
  return(out)  
}


test_crit<-function(x,limit_vec,limit_type)
{
  ## Test whether observed metric in x exceeds threshold value defined in limit_vec of type limit_type
  
  ## Necessary packages: NA
  
  ## Arguments:
  # x = vector of numerical observations of one or more metrics
  # limit_vec = vector of numerical thresholds for each metric in x, in same order as in x
  # limit_type = character vector defining each entry in limit_vec as either "max" or "min"; entries must be in same order as in limit_vec
  
  ## Value: named vector of the metrics that exceed threshold
  
  ## Function type: internal
  #  used in: rb.filter.tag
  
  ## Details:
  # metric exceeds treshold value if x>limit_vec for limit_type="max", or if x<limit_vec for limit_type="min"
  
  
  tmp<-(x>limit_vec)
  tmp[limit_type=="min"]<-(x<limit_vec)[limit_type=="min"]
  
  return(which(tmp))
}

get_pred_tag<-function(tag,dat)
{
  ### Identify whether tag was ever classified as predated from event-level classification
  
  ## Necessary packages: NA
  
  ## Arguments:
  
  # tag = character string identifying the tag
  # dat = 'Data' output of rb.filter() function
  
  ## Value: logical vector of length equal to the number of rows in dat with observations from tag,
  #         each element of the vector indicates whether the tag was ever classified as a predator (all elements are identical)
  
  dat<-dat[dat$tagID==tag,,drop=F]
  return(rep(dat$predator[nrow(dat)],nrow(dat)))
}



summarize_firstpred<-function(dlist,min_det)
{
  ### Summarize the detection counts and number of tags diagnosed as first predated at each site
  
  ## Necessary packages: NA
  
  # Arguments:
  
  # dlist = output of rb.filter() with detailed_output = TRUE
  # min_det = minimum number of detections required at site to be included in output
  
  # Value: dataframe with columns:
  #   gen_loc = general location
  #   det = number of tags detected at the general location
  #   first_pred = number of tags first classified as predated at the general location
  #   prop.fp = proportion of tags detected at the general location that were first classified as predated there
  #   pct.fp = prop.fp converted to percentage, rounded to nearest 1%
  #   print = first_pred pasted with pct.fp, suitable for printing
  
  
  # number of detections per general location
  dloc<-table(dlist$Data$gen_loc)
  
  # number of first_preds at general location
  ploc<-table(dlist$first_pred$gen_loc)
  
  # percentage of detections where tag was first classified as predated
  ddp.df<-data.frame(gen_loc=names(dloc),det=as.numeric(dloc),first_pred=NA)
  ddp.df$gen_loc<-as.character(ddp.df$gen_loc)
  ddp.df$first_pred[ddp.df$gen_loc %in% names(ploc)]<-ploc[ddp.df$gen_loc[ddp.df$gen_loc %in% names(ploc)]]
  ddp.df$first_pred[is.na(ddp.df$first_pred)]<-0
  ddp.df$prop.fp<-ddp.df$first_pred/ddp.df$det
  ddp.df$pct.fp<-round(ddp.df$prop.fp*100,0)
  ddp.df$print<-paste(ddp.df$first_pred," (",format(ddp.df$pct.fp,0),"%)",sep="")
  
  # limit to those with at least min_det detections
  if(missing(min_det)) min_det<-1
  return(ddp.df[ddp.df$det>=min_det,])
}


