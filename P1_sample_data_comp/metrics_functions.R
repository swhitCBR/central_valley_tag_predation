### functions for computing metrics

# used in function compute_metrics()


### Contents
# get_first_event                  # Add a logical (TRUE/FALSE) column indicating the first detection     
# get_det_event_ID                 # Add a integer column indicating the order of events for each tagIDs detection history, ordered by 'first_detection'
# get_time_since_release           # Calculate time since release    
# count_unique_visits              # Calculate running count of unique visits to each station
# get_time_since_last_at_station   # Calculate time since the previous visit at the station (general location)
# get_resid_time_near              # Calculate near-field residence time = duration of detection event
# assign_visit_events_MOD          # Calculate running count of unique visit events to each station, unbroken by detections at other stations
# get_resid_time_mid               # Calculate mid-field residence time = cumulative duration from first to lst detections of sequence of consecutive detection events at the same general location
# get_consec_visit_count           # Calculate running count of consecutive unbroken detection events (visits) at each station
# get_cumul_resid_time             # Calculate cumulative near field residence time at each station
# get_migration_rate               # Calculate migration rate for transition
# get_blps                         # Calculate BLPS (body lengths per second) for transition
# get_against_flow                 # Calculate logical indicator of transition directed against the direction of water flow
# get_cumul_down_up_switch_tag     # Calculate cumulative count of switches from moving downstream to moving upstream for a single tag
# get_cumul_down_up_switch         # Calculate cumulative count of switches from moving downstream to moving upstream for all tags
# get_against_flow                 # Calculate cumulative count of against-flow transitions, defined for every transition

####################################################################################################################


get_first_event <- function(dat){
  ## Description: 
  # Add a logical (TRUE/FALSE) column indicating the first detection
  
  ## Necessary packages:
  require(dplyr)
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: time_since_release
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "first_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  
  dat <- dat %>%
    arrange(tagID,first_detection) %>%
    group_by(tagID) %>%
    mutate(first_event=row_number()==1)
  
  # Author: Steve Whitlock , UW, SAFS, CBR

  ## Internal vs External function
  # External
  
}


get_det_event_ID <- function(dat){
  ## Description: 
  # Add a integer column indicating the order of events for each tagIDs detection history, ordered by 'first_detection'
  
  ## Necessary packages:
  require(dplyr)
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: time_since_release
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "first_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  
  dat <- dat %>% arrange(tagID,first_detection) %>%
    group_by(tagID) %>% mutate(det_event_ID=row_number())
  
  # Author: Steve Whitlock , UW, SAFS, CBR
  
  ## Internal vs External function
  # External
}



### Calculate time since release
get_time_since_release <- function(dat, unit = "hours")
{
  ## Description: 
  # Calculate time since release to the start of each detection event
  
  ## Necessary packages:
  require(dplyr)
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: time_since_release
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection", "rls_time"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  #    rls_time = POSIXct or POSIXlt; date and time of either actual or virtual release of tag
  
  # Assumes that separate detection events do not overlap in time
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal vs External function
  # External
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID",
            "studyID",
            "gen_loc",
            "first_detection",
            "last_detection",
            "rls_time") %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID', 
         'studyID', 'gen_loc', 'first_detection','last_detection','rls_time'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection)){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }

  if(!is.POSIXt(dat$rls_time) | !is.POSIXct(dat$rls_time)){
    stop("Error: rls_time must be of class POSIXct or POSIXt")
  }  
  
  
  ### Compute metric
  
  dat %>%
    arrange(tagID, first_detection) %>%
    mutate(time_since_release = difftime(first_detection,
                                         rls_time,
                                         unit = unit)) %>%
    return(dat)
  
}


### Calculate time between detection events
get_time_btw <- function(dat,d_mat,includeLF=FALSE)
{
  ## Description: 
  # Calculate the value for the time between detections to be used to calculate migration rate and other metrics
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies: 
  #  
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # d_mat = table of transitions and and estimated distances (distance matrix); see Details

  ## Value:
  # dat data frame with additional difftime fields: 'time_btw' and optionally 'time_btw_trans_LF'; see Details

  ## Details:
  #
  # This function assumes that a valid table of transitions (d_mat) has already been calculated
  #
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  # 
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # visit_event = sequence of one or more detection events at a given station that are unbroken by detections at other stations
  # If tag has only a single detection event (i.e., row of dat) at station (general location), then visit_event = visit = detection event
  # If tag has multiple detection events (i.e., rows of dat) at station that are not separated by detections at other stations, then
  #  visit_event includes multiple visits
  # visit = detection event = row of dat
  
  # Author: Steve Whitlock, UW, SAFS, CBR
  
  
  ## Internal or External:
  # Internal (get_time_btw)
  
  
  ##########################################################
  
  ### Error checking
  
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  ### Compute metric
  
  dat <- dat %>% 
    arrange(tagID,first_detection) %>%
    group_by(tagID) %>%
    mutate(first_event=row_number()==1,
           time_btw=case_when(first_event ~ difftime(first_detection,rls_time,units="hours"),
                              .default = difftime(first_detection,lag(first_detection),units = "hours"))) 
  
  if(includeLF){
    dat <- dat %>% 
    left_join(d_mat %>% 
                dplyr::select(transID,time_btw_trans_LF))}
  
    return(dat)

}



### Calculate running count of unique visits to each station
count_unique_visits <- function(dat)
{
  ## Description: 
  # Calculate running count of unique visits to each station
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: cumul_visit_count
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # cumul_visit_count = running total of number of detection events ("visits") to station for each tag
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal vs External function
  # External
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID', 'first_detection', 'last_detection', 'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection)  | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  
  ### Compute metric
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID, gen_loc) %>% 
    dplyr::mutate(cumul_visit_count = 1:n()) %>% 
    ungroup() %>% 
    return()
}


### Calculate time since the previous visit at the station (general location)
get_time_since_last_at_station <- function(dat, unit = "hours")
{
  ## Description: 
  # Calculate time since the previous visit at the station (general location)
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # unique_visits()
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: time_since_last_at_station
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event

  # Assumes that separate detection events do not overlap in time
  
  # produces NAs for stations with only a single detection event
  
  # time lag since last detection is measured as "first to first": 
  #   = the difference between the start of the current event and the start of the most recent event at the same station
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal vs External function
  # External  

  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            "first_detection", 
            "last_detection") %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 
          'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  ### Compute metric
  
  if(!("cumul_visit_count" %in% names(dat))){
    dat = dat %>% 
      unique_visits() %>% 
      arrange(tagID, gen_loc, cumul_visit_count) %>% 
      group_by(tagID, gen_loc) %>% 
      mutate(time_since_last_at_station = difftime(first_detection, 
                                                   lag(first_detection), 
                                                   unit = unit)) %>% 
      ungroup() %>%
      arrange(tagID, first_detection) %>%
      dplyr::select(-cumul_visit_count) 
      return(dat)
  }  
  

   dat %>% 
      arrange(tagID, gen_loc, cumul_visit_count) %>% 
      group_by(tagID, gen_loc) %>% 
      mutate(time_since_last_at_station = difftime(first_detection, 
                                                   lag(first_detection), 
                                                   unit = unit)) %>% 
      ungroup() %>%
      arrange(tagID, first_detection) %>%
      return()

}


### Calculate near-field residence time
get_resid_time_near <- function(dat, unit = "hours")
{
  ## Description: 
  # Calculate near-field residence time = duration of detection event
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: resid_time_near
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # detection event duration calculated as first to last: difference between last_detection and first_detection of the detection event
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal vs External function
  # External
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  
  ### Compute metric
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID, gen_loc) %>% 
    mutate(resid_time_near = difftime(last_detection, 
                                      first_detection, 
                                      unit = unit)) %>% 
    ungroup() %>% 
    return()  
  
  
}


### Calculate running count of unique visit events to each station, unbroken by detections at other stations
assign_visit_events_MOD <- function(dat)
{
  ## Description: 
  # Calculate running count of unique visit events to each station, unbroken by detections at other stations
  
  ## Necessary packages:
  require(dplyr)
  require(tidyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details

  ## Value:
  # dat data frame with additional field: visit_event
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # visit_event = sequence of one or more detection events at a given station that are unbroken by detections at other stations
  # If tag has only a single detection event (i.e., row of dat) at station (general location), then visit_event = visit = detection event
  # If tag has multiple detection events (i.e., rows of dat) at station that are not separated by detections at other stations, then
  #  visit_event includes multiple visits
  # visit = detection event = row of dat
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # Internal (get_resid_time_mid)
  
  
  ##########################################################
  
  ### Error checking
  
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  ### Compute metric
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID) %>% 
    mutate(diff_loc_ind = as.integer(gen_loc != 
                                       lag(gen_loc))) %>% 
    tidyr::replace_na(list(diff_loc_ind = 1)) %>% 
    ungroup() %>%
    group_by(tagID, gen_loc) %>%
    mutate(visit_event = cumsum(diff_loc_ind)) %>% 
    dplyr::select(-diff_loc_ind) %>% 
    ungroup() %>% 
    return()
  
  
}


### Calculate mid-field residence time = cumulative duration from first to lst detections of sequence of consecutive detection events at the same general location
get_resid_time_mid <- function(dat, unit = "hours")
{
  ## Description: 
  # Calculate mid-field residence time = cumulative duration from first to lst detections of sequence of consecutive detection events at the same general location
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # assign_visit_events
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  # unit = time unit for output
  
  ## Value:
  # dat data frame with additional field: resid_time_mid
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # output metric resid_time_mid = total duration of consecutive unbroken detection events at single location, including time between events;
  # has unique value for each detection event; 
  # incremental measure of mid-field residence time through the current detection event;
  # restarts from 0 for every detection event at a separate station (general location);
  # = resid_time_near for the first detection event of a visit_event (see get_resid_time_near(), assign_visit_events())
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  
  ### Compute metric
  
  if("visit_event" %in% names(dat)){
    dat = (dat %>% 
      arrange(tagID, first_detection) %>% 
      group_by(tagID) %>% 
      group_by(tagID, gen_loc, visit_event) %>%
      mutate(resid_time_mid = difftime((last_detection), # use max(last_detection) if want a single value for resid_time_mid per visit_event rather than an incremental metric
                                       min(first_detection), 
                                       unit = unit)) %>% 
      ungroup()) # %>% 
      return(dat)
  }
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID) %>% 
    assign_visit_events() %>%  
    group_by(tagID, gen_loc, visit_event) %>%
    mutate(resid_time_mid = difftime((last_detection), # use max(last_detection) if want a single value for resid_time_mid per visit_event rather than an incremental metric
                                     min(first_detection), 
                                     unit = unit)) %>% 
    ungroup() %>% 
    return()
  
  
}


### Calculate running count of consecutive unbroken detection events (visits) at each station
get_consec_visit_count <- function(dat)
{
  ## Description: 
  # Calculate running count of consecutive unbroken detection events (visits) at each station
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # assign_visit_events
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details

  ## Value:
  # dat data frame with additional field: consecutive_visits
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # output metric consecutive_visits is a running tally of consecutive visits at each station (general location);
  # tally restarts at 1 when tag moves to new station;
  # increases incrementally for each repeated detection at the same station if not detected elsewhere in between.
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  
  ##########################################################
  
  ### Error checking

  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt")
  }
  
  
  ### Compute metric
  
  if("visit_event" %in% names(dat)){
    dat = (dat %>% 
      arrange(tagID, first_detection) %>% 
      group_by(tagID) %>% 
      group_by(tagID, gen_loc, visit_event) %>%
      dplyr::mutate(consec_visit_count = 1:n()) %>%  # use = n() if want total number of consecutive vists for each detection event within the visit_event
      ungroup()) # %>% 
      return(dat)
  }
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID) %>% 
    assign_visit_events() %>%  
    group_by(tagID, gen_loc, visit_event) %>%
    dplyr::mutate(consec_visit_count = 1:n()) %>%  # use = n() if want total number of consecutive vists for each detection event within the visit_event
    ungroup() %>% 
    return()
}


### Calculate cumulative near field residence time at each station
get_cumul_resid_time <- function(dat)
{
  ## Description: 
  # Calculate cumulative near field residence time at each station
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # resid_time_near
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: cumul_resid_time
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # output metric cumul_resid_time is the cumulative near-field residence time at each station;
  # value for each station starts from 0
  # value = resid_time_near for the first detection event at a given station and increases incrementally thereafter
  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  ### Compute metric
  
  if("resid_time_near" %in% names(dat)){
    dat<- (dat %>% 
      arrange(tagID, first_detection) %>% 
      group_by(tagID, gen_loc) %>% 
      mutate(cumul_resid_time = cumsum(as.numeric(resid_time_near))) %>% 
      ungroup()) # %>% 
      return(dat)
  }
  
  dat %>% 
    arrange(tagID, first_detection) %>% 
    resid_time_near() %>% 
    group_by(tagID, gen_loc) %>% 
    mutate(cumul_resid_time = cumsum(as.numeric(resid_time_near))) %>% 
    ungroup() %>% 
    return()
  
}


### Calculate migration rate for transition
get_migration_rate <- function(dat)
{
  ## Description: 
  # Calculate migration rate for transition
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details

  ## Value:
  # dat data frame with additional field: migration_rate, mig_rate_unit
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "last_detection", "time_btw", "trans_len"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    last_detection = POSIXct or POSIXlt; date and time of last detection of detection event
  #    time_btw = difftime object; time lag between end of previous detection event and start of current detection event
  #    trans_len = numeric; distance in kilometers between general location of previous detection event and current detection event
  
  # Assumes that separate detection events do not overlap in time
  
  # units = kilometers per [time_btw unit]
  
  # Authors: Dalton Hance, USGS, WFRC; Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'last_detection') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'last_detection',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection) | 
     !is.POSIXt(dat$last_detection) | !is.POSIXct(dat$last_detection) ){
    stop("Error: first_detection/last_detection must be of class POSIXct or POSIXt ")
  }
  
  if(!all(c('time_btw','trans_len') %in% names(dat))){
    stop("Error: data must minimally contain variables 'time_btw', 'trans_len'")
  }
  
  
  ### Compute metric
  
  dat <- (dat %>% 
    arrange(tagID, first_detection) %>% 
    group_by(tagID) %>% 
    mutate(migration_rate = trans_len/as.numeric(time_btw),
           mig_rate_unit = paste0("km/", units(time_btw))) %>% 
    ungroup()) # %>% 
    return(dat)  
}


### Calculate BLPS (body lengths per second) for transition
get_blps <- function(dat)
{
  ## Description: 
  # Calculate BLPS (body lengths per second) for transition
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: BLPS
  # - fields migration_rate and mig_rate_unit will be include in output dat if not already defined in input
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "first_detection", "gen_loc", "fork_len"
  #    and either: "time_btw" and "trans_len" or "migration_rate" and "mig_rate_unit"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    gen_loc  = character; name of telemetry station (general location)
  #    fork_len = numeric; fork length of fish at tagging (mm)
  #    time_btw = [MODE]; time lag between end of previous detection event and start of current detection event
  #    trans_len = numeric [??]; distance in kilometers (??) between general location of previous detection event and current detection event
  #    migration_rate = numeric; ratio of transition distance to transition time lag (=trans_len/time_btw); 
  #    mig_rate_unit = units for migration_rate; = units(trans_len)/units(time_btw); 
  #        default value = km/hr
  #        valid distance units = m or km
  #        valid time units = secs, mins, or hours

  
  # Assumes that separate detection events do not overlap in time
  
  ## Assumes that fork_len unit is mm
  
  # Authors: Dalton Hance, USGS, WFRC; Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 
            'fork_len') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID', 'gen_loc', 'first_detection', 'fork_len")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection)){
    stop("Error: first_detection must be of class POSIXct or POSIXt ")
  }
  
  if((!all(c('time_btw','trans_len') %in% names(dat))) | (!all(c('migration_rate','mig_rate_unit') %in% names(dat))) ){
    stop("Error: data must minimally contain either variables 'time_btw' and 'trans_len' or 'migration_rate' and 'mig_rate_unit'")
  }
  
  if("mig_rate_unit" %in% names(dat))
  {
    tmp <- strsplit(as.character(unique(dat$mig_rate_unit)),"/")[[1]]
    dist_unit <- tmp[1]
    time_unit <- tmp[2]
    
    if(length(intersect(c("km","m"), c(dist_unit))) == 0) stop("distance component of mig_rate_unit must be either km or m")
    if(length(intersect(c("hours","mins","secs"), c(time_unit))) == 0) stop("time component of mig_rate_unit must be either hours, mins, or secs")
  }
  
  if(!("mig_rate_unit" %in% names(dat)) & "time_btw" %in% names(dat))
  {
    time_unit <- unique(units(dat$time_btw))
    if(!(time_unit %in% c("hours","mins","secs"))) stop("time_btw units must be either hours, mins, or secs")
  }
  
  
  ### Compute metric
  
  if(all(c("migration_rate","mig_rate_unit") %in% names(dat))){

    # have already defined dist_unit and time_unit above under 'Error checking'
    
    mm_per_dist_unit <- (1000*1)*(10^(ifelse(dist_unit=="km",3,0)))
    sec_per_time_unit <- ifelse(time_unit=="secs",1,ifelse(time_unit=="mins",60,60*60))
    
    dat <- (dat %>% 
              arrange(tagID, first_detection) %>% 
              group_by(tagID) %>% 
              mutate(BLPS = (as.numeric(migration_rate)/as.numeric(fork_len))*(mm_per_dist_unit)/(sec_per_time_unit)) %>% 
              ungroup()) # %>% 
    return(dat)  
    
    
    #BLPS = as.numeric(migration_rate)/as.numeric(fork_len)*(mm_per_dist_unit)/(sec_per_time_unit)
    
  }
  
  
  ## do not have migration_rate and mig_rate_unit defined, but do have trans_len and time_btw
  
  # define unit conversion from migration rate to BLPS
  # have already defined time_unit under 'Error checking'

  mm_per_dist_unit <- 1000*1000
  sec_per_time_unit <- ifelse(time_unit == "secs", 1, ifelse(time_unit == "mins", 60, 60*60))
  
  dat <- (dat %>% 
            arrange(tagID, first_detection) %>% 
            group_by(tagID) %>% 
            mutate(migration_rate = trans_len/as.numeric(time_btw),
                   mig_rate_unit = paste0("km/", units(time_btw))) %>% # if trans_len comes with units, use that; assume mm
            mutate(BLPS = (as.numeric(migration_rate)/as.numeric(fork_len))*(mm_per_dist_unit)/(sec_per_time_unit)) %>% 
            ungroup()) # %>% 
  
  return(dat)  
}


### Calculate against_flow
get_against_flow <- function(dat)
{
  ## Description: 
  # Calculate logical indicator of transition directed against the direction of water flow
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: against_flow
  

  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "first_detection", "gen_loc", "trans_unidir_flow", "trans_type" 

  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    gen_loc  = character; name of telemetry station (general location)
  #    trans_undir_flow = logical; indicates whether transition included only unidirectional (downstream-directed) flow (TRUE)
  #                       or if it included a region with bi-directional flow (FALSE); = NA if unknown
  #    trans_type = categorical; "downstream", "upstream", or NA (NA includes lateral transition or repeated station transition)


  # Assumes that separate detection events do not overlap in time
  
  # against_flow metric = TRUE if the transition was directed upstream in a region and time period when it is known that the river
  #                       flow was directed downstream;
  #                     = FALSE if the transition was directed downstream in a region and time period when it is known that the 
  #                       river flow was directed downstream;
  #                     = NA otherwise: 
  #                         - transition was neither downstream nor upstream (e.g., lateral or none/repeated station);
  #                         - transition was in fluvial-tidal transition zone downstream of extent of unidirectional flow 
  #                            (so the direction of river flow was not known)
  #                         - transition was in tidal zone (so the direction of river flow was not known)
  
  # Authors: Dalton Hance, USGS, WFRC; Rebecca Buchanan, UW, SAFS, CBR; Steve Whitlock, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  

  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID",
            "first_detection",
            "gen_loc",
            'trans_unidir_flow',
            'trans_type') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID', 'first_detection', 'gen_loc', 'trans_unidir_flow', 'trans_type'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection)){
    stop("Error: first_detection must be of class POSIXct or POSIXt ")
  }
  
  
  ### Compute metric
  
  dat %>% 
    arrange(tagID, first_detection) %>%
    group_by(tagID) %>%
    mutate(against_flow = case_when((trans_unidir_flow==TRUE & trans_type=="downstream") ~ FALSE,
                                    (trans_unidir_flow==TRUE & trans_type=="upstream") ~ TRUE)) %>%
    ungroup() %>%
    return()

  }


### Calculate cumulative count of down-up switches on tag level
get_cumul_down_up_switch_tag <- function(x)
{
  ## Description: 
  # Calculate cumulative count of switches from moving downstream to moving upstream for a single tag
  
  ## Necessary packages:
  #
  
  ## Dependencies:
  # 
  
  ## Arguments:
  # x = categorical vector defining transition type for each detection event: "downstream", "upstream", or NA;
  #     NA category includes every transition type that is not strictly "downstream" or "upstream", including
  #     lateral transitions and repeated detections at the same station
  
  ## Value:
  # character vector of length nrow(x) equal to the cumulative count of switches from moving downstream to moving upstream
  
  ## Details:
  # Internal function used in get_cumul_down_up_switch()
  
  
  # Author: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External: 
  # Internal (get_cumul_down_up_switch)
  
  
  ##########################################################
  
  ### Error checking
  if(any(!(x %in% c("downstream","upstream",NA))))
  {
    stop("valid entries of x must be limited to 'downstream', 'upstream', and NA")
  }
  
  
  ### Compute metric
  
  if(all(is.na(x))) return(as.numeric(x))
  
  x1 <- as.numeric(as.factor(x))
  DU.switch <- (diff(x1[!is.na(x1)])==1)
  DU.switch.full <- rep(NA,length(x))
  #na.index <- unique(c(1,which(is.na(x))))   
  na.index <- sort(c(min(which(!is.na(x))),which(is.na(x))))
  DU.switch.full[-na.index] <- DU.switch       
  cumul_down_up_switch <- rep(NA,length(x))
  cumul_down_up_switch[-na.index] <- cumsum(as.numeric(DU.switch))
  
  return(cumul_down_up_switch)
}

 
### Calculate cumul_down_up_switch (all tags)
get_cumul_down_up_switch <- function(dat)
{
  ## Description: 
  # Calculate cumulative count of switches from moving downstream to moving upstream for all tags
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # get_cumul_down_up_switch_tag()
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: cumul_down_up_switch
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "gen_loc", "first_detection", "trans_type"
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    gen_loc  = character; name of telemetry station (general location)
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    trans_type = categorical; "downstream", "upstream", or NA (NA includes lateral transition or repeated station transition)
  
  # Assumes that separate detection events do not overlap in time
  
  # wrapper function for get_cumul_down_up_switch_tag()
  
  # output metric cumul_down_up_switch is the cumulative count of switches from downstream movement to upstream movement;
  # ignores intervening transitions with trans_type defined as NA (such as lateral transitions or repeated detections at the same station);
  # detection history segments that move from trans_type="downstream" to trans_type=NA (for 1 or more rows) to 
  #           trans_type="upstream" result in increasing cumul_down_up_switch by 1

  
  # Author: Dalton Hance, USGS, WFRC
  # Edited: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  

  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID", 
            "gen_loc",
            'first_detection', 'trans_type') %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID',
          'first_detection', 'trans_type',
          'gen_loc'")
  }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection)){
    stop("Error: first_detection must be of class POSIXct or POSIXt ")
  }
  
  # Note: once called "get_cumul_down_up_switch_MOD" but replacing old version with 'group()' rather than 'group_by'
  
  ### Compute metric
  dat %>% 
    arrange(tagID, first_detection) %>%
    group_by(tagID) %>% # replaced group() with group_by
    mutate(cumul_down_up_switch = (trans_type %>% get_cumul_down_up_switch_tag())) %>%
    ungroup() %>%
    return()

}

### Calculate cumulative count of against-flow transitions, defined for every transition
get_cumul_against_flow <- function(dat)
{
  ## Description: 
  # Calculate cumulative count of against-flow transitions, defined for every transition
  
  ## Necessary packages:
  require(dplyr)
  
  ## Dependencies:
  # against_flow or trans_unidir_flow, trans_type
  
  ## Arguments:
  # dat = data frame or tibble of detection events on general location spatial scale; see Details
  
  ## Value:
  # dat data frame with additional field: cumul_against_flow
  
  ## Details:
  # Required fields in dat:
  #    "tagID", "studyID", "first_detection", and either "against_flow" or all of "gen_loc", "trans_unidir_flow", "trans_type" 
  
  #    tagID = alphanumeric code representing tag ID
  #    studyID = character; name of study ID
  #    first_detection = POSIXct or POSIXlt; date and time of first detection of detection event
  #    against_flow = logical indicator of transition directed against the direction of water flow
  #    gen_loc  = character; name of telemetry station (general location)
  #    trans_undir_flow = logical; indicates whether transition included only unidirectional (downstream-directed) flow (TRUE)
  #                       or if it included a region with bi-directional flow (FALSE); = NA if unknown
  #    trans_type = categorical; "downstream", "upstream", or NA (NA includes lateral transition or repeated station transition)
  
  
  # Assumes that separate detection events do not overlap in time
  
  # output metric cumul_against_flow is the cumulative (running) count of transitions known to be directed against the direction of river flow
  # single counting sequence per tag
  # defined for every station

  # Author: Rebecca Buchanan, UW, SAFS, CBR
  
  
  ## Internal or External:
  # External
  
  
  
  ##########################################################
  
  ### Error checking
  
  if(!all(c("tagID", 
            "studyID",
            "first_detection") %in% names(dat))){
    stop("Error: data must minimally contain variables 'tagID',
          'studyID', 'first_detection'")
  }
  
  if(!all(c("gen_loc",
            'trans_unidir_flow',
            'trans_type') %in% names(dat)) &
     !all(c("against_flow") %in% names(dat))
    ){
      stop("Error: data must contain either variable 'against_flow' or variables 'gen_loc',
           'trans_unidir_flow', 'trans_type'")
    }
  
  if(!is.POSIXt(dat$first_detection) | !is.POSIXct(dat$first_detection)){
    stop("Error: first_detection must be of class POSIXct or POSIXt ")
  }
  
  
  ### Compute metric
  
  if("against_flow" %in% names(dat)){
    dat<- (dat %>% 
             arrange(tagID, first_detection) %>% 
             group_by(tagID) %>% 
             mutate(cumul_against_flow = cumsum(coalesce(against_flow,0))) %>% 
             ungroup()) # %>% 
    return(dat)
  }
  
  dat <- (dat %>% 
    arrange(tagID, first_detection) %>% 
    get_against_flow() %>% 
    group_by(tagID) %>% 
    mutate(cumul_against_flow = cumsum(coalesce(against_flow,0))) %>% 
    ungroup()) # %>% 
  return(dat)
  
}


