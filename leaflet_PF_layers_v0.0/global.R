library(dplyr)

allzips <- readRDS("data/superzip.rds")
allzips$latitude <- jitter(allzips$latitude)
allzips$longitude <- jitter(allzips$longitude)
allzips$college <- allzips$college * 100
allzips$zipcode <- formatC(allzips$zipcode, width=5, format="d", flag="0")
row.names(allzips) <- allzips$zipcode

# cleantable <- allzips %>%
#   select(
#     City = city.x,
#     State = state.x,
#     Zipcode = zipcode,
#     Rank = rank,
#     Score = centile,
#     Superzip = superzip,
#     Population = adultpop,
#     College = college,
#     Income = income,
#     Lat = latitude,
#     Long = longitude
#   )

lmap <- readRDS("data/lmap.rds")
my_flownet <- readRDS("data/my_flownet.rds")

# my_flownet$sp

names(my_flownet$sp)

cleantable <- my_flownet$sp@data %>%
  select(
    ID,
    FLW_TYP,
    Basin = basin,
    Region = region
    # Zipcode = zipcode,
    # Rank = rank,
    # Score = centile,
    # Superzip = superzip,
    # Population = adultpop,
    # College = college,
    # Income = income,
    # Lat = latitude,
    # Long = longitude
  )

flow_leaf <- readRDS("data/flow_leaf")
reg_cols_hex <- readRDS("data/reg_cols_hex")
reg_pal <- readRDS("data/reg_pal")
bas_cols_hex <- readRDS("data/bas_cols_hex")
bas_pals <- readRDS("data/bas_pal")
TFTZ_lab_use <- readRDS("data/TFTZ_lab_use")
tftz_cols_hex <- readRDS("data/tftz_cols_hex")
tftz_pal <- readRDS("data/tftz_pal")
bas_pal <- readRDS("data/bas_pal")

flow_leaf$HRZ_fact <- factor(ifelse(flow_leaf$HRZ_ID>0,"High-risk zone",NA))
HRZ_hex_cols <- "#D13434"
hrz_pal <- leaflet::colorFactor(palette = HRZ_hex_cols,
                                domain = flow_leaf$HRZ_fact,na.color = NA)


# lmap_HRZ <- leaflet() %>%  
#   # addProviderTiles(providers$CartoDB.Positron) %>%
#   addPolylines(data=flow_leaf, color = "#000000", label=flow_leaf$ID,group = "basin",opacity = 1,weight = 1) %>%
#   addPolylines(data=flow_leaf, color = ~hrz_pal(HRZ_fact),label=flow_leaf$ID,group = "High Risk Zone",opacity = 1,weight=4)

# load("data/leaflet_env.Rdata")

