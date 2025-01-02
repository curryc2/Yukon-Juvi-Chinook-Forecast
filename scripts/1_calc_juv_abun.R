##############################################################################
##############################################################################
# SECTION 1. Juvenile Chinook salmon abundance
# Last updated: December 31, 2024
# Written by: Sabrina Garcia
##############################################################################
##############################################################################

# Load required packages
require(geosphere)
require(reshape2)
require(tidyverse)
require(readxl)
require(zoo)
require(scales)
require(ggpubr)
require(ggrepel)
require(ggpmisc)
require(tibbletime)

###############################################################################
# IMPORTANT: Annual changes needed to R script
################################################################################

#1. If station-specific MLDs are available for the most recent year, set 
# recent_year_mlds to TRUE
recent_year_mlds <- TRUE

#2. If any strata were missed in the most recent year of sampling, will need
# to manipulate final sections of code that apply a correction to the juvenile
# abundance estimates to account for missed sampling. EVentually, we could 
# write a function that can handle that part automatically

#3. Set working directory to desired folder
setwd( "C:/Users/sgarcia/Desktop/Marine Research/NBS Data/Chinook/Juvenile Forecast" )  

################################################################################
# Read in northern Bering Sea survey data
################################################################################
nbs <- read_excel( "masterNBSdata.xlsx", sheet="juvs" )

# columns are: year, region, station, date (in AKST), start and end lat/longs,
# distance (km) as calculated using the starting and ended GPS coordinates, 
# horizontal (aka wingspread, meters),horizontal_adj (wingspread adjusted for 
# incorrect measurements, m), vertical_netsounder (vertical opening as measured 
# by the netsounder if available, meters), meanfootrope_sbeR (footrope depth
# as measured by CTD sensors if available, meters), meanheadrope (depth of the 
# headrope as measured by CTD sensors if available, meters),
# headrope_qualifier (method of estimate headrope depth; either estimated 
# based on previous surveys, observed from CTD sensor, or averaged across 
# available data), areaswept (horizontal spread adjusted times distance, km2), 
# kingnum (number of juvenile kings), cpue (kingnum divided by areaswept, #/km2), 
# and mld01 (mixed layer depth based on the definition of MLD using the depth
# where density increases by 0.10 kg/m3 relative to the density at 5 m)
# NOTE: MLD needs to be calculated prior to using this code

################################################################################
# Create new column for stratum based on GPS coordinates
################################################################################

# 1 <- 60-62 N, 2 <- 62-64, 3 <- Norton Sound, 4 <- Bering Strait
nbs <- nbs %>%
     mutate( stratum = case_when ( startlat < 62.1 ~ 1, 
                                   startlong > -165.5 ~ 3,
                                   startlat > 64.1 ~ 4, 
                                   TRUE ~ 2 ) ) 

# Update 6/1/2023: One station in stratum 3 is at 67 N because it was sampled
# very far east, north of the Bering Strait (i.e., it is not in Norton Sound), 
# so it needs to be reassigned to stratum 4
nbs <- nbs %>% 
     mutate( stratum = case_when( year == 2018 & station == 49 ~ 4 , 
                                  TRUE ~ stratum ) ) 

################################################################################
# Create new column for footrope depth based on data available 
# (netsounder vs net CTD sensor)
# Ideally, footrope depth is estimated directly from the net sensor (sbe) but 
# when those data are unavailable (e.g. 2021), footrope depth is equal to 
# vertical opening + headrope depth
# Footrope depth is used to determine the MLD expansion used below
#################################################################################

nbs <- nbs %>%
     mutate( finalfootrope = case_when( is.na( meanfootrope_sbeR ) ~ 
                                             ( vertical_netsounder + meanheadrope ),
                                        TRUE ~ meanfootrope_sbeR ) )

################################################################################
# Calculate MLD expansion for each station using MLD and footrope depths
# If the trawl fishes deeper than the MLD depth, the adjustment is 1 
# (i.e., no expansion)
# If trawl is shallower than MLD depth, 
# the expansion is the ratio of the MLD depth:footrope depth
# If MLD depths are not yet available, then no expansion occurs (i.e., equal to 1)
# Calculate MLD-expanded catch for each station
################################################################################

nbs <- nbs %>%
     mutate( MLDexp = case_when( mld01 <= finalfootrope ~ 1,
                                 mld01 > finalfootrope ~ ( mld01/finalfootrope ),
                                 is.na(mld01) ~ 1 ), 
             kingexp = kingnum * MLDexp )

# Calculate overall MLD expansion by year
# The overall MLD gives an indication of overall fishing performance
# If station-specific MLDs are not available, the overall MLD is used to adjust
# the juvenile abundance later in the code
catch_table <- nbs %>%
     group_by( year ) %>%
     summarise( juvs = sum( kingnum ), expjuvs = round( sum( kingexp ), 0 ) ) %>%
     mutate( overallMLD = round( expjuvs/juvs, 2 ) ) %>%
     as.data.frame()

# Need to set recent_year_mlds to TRUE (MLDs available) or FALSE (MLDs unavailable)
# If true, next code chunk won't execute; if false, will calculate average MLD
# from the last three years
if( recent_year_mlds == FALSE ) {
     
recent_year <- catch_table[ nrow( catch_table ), "year" ]

recent_mld <- filter( nbs, year < recent_year - 1 & year > recent_year - 4 ) %>% 
     select( MLDexp ) %>% 
     summarise( mean_MLD = mean( MLDexp ) ) %>% 
     as.numeric() }

################################################################################
# Estimate CPUE by year and stratum
# Sum the MLD-expanded king catches and divide by the summed area swept 
# within each stratum 
################################################################################

cpue_table <- nbs %>%
     group_by( year, stratum ) %>%
     summarize( tot_stns = length( station ), stratum_catch = sum( kingexp ), 
                stratum_effort = sum( areaswept )  ) %>%
     mutate( stratum_cpue = stratum_catch/stratum_effort ) %>%
     as.data.frame()

################################################################################
# Estimate variance of stratum-specific CPUE
# Uses equations 1.46 and 1.47 from Quinn and Deriso 1999
################################################################################

# Select columns necessary for variance estimator
# We need station-specific data which is why we cannot use cpue.table which is a 
# stratum-level summary
var_data <- nbs %>% 
     select( year, station, stratum, kingexp, areaswept )

# Join cpue.table (stratum-level) to varNBS (station-level) to calculate variance 
# across every row of data

# Stratum 3 (NS) only has one station sampled in 2005 and 2007. 
# Because there was only one station sampled, the scalar computation results in
# Inf which then creates a NA in the variance calculation.
# According to Jim Jasper (ADF&G biometrician), we are in the realm of negative
# binomial distributed catches. Therefore, Norton Sound variance in 2005 and 2007
# is equal to (catch plus catch^2)/areaswept^2. 
# This is handled within the case_when call nested in the mutate function
var_calc <- var_data %>% 
     right_join( cpue_table, by = c( "year","stratum" ) ) %>% 
     mutate( var_numerator = ( kingexp - ( stratum_cpue * areaswept ) )^2,
             var_denominator = stratum_effort^2, scalar = tot_stns/( tot_stns-1 ) ) %>%
     group_by( year, stratum ) %>%
     mutate( var_numerator = sum( var_numerator ), 
             strat_CPUE_var = scalar * ( var_numerator/var_denominator ),
             strat_CPUE_var = case_when( is.nan( strat_CPUE_var ) ~ 
                                                 ( kingexp^2 )/( areaswept^2 ),
                                            TRUE ~ strat_CPUE_var ),
             strat_CPUE_sd = sqrt( strat_CPUE_var ) ) %>% 
     as.data.frame() 

cpue_table <- cpue_table %>%
     right_join( distinct( select( var_calc, year, stratum, strat_CPUE_var, 
                                   strat_CPUE_sd ) ), 
                 by = c( "year", "stratum" ) )

################################################################################
# Calculate survey area for each year
################################################################################

# Get average latitude for each stratum needed for longitude calculation below
avg_lat <- nbs %>% 
     group_by( year, stratum ) %>%
     summarise( lat = mean( startlat, na.rm = TRUE ) ) 

cpue_table <- cpue_table %>%
     left_join( avg_lat, by = c( "year","stratum" ) )

# Use two arbitrary longitudes, one degree apart, to get the distance of one 
# degree longitude in each strata
lon1 <- -170
lon2 <- -171

# Need to get the distance of one degree longitude for each stratum for each year by using
# the Haversine distance formula (output is in meters which is then converted to km 
# by dividing by 1000)
# One-half degree latitude equals 30 nautical miles which is converted to kilometers
# Grid area = the longitude distance (km) times the latitude distance (km)
# The final survey area =  grid area of a stratum multiplied by the number of stations sampled
# in that stratum
# Norton Sound (stratum=3) has a fixed survey area of 5,461 km^2
cpue_table <- cpue_table %>%
     mutate( lon_dist = ( distHaversine( p1 = cbind(lon1, lat), p2 = cbind(lon2, lat) )/1000 ),
             lat_grid = 30 * 1.852,
             grid_area = lon_dist * lat_grid, 
             surv_area = grid_area * tot_stns,
             surv_area = case_when( stratum == 3 ~ 5461, TRUE ~ surv_area ) ) 

################################################################################
# Calculate total survey area CPUE by weighting by area
################################################################################

nbs_area <- cpue_table %>%
     group_by ( year ) %>%
     summarise( total_area = sum( surv_area ) )

cpue_table <- cpue_table %>%  
     right_join( nbs_area, by = "year" ) %>% 
     mutate( area_prop = surv_area/total_area, weighted_cpue = area_prop * stratum_cpue )

# nbs_cpue has annual CPUE weighted by strata area
# total_nbs_cpue uses Equation 1.51
# total_nbs_var (variance for yearly weighted CPUE) uses eq. 1.52 
nbs_cpue <- cpue_table %>%
     group_by( year ) %>%
     summarise( total_nbs_cpue = sum( weighted_cpue ), 
                total_nbs_var= sum( area_prop^2 * strat_CPUE_var) ) %>%
     mutate( total_nbs_sd = sqrt( total_nbs_var ) )

################################################################################
# Calculate mean annual juvenile abundance
################################################################################

total_juvs <- nbs_cpue %>%  
     right_join( nbs_area, by = c( "year" ) ) %>% 
     mutate( juv_abun = total_area * total_nbs_cpue, 
             juv_var = total_area^2 * total_nbs_var,
             juv_sd = sqrt( juv_var ), 
             juv_CV = round( ( ( juv_sd/juv_abun )*100 ),0 ) )

################################################################################
# Account for missed sampling of Norton Sound in 2016
# If another year of Norton Sound is missed, this section will need to be updated
################################################################################

# Pull out year that has missed sampling of Norton Sound
missed_ns <- cpue_table %>% 
     filter( stratum == 3 ) %>% 
     pull( year )

sequence <- unique(nbs$year)

missing_ns_year <- sequence[!sequence %in% missed_ns]

# Calculate proportion of juveniles that are caught in Norton Sound relative to
# the whole area
ns_juvs <- cpue_table %>%
     filter( stratum == 3 ) %>% 
     select( year, stratum_cpue, surv_area ) %>%
     mutate( ns_juv_abun = stratum_cpue * surv_area )  %>%
     left_join( select( total_juvs, year, juv_abun ), by = c( "year" ) )  %>%
     mutate( ns_prop = ns_juv_abun/juv_abun )

# Determine the average Norton Sound proportion of juveniles in the two years before and
# two years after missed sampling in Norton Sound
avg_ns_prop <- ns_juvs %>%
     filter( year %in% c( (missing_ns_year-2):(missing_ns_year+2) ) )   %>%
     summarise( mean( ns_prop ) ) %>%
     unlist( use.names = F )

# Multiply 2016 juv abundance by average NS proportion to account for missed sampling
total_juvs <- total_juvs %>%
     mutate( juv_abun_adj = case_when( year == missing_ns_year ~ 
                                            juv_abun + ( juv_abun * avg_ns_prop ), 
                                     TRUE ~ juv_abun ) )

################################################################################
# Account for missed sampling of Bering Strait in 2004-2006
################################################################################

# Pull out years that have missed sampling of Bering Strait
missed_bs <- cpue_table %>% 
     filter( stratum == 4 ) %>% 
     pull( year )

missing_bs_year <- sequence[!sequence %in% missed_bs]

# Calculate Bering Strait juvenile abundance to get proportion
bs_juvs <- cpue_table %>%
     filter( stratum == 4 ) %>% 
     select( year, stratum_cpue, surv_area ) %>%
     mutate( bs_juvs = stratum_cpue * surv_area )  %>%
     left_join( select( total_juvs, year, juv_abun ), by=c( "year" ) ) %>%
     mutate( bs_prop= bs_juvs/juv_abun )

# Determine the average Norton Sound proportion of juveniles in the year before and
# after missed sampling in Bering Strait
avg_bs_prop <- bs_juvs %>%
     filter( year %in% c( (min(missing_bs_year)-1):(max(missing_bs_year)+1) ) )   %>%
     summarise( mean( bs_prop ) ) %>%
     unlist( use.names = F )

# Multiply 2004-2006 juv abundance by average BS proportion to account for missed sampling
total_juvs <- total_juvs %>%
     mutate( juv_abun_adj = case_when( year %in% c( min(missing_bs_year):max(missing_bs_year) ) ~ 
                                          juv_abun + ( juv_abun * avg_bs_prop ),
                                     TRUE ~ juv_abun_adj ) ) %>% 
     data.frame()

################################################################################
# Adjust recent year juvenile abundance if station-specific MLDs are unavailable
################################################################################

if( recent_year_mlds == FALSE ) {
     
     total_juvs <- total_juvs %>% 
     mutate( juv_abun_adj = case_when( year == recent_year ~ juv_abun * recent_mld, 
     TRUE ~ juv_abun_adj ) )     }

################################################################################
# Write csv file with the final juvenile abundance estimates
# Saves to output folder within directory
################################################################################

dir.create( paste0( getwd(),"/Output") ) # will output warning if this folder already exists
    
write.csv( total_juvs, file = paste0("Output/juvChinookAbun_",Sys.Date(),".csv") )













