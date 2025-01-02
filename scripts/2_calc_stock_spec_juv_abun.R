##############################################################################
##############################################################################
# SECTION 2. Stock-specific juvenile abundance estimates
# Last updated: December 31, 2024
# Written by: Sabrina Garcia
##############################################################################
##############################################################################

# If running scripts simultaneously, objects should already be in global 
# environment and you can run entire script

# Next line will check if object exists and will read in previous script if
# it doesn't
if( !exists( "total_juvs" ) ) source( "1_calc_juv_abun.R" )

###############################################################################
# Read in genetics data in order to calculate stock-specific juvenile abundance
# Calculate averages to fill in years where no genetic data is available
###############################################################################

gsi <- read_excel( "masterNBSdata.xlsx", sheet="GSI" )  %>%
       mutate( cdn_mean = cdn_mean/100, cdn_sd = cdn_sd/100, middle_mean = middle_mean/100,
             lower_mean = lower_mean/100, otherwak_mean = otherwak_mean/100,
             totalYukon_mean = totalYukon_mean/100, totalYukon_sd = totalYukon_sd/100 ) 

# Calculate mean genetic stock composition for years with data
# Missing genetic stock compositions will be equal to the average stock compositions
# across all years between 2003 and 2016
# Assumed 0.10 sd for Canadian-origin and Total Yukon stocks in years 2005, 2012, 2013
mean_gsi <- gsi %>%
     filter( !is.na( n ), year < 2017 ) %>%
     summarise( mean( cdn_mean, na.rm = TRUE ), 
                mean( totalYukon_mean, na.rm = TRUE ) ) %>%
     as.numeric()

gsi <- gsi %>%
     select( year, cdn_mean, cdn_sd, totalYukon_mean, totalYukon_sd ) %>% 
     mutate( cdn_mean = case_when( is.na( cdn_mean ) ~ mean_gsi[1], 
                              TRUE ~ cdn_mean ), 
             cdn_sd = case_when( is.na( cdn_sd ) ~ 0.10, 
                              TRUE ~ cdn_sd ), 
             totalYukon_mean = case_when( is.na( totalYukon_mean ) ~ mean_gsi[2], 
                              TRUE ~ totalYukon_mean ),
             totalYukon_sd = case_when( is.na( totalYukon_sd ) ~ 0.10, 
                              TRUE ~ totalYukon_sd ) )

###############################################################################
# Join genetic data to juvenile data
###############################################################################

all_juvs <- select( total_juvs, year, juv_abun_adj, juv_sd, juv_var ) %>%
     right_join( select( gsi, year, cdn_mean, cdn_sd, totalYukon_mean, totalYukon_sd ),
                 by = c( "year" ) )

# cdn.juvs has the Canadian-origin juvenile abundance with associated variance
# Remove 2005 because that year had odd survey timing and will be replaced in
# next chunk of code
all_juvs <- all_juvs %>%
     filter( !year == 2005 & !year == 2020 & !year == 2008 ) %>%         
     mutate( cdn_abun = juv_abun_adj * cdn_mean,
             cdn_joint_var = ( cdn_mean^2 * juv_sd^2 ) + ( juv_abun_adj^2 * cdn_sd^2 ) + 
                  ( 2 * juv_abun_adj * cdn_mean * cor( juv_abun_adj, cdn_mean ) * juv_sd * cdn_sd ) ,
             cdn_joint_sd = sqrt( cdn_joint_var ), cdn_joint_cv = cdn_joint_sd/cdn_abun ) %>% 
     mutate( totyuk_abun = juv_abun_adj * totalYukon_mean,
             totyuk_joint_var = ( totalYukon_mean^2 * juv_sd^2 ) + 
               ( juv_abun_adj^2 * totalYukon_sd^2 ) + 
          ( 2 * juv_abun_adj * totalYukon_mean * cor( juv_abun_adj, totalYukon_mean ) * juv_sd * totalYukon_sd ),
             totyuk_joint_sd = sqrt( totyuk_joint_var ), totyuk_joint_cv = totyuk_joint_sd/totyuk_abun )
  
################################################################################
# Hard enter 2005, 2008, and 2020 stock-specific juvenile abundance calculated 
# using Kalman smoothing by Jim Murphy (NOAA)
# These values are sent in a separate Excel file but ideally, this would happen 
# within the code to have it all in one place
################################################################################

row_2005 <- data.frame( year = 2005, cdn_abun = 2311816, totyuk_abun = 4289692 )

row_2008 <- data.frame( year = 2008, cdn_abun = 1640563, totyuk_abun = 2680333 )

row_2020 <- data.frame( year = 2020, cdn_abun = 910023, totyuk_abun = 1175648 )

all_juvs <- all_juvs %>% 
     rows_insert( rbind( row_2005, row_2008, row_2020  ) ) %>% 
     arrange( year )

################################################################################
# Write csv file with the final stock-specific juvenile abundance estimates
################################################################################

write.csv( all_juvs, paste0("Output/juvYukonChinookAbun_",Sys.Date(),".csv") )

