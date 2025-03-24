##############################################################################
##############################################################################
# SECTION 3. Prepare brood and maturity tables
# Last updated: December 31, 2024
# Written by: Sabrina Garcia
##############################################################################
##############################################################################

### NOTE: Some of the code chunks could probably be condensed into a function ###

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
require(here)

wd <- here()  
setwd(wd)

dir.data <- file.path(wd, "data")

# If running scripts simultaneously, objects should already be in global 
# environment 

# Next line will check if object exists and will read in previous script if
# it doesn't. Note script 2 will source script 1
if( !exists( "all_juvs" ) ) source( "2_calc_stock_spec_juv_abun.R" )

######################### Canadian-origin ######################################

# Read in brood tables from master file 
cdn_brood <- read_excel( file.path(dir.data, "masterNBSdata.xlsx"), sheet = "CDNbrood - Hama" )

# Create juvenile per spawner data table for completed brood years - CANADIAN
# Canadian return data is assumed to have a CV of 10%
# Survival variance calculated using taylor series approximation
cdn_all_dat <- cdn_brood %>%
     select( brood_year, juv_year, returns, spawners ) %>% 
     left_join( select( all_juvs, year, cdn_abun, cdn_joint_var ), 
                 by = c( "juv_year" = "year" ) ) %>%
     mutate( return_cv = 0.10, 
             return_sd = return_cv * returns, 
             return_var = return_sd^2,
             juv_per_spawn = cdn_abun/spawners, 
             survival = returns/cdn_abun,
             surv_var = ( returns^2 )/ ( cdn_abun^4 ) * cdn_joint_var + ( 1/( cdn_abun^2 ) ) *
                  ( return_var ),
             surv_sd = sqrt( surv_var ),
             surv_cv = surv_sd/survival  )

# Only includes brood years 2001 to present - CANADIAN      
cdn_mat <- cdn_brood %>%
     select( brood_year, three, four, five, six, seven, eight, returns) %>%
     mutate( age3 = round( three/returns, 3 ), age4 = round( four/returns, 3), 
             age5 = round( five/returns, 3 ), age6 = round( six/returns, 3),
             age7 = round( seven/returns, 3 ), age8 = round( eight/returns, 3 ) ) %>% 
     select( brood_year, age3, age4, age5, age6, age7, age8 ) %>%
     filter( brood_year > 1997 ) 

# Create average maturity tables that will be used to forecast - CANADIAN
# Average age columns are lagged so that the average ages align correctly
# E.g., the average ages aligned with brood year 1985 are the averages from 1982 - 1984
cdn_mat_forecast <- cdn_mat %>%
     mutate( age3_avg = lag( rollmean( age3, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE,
                                       fill = NA ) ),
             age4_avg = lag( rollmean( age4, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ),
             age5_avg = lag( rollmean( age5, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ),
             age6_avg = lag( rollmean( age6, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ) ) %>%
     as.data.frame() %>% 
     select( brood_year, age3_avg, age4_avg, age5_avg, age6_avg ) %>% 
     filter( brood_year > 2000 )        

############################# Total Yukon ######################################

# Read in brood tables from master file 
total_brood <- read_excel( file.path(dir.data,"masterNBSdata.xlsx"), sheet = "TOTbrood - Hama" )

# Create juvenile per spawner data table for completed brood years - TOTAL RUN
# Total Yukon return data is assumed to have a CV of 15%
total_all_dat <- total_brood %>%
     select( brood_year, juv_year, returns, spawners ) %>%
     left_join( select( all_juvs, year, totyuk_abun, totyuk_joint_var ), 
                 by = c( "juv_year" = "year" ) ) %>%
     mutate( return_cv = 0.15, 
             return_sd = return_cv * returns, 
             return_var = return_sd^2,
             juv_per_spawn = totyuk_abun/spawners, 
             survival = returns/totyuk_abun,
             surv_var = ( returns^2 )/ ( totyuk_abun^4 ) * totyuk_joint_var + 
                  ( 1/( totyuk_abun^2 ) ) * ( return_var ),
             surv_sd = sqrt( surv_var ), 
             surv_cv = surv_sd/survival  )

# Only includes brood years that are complete - TOTAL
total_mat <- total_brood %>%
     select( brood_year, three, four, five, six, seven, returns) %>%
     mutate( age3 = three/returns, age4 = four/returns, age5 = five/returns, 
             age6 = six/returns,
             age7 = seven/returns )  %>% 
     select( brood_year, age3, age4, age5, age6, age7 ) %>% 
     filter( brood_year > 1997 ) 

# Create average maturity tables that will be used to forecast - TOTAL
total_mat_forecast <- total_mat %>%
     mutate( age3_avg = lag( rollmean( age3, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ),
             age4_avg = lag( rollmean( age4, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ),
             age5_avg = lag( rollmean( age5, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ),
             age6_avg = lag( rollmean( age6, k = 3, align = "right",
                                       by.column = TRUE, partial = FALSE, 
                                       fill = NA ) ) ) %>%
     as.data.frame() %>% 
     select( brood_year, age3_avg, age4_avg, age5_avg, age6_avg ) %>% 
     filter( brood_year > 2000 )

