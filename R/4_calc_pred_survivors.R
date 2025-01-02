##############################################################################
##############################################################################
# SECTION 4. Calculate predicted survivors from juvenile abundance
# Last updated: December 20, 2024
# Written by: Sabrina Garcia
##############################################################################
##############################################################################

# If you don't run scripts 1 - 3, you will need to source the objects from them
# Next line will search for objects and source scripts if needed
if( !exists( "total_mat_forecast" ) ) source( "3_prep_brood_tables.R" )

# Function to fit linear regression and calculates predicted survivors (i.e.,
# fish that will survive marine life and return to the Yukon)
regress_func <- function( dat ) {
     
     dat <- select( dat, brood_year, juv_year, returns, matches( "abun" ) ) %>% 
          rename( juv_abun = matches("abun" ) )
     
     regress <- lm( returns ~ juv_abun, data = dat )
     
     pred_survs <- predict( regress, newdata = dat["juv_abun"],
                               interval = 'prediction', level = 0.8 ) 
     
     colnames( pred_survs ) <- c( "pred.surv.fit", "pred.surv.lwr", "pred.surv.upr" )
     
     results <- list( cbind( dat, pred_survs ), summary(regress) )

     return(results) }

# Canadian linear model fit (element 2 is summary of the lm)
cdn_model_fit <- regress_func( dat = cdn_all_dat )[[1]]

# Need to create projection data that include the linear model predictions
# along with the 3-year average maturity calculated in script 3
# mutate fills in recent brood years with the latest brood year averages
cdn_proj_data <- full_join( cdn_model_fit, cdn_mat_forecast ) %>% 
        drop_na( juv_year ) %>% 
        mutate( age3_avg = case_when( is.na( age3_avg ) ~ na.locf( age3_avg, na.rm = FALSE ),
                                      TRUE ~ age3_avg ),
                age4_avg = case_when( is.na( age4_avg ) ~ na.locf( age4_avg, na.rm = FALSE ),
                                      TRUE ~ age4_avg ),
                age5_avg = case_when( is.na( age5_avg ) ~ na.locf( age5_avg, na.rm = FALSE ),
                                      TRUE ~ age5_avg ),
                age6_avg = case_when( is.na( age6_avg ) ~ na.locf( age6_avg, na.rm = FALSE ),
                                      TRUE ~ age6_avg ) )

write.csv( cdn_proj_data, paste0("Output/cdn_proj_data_",Sys.Date(),".csv") )

# Total Yukon linear model fit
total_model_fit <- regress_func( dat = total_all_dat )[[1]]

total_proj_data <- right_join( total_model_fit, total_mat_forecast ) %>% 
        drop_na( juv_year ) %>% 
        mutate( age3_avg = case_when( is.na( age3_avg ) ~ na.locf( age3_avg, na.rm = FALSE ),
                                      TRUE ~ age3_avg ),
                age4_avg = case_when( is.na( age4_avg ) ~ na.locf( age4_avg, na.rm = FALSE ),
                                      TRUE ~ age4_avg ),
                age5_avg = case_when( is.na( age5_avg ) ~ na.locf( age5_avg, na.rm = FALSE ),
                                      TRUE ~ age5_avg ),
                age6_avg = case_when( is.na( age6_avg ) ~ na.locf( age6_avg, na.rm = FALSE ),
                                      TRUE ~ age6_avg ) )

write.csv( total_proj_data, paste0("Output/total_proj_data_",Sys.Date(),".csv") )
