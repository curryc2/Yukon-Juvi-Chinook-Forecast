##############################################################################
##############################################################################
# SECTION 5. Forecast future runs using predicted survivors and average age-
# -at-maturity
# Last updated: Janaury 2, 2025
# Written by: Sabrina Garcia
##############################################################################
##############################################################################

# If you haven't run scripts 1 - 4, you will need to source the objects from them
# Next line will search for objects and source scripts if needed
if( !exists( "total_proj_data" ) ) source( "4_calc_pred_survivors.R" )

# Function to calculate projected survivors based on prediction type
# (fit, lower 80%, upper 80%)
# This function will calculate age-specific survivors using all juvenile years
# Age compositions do not add up to 100% because we don't forecast age-7
# With juvenile years 2003 - 2023, we will have complete forecasts for 2007-2026
projSurvs <- function( dat ) {
        
        pred_list <- list() # empty list to hold results
        
        pred_cols <- select( dat, starts_with( "pred" ) ) %>% 
                data.frame() # pulls out three prediction levels
        
        for( i in 1:ncol(pred_cols) ) {
                
                pred_list[[i]] <- dat %>%
                        dplyr::select( juv_year, age3_avg, age4_avg, age5_avg, age6_avg) %>%
                        data.frame() %>% 
                        mutate( proj3 = pred_cols[,i] * age3_avg, 
                                proj4 = pred_cols[,i] * age4_avg, 
                                proj5 = pred_cols[,i] * age5_avg,
                                proj6 = pred_cols[,i] * age6_avg ) }
                names( pred_list ) <- c( colnames( pred_cols ) )                

pred_list <- lapply( pred_list, function (x) select( x, juv_year, proj3, proj4, proj5, proj6 ) ) 

pred_list <- lapply( pred_list, na.omit )

        return( pred_list ) }

# Canadian-origin projected survivors
cdn_survs <- projSurvs( cdn_proj_data )

# Need to calaculate average age-3 and age-4 components to complete last two year forecasts
# Through 2024 juveniles, we have complete forecasts through 2025 but need age-3
# component for 2026 and age-3 and age-4 component for 2027

avg_age3_cdn <- lapply( cdn_survs, function(x) mean( x[c(nrow(x):(nrow(x)-2)), "proj3"] ) )

avg_age4_cdn <- lapply( cdn_survs, function(x) mean( x[c(nrow(x):(nrow(x)-2)), "proj4"] ) )

# Total Yukon projected survivors
total_survs <- projSurvs( total_proj_data )

avg_age3_total <- lapply( total_survs, function(x) mean( x[c(nrow(x):(nrow(x)-2)), "proj3"] ) )

avg_age4_total <- lapply( total_survs, function(x) mean( x[c(nrow(x):(nrow(x)-2)), "proj4"] ) )

# Sum across each of the diagonals of the matrix to get the projections
# The earliest projection will be the age-3s in 2004, the latest projection
# will be the last juvenile year plus 4
forecast_years <- c( seq( min( total_juvs$year ) + 1, max( total_juvs$year ) + 4 ) )

calcForecast <- function(dat) {
        
        mat <- lapply( dat, function(x) as.matrix(x) )
        
        # Remove first column (juv_year)
        tidy_dat <- lapply( mat, function(x) x[,-1] )
        # Create indicator across diagonals
        diag_ind <- lapply( tidy_dat, function(x) col(x) + row(x) )
        # Split the matrix using the indicators created above
        diag_split <- lapply( tidy_dat, function(x) split( x, diag_ind[[1]] ) )
        # Sum the diagonals to create annual forecast
        forecast <- lapply( diag_split, function(x) plyr::ldply( x, sum ) ) 

final_forecast <- lapply( forecast, function(x) 
        data.frame( x ) %>% 
        select( -.id ) %>% 
        rename( forecast = V1 ) %>% 
        mutate( year = forecast_years ) )
                
        return( final_forecast )
        
}

cdn_forecast <- calcForecast( cdn_survs )

total_forecast <- calcForecast( total_survs )

# Need to add average age-3 and age-4 components to 2026 and 2027 forecasts to
# get final forecast. First melt lists into data frames for easier 
# manipulation and only include years with complete forecasts

cdn_forecast <- as.data.frame(do.call( cbind, cdn_forecast ) ) %>% 
        select( pred.surv.fit.year, pred.surv.fit.forecast, pred.surv.lwr.forecast, 
                pred.surv.upr.forecast ) %>% 
        rename( forecast_fit = pred.surv.fit.forecast,
                forecast_lwr = pred.surv.lwr.forecast,
                forecast_upr =pred.surv.upr.forecast,
                year = pred.surv.fit.year ) %>% 
        filter( year > 2006 & year < 2028 )
                
total_forecast <- as.data.frame(do.call( cbind, total_forecast ) ) %>% 
        select( pred.surv.fit.year, pred.surv.fit.forecast, pred.surv.lwr.forecast, 
                pred.surv.upr.forecast ) %>% 
        rename( forecast_fit = pred.surv.fit.forecast,
                forecast_lwr = pred.surv.lwr.forecast,
                forecast_upr =pred.surv.upr.forecast,
                year = pred.surv.fit.year ) %>% 
        filter( year > 2006 & year < 2028 )

# Add average age-3 and age-4 components to 2026 and 2027 forecasts to
# get final forecast - Canadian
adj_2026_cdn <- cdn_forecast %>%
        filter( year == 2026 ) %>%       
        mutate( forecast_fit = forecast_fit + avg_age3_cdn$pred.surv.fit,
                forecast_lwr = forecast_lwr + avg_age3_cdn$pred.surv.lwr,
                forecast_upr = forecast_upr + avg_age3_cdn$pred.surv.upr )


adj_2027_cdn <- cdn_forecast %>%
        filter( year == 2027 ) %>%       
        mutate( forecast_fit = forecast_fit + avg_age3_cdn$pred.surv.fit +
                        avg_age4_cdn$pred.surv.fit,
                forecast_lwr = forecast_lwr + avg_age3_cdn$pred.surv.lwr +
                        avg_age4_cdn$pred.surv.lwr,
                forecast_upr = forecast_upr + avg_age3_cdn$pred.surv.upr +
                        avg_age4_cdn$pred.surv.upr )

# Final Canadian forecast
cdn_final_forecast <- cdn_forecast %>% 
        filter( year < 2026 ) %>% 
        rows_insert( rbind( adj_2026_cdn, adj_2027_cdn ) )

# Write final Canadian forecast to csv file
write.csv( cdn_final_forecast, paste0("Output/CanadianForecast_",Sys.Date(),".csv") )

# Add average age-3 and age-4 components to 2026 and 2027 forecasts to
# get final forecast - Total Yukon
adj_2026_total <- total_forecast %>%
        filter( year == 2026 ) %>%       
        mutate( forecast_fit = forecast_fit + avg_age3_total$pred.surv.fit,
                forecast_lwr = forecast_lwr + avg_age3_total$pred.surv.lwr,
                forecast_upr = forecast_upr + avg_age3_total$pred.surv.upr )


adj_2027_total <- total_forecast %>%
        filter( year == 2027 ) %>%       
        mutate( forecast_fit = forecast_fit + avg_age3_total$pred.surv.fit +
                        avg_age4_total$pred.surv.fit,
                forecast_lwr = forecast_lwr + avg_age3_total$pred.surv.lwr +
                        avg_age4_total$pred.surv.lwr,
                forecast_upr = forecast_upr + avg_age3_total$pred.surv.upr +
                        avg_age4_total$pred.surv.upr )

# Final Total Yukon forecast
total_final_forecast <- total_forecast %>% 
        filter( year < 2026 ) %>% 
        rows_insert( rbind( adj_2026_total, adj_2027_total ) )

# Write final Total Yukon forecast to csv file
write.csv( total_final_forecast, paste0("Output/TotalYukonForecast_",Sys.Date(),".csv") )

