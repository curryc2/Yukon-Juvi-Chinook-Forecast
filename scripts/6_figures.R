##############################################################################
##############################################################################
# SECTION 6. Figures
# Last updated: Janaury 2, 2025
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
require(here)

wd <- here()  
setwd(wd)

dir.data <- file.path(wd, "data")
# NBS juvenile abundance
nbs_juv_plot <- ggplot( total_juvs, aes( x = year, y = juv_abun_adj/1000 ) ) + 
     geom_bar( stat = "identity", fill = "light gray" ) + 
     geom_line( data = total_juvs, 
                aes( x = year, y = mean( juv_abun_adj/1000 ) ), colour ="black", 
                linewidth = 0.75 ) +
     scale_x_continuous( name="Juvenile Year", breaks = seq( min(all_juvs$year), 
                                                             max (all_juvs$year), 2 ), 
                         labels = seq( min(all_juvs$year), 
                                      max (all_juvs$year), 2 ),
                         expand = c( 0.01,0 ) ) +
     scale_y_continuous( name ='NBS Juvenile Abundance (000s)',
                         labels = comma, limits = c( 0,7000 ), expand = c( 0,0 ) ) + 
     geom_errorbar( total_juvs, mapping = aes( ymin=( juv_abun_adj/1000 )-( juv_sd/1000 ), 
                                               ymax=( juv_abun_adj/1000 )+( juv_sd/1000 ) ), 
                    width = 0.4, position=position_dodge( 0.9 ) ) + 
     theme_classic( base_size = 14, base_family= "Times" ) +
     theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text( angle = 45, hjust = 1 ) )

# Canadian-origin juvenile abundance
cdn_juv_plot <- ggplot( all_juvs, aes( x = year, y = cdn_abun/1000 ) ) + 
     geom_bar( stat="identity", fill= "light gray" ) + 
     geom_line( data = all_juvs, aes( x = year, y = mean( cdn_abun/1000 ) ), 
                colour ="black", linewidth = 0.75 ) +
     scale_x_continuous( name = "Juvenile Year", breaks = seq( min(all_juvs$year), 
                                                               max (all_juvs$year), 2 ), 
                         labels = seq( min(all_juvs$year), 
                                       max (all_juvs$year), 2 ), expand = c( 0.01,0 ) ) +
     scale_y_continuous( name ='Canadian-origin Juvenile Abundance (000s)', 
                         labels = comma, limits=c( 0,4000 ), expand = c( 0,0 ), 
                         breaks = c( 0,500,1000,1500,2000,2500,3000,3500,4000 ) ) + 
     theme_classic( base_size = 14, base_family= "Times" ) +
     geom_errorbar( aes( ymin = (cdn_abun/1000)-(cdn_joint_sd/1000), 
                         ymax = (cdn_abun/1000)+(cdn_joint_sd/1000) ), width=.4,
                    position=position_dodge(0.9) ) 

# Total Yukon juvenile abundance
totyuk_juv_plot <- ggplot( all_juvs, aes( x = year, y = totyuk_abun/1000 ) ) + 
     geom_bar( stat="identity", fill= "light gray" ) + 
     geom_line( data = all_juvs, aes( x = year, y = mean( totyuk_abun/1000 ) ), 
                colour ="black", linewidth = 0.75 ) +
     scale_x_continuous( name = "Juvenile Year", breaks = seq( min(all_juvs$year), 
                                                               max (all_juvs$year), 2 ), 
                         labels = seq( min(all_juvs$year), 
                                       max (all_juvs$year), 2 ), expand = c( 0.01,0 ) ) +
     scale_y_continuous( name ='Total Yukon Juvenile Abundance (000s)', 
                         labels = comma, limits = c( 0,7000 ), expand = c( 0,0 ) ) + 
     theme_classic( base_size = 14, base_family = "Times" ) +
     geom_errorbar( aes( ymin = (totyuk_abun/1000)-(totyuk_joint_sd/1000), 
                         ymax = (totyuk_abun/1000)+(totyuk_joint_sd/1000) ), width = 0.4,
                    position = position_dodge(0.9) ) 

# Linear regression between CDN juveniles and CDN adult returns from completed cohorts
# Black lines are the 80% prediction interval
# Gray shaded area is the 80% confidence interval
formula <- y ~ x
cdn_lm_plot <- drop_na( cdn_model_fit, juv_abun, returns ) %>% 
     ggplot( aes( x = juv_abun/1000, y = returns/1000 ) ) +
     geom_smooth( method = lm, se = TRUE, level = 0.80, color="gray30", fill = "gray69" ) +
     stat_poly_eq( formula = y~x, aes( label = paste(..eq.label.., ..rr.label.., sep = "~~~")), label.x = "left", label.y = "top", parse = TRUE, size = 5 ) +
     geom_point( size = 3 ) +
     geom_text_repel( data = drop_na( cdn_model_fit, juv_abun, returns ), 
                      aes(x = juv_abun/1000, y = returns/1000, label= juv_year ),
                      family = "Times", point.padding = 0.5, nudge_x = -0.1, nudge_y = 0.3 ) +
     scale_y_continuous( name ='Canadian-origin Adult Returns (000s)', 
                         limits = c( 0 ,150 ), expand = c(0,0),
                         breaks = c(0, 25, 50, 75, 100, 125, 150) ) +
     scale_x_continuous( name ='Canadian-origin Juvenile Abundance (000s)', 
                         limits = c( 300,3000 ),
                         labels = comma, breaks = c(seq(300,3000,300)), expand= c( 0,20 ) ) +
     theme_classic( base_size = 14, base_family = "Times" ) +
     geom_line( aes( x = juv_abun/1000, y = pred.surv.lwr/1000 ), 
                                  color = "black", linetype = "dashed", linewidth = 1 ) +
     geom_line( aes( juv_abun/1000, y = pred.surv.upr/1000 ), color = "black", 
                linetype = "dashed", linewidth = 1 ) 

# Linear regression between total Yukon juveniles and total Yukon adult returns
# from completed cohorts
# Black lines are the 80% prediction interval
# Gray shaded area is the 80% confidence interval
tot_lm_plot <- drop_na( total_model_fit, juv_abun, returns ) %>% 
     ggplot( aes( x = juv_abun/1000, y = returns/1000 ) ) +
     geom_smooth( method = lm, se = TRUE, level = 0.80, color="gray30", fill = "gray69" ) +
     stat_poly_eq( formula = y~x, aes( label = paste(..eq.label.., ..rr.label.., sep = "~~~")), label.x = "left", label.y = "top", parse = TRUE, size = 5 ) +
     geom_point( size = 3 ) +
     geom_text_repel( data = drop_na( total_model_fit, juv_abun, returns ), 
                      aes(x = juv_abun/1000, y = returns/1000, label = juv_year ),
                      family = "Times", point.padding = 0.5, nudge_x = -0.1, nudge_y = 0.3 ) +
     scale_y_continuous( name ='Total Yukon Adult Returns (000s)', 
                         limits = c( 0, 400 ), expand = c(0,0) ) +
     scale_x_continuous( name ='Total Yukon Juvenile Abundance (000s)', limits = c( 1000,5500 ),
                         labels = comma, expand= c( 0,20 ) ) +
     theme_classic( base_size = 14, base_family = "Times" ) +
     geom_line( aes( x = juv_abun/1000, y = pred.surv.lwr/1000 ), 
                color = "black", linetype = "dashed", linewidth = 1 ) +
     geom_line( aes( juv_abun/1000, y = pred.surv.upr/1000 ), color = "black", 
                linetype = "dashed", linewidth = 1 )

# Total CDN run size and CDN adult forecast        
cdn_forecast_plot <- ggplot( data = select( cdn_brood, brood_year, total_run ) %>% 
                  filter( brood_year > 2000 & brood_year < 2027 ),
                  aes( x = brood_year, y = total_run/1000 ) ) +
     geom_bar( stat = "identity", width = 0.5, fill = "gray67"  ) +
     geom_line( data = cdn_final_forecast, aes( x = year, y = forecast_fit/1000), 
                linetype = "dashed", linewidth = 1 ) +  
     geom_errorbar( data = cdn_final_forecast,
                    aes( x = year, ymin = forecast_lwr/1000, 
                         ymax = forecast_upr/1000 ), size = 1, inherit.aes = FALSE ) +
     scale_x_continuous( name = "Run Year", labels = c( seq( 2001, 2027, 2 ) ), 
                         breaks = c( seq( 2001, 2027, 2 ) ),
                         expand = expansion( add = c( 0.10, 0.40 ) ) ) +
     scale_y_continuous( name = "Canadian-origin Chinook Salmon Run (000s)",
                         labels = comma, expand = c( 0, 0 ), limits = c( 0, 180 ) ) +
     theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = "none", 
            axis.text.x = element_text( size = 11, angle= 45, hjust = 1 ),
            axis.text.y = element_text( size = 11, angle= 45, hjust = 1 ) ) +
     theme_bw( base_size = 14, base_family = "Times")

# Total Yukon run size and Yukon adult forecast        
total_forecast_plot <- ggplot( data = select( total_brood, brood_year, total_run ) %>% 
                                  filter( brood_year > 2000 & brood_year < 2027 ),
                             aes( x = brood_year, y = total_run/1000 ) ) +
     geom_bar( stat = "identity", width = 0.5, fill = "gray67"  ) +
     geom_line( data = total_final_forecast, aes( x = year, y = forecast_fit/1000), 
                linetype = "dashed", linewidth = 1 ) +  
     geom_errorbar( data = total_final_forecast,
                    aes( x = year, ymin = forecast_lwr/1000, 
                         ymax = forecast_upr/1000 ), size = 1, inherit.aes = FALSE ) +
     scale_x_continuous( name = "Run Year", labels = c( seq( 2001, 2027, 2 ) ), 
                         breaks = c( seq( 2001, 2027, 2 ) ),
                         expand = expansion( add = c( 0.10, 0.40 ) ) ) +
     scale_y_continuous( name = "Total Yukon Chinook Salmon Run (000s)",
                         labels = comma, expand = c( 0, 0 ), limits = c( 0, 400 ) ) +
     theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = "none", 
            axis.text.x = element_text( size = 11, angle= 45, hjust = 1 ),
            axis.text.y = element_text( size = 11, angle= 45, hjust = 1 ) ) +
     theme_bw( base_size = 14, base_family = "Times")

# Canadian juveniles per spawner and Canadian spawner abundance
juv_years <- cdn_all_dat %>% 
     drop_na( juv_year ) %>% 
     select( juv_year ) 

cdn_jps_plot <- drop_na( cdn_all_dat, cdn_abun ) %>% 
     ggplot() +
     geom_bar( mapping = aes( x = juv_year, y = juv_per_spawn ), 
               stat = "identity", fill="gray67" ) +
     geom_line( mapping = aes( x = juv_year, y = spawners/1000 ),                                                   linetype = "dashed", linewidth = 1.0 ) +
     scale_x_continuous( name = "Juvenile Year", 
                         breaks = c( seq( min(juv_years), max(juv_years), 2 ) ), 
                         labels = c( seq( min(juv_years), max(juv_years), 2  ) ), 
                         expand = c( 0.01,0 ) ) +
     scale_y_continuous( name ='Canadian-origin Juveniles Per Spawner', labels = comma, 
                         limits = c( 0,100 ), expand = c( 0,0 ), 
                         breaks = c( 0,20,40,60,80,100 ), 
                         sec.axis = sec_axis( ~ . * 1, 
                         name = "Canadian-origin Spawner Abundance (000s)", 
                         breaks = c( 0, 25, 50, 75, 100 ) ) ) + 
     theme_classic( base_size = 14 , base_family = "Times" )

# Total Yukon juveniles per spawner and Total Yukon spawner abundance
total_jps_plot <- drop_na( total_all_dat, totyuk_abun ) %>% 
     ggplot() +
     geom_bar( mapping = aes( x = juv_year, y = juv_per_spawn ), 
               stat = "identity", fill = "gray67" ) +
     geom_line( mapping = aes( x = juv_year, y = spawners/5000, group = 1 ),                                                   linetype = "dashed", linewidth = 1.0 ) +
     scale_x_continuous( name = "Juvenile Year", 
                         breaks = c( seq( min(juv_years), max(juv_years), 2 ) ), 
                         labels = c( seq( min(juv_years), max(juv_years), 2  ) ), 
                         expand = c( 0.01,0 ) ) +
     scale_y_continuous( name ='Total Yukon Juveniles Per Spawner', labels = comma, 
                         limits = c( 0,60 ), expand = c( 0,0 ), 
                         breaks = c( 0,20,40,60), 
                         sec.axis = sec_axis( ~ . * 5000, 
                                              name = "Total Yukon Spawner Abundance", 
                                              breaks = c( 0, 50000, 100000, 150000, 200000, 250000, 300000 ), labels = comma ) ) +
     theme_classic( base_size = 14 , base_family= "Times" )

# Create list to store ggplots
plot.list <- list()
plot.list[[1]] = nbs_juv_plot 
plot.list[[2]] = cdn_juv_plot
plot.list[[3]] = totyuk_juv_plot
plot.list[[4]] = cdn_lm_plot
plot.list[[5]] = tot_lm_plot
plot.list[[6]] = cdn_forecast_plot
plot.list[[7]] = total_forecast_plot
plot.list[[8]] = cdn_jps_plot 
plot.list[[9]] = total_jps_plot  

# For loop will save all files to the current working directory
# Files will be names "file1", "file2", etc and will then need to be
# manually named
for( t in 1:length( plot.list ) ) {
     
     ggsave( plot = plot.list[[t]], 
             file = paste0("figs/figure_",t,"_",Sys.Date(),".jpeg") ) 
     
}
