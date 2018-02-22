# R script, JAGS models, and data used to estimate urban bat activity at multiple scales and predict activity across the Chicago metropolitan area, IL USA.
## **Gallo, T., M. Fidino, E.W. Lehrer, J. Kilgour, P. Wolff, and S. Magle. Need for multiscale planning for conservation of urban bats. _Conservation Biology_ DOI: 10.1111/cobi.13047**
# **File Descriptions:**
gallo_et_al_2018_bat_landscape_analysis.R – the only file that you need to open. Loads datasets and formats data for JAGS model. Has code for summarizing the posterior distributions and plotting model results and predictive raster maps.

gallo_et_al_20198_bat_model_lasso.R - JAGS model used to estimate bat activity at multiple scales. Read in gallo_et_al_2018_bat_landscape_analysis.R.

gallo_et_al_20198_bat_model_lasso_loocv_model.r - JAGS model using leave-one-add cross validation methods to estimate the likelihood that our original model could predict out of sample data points. This method is used to assess model fit and predictive power. Read in gallo_et_al_2018_bat_landscape_analysis.R.

gallo_et_al_20198_bat_model_lasso_loocv_null_model.r - JAGS model using leave-one-add cross validation methods to estimate the likelihood that a null model could predict out of sample data points. This method is used to assess model fit and predictive power. Read in gallo_et_al_2018_bat_landscape_analysis.R.

gallo_et_al_2018_bat_detections.RDS - Data on the number of days big brown bat, eastern red bat, and silver-haired bat were detected at each sampling site during each sampling season.

gallo_et_al_2018_j_matrix.RDS - Data on the number of days an acoustic bat monitor was active at each site during each season. If a site was not sampled a zero is reported.

gallo_et_al_2018_landcover_covariates.RDS - The proportion of tree cover, open vegetation cover, and impervious cover at each site at 100-m, 500-m, and 1000-m scales. Data was extracted using the 2010_HRLC_North_Central_Merge.tiff data and a R function located at https://github.com/TravisGallo/Functions. For the sake of landowner privacy we have decided not to include the coordinates of individual study sites in this analysis.

gallo_et_al_2018_bat_detection_FullData.RDS - The number of days each species (all species detected) was detected and is used to caculate summary statistics about species richness and compare relative activity.

2010_HRLC_North_Central_Merge.tiff - THIS DATA IS NOT LOCATED IN THE REPOSITORY DUE TO SIZE LIMITS. AVALIABLE UPON REQUEST. A merged raster of the 2010 High-Resolution Land Cover, NE Illinois and NW Indiana data used to extract landcover data at multiple scales. Origninal files can  be found at https://datahub.cmap.illinois.gov/dataset/high-resolution-land-cover-ne-illinois-and-nw-indiana-2010. Associated file is 2010_HRLC_North_Central_Merge.tif.aux.xml. CRS:26916

gallo_et_al_2018_grid_100m_reduced.shp (and associated files) - Grid of spatial points spaced 100m apart placed across the greater study region. These points were used to extract landcover data and predict bat activity across the greater Chicago metropolitan area. CRS:26916

gallo_et_al_2018_urban_points.shp (and associated files) - Spatial points associated with gallo_et_al_2018_grid_100m_reduced.shp that indicated which points within gallo_et_al_2018_grid_100m_reduced.shp are "urban" sites. CRS:26916

ChicagoRiver_clip2.shp (and associated files) - Spatial polygon of the Chicago river to plot with predictive raster data. CRS:26916

FoxRiver_Clip.shp (and associated files) - Spatial polygon of the Fox river to plot with predictive raster data. CRS:26916

lake_michigan_shoreline.shp (and associated files) - Spatial polygon of Lake Michigan to plot with predictive raster data. CRS:26916

**Note:** All of these files must be within your working directory for the analysis to work.  Our analysis was done in parallel and used all but two of the cores in a computer. Therefore, if you have two or less cores on your computer you will need to adjust the settings annotated within gallo_et_al_2018_bat_landscape_analysis.R.
