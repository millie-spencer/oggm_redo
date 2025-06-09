
rm(list = ls())

library(ncdf4)
library(raster)

#
#  LAMADO CR2MET y DEM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

#abrir CR2met
vtmean = nc_open("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/CR2MET_tmean_v2.5_mon_1960_2021_005deg.nc")
vp = nc_open    ("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/CR2MET_pr_v2.5_month_1960_2021_005deg.nc")
dem = raster("E:/HBV/Elevation_TC/actualizar_con_nuevos_glaciares/input/raster/srtm_100m_geo.tif")
print(vtmean)
print(vp)
print(dem)


#
#  resample DEM        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

# Raster dem to same climate resolution (misma resolucion espacial y extensi�n)

  # lectura 1er raster desde nc tem
t_grid = brick("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/CR2MET_tmean_v2.5_mon_1960_2021_005deg.nc", varname="tmean")
#plot(t_grid[[1]])
t_grid_1 = t_grid[[1]]
#plot(t_grid_1)

  # resample de DEM and clip extension usando t_grid_1
dem_resam = resample(dem,t_grid_1,method='ngb')

t_grid[[1]]       # ac� # 2.5 resolution : 0.05, 0.05  (x, y)                     dimensions : 800, 220, 176000  (nrow, ncol, ncell)

print(dem)        # resolution : 0.001225577, 0.001225577  (x, y)   dimensions : 60372, 31054, 1874792088  (nrow, ncol, ncell)
print(dem_resam)  # resolution : 0.05, 0.05  (x, y)          dimensions : 800, 220, 176000  (nrow, ncol, ncell)

plot(dem_resam) 
plot(t_grid_1) 


# 1er archivo. Temp mean
# NUEVO NETCDF ++++++++++++++++++++++++++++++++++++++++++++++

      # First, create the netCDF filename:
ncpath <- "C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/"
ncname <- "CR2met_t2m_hgt_2022_1960_dic_2021_2.5"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

      # Then define the contents of the file:
 
         # extracion te mean
tmp_mean <- ncvar_get(vtmean,"tmean")

         # extracci�n desde original
lon <- ncvar_get(vtmean,"lon")
lat <- ncvar_get(vtmean,"lat")
#time <- ncvar_get(vt,"time")
time = as.numeric(c(1:744))
tunits <- ncatt_get(vtmean,"time","units")
tunits =  "months since 1960-01-01"
         # create and write the netCDF file -- ncdf4 version
         # define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
latdim <- ncdim_def("lat","degrees_north",as.double(lat)) 
timedim <- ncdim_def("time",tunits,as.double(time))

          # define variables
fillvalue <- 1e32
dlname <- "air temperature month"
tmp_def <- ncvar_def("temp","degC",list(londim,latdim,timedim),fillvalue,dlname,prec="single")
dlname2 <- "elevation"
ele.def <- ncvar_def("hgt","m",list(londim,latdim),fillvalue,dlname2,prec="single")

      # Next, create the file, and put the variables into it, 
      # along with additional variable and "global" attributes (those that apply to the whole file). 
      # Note that the attributes are of key importance to the self-documenting properties of netCDF files.

           # create netCDF file and put arrays
ncout <- nc_create(ncfname,list(tmp_def,ele.def),force_v4=TRUE)

           # put variables
               # to array
tmp_array <- array(tmp_mean, dim=c(NROW(lon),NROW(lat),NROW(time)))
#ele_array <- array(dem_resam, dim=c(NROW(lon),NROW(lat)))
#ele_array = rev(ele_array)
library(oceanmap)
ele_array = raster2array(dem_resam)

           # env�o de variables 
ncvar_put(ncout,tmp_def,tmp_array)
ncvar_put(ncout,ele.def,ele_array)


      # put additional attributes into dimension and data variables
#ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
#ncatt_put(ncout,"lat","axis","Y")
#ncatt_put(ncout,"time","axis","T")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)

# revision visual

pp = brick("/Users/milliespencer/Desktop/CR2_OGGM_Paper/files_chile_OGGM_climate_comparison/CR2met_t2m_hgt_2022_1960_dic_2021_2.5.nc", varname="temp")
plot(pp[[744]]) # ene 1960 a dic 2021





# 2 archivo. Solo Precipitacion 
# NUEVO NETCDF ++++++++++++++++++++++++++++++++++++++++++++++

# First, create the netCDF filename:
ncpath <- "C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/"
ncname <- "CR2met_pr_2022_1960_dic_2021_2.5"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

# Then define the contents of the file:

# extracion variables
pre <- ncvar_get(vp,"pr_month")

# extracci�n desde original
lon <- ncvar_get(vp,"lon")
lat <- ncvar_get(vp,"lat")
#time <- ncvar_get(vp,"time")
time = as.numeric(c(1:744))
tunits <- ncatt_get(vp,"time","units")
#tunits = as.character(tunits[2])
tunits =  "months since 1960-01-01"

# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
latdim <- ncdim_def("lat","degrees_north",as.double(lat)) 
timedim <- ncdim_def("time",tunits,as.double(time))

# define variables
fillvalue <- 1e32
dlname <- "precipitation month"
pr.def <- ncvar_def("prcp","mm",list(londim,latdim,timedim),fillvalue,dlname,prec="single")

# Next, create the file, and put the variables into it, 
# along with additional variable and "global" attributes (those that apply to the whole file). 
# Note that the attributes are of key importance to the self-documenting properties of netCDF files.

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(pr.def),force_v4=TRUE)

# put variables
# to array
pr_array <- array(pre, dim=c(NROW(lon),NROW(lat),NROW(time)))

# env�o de variables 
ncvar_put(ncout,pr.def,pr_array)


# put additional attributes into dimension and data variables
#ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
#ncatt_put(ncout,"lat","axis","Y")
#ncatt_put(ncout,"time","axis","T")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)


# revision visual

pp = brick("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/CR2met_pr_2022_1960_dic_2021_2.5.nc", varname="prcp")
plot(pp[[744]]) 

#mean(pp[[1]])
#mean(pp[[13]])
#mean(pp[[25]])
#mean(pp[[2]])
#mean(pp[[14]])







# comparaci�n P 2.5 y P 2.0

p25x = brick("/Users/milliespencer/Desktop/CR2_OGGM_Paper/files_chile_OGGM_climate_comparison/CR2met_pr_2022_1960_dic_2021_2.5.nc", varname="prcp")
p20 = brick("C:/Users/caropara/Desktop/OGGM/cr2met/CR2met_pr_2022_1950_dic_2020.nc", varname="prcp")

library(raster)
library(rgdal)
mask = shapefile("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/shp/mask.shp")
plot(mask)

P25_2 = crop(p25x,mask)
P20_2 = crop(p20,mask)


P25_2 = cellStats(P25_2, stat='mean', na.rm=TRUE)
P20 = cellStats(P20_2, stat='mean', na.rm=TRUE)
P25 = as.data.frame(as.numeric(P25_2))
P20 = as.data.frame(as.numeric(P20))


P20$date = seq(as.Date("1950/1/1"), by = "month", length.out = 852)
P25$date = seq(as.Date("1960/1/1"), by = "month", length.out = 744)

P = merge(P20,P25,by="date")
P$year = lubridate::year(P$date)

plot(P$date,P$`as.numeric(P20)`, type="l")
lines(P$date,P$`as.numeric(P25_2)`, type="l", col="red")

mean(P$`as.numeric(P20)`)
mean(P$`as.numeric(P25_2)`)


library(stats)
P_mean = aggregate(. ~ year, data = P, FUN = sum) # promedio por grupo=year

mean_25 = round(mean(P_mean$`as.numeric(P25_2)` ),0)
mean_20 = round(mean(P_mean$`as.numeric(P20)`),0)  


f1=
plot(P_mean$year, P_mean$`as.numeric(P25_2)`, type='l', col = 'red',ylim = c(0, 1700) )
lines(P_mean$year, P_mean$`as.numeric(P20)`, type='l', col = 'blue')
legend("topleft",
       c(paste("prec 2.5   ",mean_25,' mm/yr',sep = ''),paste("prec 2.0   ",mean_20,' mm/yr',sep = '')),
       fill=c("red","blue")
)



# comparaci�n T 2.5 y T 2.0

p25 = brick("C:/Users/caropara/Desktop/OGGM/cr2met/cr2met_2.5_octubre/CR2met_t2m_hgt_2022_1960_dic_2021_2.5.nc", varname="temp")
p20 = brick("C:/Users/caropara/Desktop/OGGM/cr2met/CR2met_t2m_hgt_2022_1950_dic_2020.nc", varname="temp")


P25_2 = crop(p25,mask)
P20_2 = crop(p20,mask)


P25_2 = cellStats(P25_2, stat='mean', na.rm=TRUE)
P20 = cellStats(P20_2, stat='mean', na.rm=TRUE)
P25 = as.data.frame(as.numeric(P25_2))
P20 = as.data.frame(as.numeric(P20))


P20$date = seq(as.Date("1950/1/1"), by = "month", length.out = 852)
P25$date = seq(as.Date("1960/1/1"), by = "month", length.out = 744)


P2 = merge(P20,P25,by="date")
P2$year = lubridate::year(P2$date)

plot(P2$date,P$`as.numeric(P20)`, type="l")
lines(P2$date,P$`as.numeric(P25_2)`, type="l", col="red")

mean(P2$`as.numeric(P20)`)
mean(P2$`as.numeric(P25_2)`)


library(stats)
P_mean = aggregate(. ~ year, data = P2, FUN = mean) # promedio por grupo=year

mean_25 = round(mean(P_mean$`as.numeric(P25_2)` ),1)
mean_20 = round(mean(P_mean$`as.numeric(P20)`),1)  


f2=
  plot(P_mean$year, P_mean$`as.numeric(P25_2)`, type='l', col = 'red',ylim = c(7, 11) )
lines(P_mean$year, P_mean$`as.numeric(P20)`, type='l', col = 'blue')
legend("topleft",
       c(paste("temp 2.5   ",mean_25,' mm/yr',sep = ''),paste("temp 2.0   ",mean_20,' mm/yr',sep = '')),
       fill=c("red","blue")
)




