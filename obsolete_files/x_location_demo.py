import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pickle
import matplotlib.patches as patches

fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(10,8),subplot_kw={'projection':ccrs.PlateCarree()})
gl=axs.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
axs.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
axs.add_feature(cfeature.STATES,edgecolor='black',linewidth=0.5,facecolor='none',linestyle='--')

esa_dict=pickle.load(open('ESA.pkl','rb'))
latCell=esa_dict['lat_1d']
lonCell=esa_dict['lon_1d']
great_plains_mask=(latCell>=37)*(latCell<=44)*(lonCell>=360-104)*(lonCell<=360-97)
california_mask=(latCell>=32.5)*(latCell<=42)*(lonCell>=360-124.5)*(lonCell<=360-114.1)
florida_mask=(latCell>=24.4)*(latCell<=31)*(lonCell>=360-87.6)*(lonCell<=360-80)
russia_mask=(latCell>=70)*(latCell<=90)*(lonCell>=90)*(lonCell<=120)
Wcasp_mask=(latCell>=45)*(latCell<=50)*(lonCell>=30)*(lonCell<=40)
SEcasp_mask=(latCell>=30)*(latCell<=40)*(lonCell>=55)*(lonCell<=65)
Waus_mask=(latCell>=-25)*(latCell<=-15)*(lonCell>=95)*(lonCell<=105)

lon_min1=360-104
lon_min2=360-124.5
lon_min3=360-87.6
lon_min4=90
lon_min5=30
lon_min6=55
lon_min7=95

lon_max1=360-97
lon_max2=360-114.1
lon_max3=360-80
lon_max4=120
lon_max5=40
lon_max6=65
lon_max7=105

lat_min1=37
lat_min2=32.5
lat_min3=24.4
lat_min4=70
lat_min5=45
lat_min6=30
lat_min7=-25

lat_max1=44
lat_max2=42
lat_max3=31
lat_max4=90
lat_max5=50
lat_max6=40
lat_max7=-15

ranges = [
    (lon_min1, lon_max1, lat_min1, lat_max1),
    (lon_min2, lon_max2, lat_min2, lat_max2),
    (lon_min3, lon_max3, lat_min3, lat_max3),
    (lon_min4, lon_max4, lat_min4, lat_max4),
    (lon_min5, lon_max5, lat_min5, lat_max5),
    (lon_min6, lon_max6, lat_min6, lat_max6),
    (lon_min7, lon_max7, lat_min7, lat_max7),
]

for lon_min, lon_max, lat_min, lat_max in ranges:
    lon_min = (lon_min + 180) % 360 - 180
    lon_max = (lon_max + 180) % 360 - 180

    rect = patches.Rectangle(
        (lon_min, lat_min),
        lon_max - lon_min,
        lat_max - lat_min,
        linewidth=2,
        edgecolor='red',
        facecolor='none',
        transform=ccrs.PlateCarree()
    )
    axs.add_patch(rect)

axs.set_global()
plt.savefig('globe_demo.png')
plt.close()