# -*- coding: utf-8 -*-
"""
To iterate calculation of Fresnel Zone and show reflexion points for GNSS-R

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede, NL), Aug 26, 2022
"""


# Usefull librairies
from fresnelzone import *
from getorbits import *
from plotfresnel import *
from geod import *
import geopandas as gpd

# Usefull constants
c = 299792458                   # m.s-1 Speed of light
L1_GPS = 1575.42e6              # Hz L1 frequency for GPS
L2_GPS = 1227.60e6              # Hz L2 frequency for GPS
L1_Glo = 1602.0e6               # Hz L1 frequency for GLONASS
L2_Glo = 1246.0e6               # Hz 21 frequency for GLONASS
lambda_L1_GPS = (c/L1_GPS)      # m wavelenght for L1 GPS
lambda_L2_GPS = (c/L2_GPS)      # m wavelenght for L2 GPS
lambda_L1_Glo = (c/L1_Glo)      # m wavelenght for L1 Glo
lambda_L2_Glo = (c/L2_Glo)      # m wavelenght for L2 Glo

###############################################################################

def reflZone(h, elev_list, lon, lat, output='gpkg', dirName=None, df_sp3=None):
    """
    Iterate the calculation of the first Fresnel zone by taking into account the visible satellites.
    Possibily to choose output format, default is gpkg but can be set to shp or kml.
    
    Parameters
    ----------
    h: float
        Hight of the receiver.
    elev: list
        List of elevation angle in degrees (should not exceed 4 values for vizualisation).
    lon,lat: float
        Position of the receiver in geographical coordinates (degrees).
    h: float
        Hight of the receiver in meters. 
    output: str
        Output format of file, can be either shp, gpkg or kml. Default is gpkg.
    dirName: str
        Directory on which files will be stored.
    df_sp3: DataFrame
        If you wantto give a modified DataFrame of avaible satellites. Default
        uses latest orbits file with no modification on it.
    
    Return
    ------
    gdf: GeoDataFrame
        GeoDataFrame that is exported to the desired format.
    """
    if df_sp3==None:
        # Get orbits file
        orb = retrieve_orbits()
        df = read_sp3(orb)

        # Calculate elevation and azimuth
        df_sp3 = elevazim(df,lon,lat,h)

        # Remove all values bellow 0° elevation, i.e satellites that are not visible
        df_sp3 = df_sp3[df_sp3['elevation'] > 0]
    
    # Empty list to store Polygon and Polygons area
    pol = []    
    el = []
    larea = []
    az = []
    # Create directory
    if dirName==None:
        dirName = 'FresnelZones'
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print(f"Files will be stored in directory {dirName}")
            
    # Check if elevation is a list or a single value
    if isinstance(elev_list, list):
        # Calcule First Fresnel Surface for each elevation angles
        c = 0
        for elev in elev_list:
            df_int = df_sp3.loc[(df_sp3['elevation'] == elev)]
            azim = df_int["azimuth"]
            azim = azim.sort_values()
            # Prevent to calculate several time the same azimuth
            list_azim = []
            for val in azim:
                if val not in list_azim:
                    list_azim.append(val)
            for angle in list_azim:
                a,b,R = firstFresnelZone(L1_GPS, h, elev)
                ellr,area=plotEllipse(a, b, R, lon, lat, h, angle)
                pol.append(ellr)
                el.append(elev)
                if angle>360:
                    angle=angle-360
                az.append(angle)
                larea.append(area)

            c+=1
    
    # Calcule First Fresnel Surface for the only elevation angle
    else:
        elev = elev_list
        df_sp3 = df_sp3.loc[(df_sp3['elevation'] == elev)]
        azim = df_sp3["azimuth"]
        azim = azim.sort_values()
        # Prevent to calculate several time the same azimuth
        list_azim = []
        for val in azim:
            if val not in list_azim:
                list_azim.append(val)
        # Calcule First Fresnel Surface for each elevation angle
        for angle in list_azim:
            a,b,R = firstFresnelZone(L1_GPS, h, elev)
            ellr, area=plotEllipse(a, b, R, lon, lat, h, angle)
            pol.append(ellr)
            el.append(elev)
            if angle>360:
                angle=angle-360
            az.append(angle)
            larea.append(area)
     
    # Create output file with GeoPandas
    d = {'azimuth':az,'elevation':el, 'area':larea, 'geometry':pol}
    gdf = gpd.GeoDataFrame(d)
    gdf.insert(0, 'ID', range(len(gdf)))
    gdf.set_crs(epsg=4326,inplace=True)

    # Check format for output file
    if output=='shp':
        gdf.to_file(f"{dirName}/Fresnel{elev_list}_{h}.shp", driver='ESRI Shapefile')
        if os.path.exists(f'{dirName}/Fresnel{elev_list}_{h}.shp') is True:
            print(f"File Fresnel{elev_list}_{h}.shp successfully created in {dirName}")
        else:
            print("Failed to create file")
    elif output=='gpkg':
        gdf.to_file(f"{dirName}/Fresnel{elev_list}_{h}.gpkg", driver="GPKG")
        if os.path.exists(f'{dirName}/Fresnel{elev_list}_{h}.gpkg') is True:
            print(f"File Fresnel{elev_list}_{h}.gpkg successfully created in {dirName}")
        else:
            print("Failed to create file")
    elif output=='kml':
        gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
        gdf.to_file(f"{dirName}/Fresnel{elev_list}_{h}.kml", driver="KML")
        if os.path.exists(f'{dirName}/Fresnel{elev_list}_{h}.kml') is True:
            print(f"File Fresnel{elev_list}_{h}.kml successfully created in {dirName}")
        else:
            print("Failed to create file")
    else:
        raise Exception("Wrong output format specified")    
        
    return gdf
