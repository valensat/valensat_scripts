
from osgeo import gdal, osr
import netCDF4 as nc
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from sklearn.linear_model import LinearRegression
from scipy import interpolate
from pyproj import Proj

varDir = 'C:\\Curro\\Proyectos\\2021 - Morera\\Python_SSEBI_Var\\'
varName = ['ALBEDO', 'EMISS', 'LST', 'NDVI']

era5path = 'C:\\Curro\\Proyectos\\2021 - Morera\\Python_SSEBI\\'

albedoRange = [0.003, 0.5]
lstRange = [275., 330.]
percLE = [0.0001, 0.0002, 0.0005, 0.001]
percH = [0.001, 0.005, 0.001, 0.015]
percA = 0.999
sigma = 5.670367e-8

SSEBInpts = 20

verbose = True
graph = False # set to True to view lst = f(alb) plot
nptGraph = 100 # number of points for graph plot

class MoreraEvapoTranspiration:

    def __init__(self, varDir, era5path, verbose):
        self.varDir = varDir
        self.era5path = era5path
        self.verbose = verbose

    def run(self):
        overpassTime = self.overpass_hour('RPAS_' + varName[0] + '_140711_1201')
        print('{0} -> Processing Evapotranspiration for {1}...'. format(self.now(), overpassTime))

        NDV, xsize, ysize, GeoT, Projection, DataType = self.GetGeoInfo(self.varDir + 'RPAS_' + varName[0] + '_140711_1201' + '.tif')
        # print(GeoT, Projection)
        # https://gdal.org/tutorials/raster_api_tut.html
        x0 = GeoT[0]
        y0 = GeoT[3]
        xstep = GeoT[1]
        ystep = GeoT[5]

        alb = self.read_variable('RPAS_' + varName[0] + '_140711_1201')
        # print(alb.shape)
        emi = self.read_variable('RPAS_' + varName[1] + '_140711_1201')
        lst = self.read_variable('RPAS_' + varName[2] + '_140711_1201')
        ndvi = self.read_variable('RPAS_' + varName[3] + '_140711_1201')
        cloudMask = self.mask_clouds(alb, lst)
        Rup, Rsol, Ldown, time, lonsurf, latsurf = self.retrieve_rai(era5path, overpassTime, Projection)
        LEmod, Hmod = self.SSEBIparams(alb, lst, cloudMask, percH, percLE, percA, SSEBInpts, graph)
        ETd, Lfnet, Rnet, G, ETi, Hi = self.estimate_ET(alb, lst, ndvi, emi, cloudMask, Rup, Rsol, Ldown, time, lonsurf, latsurf, LEmod, Hmod, overpassTime, xstep, ystep, x0, y0)
        ETdName = self.CreateGeoTiff(self.varDir+'RPAS_ETd_140711_1201', ETd, gdal.GetDriverByName('GTiff'), 0.,
                      xsize, ysize, GeoT, Projection, DataType)
#        LfnetdName = self.CreateGeoTiff(self.varDir+'RPAS_Lfnet_140711_1201', Lfnet, gdal.GetDriverByName('GTiff'), 0.,
#                      xsize, ysize, GeoT, Projection, DataType)
#        RnetName = self.CreateGeoTiff(self.varDir+'RPAS_Rnet_140711_1201', Rnet, gdal.GetDriverByName('GTiff'), 0.,
#                      xsize, ysize, GeoT, Projection, DataType)
#        GName = self.CreateGeoTiff(self.varDir+'RPAS_G_140711_1201', G, gdal.GetDriverByName('GTiff'), 0.,
#                      xsize, ysize, GeoT, Projection, DataType)
#        ETiName = self.CreateGeoTiff(self.varDir+'RPAS_ETi_140711_1201', ETi, gdal.GetDriverByName('GTiff'), 0.,
#                      xsize, ysize, GeoT, Projection, DataType)
#        HiName = self.CreateGeoTiff(self.varDir+'RPAS_Hi_140711_1201', Hi, gdal.GetDriverByName('GTiff'), 0.,
#                      xsize, ysize, GeoT, Projection, DataType)
        if self.verbose: print('     {0} -> Evapotranspiration estimation competed.'.format(self.now()))

    def now(self):
        now = datetime.datetime.now()
        return now.strftime("%Y-%m-%d %H:%M:%S")

    def read_variable(self, varNameV):
        if self.verbose: print('     {0} -> Reading {1}...'.format(self.now(), varNameV))
        data = gdal.Open(r''+self.varDir+varNameV+'.tif')
        array = data.ReadAsArray()
        array[array < 0] = 0.
        data = None
        varNameV = None
        return np.float64(array)

    def GetGeoInfo(self, FileName):
        SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
        NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
        xsize = SourceDS.RasterXSize
        ysize = SourceDS.RasterYSize
        GeoT = SourceDS.GetGeoTransform()
        Projection = osr.SpatialReference()
        Projection.ImportFromWkt(SourceDS.GetProjectionRef())
        DataType = SourceDS.GetRasterBand(1).DataType
        DataType = gdal.GetDataTypeName(DataType)
        return NDV, xsize, ysize, GeoT, Projection, DataType

    def CreateGeoTiff(self, Name, Array, driver, NDV,
                      xsize, ysize, GeoT, Projection, DataType):
        if self.verbose: print('     {0} -> Saving {1} ...'.format(self.now(), Name+'.tif'))
        if DataType == 'Float32':
            DataType = gdal.GDT_Float32
        NewFileName = Name + '.tif'
        # Set nans to the original No Data Value
        Array[np.isnan(Array)] = NDV
        # Set up the dataset
        DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
        # the '1' is for band 1.
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projection.ExportToWkt())
        # Write the array
        DataSet.GetRasterBand(1).SetNoDataValue(NDV)
        DataSet.GetRasterBand(1).WriteArray(Array)
        return NewFileName

    def mask_clouds(self, alb, lst):
        if self.verbose: print('     {0} -> Creating cloud masks ...'.format(self.now()))
        mask1 = (alb >= albedoRange[0]) & (alb <= albedoRange[1]) & (lst >= lstRange[0]) & (lst <= lstRange[1])
        mask2 = ((alb <=.2) & (lst >= lstRange[0])) | ((alb >=.2) & (lst >= lstRange[0]+6.))
        mask = mask1 & mask2
        mask1 = None
        mask2 = None
        return mask

    def overpass_hour(self, varStr):
        filds = varStr.split('_')
        overpass = datetime.datetime.strptime('20'+'_'.join(filds[-2:]), '%Y%m%d_%H%M')
        return overpass

    def convert_accumulated_to_instantaneous(self, var):
        var2 = np.zeros((48, 2, 2), dtype=float)
        var2[1, :, :] = var[1, :, :]
        var2[2:24, :, :] = var[2:24, :, :] - var[1:23, :, :]
        var2[25, :, :] = var[25, :, :]
        var2[26:47, :, :] = var[26:47, :, :] - var[25:46, :, :]
        return var2

    def retrieve_rai(self, era5path, overpass, Projection):
        if self.verbose: print('     {0} -> Retrieving radiation from ERA5...'.format(self.now()))
        dayStr = overpass.strftime('%Y%m%d')[2:]
        """
        dataset = gdal.Open("NETCDF:{0}".format(era5path+dayStr+'.nc'))
        gdal.Info(dataset)
        1: fal: [48x2x2] fal (16-bit integer)
        2: slfh: [48x2x2] surface_upward_latent_heat_flux (16-bit integer)
        3: ssr: [48x2x2] surface_net_downward_shortwave_flux (16-bit integer)
        4: str: [48x2x2] surface_net_upward_longwave_flux (16-bit integer)
        5: sshf: [48x2x2] surface_upward_sensible_heat_flux (16-bit integer)
        6: ssrd: [48x2x2] surface_downwelling_shortwave_flux_in_air (16-bit integer)
        7: strd: [48x2x2] strd (16-bit integer)
        dimensions: 'longitude' (2), 'latitude' (2), 'time' (48)
        """

        # https://towardsdatascience.com/read-netcdf-data-with-python-901f7ff61648
        ds = nc.Dataset(era5path+dayStr+'.nc')
        time = ds['time']
        latsurf = ds['latitude']
        lonsurf = ds['longitude']
        Ru = ds['str']
        Rso = ds['ssrd']
        Ldow = ds['strd']

        myProj = Proj(Projection.ExportToProj4())
        UTMx, UTMy = myProj(lonsurf, latsurf)

        Rup0 = (self.convert_accumulated_to_instantaneous(Ru))/3600.
        Rsol = (self.convert_accumulated_to_instantaneous(Rso))/3600.
        Ldown = (self.convert_accumulated_to_instantaneous(Ldow))/3600.
        Rup = Ldown - Rup0

        return Rup, Rsol, Ldown, time, UTMx, UTMy

    def SSEBIparams(self, alb, lst, mask, percH, percLE, percA, SSEBInpts, graph):
        if self.verbose: print('     {0} -> Estimating S-SEBI parameters...'.format(self.now()))
        npx = len(lst[mask])
        vLEalb = np.zeros((SSEBInpts, len(percH)), dtype=float)
        vHalb = np.zeros((SSEBInpts, len(percH)), dtype=float)
        vLElst = np.zeros((SSEBInpts, len(percH)), dtype=float)
        vHlst = np.zeros((SSEBInpts, len(percH)), dtype=float)

        sortedLst = np.sort(lst[mask])
        for p in range(len(percH)):
            minLst = sortedLst[round(percLE[p]*npx)]
            maxLst = sortedLst[round((1.-percH[p])*npx)]

            lowAlb = np.sort(alb[(np.abs(lst-minLst) < 1.) & mask])
            hiAlb = np.sort(alb[(np.abs(lst-maxLst) < 1.) & mask])
            sortedAlb = np.sort(alb[mask])

            minLEalb = lowAlb[round(percLE[p]*len(lowAlb))]
            meanAlb = np.mean(alb[mask])
            maxAlb = sortedAlb[round(percA*npx)]

            stepLE = (maxAlb-minLEalb)/SSEBInpts
            vLEalbP = np.arange(minLEalb, maxAlb, stepLE) + stepLE
            stepH = (maxAlb-meanAlb)/SSEBInpts
            vHalbP = np.arange(meanAlb, maxAlb, stepH) + stepH
            vLEalb[:,p] = vLEalbP
            vHalb[:,p] = vHalbP

            for a in range(SSEBInpts-1):
                lstLE = np.sort(lst[(np.abs(alb-vLEalbP[a]) <.001) & mask])
                lstH = np.sort(lst[(np.abs(alb-vHalbP[a]) <.001) & mask])
                vLElst[a, p] = lstLE[min(round(len(lstLE)*percLE[p]), len(lstLE)-1)]
                vHlst[a, p] = lstH[min(round(len(lstH)*(1.-percH[p])), len(lstH)-1)]
            vLElst[SSEBInpts-1, p] = np.mean(lst[(np.abs(alb - vLEalbP[SSEBInpts-1]) < .001) & mask])
            vHlst[SSEBInpts - 1, p] = vLElst[SSEBInpts-1, p]

        albVal = np.arange(0, np.max(alb[mask]), .01).reshape((-1, 1))

        # https://realpython.com/linear-regression-in-python/"
        Hmod = np.zeros((2, len(percH)), dtype=float)
        LEmod = np.zeros((2, len(percH)), dtype=float)
        for p in range(len(percH)):
            modelLE = LinearRegression().fit(vLEalb[:,p].reshape((-1, 1)), vLElst[:,p].reshape((SSEBInpts)))
            modelH = LinearRegression().fit(vHalb[:,p].reshape((-1, 1)), vHlst[:,p].reshape((SSEBInpts)))
            LEmod[:,p] = [modelLE.coef_[0], modelLE.intercept_]
            Hmod[:,p] = [modelH.coef_[0], modelH.intercept_]

        pLE = np.argmax(LEmod[0,:])
        pH = np.argmin(Hmod[0,:])
        modelLE = LinearRegression().fit(vLEalb[:, pLE].reshape((-1, 1)), vLElst[:, pLE].reshape((SSEBInpts)))
        modelH = LinearRegression().fit(vHalb[:, pH].reshape((-1, 1)), vHlst[:, pH].reshape((SSEBInpts)))
        predLE = modelLE.predict(albVal)
        predH = modelH.predict(albVal)
        LE = modelLE.predict(meanAlb.reshape((-1, 1)))
        H = modelH.predict(meanAlb.reshape((-1, 1)))

        if graph:
            if self.verbose: print('     {0} -> Preparing S-SEBI graph...'.format(self.now()))
            #prepare density plot
            x = alb[mask]
            y = lst[mask]
            fig, ax = plt.subplots()
            data, x_e, y_e = np.histogram2d(x, y, bins=[nptGraph,nptGraph], density=True)
            z = interpn((0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])), data, np.vstack([x, y]).T,
                        method="splinef2d", bounds_error=False)
            z[np.where(np.isnan(z))] = 0.0
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]

            # normalize color density
            norm = Normalize(vmin=np.min(z), vmax=np.max(z))
            cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='jet'), ax=ax)
            cbar.ax.set_ylabel('Density')

            # actual plot
            plt.scatter(x, y, 1, c=z, cmap='jet')
            """for i in range(len(percH)):
                plt.plot(vLEalb[:,i], vLElst[:,i], 'go')
                plt.plot(vHalb[:,i], vHlst[:,i], 'rx')"""
            plt.plot(albVal, predLE, 'g-')
            plt.plot(albVal, predH, 'r-')
            plt.plot((meanAlb, meanAlb), (H, LE), 'k-')
            plt.xlabel('Albedo')
            plt.ylabel('LST (in K)')
            plt.savefig(varDir + "\\LST-ALB-graph.jpg", dpi=300, transparent=True, bbox_inches='tight')
            plt.show()
            plt.close()


        return LEmod[:,pLE].reshape(2), Hmod[:,pH].reshape(2)

    def interpolate_image(self, xin, yin, zin, xstep, ystep, x0, y0, dimIm):
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html
        f = interpolate.interp2d(xin, yin, zin[:, ::-1], kind='linear')
        xout = np.arange(x0, x0+dimIm[1]*xstep, xstep)
        yout = np.arange(y0+dimIm[0]*ystep, y0, np.abs(ystep))

        zout = f(xout, yout)

        return zout

    def estimate_ET(self, alb, lst, ndvi, emi, cloudMask, Rup, Rsol, Ldown, time, lonsurf, latsurf, LEmod, Hmod,
                    overpass, xstep, ystep, x0, y0):
        if self.verbose: print('     {0} -> Estimating daily evapotranspiration ...'.format(self.now()))
        Lup = sigma * (lst**4.)
        dimIm = lst.shape

        timeg = []
        timed = []
        for t in time[:]:
            timeg.append(datetime.datetime(1900, 1, 1, 0, 0, 0) + datetime.timedelta(seconds=t*3600.))
            timed.append((datetime.datetime(1900, 1, 1, 0, 0, 0) + datetime.timedelta(seconds=t*3600.) - overpass).total_seconds())

        t0 = np.argmax([timed[i] for i in range(len(timed)) if timed[i] < 0])
        t1 = t0 + 1

        dayStart = 1 if t0 <= 24 else 25
        Rsold = np.mean(Rsol[dayStart:dayStart+23, :, :], axis=0)
        Ldownd = np.mean(Ldown[dayStart:dayStart+23, :, :], axis=0)

        Rupb = Rup[t0]+(overpass-timeg[t0]).total_seconds()*(Rup[t1]-Rup[t0])/(timeg[t1]-timeg[t0]).total_seconds()
        Rsolb = Rsol[t0]+(overpass-timeg[t0]).total_seconds()*(Rsol[t1]-Rsol[t0])/(timeg[t1]-timeg[t0]).total_seconds()
        Ldownb = Ldown[t0]+(overpass-timeg[t0]).total_seconds()*(Ldown[t1]-Ldown[t0])/(timeg[t1]-timeg[t0]).total_seconds()

        RsoldIm = self.interpolate_image(lonsurf, latsurf, Rsold, xstep, ystep, x0, y0, dimIm)
        LdowndIm = self.interpolate_image(lonsurf, latsurf, Ldownd, xstep, ystep, x0, y0, dimIm)
        RupbIm = self.interpolate_image(lonsurf, latsurf, Rupb, xstep, ystep, x0, y0, dimIm)
        RsolbIm = self.interpolate_image(lonsurf, latsurf, Rsolb, xstep, ystep, x0, y0, dimIm)
        LdownbIm = self.interpolate_image(lonsurf, latsurf, Ldownb, xstep, ystep, x0, y0, dimIm)

        Cdup = Lup/RupbIm
        RupdIm = np.zeros(dimIm, dtype=float)
        for i in range(dayStart, dayStart+23):
            Rupd = Rup[i, :, :]
            RupIm = self.interpolate_image(lonsurf, latsurf, Rupd, xstep, ystep, x0, y0, dimIm)
            RupIm = RupIm*Cdup
            RupdIm += RupIm
        RupdIm = RupdIm/23.

        Rnet = np.multiply((1. - alb), RsolbIm) - np.multiply(Lup, emi) + np.multiply(LdownbIm, emi)
        G = np.multiply(np.multiply((lst - 273.15), (0.0038 + 0.0074 * alb)), np.multiply((1. - 0.98 * (ndvi**4.)), Rnet))
        Lfraccion = np.divide(Hmod[0] * alb + Hmod[1] - lst, Hmod[0] * alb + Hmod[1] - (LEmod[0] * alb + LEmod[1]))
        Lfraccion[Lfraccion > 1.] = 1.
        Lfraccion[Lfraccion < 0.] = 0.
        Evapo = np.multiply(Lfraccion, (Rnet - G))
        Hlat = Rnet - G - Evapo

        bad = Evapo > 900 | np.isnan(Evapo) | np.isinf(Evapo)
        Rnet[bad] = 0.
        G[bad] = 0.
        Evapo[bad] = 0.
        Hlat[bad] = 0.
        Lfraccion[bad] = 0.

        Rnetd = np.multiply((1-alb), RsoldIm) - np.multiply(RupdIm, emi) + np.multiply(LdowndIm, emi)
        Cdi = np.divide(Rnetd, Rnet)
        Evapod = np.multiply(np.multiply(Cdi, Lfraccion), Rnet)*0.035265
        Evapod[bad] = 0.

        fig = plt.figure()
        plt.title('ETd')
        im = plt.imshow(Evapod)
        pos = fig.add_axes([0.85, 0.1, 0.02, 0.78])
        fig.colorbar(im, cax=pos)
        plt.show()
        fig.savefig(self.varDir + "\\ETd.jpg", dpi=300, transparent=True, bbox_inches='tight')
        plt.close

        return Evapod, Lfraccion, Rnet, G, Evapo, Hlat


test = MoreraEvapoTranspiration(varDir, era5path, verbose)
test.run()