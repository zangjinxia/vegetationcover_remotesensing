'''
实现GF1/2/6和HY1C数据的植被覆盖度的计算
基本流程为根据shp进行海陆掩模-去除异常值-NDVI计算-频率分别为5%和95%的NDVI值分别为NDVI_soil和NDVI_veg
-利用植被覆盖度计算公式进行计算
本版本为去除裁剪后的版本，需要利用shp文件进行海陆掩模，在浙江二期项目中使用的版本
'''
from osgeo import ogr,gdal,gdalconst
import numpy as np
import scipy.stats as st
import gdalconst
import os



def read_img(filename):
    dataset = gdal.Open(filename)
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    band = dataset.RasterCount
    im_data = dataset.ReadAsArray(0, 0, width, height)

    geotrans = dataset.GetGeoTransform()
    proj = dataset.GetProjection()
    # data = np.zeros([width, height, band])

    return im_data, proj, geotrans, band, width, height


def write_tiff(filename, proj, geotrans, data):
    # gdal数据类型包括
    # gdal.GDT_Byte,
    # gdal .GDT_UInt16, gdal.GDT_Int16, gdal.GDT_UInt32, gdal.GDT_Int32,
    # gdal.GDT_Float32, gdal.GDT_Float64
    # 判断栅格数据的数据类型
    if 'int8' in data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    # 判读数组维数
    if len(data.shape) == 3:
        bands, height, width = data.shape
    else:
        bands = 1
        height, width = data.shape
    # 创建文件
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(filename, width, height, bands, datatype)

    dataset.SetGeoTransform(geotrans)
    dataset.SetProjection(proj)

    if bands == 1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for i in range(bands):
            dataset.GetRasterBand(i + 1).WriteArray(data[i])
    del dataset

def NDVI(B1,B2):
    """
    计算NDVI
    :param B1: 近红波段
    :param B2: 红波段
    :return: NDVI矩阵
    """
    result1 = (B1 - B2) / (B1 + B2)
    return result1

def shp2Raster(shp, templatePic, output):
    """
    shp:字符串，一个矢量，从0开始计数，整数
    templatePic:字符串，模板栅格，一个tif，地理变换信息从这里读，栅格大小与该栅格一致
    output:字符串，输出栅格，一个tif
    field:字符串，栅格值的字段
    nodata:整型或浮点型，矢量空白区转换后的值
    """
    # 读取栅格模板的大小形状和投影
    ndsm = templatePic
    data = gdal.Open(ndsm, gdalconst.GA_ReadOnly)
    geo_transform = data.GetGeoTransform()
    proj = data.GetProjection()
    # source_layer = data.GetLayer()
    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_max = x_min + geo_transform[1] * data.RasterXSize
    y_min = y_max + geo_transform[5] * data.RasterYSize
    x_res = data.RasterXSize
    y_res = data.RasterYSize

    # 打开shp文件，读取信息
    mb_v = ogr.Open(shp)
    mb_l = mb_v.GetLayer()
    pixel_width = geo_transform[1]

    # 输出影像为16位整型
    target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 1, gdal.GDT_Int16)

    target_ds.SetGeoTransform(geo_transform)  # %设置输出文件的仿射变换六参数
    target_ds.SetProjection(proj)
    band = target_ds.GetRasterBand(1)
    NoData_value = 0
    band.SetNoDataValue(NoData_value)
    band.FlushCache()

    gdal.RasterizeLayer(target_ds, [1], mb_l, burn_values=[1],
                        options=['ALL_TOUCHED=TRUE'])  # ATTRIBUTE：指定矢量的字段值写入栅格数据中
    # 输入的栅格数据，栅格数据要更新的波段列表，矢量图层数组，输出栅格图像的值，ALL_TOUCHED=TRUE表示将所有矢量转为栅格，否则只是矢量的中心某区域转化
    target_ds = None


if __name__ == '__main__':
    #输入参数，包括掩模shp文件，需计算植被覆盖度的文件，植被覆盖度生成图
    shp = 'D:/AAdata/testdata/shp/islandshp/GF63704_islandpoly1.shp'
    templatePic = 'D:/AAdata/testdata/GF6_PMS_E121.4_N27.6_20200826_L1A1120030364/GF6_PMS_E121.4_N27.6_20200826_L1A1120030364-MUX-orthp-pipei8_caijian.tiff'
    # output = 'D:/AAdata/testData/mask.tif'
    output = 'D:/AAdata/testData/vegetationCover10.tif'

    #生成植被覆盖度的路径和文件名
    mask = os.path.dirname(output)+'/'+'mask.tif'

    #制作掩模文件
    shp2Raster(shp,templatePic,mask)

    # 利用掩模文件进行掩模
    im_data, proj, geotrans, band, width, height = read_img(templatePic)
    templatePic_data, tp_proj, tp_geo, tp_band, tp_width, tp_height = read_img(mask)
    img = np.array(im_data)
    tem = np.array(templatePic_data)
    out_data = np.empty(shape=[band, height, width])
    for i in range(band):
        out_data[i] = img[i] * tem  # 进行掩模
    print("掩模完成")

    #计算NDVI
    ndvi = NDVI(out_data[3],out_data[2])
    # print(ndvi.shape)

    #计算植被覆盖度
    mean = np.nanmean(ndvi)
    std = np.nanstd(ndvi)
    conf_intveral = st.norm.interval(0.95, loc=mean, scale=std)
    # print(conf_intveral)


    for i in range(height):
        for j in range(width):
            if ndvi[i,j] < conf_intveral[0]:
                ndvi[i,j] = 0
            elif (ndvi[i,j]>=conf_intveral[0])&(ndvi[i,j]<=conf_intveral[1]):
                ndvi[i,j] = (ndvi[i,j] - conf_intveral[0]) / (conf_intveral[1] - conf_intveral[0])

    write_tiff(output, proj, geotrans, ndvi)

    #删除生成的掩模文件mask
    os.chdir(os.path.dirname(output))
    os.remove('mask.tif')


    print('植被覆盖度反演完成！')




