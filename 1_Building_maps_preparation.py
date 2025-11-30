# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:36:40 2022

1. Calculate building footprints 

@author: Janneke van Oorschot, Mike Slootweg

这份脚本的主要作用：
- 以荷兰各个市镇（municipality）为单位
- 把 BAG（建筑基础数据）和 PBL 场景（建设/拆除）数据
  叠加到 100m fishnet 网格上
- 计算每个网格内建筑占地面积（Ground Area）
- 按建筑占地规模分级输出（0–1000, 1000–4000, >4000 m2）
"""
#%% 
print("1. Package installation & import")

# ArcGIS 的 Python 库，用于各种 GIS 操作（矢量、栅格、空间分析等
import arcpy, arcinfo
# pandas / numpy 用于表格和数值运算
import pandas as pd
import numpy as np
# glob / os 用于文件系统操作（遍历文件、路径等）
import glob, os

from arcpy import env
from arcpy.sa import * # Spatial Analyst 模块（栅格计算等）
from osgeo import gdal  # GDAL，用于栅格处理
from os.path import exists

# geopandas 用于操作 shapefile 等矢量数据（与 pandas 很像）
import geopandas as gpd
import pandas as pd

# 允许覆盖已有输出
arcpy.env.overwriteOutput = True
#Enable ArcGIS extensions# 启用 ArcGIS 的 Spatial Analyst 扩展
arcpy.CheckOutExtension("Spatial")

#%% 
print("2. Load data and parameters")

# Set directory # 设置主数据目录（所有路径基于这个）
directory = "E:/Paper IV/Data/"


# ============================================================
# 2.1 读取市镇边界，并清洗字段名称
# ============================================================
# Creating a list of municipality names --------------------------------------------------------------------------------
# 读取 CBS 的市镇 shapefile
muniii = gpd.read_file(directory + "GIS_data/CBS_gemeente/CBS_gemeente.shp")
# 对 'statnaam' 字段做字符替换，避免后面生成文件名时出问题
muniii['statnaam'] = muniii['statnaam'].str.replace('-', '_')
muniii['statnaam'] = muniii['statnaam'].str.replace('(', '_')
muniii['statnaam'] = muniii['statnaam'].str.replace(')', '_')
muniii['statnaam'] = muniii['statnaam'].str.replace('"', '_')
muniii['statnaam'] = muniii['statnaam'].str.replace('.', '_')
muniii['statnaam'] = muniii['statnaam'].str.replace(',', '_')
# 手动修正几个特殊市镇名称（荷兰城市名）
muniii.at[292,'statnaam'] = '_s_Gravenhage'
muniii.at[187,'statnaam'] = '_s_Hertogenbosch'
muniii.at[243,'statnaam'] = 'Bergen _L_'
muniii.at[188,'statnaam'] = 'Bergen _NH_'

# 保存清洗后的市镇 shapefile，为后续 ArcPy 使用
muniii.to_file(directory+'GIS_data/CBS_Gemeente/Municipalities.shp')

# 创建一个排序后的市镇名称列表（用于循环）
munlist_unsorted = list(muniii['statnaam'])
munlist = sorted(munlist_unsorted)

# 再拷贝一份，专门把空格替换为下划线，用于生成中间输出名
muniii_copy = muniii.copy()
muniii_copy['statnaam'] = muniii_copy['statnaam'].str.replace(' ', '_')
munlist_copy_unsorted = list(muniii_copy['statnaam'])
munlist_adj = sorted(munlist_copy_unsorted)

# 如果模型中途在某个市镇停了，可以通过切片从中间继续跑 # Model gestopt bij Zeist, bij hervatten beginnen bij. 
#munlist = munlist[333:]
#munlist_adj = munlist_adj[333:]

# ============================================================
# 2.2 基础 GIS 数据路径（fishnet 网格 & 市镇面）
# ============================================================
# 100m fishnet 网格（全国统一网格）
fishnet = directory + "GIS_data/Fishnet/Fishnet_100m.shp"
# 市镇多边形（上一步保存的）
municipality = directory + "GIS_data/CBS_gemeente/Municipalities.shp"


# ============================================================
# 2.3 读取 Excel 参数：不同建筑类型的占地系数等
# ============================================================
#read csv that contains the factors for the scenarios
#factors = pd.read_excel(directory + "Documents/Data/Grondoppervlak_scenarios.xlsx", sheet_name = 'simplified')
factors = pd.read_excel(directory+"Excel_data/BGA.xlsx", sheet_name = 'Totals')

# 把某一列设为 index，方便后面用 .loc 按行名取值
factors = factors.set_index('Unnamed: 0')
# header = pd.MultiIndex.from_product([['low_rise', 'mid_rise', 'high_rise'], ['BAU','small','large']], names = ['height','UFA'])

# ============================================================
# 2.4 预设场景参数
# ============================================================
# Parameters -------------------------------------------------------------------------------------------------
# WLO 场景（这里只用 'hoog'，即高情景）
WLO = ['hoog']
# PBL 情景：靠近（DichtBij）、宽松（Ruim）
scenario = ['_DichtBij', '_Ruim']
# 活动类型：建设、拆除（暂时没直接用到）
activiteit = ['Bouw', 'Sloop']

# ============================================================
# 2.5 建筑类型列表（荷兰 BAG / PBL 分类）
# ============================================================
# Classifciations used in PBL + BAG analysis
#gebouwtypen = ['appartement', 'Detailhandel', 'Kantoor', 'losstaand', 'NijverheidEnLogistiek', 'Overheid_kw_diensten',  'rijtjeswoning']
# BAG 
# BAG 建筑类型
gebouwtypen = ['appartement', 'bedrijfshal', 'distributie', 'kantoor_groot', 'kantoor_klein', 'onderwijs', 'serieel', 'vrijstaand', 'winkel', 'woonflat']


# ============================================================
# 2.6 UFA（使用面积）情景和建筑高度情景
# ============================================================
# UFA scenarios ---------------------------------------------------------------------------------------------------
# Ground area per residential building type (non-residential building data is already provided in m2 footprint) 
# Ground area per residential building type (住宅的建筑占地系数)
# 非住宅建筑（比如工业、办公）占地信息已经直接给到 m² footprint
UFA_scenarios = ['BAU']#,'small', 'large']
# 建筑高度情景，目前只算 'low_rise'，可以扩展 'mid_rise', 'high_rise'
height_scenarios = ['low_rise']#,'high_rise']# 'mid_rise', 'high_rise'] # only applies for apartments & offices# 只对公寓和办公有效
        
#%% #%% 主循环开始
print('2. Loop trough each municipality, and add BAG & PBL info to fishnet')      

# 这里通过 zip 把原始市镇名（munlist）和替换空格后的（munlist_adj）配对for m, m_adj in zip(munlist, munlist_adj):
    print('Step 2, starting with municipality:')# 当前处理哪个市镇
    print(m)
    municipality_name = m
     # 当前市镇对应的 BAG 数据（已经预处理过，按市镇分好）
    BAG_municipality = directory + "GIS_data/GIS_building_data_output/BAG_per_municipality_reclassified_materialized/" + municipality_name + "_voorraden"
    
    # Define the municipality and make selection according to given municipality name
    # after this, geodatabase files are created using the name of the municipality.
    # 'inter_gdb' is used for intermediate data and steps
    # 'results_gdb' is used for final datasets
    # ============================================================
    # 2.x 为当前市镇创建中间 & 结果 geodatabase
    # ============================================================

    # 用于 Select 的查询语句（按市镇名称筛选）
    search_name = "statnaam = " + "'" + municipality_name + "'" 
    # The exists tool checks if the directory exists. 
    # With the if-statement we check whether directory,if not, it creates a new one. # 判断中间 geodatabase 是否存在，不存在则创建
    inter_gdb_exists = exists(directory +"GIS_data/footprint_municipality/" + municipality_name + "_intermediate.gdb")
    results_gdb_exists = exists(directory +"GIS_data/footprint_municipality/"+ municipality_name + "_results.gdb")
    
    if inter_gdb_exists == False:
        arcpy.CreateFileGDB_management(directory + "GIS_data/footprint_municipality/", municipality_name + "_intermediate.gdb")
    
    if results_gdb_exists == False:
        arcpy.CreateFileGDB_management(directory + "GIS_data/footprint_municipality/", municipality_name + "_results.gdb")
    # 中间数据 gdb 和结果 gdb 的路径（字符串末尾加 / 便于后面拼路径）
    # Creating geodatabase files for intermediate output and for final results, for municipality x
    inter_gdb = directory + "GIS_data/footprint_municipality/" + municipality_name + "_intermediate.gdb/"
    results_gdb = directory + "GIS_data/footprint_municipality/" + municipality_name + "_results.gdb/"
    
    # 为当前市镇创建一个 shapefile 输出文件夹（存最终 shp） # Folder for shapefiles
    arcpy.management.CreateFolder(directory+'GIS_data/footprint_municipality/', municipality_name)
    
    # ============================================================
    # 2.x 从全国市镇图层中，选出当前市镇的 polygon
    # ============================================================ # Search municipality field according to the given municipality name    
    municipality_selection = arcpy.analysis.Select(municipality, directory+ 'GIS_data/footprint_municipality/'+municipality_name+'_intermediate.gdb/'+ m_adj+'selection' , "statnaam = '{}'".format(municipality_name))
                                                      
    # Make fishnet selection of the new municipality shapefile
    # Clipping the fishnet cuts the squares according to the shape of the municipality
    # Selecting by location is applied to keep all the squares fully around the edges.
    fishnet_selection = arcpy.management.SelectLayerByLocation(fishnet, "HAVE_THEIR_CENTER_IN", municipality_selection,None, "NEW_SELECTION", "NOT_INVERT")
    arcpy.CopyFeatures_management(fishnet_selection, inter_gdb + "fishnet_selection")
     # ============================================================
    # 2.x 用当前市镇 polygon 去选取 fishnet 网格
    # ============================================================ # 选择 fishnet 中中心点落在市镇范围内的格子
    fishnet_BAG = arcpy.CopyFeatures_management(fishnet_selection, inter_gdb + "fishnet_BAG")

    # Prepare BAG dataset of municipality 
    arcpy.analysis.TabulateIntersection(fishnet_BAG, 'OBJECTID' , BAG_municipality, inter_gdb + "BAG_zone_selection", "Class")
    arcpy.management.PivotTable(inter_gdb + "BAG_zone_selection", "OBJECTID_1", "CLass", "AREA", inter_gdb + "BAG")  
    
    print('all parameters are set')

# -----------------------------------------------------------------------------------------------------------------------------------         

    print("3. Data editing and sorting: PBL scenario data to fishnet, spatially joined per scenario and spatial join with BAG")
    
    targetfeatures = fishnet_selection
    
    for i in range(len(WLO)):
        print('start with: ')
        for j in range(len(scenario)):
            os.chdir(directory+"/GIS_data/PBL_data/" + WLO[i] + scenario[j])
            targetfeatures = fishnet_selection
            for file1 in glob.glob("*.tif"):
                
                filename = os.path.splitext(file1)[0] 
                print(filename)
                
                #Create abbreviation of filename for field names in spatial join
                name_split = filename.split("_")
                abr = name_split[3:]
                copy_abr = abr[:]
                insert_at = 1
                copy_abr[insert_at:insert_at] = ["_"]
                abr_name = ''.join(copy_abr)
                abr_filename = abr_name[:8]            
                
                #3.1. Define projection ---------------------------------------------------------------------------------------------
                
                arcpy.management.DefineProjection(filename+".tif",'PROJCS["RD_New",GEOGCS["GCS_Amersfoort",DATUM["D_Amersfoort",SPHEROID["Bessel_1841",6377397.155,299.1528128]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Double_Stereographic"],PARAMETER["False_Easting",155000.0],PARAMETER["False_Northing",463000.0],PARAMETER["Central_Meridian",5.38763888888889],PARAMETER["Scale_Factor",0.9999079],PARAMETER["Latitude_Of_Origin",52.15616055555555],UNIT["Meter",1.0]]')
                
                #3.2. Round to nearest integer: raster+0.5 & Int ---------------------------------------------------------------------------------------
                # something is going wrong here, some cells with value 0 get a value of 1, which is strange. 
                # question for Bart:  1% of cells (checked for Amsterdam only) do only have demolition and no construction, I thought the model assumption was only demolition in case of densification (i.e. construction)
                raster_name = directory+"/GIS_data/PBL_data/" + WLO[i] + scenario[j] + "/" + filename + ".tif"
                out_rc = RasterCalculator([raster_name], ["x"], "x+0.5", "FirstOf")
                out_rc_Int = Int(out_rc)
                out_rc_Int.save(directory+"/GIS_data/PBL_data_roundoff/" + filename + '_roundoff.tif')
                
                #3.3. Extract by mask --------------------------------------------------------------------------------------------
                outExtractByMask = ExtractByMask(directory +"/GIS_data/PBL_data_roundoff/"+ filename + "_roundoff.tif",  municipality_selection)
                outExtractByMask.save(inter_gdb + filename + "_selection")
                int_temp = Int(inter_gdb + filename + "_selection")
                
                arcpy.management.Delete(directory + filename + "_roundoff.tif")
                
                #3.4. Raster to Polygon --------------------------------------------------------------------------------------------
                arcpy.RasterToPolygon_conversion(int_temp, inter_gdb + filename + "_selection_pol", "NO_SIMPLIFY", "Value")
                
                #3.5. Spatial join ------------------------------------------------------------------------------------------------- 
                
                joinfeatures = inter_gdb + filename + "_selection_pol"
            
                fieldmappings = arcpy.FieldMappings()
                fieldmappings.addTable(targetfeatures)
                fieldmappings.addTable(joinfeatures)
    
                gridcode_index = fieldmappings.findFieldMapIndex("gridcode")
                fieldmap = fieldmappings.getFieldMap(gridcode_index)
    
                field = fieldmap.outputField
    
                field.name = abr_filename
                field.aliasName = abr_filename
                fieldmap.outputField = field
                fieldmappings.replaceFieldMap(gridcode_index, fieldmap)
    
                spatialjoin_temp = arcpy.SpatialJoin_analysis(targetfeatures, joinfeatures, inter_gdb + filename + "_sj", "#", "#", fieldmappings,"HAVE_THEIR_CENTER_IN")
                targetfeatures = spatialjoin_temp
            
            
            arcpy.CopyFeatures_management(inter_gdb+'WLO_'+WLO[i]+scenario[j]+'_Sloop_rijtjeswoning_sj', inter_gdb+'WLO_'+WLO[i]+scenario[j]+'_PBL_sj')
            
            print('done with SJ of PBL maps for scenario')        
            
            #3.6 join with BAG dataset ----------------------------------------------------------------------------------------------
            
            arcpy.MakeFeatureLayer_management(inter_gdb+'WLO_'+WLO[i]+scenario[j]+'_PBL_sj', inter_gdb + WLO[i] + scenario[j] + "_sj_lyr")
            arcpy.management.AddJoin(inter_gdb + WLO[i] + scenario[j] + "_sj_lyr", "OBJECTID", inter_gdb + "BAG", "OBJECTID_1")
            arcpy.CopyFeatures_management(inter_gdb + WLO[i] + scenario[j] + "_sj_lyr", inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG")
            
            path = inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG"
            fieldObs = arcpy.ListFields(path)  
            fieldNames = []  
            for field in fieldObs:  
                fieldNames.append(field.name)  
            del fieldObs  
            fieldCount = len(fieldNames)  
    
            with arcpy.da.UpdateCursor(path, fieldNames) as curU:  
                for row in curU:  
                    rowU = row  
                    for field in range(fieldCount):  
                        if rowU[field] == None:  
                            rowU[field] = "0"  
          
          
                    curU.updateRow(rowU)
    
            del curU
            print('done SJ with BAG')
            
            #3.7 Delete repititive and useless columns ----------------------------------------------------------
            
            targetfid = []
            joincount = []
            target_fid_list = arcpy.ListFields(inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG", "*TARGET_FID_*")
            join_count_list = arcpy.ListFields(inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG", "*Join_count*")
            
            for field in target_fid_list:
                targetfid.append(field.name)
    
            for field in join_count_list:
                joincount.append(field.name)
            
            if targetfid: 
                arcpy.management.DeleteField(inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG", targetfid, "DELETE_FIELDS")
            if joincount:
                arcpy.management.DeleteField(inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG", joincount, "DELETE_FIELDS")
            
            
            print('completely done with this scenario')
        
        
# -----------------------------------------------------------------------------------------------------------------------------------         
    
    print("4. Calculate constructed ground area per scenario")
    
    # List of field names of which PBL file & BAG prefix need to be removed
    PBL_names = ['Bouw_app','Bouw_Det', 'Bouw_Kan', 'Bouw_los', 'Bouw_Nij', 'Bouw_Ove', 'Bouw_rij', 'Sloop_ap', 'Sloop_De', 'Sloop_Ka', 'Sloop_lo', 'Sloop_Ni', 'Sloop_Ov', 'Sloop_ri']
    BAG_names = [ 'appartement', 'bedrijfshal', 'distributie', 'kantoor_groot', 'kantoor_klein', 'onderwijs', 'serieel', 'vrijstaand', 'winkel', 'woonflat', 'zorg']
    
    Sloop_list = ['Sloop_ap', 'Sloop_De', 'Sloop_Ka', 'Sloop_lo', 'Sloop_Ni', 'Sloop_Ov', 'Sloop_ri']
    Bouw_list = ['det_A_c','row_A_c','ret_A_c', 'ind_A_c', 'off_A_c', 'ap_A_c', 'ser_A_c']
    BAG_names_short = [ 'appartemen', 'bedrijfsha', 'distributi', 'kantoor_gr', 'kantoor_kl', 'onderwijs', 'serieel', 'vrijstaand', 'winkel', 'woonflat', 'zorg']
    
    results_shp_mun = directory+'GIS_data/footprint_municipality/' + municipality_name+ '/'
    
    for i in range(len(WLO)):
        for j in range(len(scenario)):
            for k in range(len(UFA_scenarios)):
                for l in range(len(height_scenarios)):
                    
                    arcpy.CopyFeatures_management(inter_gdb + WLO[i] + scenario[j] + "_sj_with_BAG", inter_gdb + WLO[i] + scenario[j] + "_sj_copy")
                    
                    # Shorten fieldnames (required for feature to shapefile)
                    for name in PBL_names:                      
                        arcpy.management.AlterField(inter_gdb + WLO[i] + scenario[j] + "_sj_copy", 'WLO_'+WLO[i]+scenario[j]+'_PBL_sj_'+name, name, name, "LONG", 4, "NULLABLE", "DO_NOT_CLEAR")                    
                    for name1 in BAG_names:
                        # Shorten BAG names if 'BAG_' is included in field name
                        featureclass = inter_gdb + WLO[i] + scenario[j] + "_sj_copy"
                        field_names = [f.name for f in arcpy.ListFields(featureclass)]
                        list_BAG = []
                        for fn in field_names:
                            if 'BAG' in fn:
                                BG = fn[4:]
                                list_BAG.append(BG)
                        if name1 in list_BAG:
                            arcpy.management.AlterField(inter_gdb + WLO[i] + scenario[j] + "_sj_copy",'BAG_'+name1, name1[:10], '', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
                    
                    # Create copy of fishnet for height & UFA scenario:     
                    arcpy.CopyFeatures_management(inter_gdb + WLO[i] + scenario[j] + "_sj_copy", results_shp_mun + "fishnet_" + WLO[i] + "_" + scenario[j] + "_" + UFA_scenarios[k] + "_"+  height_scenarios[l])
                    
                    # Loading the copy as a geopandas dataframe                
                    df = gpd.read_file(results_shp_mun + "fishnet_" + WLO[i] + "_" + scenario[j] + "_" + UFA_scenarios[k] + "_"+  height_scenarios[l]+'.shp')
                    filename = "fishnet_" + WLO[i] + "_" + scenario[j] + "_" + UFA_scenarios[k] + "_"+  height_scenarios[l]
                    print('Calculating BGA for scenario:')
                    print(filename)
                    
                    # Calculate ground area of buildings according to the factors in Excel file                
                    df['det_A_c'] = df['Bouw_los'] * factors.loc['detachedHouse', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    df['row_A_c'] = df['Bouw_rij'] * factors.loc['rowHouse', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    
                    df['ret_A_c'] = df['Bouw_Det'] * factors.loc['retail', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    df['ind_A_c'] = df['Bouw_Nij'] * factors.loc['industry', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    df['off_A_c'] = df['Bouw_Kan'] * factors.loc['office', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    df['ser_A_c']= df['Bouw_Ove'] * factors.loc['services', [UFA_scenarios[k]+'_'+height_scenarios[l]]]
                    
                    # select appartment building ground area based on number of apartments & scenario
                    s1 = factors.start.values # minimum
                    s2 = factors.stop.values # maximum
                    s = df.Bouw_app.values[:,None]
                    scenario_name = UFA_scenarios[k]+'_'+height_scenarios[l]
                    df['ap_A_c'] = np.dot((s>=s1)&(s<=s2),factors[scenario_name])
                    
                    # Calculate building ground area per cell (2050)
                    # some municipalities do not have all building types (e.g. zorg), so make a list of building types per municipality
                    BAG_names_short_adj = []
                    for fn in list(df):
                        if fn in BAG_names_short:
                            BAG_names_short_adj.append(fn)
                    
                    df['Sl_check'] = df[Sloop_list].sum(axis=1)
                    df['Gr_A'] = np.nan
                    df.loc[df['Sl_check'] != 0, 'Gr_A'] = df[Bouw_list].sum(axis=1) # in case of demolition,  everything will be demolished within gridcell: Gr_A = sum of construction activities
                    df.loc[df['Sl_check']==0, 'Gr_A'] = df[BAG_names_short_adj].sum(axis=1)+df[Bouw_list].sum(axis=1) # in case of no demolition, Gr_A = sum of construction activities + building stock 2018              
                    # df.loc[df['Gr_A']>10000, 'Gr_A']= 10000 # ground area scenarios will in some cells result in a building ground area larger than the cell size, needs to be corrected. 
                    df['activity'] =df[PBL_names].sum(axis=1) # if Activity != 0, demolition or construction
                    
                    # Export file
                    df_small = df[['Gr_A','geometry', 'activity']].copy() # only column relevant for wallpapering
                    
                    df_small.to_file(results_shp_mun+"Gr_A"+filename+".shp")
                    df.to_file(results_shp_mun+"buildings_"+filename+'.shp')
                                  
                    #arcpy.conversion.FeatureClassToFeatureClass(results_shp_mun+"fishnet_buildings_"+filename+'.shp', results_gdb, filename)
                    
                    # split data based on building ground area
                    
# -----------------------------------------------------------------------------------------------------------------------------------         
                
                    print('5. Split files based on building footprint')         
              
                    filename = "fishnet_" + WLO[i] + "_" + scenario[j] + "_" + UFA_scenarios[k] + "_"+  height_scenarios[l]
                    #splitdata = gpd.read_file(results_shp_mun+"buildings_"+filename+'.shp')
                    splitdata = gpd.read_file(results_shp_mun+"Gr_A"+filename+".shp")
                    
                    split_10 = splitdata.loc[(splitdata['Gr_A']>0) & (splitdata['Gr_A']<= 1000) & (splitdata['activity']!=0)]
                    split_40 = splitdata.loc[(splitdata['Gr_A']>1000) & (splitdata['Gr_A']<=4000) & (splitdata['activity']!=0)]
                    split_100 = splitdata.loc[(splitdata['Gr_A']>4000) & (splitdata['activity']!=0)]
                    
                    
                    if split_10.empty == False: 
                        split_10.to_file(results_shp_mun+"buildings_"+filename+'_10.shp')
                    if split_40.empty == False:
                        split_40.to_file(results_shp_mun+"buildings_"+filename+'_40.shp')
                    if split_100.empty == False:
                        split_100.to_file(results_shp_mun+"buildings_"+filename+'_100.shp')
                        
#%%
print('6. Loop trough each municipality, and delete geodatabases')      
for m, m_adj in zip(munlist, munlist_adj):
    print('Step 2, starting with municipality:')
    print(m)
    municipality_name = m
    os.remove(directory +"GIS_data/footprint_municipality/" + municipality_name + "_intermediate.gdb")
    os.remove(directory +"GIS_data/footprint_municipality/" + municipality_name + "_results.gdb")
#%% 
print('7. Merge municipality data, specifically per building footprint classification ')       
directory = "E:/Paper IV/Data/"
classes = ['_10.shp', '_40.shp', '_100.shp']

for c in classes:
    for i in range(len(WLO)):
        for j in range(len(scenario)):
            for k in range(len(UFA_scenarios)):
                for l in range(len(height_scenarios)):
                    filename = "fishnet_" + WLO[i] + "_" + scenario[j] + "_" + UFA_scenarios[k] + "_"+  height_scenarios[l]
                    bf_list = []    
                    for m in munlist:
                        dir_mun = directory+'GIS_data/footprint_municipality/' + m + '/'+"buildings_"+filename+c
                        if exists(dir_mun):
                            bf_list.append(dir_mun)
                    arcpy.management.Merge(bf_list, directory+'GIS_data/footprint_NL/'+filename+c, "", "ADD_SOURCE_INFO")

#%% BAG only (run separetely)
'''
print('7. Merge municipality data, specifically per building footprint classification, BAG ')       

BAG_names_short = [ 'appartemen', 'bedrijfsha', 'distributi', 'kantoor_gr', 'kantoor_kl', 'onderwijs', 'serieel', 'vrijstaand', 'winkel', 'woonflat', 'zorg']
    
          
list_BAG_merge = []
for m in munlist:
    filename = "buildings_fishnet_hoog__DichtBij_BAU_high_rise.shp" #can be any scenario
    featureclass = gpd.read_file(directory+'GIS_data/footprint_municipality/' + m+ '/'+filename)
    feature = directory+'GIS_data/footprint_municipality/' + m+ '/'+filename
    field_names = [f.name for f in arcpy.ListFields(feature)]
    list_BAG = ['geometry']
    for B in BAG_names_short:
        if B in field_names:
            list_BAG.append(B)
    BAG_ftrclss = featureclass[list_BAG].copy()
    BAG_ftrclss.to_file(directory+'GIS_data/footprint_municipality/' + m+ '/'+'BAG_'+filename)
    list_BAG_merge.append(directory+'GIS_data/footprint_municipality/' + m+ '/'+'BAG_'+filename)    
print('Start merging BAG of municipalities')
arcpy.management.Merge(list_BAG_merge, directory+'GIS_data/footprint_municipality/footprint_NL/'+'BAG_'+filename, "", "ADD_SOURCE_INFO")

'''

#%%
# lstFields = arcpy.ListFields(inter_gdb + WLO[i] + scenario[j] + "_sj_copy")

# x = False

# for field in lstFields:
#     if field.name == "serieel":
#         print("Field exists")
#         x = True  
# if x == False:

#     print("Field does not exist")
