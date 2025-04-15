# -*- coding: utf-8 -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import math
import job
import numpy as np
import os
import csv
import re

working_directory=r"D:/TEMP"         
#os.makedirs(working_directory)
os.chdir(working_directory)

def read_data(filename):
    data = {}
    
    with open(filename, 'r') as file:
        for line in file:
            # 去除行首尾的空白字符，并跳过空行或注释行
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # 按等号分割每行数据
            key, value = line.split(':', 1)
            key = key.strip()  # 去除键名两端的空格
            
            value = value.strip()  # 去除值两端的空格
            
            # 检查是否是列表（即含有逗号）
            if ',' in value:
                # 处理为列表类型，使用逗号分割，并转换为整数
                value = [int(v.strip()) for v in value.split(',')]
            else:
                # 处理为单一数值类型
                if '.' in value:
                    value = float(value)  # 转换为浮点数
                else:
                    value = int(value)  # 转换为整数           
            data[key] = value
    return data

# 读取数据

with open(r"D:\TEMP\data.txt", "r") as file:
    data = file.read()

#--------------新建模型-------------------------------
mdb.Model(name='Model-1', modelType=STANDARD_EXPLICIT) 

#创建基准周期球壳
radius=data['radius']
angle_period=data['angle_period']
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, radius), point2=(radius, 0.0), direction=CLOCKWISE)
p = mdb.models['Model-1'].Part(name='Part-shell', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-shell']
p.BaseShellRevolve(sketch=s, angle=angle_period, flipRevolveDirection=OFF)
del s

# 主级纬向筋条划分，在网格宽度减半时经向筋减半，分为三个区域
longitude_hh_main=data['longitude_hh']
longitude_number_1= int(720.0 / angle_period )
longitude_number=int(longitude_number_1 * 1.5)
longitude_number_main= int(np.ceil((longitude_number - 2) / 4.0 ))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.Line(point1=(radius, 0.0), point2=(radius+longitude_hh_main, 0.0))
s.radialPattern(geomList=(g[3], ), vertexList=(), number=longitude_number_main, totalAngle=90.0/longitude_number*(longitude_number_main - 1) * 4.0 , centerPoint=(0.0, 0.0))#阵列主级筋条
s.radialPattern(geomList=(g[3], ), vertexList=(), number=2, totalAngle=90.0/longitude_number*(longitude_number - 2)  , centerPoint=(0.0, 0.0))#陈列顶部主级筋条
p = mdb.models['Model-1'].Part(name='Part-longitude-main', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-longitude-main']
p.BaseShellRevolve(sketch=s, angle=angle_period, flipRevolveDirection=OFF)
del s

# 次级纬向筋条划分，在网格宽度减半时经向筋减半，分为三个区域
longitude_hh_secondary=longitude_hh_main / 2
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.Line(point1=(radius, 0.0), point2=(radius+longitude_hh_secondary, 0.0))
s.radialPattern(geomList=(g[3], ), vertexList=(), number=longitude_number - 1, totalAngle=90.0/longitude_number*(longitude_number - 2), centerPoint=(0.0, 0.0))
start_idx = 3  # 主级筋条的起始编号
step = 4  # 每次跳过的索引间隔
delete_indices = [start_idx + i * step for i in range(longitude_number_main)]# 生成待删除的筋条的索引列表
s.delete(objectList=[g[i] for i in delete_indices])# 删除对应的主级筋条
s.delete(objectList=(g[start_idx + longitude_number - 2],))# 删除顶部最后一条的主级筋条
p = mdb.models['Model-1'].Part(name='Part-longitude-secondary', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-longitude-secondary']
p.BaseShellRevolve(sketch=s, angle=angle_period, flipRevolveDirection=OFF)
del s


# 主级经向筋条
latitude_hh_main=longitude_hh_main
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.Line(point1=(radius, 0.0), point2=(radius+latitude_hh_main, 0.0))
p = mdb.models['Model-1'].Part(name='Part-latitude-main', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-latitude-main']
p.BaseShellRevolve(sketch=s, angle=90.0/longitude_number*(longitude_number - 2), flipRevolveDirection=OFF)
del s

# 次级经向筋条-1
# 根据给定的 longitude_number_1 计算 longitude_number_2 
latitude_hh_secondary=latitude_hh_main /2
longitude_number_2 = int(round(longitude_number_1 * 0.26))  # 4 的数量，四舍五入后强制转换为整数
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.Line(point1=(radius, 0.0), point2=(radius+latitude_hh_secondary, 0.0))
p = mdb.models['Model-1'].Part(name='Part-latitude-secondary1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-latitude-secondary1']
p.BaseShellRevolve(sketch=s, angle=90.0/longitude_number*(longitude_number_1 + longitude_number_2), flipRevolveDirection=OFF)
del s

# 次级经向筋条-2
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.Line(point1=(radius, 0.0), point2=(radius+latitude_hh_secondary, 0.0))
p = mdb.models['Model-1'].Part(name='Part-latitude-secondary2', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-latitude-secondary2']
p.BaseShellRevolve(sketch=s, angle=90.0/longitude_number*(longitude_number_1), flipRevolveDirection=OFF)
del s


#装配
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-shell']
a.Instance(name='Part-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-longitude-main']
a.Instance(name='Part-2', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-longitude-secondary']
a.Instance(name='Part-3', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-latitude-main']
a.Instance(name='Part-4', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-latitude-secondary1']
a.Instance(name='Part-5', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-latitude-secondary2']
a.Instance(name='Part-6', part=p, dependent=ON)
a.rotate(instanceList=('Part-4', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
a.RadialInstancePattern(instanceList=('Part-4', ), axis=(0.0, -1.0, 0.0), number=2, totalAngle=angle_period / 2.0)#主经向加强筋旋转生成两条
a.rotate(instanceList=('Part-5', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
a.rotate(instanceList=('Part-5', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=-angle_period / 4.0)
a.RadialInstancePattern(instanceList=('Part-5', ), axis=(0.0, -1.0, 0.0), number=2, totalAngle=angle_period / 2.0)#次级经向加强筋1旋转生成两条
a.rotate(instanceList=('Part-6', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
a.rotate(instanceList=('Part-6', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=-angle_period / 8.0)
a.RadialInstancePattern(instanceList=('Part-6', ), axis=(0.0, -1.0, 0.0), number=4, totalAngle=angle_period / 4.0 * 3.0)#次级经向加强筋2旋转生成四条


#材料属性
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Density(table=((2.78e-09, ), ))
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((73000.0, 0.33), ))
mdb.models['Model-1'].materials['Material-1'].Plastic(table=(
    (290.0, 0.0), (300.0, 0.00789), (310.0, 0.01121), (320.0, 0.01614), (330.0, 0.02163), (340.0, 0.02586), 
    (350.0, 0.02881), (370.0, 0.036036), (380.0, 0.04172), (390.0, 0.049), (400.0, 0.05802),(410.0, 0.06895), 
    (420.0, 0.08198), (430.0, 0.09738), (440.0, 0.11547), (450.0, 0.13659)))

#截面属性
thickness_stiffer_main = data['thickness_stiffer']
thickness_stiffer_secondary = thickness_stiffer_main
thickness_shell = data['thickness_shell']
thickness_top_shell = 3 * thickness_shell
mdb.models['Model-1'].HomogeneousShellSection(name='Section-stiffer-main', material='Material-1', thickness = thickness_stiffer_main)
mdb.models['Model-1'].HomogeneousShellSection(name='Section-stiffer-secondary', material='Material-1', thickness = thickness_stiffer_secondary)
mdb.models['Model-1'].HomogeneousShellSection(name='Section-shell', material='Material-1', thickness = thickness_shell)
mdb.models['Model-1'].HomogeneousShellSection(name='Section-top_shell', material='Material-1', thickness = thickness_top_shell)

#截面赋予
#主级纬向筋
p = mdb.models['Model-1'].parts['Part-longitude-main']#创建纬向筋集合
p.Set(faces=p.faces[:], name='Set-wei-main')
region = p.sets['Set-wei-main']
p.SectionAssignment(region=region, sectionName='Section-stiffer-main', offset=0.0, offsetType=MIDDLE_SURFACE, )#纬向筋截面赋予
#次级纬向筋
p = mdb.models['Model-1'].parts['Part-longitude-secondary']#创建纬向筋集合
p.Set(faces=p.faces[:], name='Set-wei-secondary')
region = p.sets['Set-wei-secondary']
p.SectionAssignment(region=region, sectionName='Section-stiffer-secondary', offset=0.0, offsetType=MIDDLE_SURFACE, )#纬向筋截面赋予
#主级经向筋
p = mdb.models['Model-1'].parts['Part-latitude-main']#创建经向筋集合
p.Set(faces=p.faces[:], name='Set-jing-main')
region = p.sets['Set-jing-main']
p.SectionAssignment(region=region, sectionName='Section-stiffer-main', offset=0.0, offsetType=MIDDLE_SURFACE, )#经向筋截面赋予
#次级经向筋1
p = mdb.models['Model-1'].parts['Part-latitude-secondary1']#创建纬向筋集合
p.Set(faces=p.faces[:], name='Set-jing-secondary1')
region = p.sets['Set-jing-secondary1']
p.SectionAssignment(region=region, sectionName='Section-stiffer-secondary', offset=0.0, offsetType=MIDDLE_SURFACE, )#纬向筋截面赋予
#次级经向筋2
p = mdb.models['Model-1'].parts['Part-latitude-secondary2']#创建纬向筋集合
p.Set(faces=p.faces[:], name='Set-jing-secondary2')
region = p.sets['Set-jing-secondary2']
p.SectionAssignment(region=region, sectionName='Section-stiffer-secondary', offset=0.0, offsetType=MIDDLE_SURFACE, )#纬向筋截面赋予


#创建全模型
#单周期模型1，向上偏移4000得模型2
a = mdb.models['Model-1'].rootAssembly
instance_names = a.instances.keys()# 获取所有实例名称
part_names_1 = [name for name in instance_names if re.match(r"^Part-\d+", name)]
a.LinearInstancePattern(instanceList=(part_names_1), direction1=(0.0, 1.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=2, number2=1, spacing1=4000.0, spacing2=4000.0)
#选中模型2所有实例
instance_names = a.instances.keys()
instances_2_select = [name for name in instance_names if 'lin-' in name]# 筛选出所有包含 'lin-' 的实例名称，不管其他部分是什么
part_names_2 = [name for name in instances_2_select if re.match(r"^Part-\d+", name)]
a.rotate(instanceList=(part_names_1), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=180.0)#模型1旋转180
a.rotate(instanceList=(part_names_2), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle = angle_period)#模型2旋转45对齐
a.translate(instanceList=(part_names_2), vector=(0.0, -4000.0, 0.0))#模型2下移4000
a.RadialInstancePattern(instanceList=(part_names_1 ), point=(0.0, -1.0, 0.0), axis=(0.0, 1.0, 0.0), number=int(360.0 / angle_period), totalAngle=360.0)#阵列全模型
a.RadialInstancePattern(instanceList=(part_names_2 ), point=(0.0, -1.0, 0.0), axis=(0.0, 1.0, 0.0), number=int(360.0 / angle_period), totalAngle=360.0)#阵列全模型
# 执行布尔合并操作
instance_names = a.instances.keys()# 获取所有实例名称
instances_to_merge = tuple(a.instances[name] for name in instance_names)# 使用实例名称创建实例列表
a.InstanceFromBooleanMerge(name='Part-4', instances=instances_to_merge, originalInstances=DELETE, domain=GEOMETRY)

##创建球壳外表面set
# 每层经向数量
def generate_latitude_sequence(longitude_number_1):
    # 根据给定的 longitude_number_1 计算 longitude_number_2 和 longitude_number_3
    longitude_number_2 = int(round(longitude_number_1 * 0.26))  # 4 的数量，四舍五入后强制转换为整数
    longitude_number_3 = int(round(longitude_number_1 * 0.24))  # 0 和 2 的总数量，四舍五入后强制转换为整数
    print(longitude_number_3)
    # 如果 longitude_number_3 小于 2，则只取 0，不取 2
    if longitude_number_3 <= 2:
        count_0 = longitude_number_3  # 0 和 2 的数量各自是 longitude_number_3 的一半
        count_2 = 0  # 2 的数量与 0 的数量相等
    else:
        count_0 = 2  # 只取 0
        count_2 = longitude_number_3 - 2  # 2 的数量为 0
     # 计算 8 的数量
    count_8 = longitude_number_1  # 8 的数量是 longitude_number_1   
    # 创建纬度数列
    latitude_sequence = [0] * count_0 + [2] * count_2  # 0 和 2 各自的数量是 count_0_2
    latitude_sequence += [4] * longitude_number_2  # 4 的数量是 longitude_number_2
    latitude_sequence += [8] * count_8  # 8 的数量是 longitude_number_1
    
    # 返回最终的纬度数列
    return latitude_sequence

# 生成对应的纬度数列
latitude_sequence = generate_latitude_sequence(longitude_number_1)

# 初始化球面点列表
out_surface_sphere = []
# 计算经度分隔角度
angle_longitude_period = 90.0 / longitude_number  # 等分90度范围
for i in range(2):
    if i == 0:
        angle_1_0 = 90 - angle_longitude_period / 2  # 起始纬度角
    else:
        angle_1_0 = -90 + angle_longitude_period / 2
    for num_out_surface_sphere in latitude_sequence:# 遍历每个纬度圈
        if num_out_surface_sphere == 0:                    
            if i == 0:
                angle_1_0 -= angle_longitude_period# 顶点特殊处理（仅一个点）
                out_surface_sphere.append(((0, radius, 0), ))  # 添加上极点
                continue
            else:
                angle_1_0 += angle_longitude_period# 顶点特殊处理（仅一个点）
                out_surface_sphere.append(((0, -radius, 0), ))  # 添加下极点                  
                continue
        else:        
            angle_latitude_period = angle_period / num_out_surface_sphere# 计算当前纬度圈上的经度分隔角度
            angle_2_0 = angle_latitude_period / 2  # 起始经度角
            num_out_surface_sphere_all = num_out_surface_sphere * 360.0 / angle_period
            for j in range(int(num_out_surface_sphere_all)):
                # 球面点的参数化计算
                x = radius * math.cos(math.radians(angle_1_0)) * math.cos(math.radians(angle_2_0))
                y = radius * math.sin(math.radians(angle_1_0))
                z = radius * math.cos(math.radians(angle_1_0)) * math.sin(math.radians(angle_2_0))
                out_surface_sphere.append(((x, y, z),))  # 添加点坐标
                angle_2_0 += angle_latitude_period
            if i == 0:
                angle_1_0 -= angle_longitude_period# 更新到下一纬度
            else:
                angle_1_0 += angle_longitude_period# 更新到下一纬度


#赋予截面-球壳外表面除去顶部区域
out_surface_sphere_1 = [point for point in out_surface_sphere if point != (((0, radius, 0), ))]#删除上顶面
out_surface_sphere_1 = [point for point in out_surface_sphere_1 if point != (((0, -radius, 0), ))]#删除下底面
p = mdb.models['Model-1'].parts['Part-4']
faces = p.faces.findAt(*out_surface_sphere_1)  # 传入点列表（解包）
p.Set(faces=faces, name='Set-shell')
region = p.sets['Set-shell']
p.SectionAssignment(region=region, sectionName='Section-shell', offset=0.0, offsetType=BOTTOM_SURFACE)

#赋予截面-顶部区域无加强筋且截面加厚
p = mdb.models['Model-1'].parts['Part-4']#创建球面顶部集合
faces = p.faces.findAt(((0, radius, 0), ), ((0, -radius, 0), ))
p.Set(faces=faces, name='Set-top_shell')
region = p.sets['Set-top_shell']
p.SectionAssignment(region=region, sectionName='Section-top_shell', offset=0.0, offsetType=TOP_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)#球面顶部截面赋予

# 创建实心支撑腿
leg_radius = data['leg_radius']
leg_height = data['leg_height']
leg_count = data['leg_count']

# 创建支撑腿部件
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(leg_radius, 0.0))
p = mdb.models['Model-1'].Part(name='Part-leg', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-leg']
p.BaseSolidExtrude(sketch=s, depth=leg_height)
del s

# 装配支撑腿
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-leg']

# 计算球壳赤道半径(y=0处)
equator_radius = radius

# 计算支撑腿的角度间隔
angle_step = 360.0 / leg_count

# 装配支撑腿
for i in range(leg_count):
    # 计算当前支撑腿的角度位置 
    angle = i * angle_step
    
    # 计算支撑腿在赤道处的位置
    x = equator_radius * math.cos(math.radians(angle))
    z = equator_radius * math.sin(math.radians(angle))
    
    # 创建支撑腿实例
    leg_instance = a.Instance(name=f'Leg-{i+1}', part=p, dependent=ON)
    
    # 移动支撑腿到赤道位置
    a.translate(instanceList=(leg_instance.name,), vector=(x, -leg_height, z))

# 将支撑腿与球壳相切并删除内部部分
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-4'] 

for i in range(leg_count):
    # 布尔切割,保留外部
    a.InstanceFromBooleanCut(name=f'Leg-{i+1}-cut', 
                            instanceToBeCut=a.instances[f'Leg-{i+1}'],
                            cuttingInstances=(a.instances['Part-4-1'], ),
                            originalInstances=SUPPRESS)

# 施加截面和材料属性
p = mdb.models['Model-1'].parts['Part-leg']
p.Set(cells=p.cells[:], name='Set-leg')
region = p.sets['Set-leg']
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-leg', material='Material-1')
p.SectionAssignment(region=region, sectionName='Section-leg', offset=0.0)

##创建屈曲分析步
mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', numEigen=1, vectors=18, maxIterations=3000)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('MISES', 'PE', 'PEEQ', 'U')) #设置输出文件   

##创建球壳外表面surface               
a = mdb.models['Model-1'].rootAssembly
side1Faces1 = a.instances['Part-4-1'].faces.findAt(*out_surface_sphere)
a.Surface(side2Faces=side1Faces1, name='Surf-shell')

#均布外压载荷
a = mdb.models['Model-1'].rootAssembly
region = a.surfaces['Surf-shell']
mdb.models['Model-1'].Pressure(name='pressure-1', createStepName='Step-1', region=region, distributionType=UNIFORM, field='', magnitude=1)

##创建位移约束
#创建拆分，进而选择固定约束边界
p = mdb.models['Model-1'].parts['Part-4']#顶部表面通过基准面切割
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
pickedFaces = p.faces.findAt(((0, radius, 0), ),((0, -radius, 0), ))
p.PartitionFaceByDatumPlane(datumPlane=p.datums[9], faces=pickedFaces)#第一次拆分，一分为二
pickedFaces = p.faces.findAt(((0, radius * math.cos( math.radians(0.1)), radius * math.sin( math.radians(0.1))), ), ((0, radius * math.cos( math.radians(0.1)), -radius * math.sin( math.radians(0.1))), ),
    ((0, -radius * math.cos( math.radians(0.1)), radius * math.sin( math.radians(0.1))), ), ((0, -radius * math.cos( math.radians(0.1)), -radius * math.sin( math.radians(0.1))), ))
p.PartitionFaceByDatumPlane(datumPlane=p.datums[10], faces=pickedFaces)#第二次拆分，二分为四

#创建集合施加位移边界条件
p = mdb.models['Model-1'].parts['Part-4']
verts = p.vertices.findAt(((0.0, radius, 0.0), ))
p.Set(vertices=verts, name='Set-top-point')
verts = p.vertices.findAt(((0.0, -radius, 0.0), ))
p.Set(vertices=verts, name='Set-up-point')

#位移约束
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Part-4-1'].sets['Set-top-point']#顶点只允许y方向位移
mdb.models['Model-1'].DisplacementBC(name='BC-top-point', createStepName='Step-1', region=region, u1=0.0, u2=UNSET, u3=0.0,ur1=0.0, ur2=0.0, ur3=0.0, )
region = a.instances['Part-4-1'].sets['Set-up-point']#下顶点固定约束
mdb.models['Model-1'].EncastreBC(name='BC-fix', createStepName='Step-1', region=region, localCsys=None, buckleCase=PERTURBATION_AND_BUCKLING)

# 创建支撑腿底部完全约束
a = mdb.models['Model-1'].rootAssembly
for i in range(leg_count):
    region = a.instances[f'Leg-{i+1}'].faces.findAt(((0.0, -leg_height, 0.0),))  # 找到支撑腿底部面
    mdb.models['Model-1'].EncastreBC(name=f'BC-Leg-{i+1}', createStepName='Step-1', region=region, localCsys=None)

##网格划分

#经向加强筋区域面分割
faces_jing_main = p.sets['Set-jing-main'].faces# 获取经向筋集合中的面
centroids = []# 创建一个空列表来存储所有面心坐标

for face in faces_jing_main:# 遍历主经向筋并获取面心坐标
    coord = face.pointOn  # 获取面心坐标
    centroids.append(coord[0])  # 将面心坐标添加到列表中

for coord in centroids:# 遍历面心坐标并查找对应的面，进行拆分操作
    try:
        face = p.faces.findAt(coord,)  # 使用面心坐标查找面
        p.PartitionFaceByAuto(face=face)               
    except Exception:#避免拆分失败跳出循环
        pass

#纬向加强筋区域面分割
faces_wei_main = p.sets['Set-wei-main'].faces
centroids1 = []# 创建一个空列表来存储所有面心坐标

for face in faces_wei_main:# 遍历主纬向筋并获取面心坐标
    coord = face.pointOn  # 获取面心坐标
    centroids1.append(coord[0])  # 将面心坐标添加到列表中

for coord in centroids1:# 遍历面心坐标并查找对应的面，进行拆分操作
    try:
        face = p.faces.findAt(coord,)  # 使用面心坐标查找面
        p.PartitionFaceByAuto(face=face)               
    except Exception:#避免拆分失败跳出循环
        pass


#全局布种
p = mdb.models['Model-1'].parts['Part-4']
mesh_size=data['mesh_size']
size1=angle_period / 360 * math.pi * radius / 4 / mesh_size  #全局种子尺寸
p.seedPart(size=size1, deviationFactor=0.1, minSizeFactor=0.1)

#纬线局部布种
point_mesh = []# 初始化点列表
number_longitude_period = 360.0 / angle_period * 8  # 单个网格的宽度； 计算经度分隔角度
for i in range(int(number_longitude_period)):
    angle_longitude = angle_period / 8 * (i +0.5)    
    x = radius * math.cos(math.radians(angle_longitude)) # 球面点的参数化计算
    y = 0
    z = radius * math.sin(math.radians(angle_longitude))
    point_mesh.append(((x, y, z),))  # 添加点坐标

pickedEdges = p.edges.findAt(*point_mesh)
p.seedEdgeBySize(edges=pickedEdges, size=size1, deviationFactor=0.1, constraint=FINER) 


#经线局部布种
target_lengths = [latitude_hh_secondary, latitude_hh_main]#加强筋尺寸
tolerance = 1e-6
selected_edges = []#创建垂直于球面边集合,划分局部网格

for edge in p.edges:#选择所有边长为加强筋尺寸的边
    length = edge.getSize()        
    if any(abs(length - target) < tolerance for target in target_lengths):# 检查是否匹配目标长度
        selected_edges.append(edge)

p.seedEdgeBySize(edges=selected_edges, size=latitude_hh_main / 6.0 , deviationFactor=0.1, constraint=FINER)#加强筋高度方向网格为主加强筋/6


#设置网格属性，划分网格
all_faces = p.faces[:]  # 获取部件中的所有面
top_surface_sphere = []#球面顶部区域
p = mdb.models['Model-1'].parts['Part-4']
for point in [
    (((radius * math.sin(math.radians(0.1))) / math.sqrt(2), radius * math.cos(math.radians(0.1)), radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((-radius * math.sin(math.radians(0.1))) / math.sqrt(2), radius * math.cos(math.radians(0.1)), radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((radius * math.sin(math.radians(0.1))) / math.sqrt(2), radius * math.cos(math.radians(0.1)), -radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((-radius * math.sin(math.radians(0.1))) / math.sqrt(2), radius * math.cos(math.radians(0.1)), -radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((radius * math.sin(math.radians(0.1))) / math.sqrt(2), -radius * math.cos(math.radians(0.1)), radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((-radius * math.sin(math.radians(0.1))) / math.sqrt(2), -radius * math.cos(math.radians(0.1)), radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((radius * math.sin(math.radians(0.1))) / math.sqrt(2), -radius * math.cos(math.radians(0.1)), -radius * math.sin(math.radians(0.1)) / math.sqrt(2)),),
    (((-radius * math.sin(math.radians(0.1))) / math.sqrt(2), -radius * math.cos(math.radians(0.1)), -radius * math.sin(math.radians(0.1)) / math.sqrt(2)),)]:
    top_surface_sphere.append(point)

# if angle_period > 30.0:
    # pickedRegions = p.faces.findAt(*top_surface_sphere)  #  选中球面顶部区域
    # pickedRegions_set = {face.index for face in pickedRegions}  # 将 pickedRegions 转为面对象的索引集合存储面对象的索引
    # remaining_faces = [face for face in all_faces if face.index not in pickedRegions_set] # 过滤掉球面顶部区域，获得剩余面
    # p.setMeshControls(regions=remaining_faces, elemShape=QUAD, technique=SWEEP)#设置为扫掠网格
    # p.generateMesh()
    # p.setMeshControls(regions=pickedRegions, elemShape=QUAD_DOMINATED, technique=STRUCTURED) #设置网格属性结构化网格
    # p.generateMesh(regions=pickedRegions)
# else:
p.setMeshControls(regions=all_faces, elemShape=QUAD_DOMINATED, technique=STRUCTURED)
p.generateMesh()

# 支撑腿网格划分
for i in range(leg_count):
    p = mdb.models['Model-1'].parts['Part-leg']
    # 全局网格种子
    p.seedPart(size=leg_radius / 4.0, deviationFactor=0.1, minSizeFactor=0.1)
    
    # 加密与球壳连接部分的网格
    edges_to_refine = p.edges.findAt(((0.0, 0.0, 0.0),))  # 找到靠近球壳的边
    p.seedEdgeBySize(edges=edges_to_refine, size=leg_radius / 8.0, deviationFactor=0.1, constraint=FINER)
    
    # 设置网格控制
    p.setMeshControls(regions=p.cells[:], elemShape=TET, technique=FREE)
    
    # 网格划分
    p.generateMesh()

# 确保支撑腿与球壳连接部分的网格一致性
a = mdb.models['Model-1'].rootAssembly
for i in range(leg_count):
    leg_instance = a.instances[f'Leg-{i+1}-cut']
    shell_instance = a.instances['Part-4-1']
    a.generateMesh(regions=(leg_instance, shell_instance))
##创建作业
mdb.Job(name='buckle', model='Model-1', numCpus=10, numDomains=10)


#提取位移场
mdb.models['Model-1'].keywordBlock.synchVersions(storeNodesAndElements=False)
keyword_block = mdb.models['Model-1'].keywordBlock.sieBlocks# 获取关键字块内容

for i, block in enumerate(keyword_block):# 动态查找插入位置
    if '** OUTPUT REQUESTS' in block:
        insert_position = i   # 在 "OUTPUT REQUESTS" 之后插入
        break
  
mdb.models['Model-1'].keywordBlock.insert(insert_position, """
*node file
u
""")

mdb.saveAs(pathName='1')

#-------------------提交作业------------------
mdb.jobs['buckle'].submit()
mdb.jobs['buckle'].waitForCompletion()


#-------------------后屈曲设置------------------
mdb = openMdb(pathName='D:/加筋球壳/3.28/1.cae')#打开cae
mdb.Model(name='Model-1-post', objectToCopy=mdb.models['Model-1'])#复制模型
maxNum_Inc = data['maxNum_Inc']
mdb.models['Model-1-post'].StaticRiksStep(name='Step-1', previous='Initial', maintainAttributes=True, maxNumInc=maxNum_Inc, initialArcInc=0.005, maxArcInc=0.1, nlgeom=ON)#更改为弧长法
mdb.models['Model-1-post'].keywordBlock.setValues(edited = 0)#删除关键字原有编辑
#引入特征值缺陷
scale=data['scale']
mdb.models['Model-1-post'].keywordBlock.synchVersions(storeNodesAndElements=False)
keyword_block = mdb.models['Model-1-post'].keywordBlock.sieBlocks# 获取关键字块内容

for i, block in enumerate(keyword_block):# 动态查找插入位置
    if '** STEP: Step-1' in block:
        insert_position = i   # 在 "OUTPUT REQUESTS" 之后插入
        break
 
# 使用f-string动态替换比例因子
imperfection_block = """
*IMPERFECTION, FILE=buckle, STEP=1
1,{}
""".format(scale)

# 插入关键字块
mdb.models['Model-1-post'].keywordBlock.insert(
    insert_position, 
    imperfection_block
)

mdb.Job(name='post-buckle', model='Model-1-post', numCpus=10, numDomains=10)#创建作业
mdb.save()# 保存cae文件

#-------------------提交作业------------------
mdb.jobs['post-buckle'].submit()
mdb.jobs['post-buckle'].waitForCompletion()


#-------------------读取装配体总质量-----------
assembly = mdb.models['Model-1'].rootAssembly
mass_properties = assembly.getMassProperties()
mass_in_tons = mass_properties['mass']# 获取质量值并转换为千克 (1 吨 = 1000 千克)
mass_in_kg = mass_in_tons * 1000  # 将质量从吨转换为千克
total_mass = ["mass =", str(mass_in_kg)]  # 转换为千克后保存
with open('D:/TEMP/328/opt.txt', 'w') as file:
    file.write("mass = %s\n" % (mass_in_kg))    #使用追加模式 ('a')


#-------------------读取结果中特征值-----------
myOdb = session.openOdb(name='D:/加筋球壳/3.28/buckle.odb', readOnly=False)# 打开 ODB 文件
myFrames = myOdb.steps["Step-1"].frames# 获取 Step-1 中的帧
mydes = myFrames[1].description# 获取第二帧的描述信息
with open('D:/TEMP/328/opt.txt', 'a') as file:
    file.write("%s\n" % (mydes))


#-------------------读取弧长第一个峰值为临界载荷-----------
myOdb = session.openOdb(name='D:/加筋球壳/328/post-buckle.odb', readOnly=True)  # 以只读方式打开# 1. 打开 odb 文件
step = myOdb.steps['Step-1']# 2. 获取分析步骤和输出
history_region = step.historyRegions['Assembly ASSEMBLY']  # 获取历史区域# 请替换为你的区域名
lpf_output = history_region.historyOutputs['LPF'].data  # 获取 LPF 输出数据
peak_risk, peak_value = None, None
max_idx = max(range(len(lpf_output)), key=lambda i: lpf_output[i][1])
peak_risk, peak_value = lpf_output[max_idx]## 查找数据中最大值对应的索引并赋值

with open('D:/TEMP/328/opt.txt', 'a') as file:
    file.write("risk-%s=%s\n" % (peak_risk, peak_value))

#-------------------记录特征值运行时间-----------
with open('D:/TEMP/328/buckle.msg', 'r') as file:# 读取 .msg 文件的最后一行
    last_line = file.readlines()[-1].strip() 

with open('D:/TEMP/328/opt.txt', 'a') as output_file:# 直接将最后一行写入到 opt.txt 文件
    output_file.write(last_line + "\n")  # 写入最后一行


#-------------------记录后屈曲运行时间-----------
with open('D:/TEMP/328/post-buckle.msg', 'r') as file:# 读取 .msg 文件的最后一行
    last_line = file.readlines()[-1].strip() 

with open('D:/TEMP/328/opt.txt', 'a') as output_file:# 直接将最后一行写入到 opt.txt 文件
    output_file.write(last_line + "\n")  # 写入最后一行