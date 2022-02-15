import os
import sys
import numpy as np

file1 = np.load('rectangle_v2.npz')
# V = file['V']
EToV = file1['SEToV'];
phys_tag = file1['SE_phys_tag'];
geom_tag = file1['SE_geom_tag'];
part_tag = file1['SE_part_tag'];
Etype = file1['SE_Etype'];

print(EToV)
print(phys_tag)
print(geom_tag)
print(part_tag)
print(Etype)

file2 = np.load('cuboid_v2.npz')
# V = file['V']
EToV = file2['VEToV']
phys_tag = file2['VE_phys_tag']
geom_tag = file2['VE_geom_tag']
part_tag = file2['VE_part_tag']
Etype = file2['VE_Etype']

print(EToV)
print(phys_tag)
print(geom_tag)
print(part_tag)
print(Etype)

file3 = np.load('rectangle_v4.npz')
# V = file['V']
EToV = file3['SEToV']
phys_tag = file3['SE_phys_tag']
geom_tag = file3['SE_geom_tag']
part_tag = file3['SE_part_tag']
Etype = file3['SE_Etype']

print(EToV)
print(phys_tag)
print(geom_tag)
print(part_tag)
print(Etype)

file4 = np.load('cuboid_v4.npz')
# V = file['V']
EToV = file4['VEToV']
phys_tag = file4['VE_phys_tag']
geom_tag = file4['VE_geom_tag']
part_tag = file4['VE_part_tag']
Etype = file4['VE_Etype']

print(EToV)
print(phys_tag)
print(geom_tag)
print(part_tag)
print(Etype)