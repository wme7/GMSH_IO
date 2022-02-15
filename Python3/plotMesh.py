import os
import sys
import numpy as np

file = np.load('FEMmeshV2.npz')
V = file['V']
points = file['PEToV']
curves = file['LEToV']
surfaces= file['SEToV']
volumes = file['VEToV']