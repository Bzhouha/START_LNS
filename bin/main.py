# import numpy as np
import pandas as pd

# HLNS
# Wave numbers and frequency
hlns_info = pd.read_csv('out/hlns_info.csv', dtype=float)
data = hlns_info.values
In = int(data[0][0])
Jn = int(data[0][1])
Kn = int(data[0][2])
hlns_a = complex(data[0][3], data[0][4])
hlns_b = complex(data[0][5], data[0][6])
hlns_o = complex(data[0][7], data[0][8])
del hlns_info
del data
# Shape function
hlns_res = pd.read_csv('out/turtle.csv', dtype=float)
data = hlns_res.values
hlns_rho_r = data[:, 0].reshape(In, Jn, Kn)
hlns_rho_i = data[:, 1].reshape(In, Jn, Kn)
hlns_rho_m = data[:, 2].reshape(In, Jn, Kn)
hlns_u_r = data[:, 3].reshape(In, Jn, Kn)
hlns_u_m = data[:, 4].reshape(In, Jn, Kn)
hlns_u_i = data[:, 5].reshape(In, Jn, Kn)
hlns_v_r = data[:, 6].reshape(In, Jn, Kn)
hlns_v_i = data[:, 7].reshape(In, Jn, Kn)
hlns_v_m = data[:, 8].reshape(In, Jn, Kn)
hlns_w_r = data[:, 9].reshape(In, Jn, Kn)
hlns_w_i = data[:, 10].reshape(In, Jn, Kn)
hlns_w_m = data[:, 11].reshape(In, Jn, Kn)
hlns_T_r = data[:, 12].reshape(In, Jn, Kn)
hlns_T_i = data[:, 13].reshape(In, Jn, Kn)
hlns_T_m = data[:, 14].reshape(In, Jn, Kn)
del hlns_res
del data

# 读取网格坐标
grid = pd.read_csv('out/grid.csv', dtype=float)
data = grid.values
xx = data[:, 0].reshape(In, Jn, Kn)
yy = data[:, 1].reshape(In, Jn, Kn)
zz = data[:, 2].reshape(In, Jn, Kn)
del data
del grid

# 读取基本流流场
flow = pd.read_csv('out/flow.csv', dtype=float)
data = flow.values
rho0 = data[:, 0].reshape(In, Jn, Kn)
u0 = data[:, 1].reshape(In, Jn, Kn)
v0 = data[:, 2].reshape(In, Jn, Kn)
w0 = data[:, 3].reshape(In, Jn, Kn)
T0 = data[:, 4].reshape(In, Jn, Kn)
del data
del flow

# LPSE
# Wave numbers and frequency
lpse_info = pd.read_csv('out/lpse_info.plt', dtype=float, sep='\s+')
data = lpse_info.values
lpse_Ma = data[:, 2].reshape(In)
lpse_Re = data[:, 3].reshape(In)
lpse_Te = data[:, 4].reshape(In)
lpse_o_r = data[:, 5].reshape(In)
lpse_o_i = data[:, 6].reshape(In)
lpse_a_r = data[:, 7].reshape(In)
lpse_a_i = data[:, 8].reshape(In)
lpse_b_r = data[:, 9].reshape(In)
lpse_b_i = data[:, 10].reshape(In)
lpse_c_r = data[:, 11].reshape(In)
lpse_c_i = data[:, 12].reshape(In)
del lpse_info
del data

# Shape function
lpse_res = pd.read_csv('out/lpse.csv', dtype=float)
data = lpse_res.values
lpse_rho_r = data[:, 0].reshape(In, Jn, Kn)
lpse_rho_i = data[:, 1].reshape(In, Jn, Kn)
lpse_rho_m = data[:, 2].reshape(In, Jn, Kn)
lpse_u_r = data[:, 3].reshape(In, Jn, Kn)
lpse_u_m = data[:, 4].reshape(In, Jn, Kn)
lpse_u_i = data[:, 5].reshape(In, Jn, Kn)
lpse_v_r = data[:, 6].reshape(In, Jn, Kn)
lpse_v_i = data[:, 7].reshape(In, Jn, Kn)
lpse_v_m = data[:, 8].reshape(In, Jn, Kn)
lpse_w_r = data[:, 9].reshape(In, Jn, Kn)
lpse_w_i = data[:, 10].reshape(In, Jn, Kn)
lpse_w_m = data[:, 11].reshape(In, Jn, Kn)
lpse_T_r = data[:, 12].reshape(In, Jn, Kn)
lpse_T_i = data[:, 13].reshape(In, Jn, Kn)
lpse_T_m = data[:, 14].reshape(In, Jn, Kn)
del lpse_res
del data
