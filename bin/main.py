import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.set_printoptions(precision=16)

# 1.读取信息
# 1.1 读取HLNS信息
# 1.1.1 Wave numbers and frequency
hlns_info = pd.read_csv('out/hlns_info.csv', dtype=np.float64)
data = hlns_info.values
In = int(data[0][0])
Jn = int(data[0][1])
Kn = int(data[0][2])
hlns_a_r = data[0][3]
hlns_a_i = data[0][4]
hlns_b_r = data[0][5]
hlns_b_i = data[0][6]
hlns_o_r = data[0][7]
hlns_o_i = data[0][8]
del hlns_info
del data

# 1.1.2 Shape function
hlns_res = pd.read_csv('out/turtle.csv', dtype=np.float64)
data = hlns_res.values
hlns_rho_r = data[:, 0].reshape(In, Jn, Kn)
hlns_rho_i = data[:, 1].reshape(In, Jn, Kn)
hlns_rho_m = data[:, 2].reshape(In, Jn, Kn)
hlns_u_r = data[:, 3].reshape(In, Jn, Kn)
hlns_u_i = data[:, 4].reshape(In, Jn, Kn)
hlns_u_m = data[:, 5].reshape(In, Jn, Kn)
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

# 1.1.3 HLNS数据
hlns_a = hlns_a_r + 1j * hlns_a_i
hlns_b = hlns_b_r + 1j * hlns_b_i
hlns_o = hlns_o_r + 1j * hlns_o_i
hlns_rho = hlns_rho_r + 1j * hlns_rho_i
hlns_u = hlns_u_r + 1j * hlns_u_i
hlns_v = hlns_v_r + 1j * hlns_v_i
hlns_w = hlns_w_r + 1j * hlns_w_i
hlns_T = hlns_T_r + 1j * hlns_T_i
hlns = np.zeros((5, In, Jn, Kn), np.complex128)
for i in range(In):
    for j in range(Jn):
        for k in range(Kn):
            hlns[0, i, j, k] = hlns_rho[i, j, k]
            hlns[1, i, j, k] = hlns_u[i, j, k]
            hlns[2, i, j, k] = hlns_v[i, j, k]
            hlns[3, i, j, k] = hlns_w[i, j, k]
            hlns[4, i, j, k] = hlns_T[i, j, k]
del hlns_a_r, hlns_b_r, hlns_o_r
del hlns_a_i, hlns_b_i, hlns_o_i
del hlns_rho_r, hlns_u_r, hlns_v_r, hlns_w_r, hlns_T_r
del hlns_rho_i, hlns_u_i, hlns_v_i, hlns_w_i, hlns_T_i
del hlns_rho_m, hlns_u_m, hlns_v_m, hlns_w_m, hlns_T_m
del hlns_rho, hlns_u, hlns_v, hlns_w, hlns_T

# 1.2 读取LPSE信息
# 1.2.1 Wave numbers and frequency
lpse_info = pd.read_csv('out/lpse_info.plt', dtype=np.float64, sep='\s+')
data = lpse_info.values
lpse_x = data[:, 1].reshape(In)
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

# 1.2.2 Shape function
lpse_res = pd.read_csv('out/lpse.csv', dtype=np.float64)
data = lpse_res.values
lpse_rho_r = data[:, 0].reshape(In, Jn, Kn)
lpse_rho_i = data[:, 1].reshape(In, Jn, Kn)
lpse_rho_m = data[:, 2].reshape(In, Jn, Kn)
lpse_u_r = data[:, 3].reshape(In, Jn, Kn)
lpse_u_i = data[:, 4].reshape(In, Jn, Kn)
lpse_u_m = data[:, 5].reshape(In, Jn, Kn)
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

# 1.2.3 LPSE数据
lpse_a = lpse_a_r + 1j * lpse_a_i
lpse_b = lpse_b_r + 1j * lpse_b_i
lpse_o = lpse_o_r + 1j * lpse_o_i
lpse_rho = lpse_rho_r + 1j * lpse_rho_i
lpse_u = lpse_u_r + 1j * lpse_u_i
lpse_v = lpse_v_r + 1j * lpse_v_i
lpse_w = lpse_w_r + 1j * lpse_w_i
lpse_T = lpse_T_r + 1j * lpse_T_i
lpse = np.zeros((5, In, Jn, Kn), np.complex128)
for i in range(In):
    for j in range(Jn):
        for k in range(Kn):
            lpse[0, i, j, k] = lpse_rho[i, j, k]
            lpse[1, i, j, k] = lpse_u[i, j, k]
            lpse[2, i, j, k] = lpse_v[i, j, k]
            lpse[3, i, j, k] = lpse_w[i, j, k]
            lpse[4, i, j, k] = lpse_T[i, j, k]
del lpse_a_r, lpse_b_r, lpse_o_r, lpse_c_r
del lpse_a_i, lpse_b_i, lpse_o_i, lpse_c_i
del lpse_rho_r, lpse_u_r, lpse_v_r, lpse_w_r, lpse_T_r
del lpse_rho_i, lpse_u_i, lpse_v_i, lpse_w_i, lpse_T_i
del lpse_rho_m, lpse_u_m, lpse_v_m, lpse_w_m, lpse_T_m
del lpse_rho, lpse_u, lpse_v, lpse_w, lpse_T

# 1.3 读取基本流
# 1.3.1 读取网格坐标
grid = pd.read_csv('out/grid.csv', dtype=np.float64)
data = grid.values
xx = data[:, 0].reshape(In, Jn, Kn)
yy = data[:, 1].reshape(In, Jn, Kn)
zz = data[:, 2].reshape(In, Jn, Kn)
del grid
del data

# 1.3.2 读取基本流流场
flow = pd.read_csv('out/flow.csv', dtype=np.float64)
data = flow.values
rho0 = data[:, 0].reshape(In, Jn, Kn)
u0 = data[:, 1].reshape(In, Jn, Kn)
v0 = data[:, 2].reshape(In, Jn, Kn)
w0 = data[:, 3].reshape(In, Jn, Kn)
T0 = data[:, 4].reshape(In, Jn, Kn)
del flow
del data

# 1.3.3 基本流数据
loc = xx[:, 0, 0]
base = np.zeros((5, In, Jn, Kn), np.float64)
for i in range(In):
    for j in range(Jn):
        for k in range(Kn):
            base[0, i, j, k] = rho0[i, j, k]
            base[1, i, j, k] = u0[i, j, k]
            base[2, i, j, k] = v0[i, j, k]
            base[3, i, j, k] = w0[i, j, k]
            base[4, i, j, k] = T0[i, j, k]
del rho0, u0, v0, w0, T0

# 2.处理数据

# 2.1 HLNS
ln = np.exp(1j * hlns_a * loc)
hlns_dist = np.zeros((5, In, Jn, Kn), np.complex128)
for l in range(5):
    for i in range(In):
        for j in range(Jn):
            for k in range(Kn):
                hlns_dist[l, i, j, k] = hlns[l, i, j, k] * ln[i]
del ln
del hlns
hlns = np.zeros((5, In, Jn, Kn), np.float64)
for l in range(5):
    for i in range(In):
        for j in range(Jn):
            for k in range(Kn):
                hlns[l, i, j, k] = abs(hlns_dist[l, i, j, k])
del hlns_dist

# 2.2 LPSE
# 300个小边长
small_edges = np.zeros(In - 1, np.float64)
for i in range(In - 1):
    small_edges[i] = loc[i + 1] - loc[i]
# 300个小梯形的面积
small_area = np.zeros(In - 1, np.complex128)
for i in range(In - 1):
    small_area[i] = 0.5 * (lpse_a[i] + lpse_a[i + 1]) * small_edges[i]
# 依次摞加得到积分
jifen = np.zeros(In, np.complex128)
for i in range(1, In):
    jifen[i] = jifen[0] + small_area[i - 1]
# e底积分
ln = np.zeros(In, np.complex128)
for i in range(In):
    ln[i] = np.exp(1j * jifen[i])
lpse_dist = np.zeros((5, In, Jn, Kn), np.complex128)
for l in range(5):
    for i in range(In):
        for j in range(Jn):
            for k in range(Kn):
                lpse_dist[l, i, j, k] = lpse[l, i, j, k] * ln[i]
del ln
del lpse
lpse = np.zeros((5, In, Jn, Kn), np.float64)
for l in range(5):
    for i in range(In):
        for j in range(Jn):
            for k in range(Kn):
                lpse[l, i, j, k] = abs(lpse_dist[l, i, j, k])
del lpse_dist

# 3.绘图
# 3.1 流向位置
iloc = 275

# 3.2 准备要使用的数据
# 3.2.1 HLNS
drt1 = np.zeros(In, np.float64)
for i in range(In):
    drt1[i] = np.max(hlns[1, i, :, 0])
# 3.2.2 LPSE
drt2 = np.zeros(In, np.float64)
for i in range(In):
    drt2[i] = np.max(lpse[1, i, :, 0])

# 3.3 绘图-最大值
plt.figure(figsize=(9, 9), dpi=200)
plt.subplot(211)
plt.plot(loc, drt1, label='HLNS', linewidth=0.8, linestyle="-")
plt.title("HLNS-Sample:x=" + str(iloc))
plt.xlabel("iloc")
plt.ylabel("$|u_{max}|$")
plt.legend()
plt.subplot(212)
plt.scatter(loc, drt2, label='LPSE', s=0.3)
plt.title("LPSE-Sample:x=" + str(iloc))
plt.xlabel("iloc")
plt.ylabel("$|u_{max}|$")
plt.legend()
plt.show()
