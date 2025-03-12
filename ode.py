#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:50:56 2025

@author: heyubo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Wilson-Cowan 典型的发放率函数
def f(x, beta=60):
    return 1 / (1 + np.exp(-beta * x))

# 无时滞 Wilson-Cowan 方程
def wilson_cowan_no_delay(t, Y, c1, c2, c3, c4, P, Q, beta):
    u, v = Y  # 当前时间点的 u 和 v
    du_dt = -u + f(c1 * u + c2 * v + P, beta)
    dv_dt = -v + f(c3 * u + c4 * v + Q, beta)
    return [du_dt, dv_dt]

# 参数设置
beta = 1
P, Q = -1, -4
c1, c2, c3, c4 = 10, -10, 10, 2

# 初始条件
Y0 = [0.1, 0.1]

# 时间范围
t_values = np.linspace(0,50, 500)

# 求解无时滞的 Wilson-Cowan 方程
sol_no_delay = solve_ivp(
    wilson_cowan_no_delay, 
    [t_values[0], t_values[-1]], 
    Y0, 
    t_eval=t_values, 
    args=(c1, c2, c3, c4, P, Q, beta)
)

# 提取结果
u_no_delay, v_no_delay = sol_no_delay.y

# 绘制无时滞的 Wilson-Cowan 结果
plt.figure(figsize=(10, 5))
plt.plot(t_values, u_no_delay, label="Excitatory (u) - No Delay", color="blue")
plt.plot(t_values, v_no_delay, label="Inhibitory (v) - No Delay", color="red")
plt.xlabel("Time")
plt.ylabel("Activity")
plt.legend()
plt.title("Wilson-Cowan Model without Delay")
plt.grid()
plt.show()
