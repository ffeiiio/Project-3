#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:33:09 2025

@author: heyubo
"""

import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint  # 时滞微分方程（DDE）求解库

# 定义 Wilson-Cowan 发放率函数
def f(x, beta=60):
    return 1 / (1 + np.exp(-beta * x))

# 单时滞 Wilson-Cowan 方程
def wilson_cowan_single_delay(Y, t, c1, c2, c3, c4, P, Q, tau, beta):
    u, v = Y(t)  # 当前位置的变量值
    u_tau, v_tau = Y(t - tau) if t >= tau else Y(0)  # 过去 tau 时间的变量值

    du_dt = -u + f(c1 * u_tau + c2 * v_tau + P, beta)
    dv_dt = -v + f(c3 * u_tau + c4 * v_tau + Q, beta)
    return np.array([du_dt, dv_dt])

# 多重时滞 Wilson-Cowan 方程
def wilson_cowan_multiple_delay(Y, t, c1, c2, c3, c4, P, Q, tau1, tau2, beta):
    u, v = Y(t)  # 当前位置的变量值
    u_tau1, v_tau1 = Y(t - tau1) if t >= tau1 else Y(0)  # 过去 tau1 时间的变量值
    u_tau2, v_tau2 = Y(t - tau2) if t >= tau2 else Y(0)  # 过去 tau2 时间的变量值

    du_dt = -u + f(c1 * u_tau1 + c2 * v_tau2 + P, beta)
    dv_dt = -v + f(c3 * u_tau2 + c4 * v_tau1 + Q, beta)
    return np.array([du_dt, dv_dt])

# 参数设置
beta = 60
P, Q = 0.65, 0.5
c1, c2, c3, c4 = -1, -0.4, -1, 0
tau_single = 0.2  # 单时滞
tau1_multiple, tau2_multiple = 0.3, 0.1  # 多时滞

# 初始条件（历史函数）
def history_function(t):
    return np.array([0.1, 0.1])  # 初始时刻之前，u 和 v 都设定为 0.1

# 时间范围
t_values = np.linspace(0, 10, 500)

# 计算单时滞数值解
sol_single = ddeint(wilson_cowan_single_delay, history_function, t_values, fargs=(c1, c2, c3, c4, P, Q, tau_single, beta))

# 计算多重时滞数值解
sol_multiple = ddeint(wilson_cowan_multiple_delay, history_function, t_values, fargs=(c1, c2, c3, c4, P, Q, tau1_multiple, tau2_multiple, beta))

# 提取 u 和 v 变量
u_single, v_single = sol_single[:, 0], sol_single[:, 1]
u_multiple, v_multiple = sol_multiple[:, 0], sol_multiple[:, 1]

# 绘制单时滞结果
plt.figure(figsize=(10, 5))
plt.plot(t_values, u_single, label="Excitatory (u) - Single Delay", color="blue")
plt.plot(t_values, v_single, label="Inhibitory (v) - Single Delay", color="red")
plt.xlabel("Time")
plt.ylabel("Activity")
plt.legend()
plt.title("Wilson-Cowan Model with Single Delay (τ=0.2)")
plt.grid()
plt.show()

# 绘制多时滞结果
plt.figure(figsize=(10, 5))
plt.plot(t_values, u_multiple, label="Excitatory (u) - Multiple Delay", color="blue")
plt.plot(t_values, v_multiple, label="Inhibitory (v) - Multiple Delay", color="red")
plt.xlabel("Time")
plt.ylabel("Activity")
plt.legend()
plt.title("Wilson-Cowan Model with Multiple Delays (τ1=0.3, τ2=0.1)")
plt.grid()
plt.show()
