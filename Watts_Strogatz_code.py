#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import networkx as nx
import numpy as np
from ddeint import ddeint
import matplotlib.pyplot as plt


N = 50  
K = 10   
p = 0.3  
num_realisations = 
t = np.linspace(0, 20, 1000)  
tau = 3  


c1 = 10
c2 = -10
c3 = 10
c4 = 2
P = -4
Q = -2.4
epsilon = 5  
beta = 1  


def f(x):
    return 1 / (1 + np.exp(-beta * x))


def history(t):
    return np.zeros(2 * N)  

all_d_solutions = []  


for realisation in range(num_realisations):
    
    G = nx.watts_strogatz_graph(N, K, p)
    adj_matrix = nx.to_numpy_array(G)
    adj_matrix = np.where(adj_matrix > 0, 1, 0)  

    
    def w_c(Y, t, adj_matrix):
        ud = Y(t)[:N]  
        vd = Y(t)[N:]  
        u_tau = Y(t - tau)[:N]  

        dud = np.zeros(N)  
        dvd = np.zeros(N)  

        for i in range(N):
            sum_u_tau = np.sum(adj_matrix[i] * u_tau)  
            dud[i] = -ud[i] + f(c1 * ud[i] + c2 * vd[i] + P + epsilon * sum_u_tau)
            dvd[i] = -vd[i] + f(c3 * ud[i] + c4 * vd[i] + Q)

        return np.concatenate([dud, dvd])

    
    del_sol = ddeint(w_c, history, t, fargs=(adj_matrix,))
    all_d_solutions.append(del_sol)


all_d_solutions_array = np.array(all_d_solutions)
avg_sol_d = np.mean(all_d_solutions_array, axis=0)


ud_sol = avg_sol_d[:, :N]  


plt.figure(figsize=(10, 6))
plt.imshow(
    ud_sol.T,  
    aspect='auto',
    cmap='viridis',
    extent=[t[0], t[-1], 0, N]
)
plt.colorbar(label='Excitatory Activity (u)')
plt.xlabel('Time')
plt.ylabel('Node Index')
plt.title('Heatmap of Neural Activity (Delayed System, Epsilon = 5)')
plt.savefig("delayed_epsilon_5_heatmap.png", dpi=300)
plt.show()

