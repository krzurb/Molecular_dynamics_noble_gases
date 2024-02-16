from numba import jit
import numpy as np
import time
import sys

#constant
k_B=8.31e-3
#creating initial state of the crystal
@jit(nopython=True)
def initial_state(n,a,N):
    b0 = np.array([a, 0, 0])
    b1 = np.array([a/2, (a * np.sqrt(3))/2, 0])
    b2 = np.array([a/2, (a * np.sqrt(3))/6, (a * np.sqrt((2/3)))])
    r_init=np.zeros((N,3))
    for i2 in range(n):
        for i1 in range(n):
            for i0 in range(n):
                i=i0+i1*n+i2*n**2
                r_init[i]=(i0-(n-1)/2)*b0+(i1-(n-1)/2)*b1+(i2-(n-1)/2)*b2
    return r_init
#intial momemntum and kinetic energy
@jit(nopython=True)
def initial_momentum(r_x,r_y,r_z,r1,r2,r3,N,k_B,T_0,m):
    E_kin_x=-1/2 * k_B * T_0 * np.log(r_x)
    E_kin_y=-1/2 * k_B * T_0 * np.log(r_y)
    E_kin_z=-1/2 * k_B * T_0 * np.log(r_z)
    p_init = np.zeros((N, 3),dtype=np.float64)
    p_init[:, 0] = r1 * np.sqrt(2 * m * E_kin_x)
    p_init[:, 1] = r2 * np.sqrt(2 * m * E_kin_y)
    p_init[:, 2] = r3 * np.sqrt(2 * m * E_kin_z)
    P=np.sum(p_init,axis=0)
    p_res=p_init-P/N
    return p_res
#potential and forces
@jit(nopython=True)
def V_P(ri,rj,e,R):
    r_ij = (ri-rj)
    r_mod=np.sqrt(r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2])
    RR=R/r_mod
    V=e*(RR*RR*RR*RR*RR*RR*RR*RR*RR*RR*RR*RR-2*RR*RR*RR*RR*RR*RR)
    return V
@jit(nopython=True)
def V_S(r, L,f):
    r_abs=np.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
    condition = r_abs < L
    if condition:
        return 0.0
    else:
        return 0.5 * f * (r_abs - L) *(r_abs-L)
@jit(nopython=True)
def F_P(ri,rj,e,R):
    r_ij = (ri-rj)
    r_mod=np.sqrt(r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2])
    RR=R/r_mod
    return 12*e*(RR*RR*RR*RR*RR*RR*RR*RR*RR*RR*RR*RR-RR*RR*RR*RR*RR*RR)*((r_ij)/(r_mod*r_mod))
@jit(nopython=True)
def F_S(r, L,f):
    eps=1e-12
    r_abs = np.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
    V0=np.array([0.0,0.0,0.0],dtype=np.float64)
    condition = r_abs < L
    if condition:
        return V0
    else:
        return f * (L - r_abs) * (r / (r_abs+eps))
@jit(nopython=True)
def force(N, r, e, R, L, f):
    F = np.zeros((N, 3),dtype=np.float64)
    P = 0.0
    constP = 1 / (4 * np.pi * L * L)
    for i in range(N):
        FS = F_S(r[i], L, f)
        for j in range(i + 1, N):
            force_ij = F_P(r[i], r[j], e, R)
            F[i] += force_ij
            F[j] -= force_ij
        F[i] += FS
        P += constP * np.sqrt(FS[0]*FS[0]+FS[1]*FS[1]+FS[2]*FS[2])
    return F, P
@jit(nopython=True)
def potential_energy(N, r, e, R, L, f):
    V_Sum = 0
    for i in range(N):
        V_Sum = V_Sum + V_S(r[i], L, f)
        for j in range(i + 1, N):
            V_Sum = V_Sum + V_P(r[i], r[j], e, R)
    return V_Sum
#########################################################################################          main function
@jit(nopython=True)
def time_simulation(Ek_init, V_init, Pre_init, T_0, p_init, F_init, S_o, S_d, S_out, S_xyz, N, tau, r_init, e, R, L, f, m):
    r = np.zeros((S_o + S_d, N, 3),dtype=np.float64)
    r[0] = r_init
    p = np.zeros((S_o + S_d, N, 3),dtype=np.float64)
    p[0] = p_init
    Force = F_init
    r_xyz = np.zeros((int((S_o + S_d) / S_xyz), N, 3),dtype=np.float64)
    Ekin = np.zeros((int((S_o + S_d) / S_xyz), N),dtype=np.float64)
    print(0, Ek_init + V_init, V_init, T_0, Pre_init)
    half_m=0.5/m
    Tem_const=2 / (3 * N * k_B)
    count=0
    H_bar=0.0
    T_bar=0.0
    P_bar=0.0
    for t in range(1, S_o + S_d):
        p_half_tau = p[t - 1] + 0.5 * Force * tau
        r[t] = r[t - 1] + (p_half_tau / m) * tau
        Force, Pre = force(N, r[t], e, R, L, f)
        p[t] = p_half_tau + 0.5 * Force * tau
        if t % S_out == 0:
            EK = half_m * np.sum(p[t] * p[t])
            Tempe=Tem_const*EK
            V = potential_energy(N, r[t], e, R, L, f)
            print(t * tau, EK+V, V, Tempe, Pre)
            if t>S_o:
                H_bar+=EK+V
                T_bar+=Tempe
                P_bar+=Pre
                count=count+1
        if t % S_xyz == 0:
            index = int(t / S_xyz)
            r_xyz[index] = r[t]
            Ekin[index] = half_m * np.sum(p[t] * p[t], axis=1)
    H_bar=H_bar/count
    T_bar=T_bar/count
    P_bar=P_bar/count
    print(H_bar,T_bar,P_bar,0,0)
    return r_xyz, Ekin
#############################################################################################################
start_time = time.time()
#read from parameters.txt
data = np.loadtxt('parameters.txt', dtype={'names': ('key', 'value'), 'formats': ('S5', 'f8')})
n = int(data[data['key'] == b'n']['value'][0])
N=n**3
m = data[data['key'] == b'm']['value'][0]
a = data[data['key'] == b'a']['value'][0]
T_0 = data[data['key'] == b'T_0']['value'][0]
R = data[data['key'] == b'R']['value'][0]
e = data[data['key'] == b'e']['value'][0]
f = data[data['key'] == b'f']['value'][0]
tau = data[data['key'] == b'tau']['value'][0]
S_o = int(data[data['key'] == b'S_o']['value'][0])
S_d = int(data[data['key'] == b'S_d']['value'][0])
S_out = int(data[data['key'] == b'S_out']['value'][0])
S_xyz = int(data[data['key'] == b'S_xyz']['value'][0])
#execute functions
r_init=initial_state(n,a,N) #initial state for positions
r_x = np.random.uniform(1e-10, 1, N)
r_y = np.random.uniform(1e-10, 1, N)
r_z = np.random.uniform(1e-10, 1, N)
r1 = np.random.choice([1, -1], N)
r2 = np.random.choice([1, -1], N)
r3 = np.random.choice([1, -1], N)
p_init=initial_momentum(r_x,r_y,r_z,r1,r2,r3,N,k_B,T_0,m) #initial kinetic energy and momentum
L=1.22*a*(n-1)*1.05
V_init=potential_energy(N,r_init,e,R,L,f)
E_kin_init = 1/(2*m)*np.sum((p_init[:, 0]*p_init[:, 0] + p_init[:, 1]*p_init[:, 1] + p_init[:, 2]*p_init[:, 2]))
print('N: ',N)
print('T_0: ',T_0)
print('Initial kinetic energy: ',E_kin_init)
print('Initial potential energy: ',V_init)
print('Total initial energy: ',V_init+E_kin_init)
F_init,Pre_init=force(N,r_init,e,R,L,f)#forces and pressure
original_stdout = sys.stdout
file = open("result.txt", "w")
sys.stdout = file
#execute main function
r_xyz,Ekin_i=time_simulation(E_kin_init,V_init,Pre_init,T_0,p_init,F_init,S_o,S_d,S_out,S_xyz,N,tau, r_init, e, R, L, f,m)
sys.stdout = original_stdout
file.close()
#print x,y,z positions
Txyz, Nxyz, _ = r_xyz.shape
with open('avs.xyz', 'w') as f3:
    for t in range(0,Txyz):
        f3.write(f"{len(r_init)}\n")
        f3.write(f"Argon\n")
        for i in range(N):
            if t==0:
                f3.write(f"Ar {r_init[i,0]:.3f} {r_init[i,1]:.3f} {r_init[i,2]:.3f}\n")
            else:
                f3.write(f"Ar {r_xyz[t,i,0]:.3f} {r_xyz[t,i,1]:.3f} {r_xyz[t,i,2]:.3f}\n")
end_time = time.time()
print('Code executed in ',end_time-start_time)
