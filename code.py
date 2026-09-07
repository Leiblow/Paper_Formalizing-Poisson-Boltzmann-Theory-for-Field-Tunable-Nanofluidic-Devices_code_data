import numpy as np
from scipy.integrate import solve_ivp
import time

def EDL_calculation_nanoslit(h,l_D,l_GC,N_solve=2000,N_con=2000,alpha=2,phi_m_inital=1e-9):#in nanometer
    # physical constants
    e = 1.602e-19           
    kB = 1.381e-23          
    T = 298                 
    epsilon_r = 78.5        
    epsilon_0 = 8.854e-12   
    epsilon = epsilon_r * epsilon_0
    
    #converting from nm to m
    h,l_D,l_GC=h*1e-9,l_D*1e-9,l_GC*1e-9

    z_span=(h/2,h)
    z_span_1 = (h/2,h-10*l_D)
    z_span_2=(h-10*l_D,h)

    Sigma = 2*epsilon*kB*T/e**2/l_GC # [m^-2]
    # target surface field
    target_phi_prime_h = e**2 * Sigma / (epsilon * kB * T)

    # Poisson–Boltzmann equations
    def pb_eq(z, y):
        return [y[1], np.sinh(y[0]) / l_D**2]

    def solve_phi(phi_m):
        if 10*l_D<h/2:
            y0 = [phi_m, 0.0]  # φ(h/2), φ′ = 0
            sol_1 = solve_ivp(pb_eq, z_span_1, y0, method='RK45', t_eval=np.linspace(h/2,h-10*l_D, N_solve//2))
            y_1=[sol_1.y[0][-1],sol_1.y[1][-1]]
            sol_2 = solve_ivp(pb_eq, z_span_2, y_1, method='RK45', t_eval=np.linspace(h-10*l_D,h, N_solve//2))
            return np.concatenate([sol_1.t,sol_2.t]),np.concatenate([sol_1.y[0],sol_2.y[0]]),np.concatenate([sol_1.y[1],sol_2.y[1]])  # z, φ, φ′
        else:
            y0 = [phi_m, 0.0]  # φ(h/2), φ′ = 0
            sol = solve_ivp(pb_eq, z_span, y0, method='RK45', t_eval=np.linspace(h/2, h, N_solve))
            return sol.t, sol.y[0], sol.y[1]  # z, φ, φ′

    time_cost=time.time()
    phi_m_temp=phi_m_inital+0

    #The first round of iteration
    _, _, phi_prime = solve_phi(phi_m_temp)
    delta_log_phi_prime_rate=(-np.log(phi_prime[-1])+np.log(target_phi_prime_h))/max(abs(np.log(target_phi_prime_h))+1e-9,abs(np.log(phi_prime[-1]))+1e-9)
    sign_old=np.sign(delta_log_phi_prime_rate)
    sign_new=np.sign(delta_log_phi_prime_rate)
    delta_log_phi=alpha*np.tanh(delta_log_phi_prime_rate)
    phi_m_temp=phi_m_temp*np.exp(delta_log_phi)
    

    for n in range(1,N_con):
        _, _, phi_prime = solve_phi(phi_m_temp)
        #adjust according to the relative error
        delta_log_phi_prime_rate=(-np.log(phi_prime[-1])+np.log(target_phi_prime_h)+1e-15)/max(abs(np.log(target_phi_prime_h))+1e-9,abs(np.log(phi_prime[-1]))+1e-9)
        sign_old=sign_new+0
        if abs((phi_prime[-1]-target_phi_prime_h)/target_phi_prime_h) < 1e-6:
            phi_m_solution = phi_m_temp+0
            print(f"φ_m = {phi_m_solution:.4f}")
            break
        else:
            sign_new=np.sign(delta_log_phi_prime_rate)
            alpha=alpha*2**((int(sign_new*sign_old)-1)/2)
            delta_log_phi=np.tanh(delta_log_phi_prime_rate)*alpha
            delta_log_phi=8*np.tanh(delta_log_phi/8)
            phi_m_temp=phi_m_temp*np.exp(delta_log_phi)
            
    else:
        phi_m_solution = phi_m_temp+0
        print(f"φ_m = {phi_m_solution:.4f}")
    time_cost=time.time()-time_cost
    print(f"Iteration rounds: {n:.0f}")
    delta_numerical=100*abs((phi_prime[-1]-target_phi_prime_h)/target_phi_prime_h)
    print(f"The relateive error for φ_m={delta_numerical:.8f}")
    # semi-channel solution
    z_half, phi_half, phi_prime_half = solve_phi(phi_m_solution)
    print('Time (s):',time_cost)
    # full-channel solution
    z_full = np.concatenate([h - z_half[::-1],z_half])
    phi_full = np.concatenate([phi_half[::-1], phi_half])
    phi_prime_full = np.concatenate([-phi_prime_half[::-1], phi_prime_half])  

    return {'z':z_full,'phi':phi_full,'E':phi_prime_full}