#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import random

num_vortices = 2
R = 1
dt = 0.00001
total_time = 5

# random point or custom input
random_init = True
position_init = []
circulations = []


# numerical method: RK4 or Euler
num_method = 'RK4'

plot_trajectory = True
generate_poincare = False
calculate_lyapunov = False

def Hamiltonian(R, circulations):
    x_vars = [sp.symbols(f'x{i}') for i in range(1,len(circulations)+1)]
    y_vars = [sp.symbols(f'y{i}') for i in range(1,len(circulations)+1)]
    
    omega_0 = 1/(R**2)
    vortex_interaction = 0
    boundary_interaction = 0
    
    
    
    omega_int = 0.1 * omega_0
    
    # computing vortex-vortex interaction term via summation
    for i in range(len(circulations)):
        for j in range(i+1, len(circulations)):
            
            dx = x_vars[i] - x_vars[j]
            dy = y_vars[i] - y_vars[j]
            
            vortex_interaction += circulations[i] * circulations[j] * sp.log((dx**2 + dy**2)/R**2)
           
    # multiply by appropriate factor
    vortex_interaction *= (-(R**2) * omega_int)/2
    
    # compute vortex-boundary interaction term via summation
    
    for i in range(len(circulations)):
        boundary_interaction += (circulations[i]**2) * sp.log(1-((x_vars[i]**2 + y_vars[i]**2)/R**2))
        
    # multiply by appropriate factor
    boundary_interaction *= ((R**2) * omega_0)/2
    
    
    return vortex_interaction + boundary_interaction, x_vars, y_vars
    
def velocity_coords(R, circulations):
    
    H, x_vars, y_vars = Hamiltonian(R, circulations)
    
    coords_change = []
    
    for i in range(len(circulations)):
        
        vx = (1/circulations[i]) * sp.diff(H, y_vars[i])
        vy = (1/circulations[i]) * -sp.diff(H, x_vars[i])
        
        vx_func = sp.lambdify((*x_vars, *y_vars), vx)
        vy_func = sp.lambdify((*x_vars, *y_vars), vy)
        
        coords_change.append((vx_func, vy_func))
        
    return coords_change, x_vars, y_vars

def rk4_step(x, y, dt, coords_change):
    
    
    new_x, new_y = np.zeros(len(x)), np.zeros(len(x))
    
    for i in range(len(x)):
        
        vx, vy = coords_change[i]
        
        k1_x = vx(*x, *y)
        k1_y = vy(*x, *y)
        
        k2_x = vx(*(x + 0.5 * dt * k1_x), *(y + 0.5 * dt * k1_y))
        k2_y = vy(*(x + 0.5 * dt * k1_x), *(y + 0.5 * dt * k1_y))
        
        k3_x = vx(*(x + 0.5 * dt * k2_x), *(y + 0.5 * dt * k2_y))
        k3_y = vy(*(x + 0.5 * dt * k2_x), *(y + 0.5 * dt * k2_y))
        
        k4_x = vx(*(x + dt * k3_x), *(y + dt * k3_y))
        k4_y = vy(*(x + dt * k3_x), *(y + dt * k3_y))
        
        new_x[i] = x[i] + (dt/6)*(k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        new_y[i] = y[i] + (dt/6)*(k1_y + 2 * k2_y + 2 * k3_y + k4_y)
        
    return new_x, new_y
    
def euler_step(x, y, dt, coords_change):
    
    new_x, new_y = np.zeros(len(x)), np.zeros(len(y))
    
    for i in range(len(x)):
        
        vx, vy = coords_change[i]
        
        new_x[i] = x[i] + dt*vx(*x, *y)
        new_y[i] = y[i] + dt*vy(*x, *y)
        
    return new_x, new_y

def trajectory(x_init, y_init, dt, total_time, coords_change):
    
    x_trajectories = [np.array(x_init)]
    y_trajectories = [np.array(y_init)]
    
    x, y = np.array(x_init), np.array(y_init)
    
    for _ in range(int(total_time/dt)):
        
        x, y = rk4_step(x, y, dt, coords_change)
        x_trajectories.append(x)
        y_trajectories.append(y)
        
    return np.array(x_trajectories), np.array(y_trajectories)

def trajectory_euler(x_init, y_init, dt, total_time, coords_change):
    
    x_trajectories = [np.array(x_init)]
    y_trajectories = [np.array(y_init)]
    
    x, y = np.array(x_init), np.array(y_init)
    
    for _ in range(int(total_time/dt)):
        
        x, y = euler_step(x, y, dt, coords_change)
        x_trajectories.append(x)
        y_trajectories.append(y)
        
    return np.array(x_trajectories), np.array(y_trajectories)

def random_point(R = 1.0):
    
    theta = random.uniform(0,2*np.pi)
    r = random.uniform(0,0.9)
    x = (R*np.sqrt(r)) * np.cos(theta)
    y = (R*np.sqrt(r)) * np.sin(theta)
    
    return (round(x,2),round(y,2))

def poincare_section(x_traj, y_traj, velocities, dt):
   
    
    poincare_y = []
    poincare_vy = []

    for i in range(len(x_traj[0])):  # Iterate over each vortex
        for t in range(1, len(x_traj)):  # Iterate over each time step
            if x_traj[t-1][i] * x_traj[t][i] < 0:  # Check if x crosses zero
                # Interpolate y-coordinate at the crossing point
                y_cross = y_traj[t-1][i] + (y_traj[t][i] - y_traj[t-1][i]) * abs(x_traj[t-1][i]) / (abs(x_traj[t-1][i]) + abs(x_traj[t][i]))
                
                # Calculate the y-velocity at this point
                vx_func, vy_func = velocities[i]
                vy_cross = vy_func(*x_traj[t], *y_traj[t])  # Use the current velocity field
                
                # Record the y-coordinate and y-velocity
                
                
                poincare_y.append(y_cross)
                poincare_vy.append(vy_cross)

    return poincare_y, poincare_vy
    
def lyapunov_exponent(x_init, y_init, dt, total_time, coords_change, perturbation = 1e-5):
    
    perturbed_x_init = x_init.copy()
    perturbed_x_init[0] += perturbation
    perturbed_y_init = y_init.copy
    
    distances = []
    
    x, y = np.array(x_init), np.array(y_init)
    x_perturbed, y_perturbed = np.array(perturbed_x_init), np.array(perturbed_y_init)
    
    for _ in range(int(total_time/dt) + 1):
        
        x, y = rk4_step(x, y, dt, coords_change)
        x_perturbed, y_perturbed = rk4_step(x, y, dt, coords_change)
        
        dist = np.sqrt(np.sum((x-x_perturbed)**2 + (y-y_perturbed)**2))
        distances.append(dist)
        
        if dist > perturbation*1e5:
            scale_factor = perturbation/dist
            x_perturbed = x + scale_factor * (x_perturbed - x)
            y_perturbed = y + scale_factor * (y_perturbed - y)
            
    distances = np.array(distances)
    times = np.arange(0,total_time,dt)
    return times, distances


    
if random_init:
    
    if not circulations:
        circulations = [random.choice([5.0,-5.0]) for i in range(num_vortices)]
    
    if not position_init:
        position_init = [random_point(R) for i in range(num_vortices)]
        
print(circulations)
print(position_init)        
        
x_init = [position_init[i][0] for i in range(len(position_init))]
y_init = [position_init[i][1] for i in range(len(position_init))]

velocities = velocity_coords(R, circulations)

if num_method == 'RK4':
    x_traj, y_traj = trajectory(x_init, y_init, dt, total_time, velocities[0])
    
if num_method == 'Euler':
    x_traj, y_traj = trajectory_euler(x_init, y_init, dt, total_time, velocities[0])

if plot_trajectory:
    
    for i in range(len(circulations)):
        line, = plt.plot(x_traj[:, i], y_traj[:, i], label = f'Vortex {i+1}', zorder = 1)
        plt.plot(x_traj[0][i], y_traj[0][i], marker = 's', markersize=6, markeredgecolor = 'black', markerfacecolor = line.get_color(), zorder = 2)
        plt.plot(x_traj[-1][i], y_traj[-1][i], marker = 'o', markersize = 6, markeredgecolor = 'black', markerfacecolor = line.get_color(), zorder = 2)

    plt.plot(R * np.cos(np.linspace(0, 2 * np.pi, 100)), R * np.sin(np.linspace(0, 2 * np.pi, 100)), color = 'black')

    plt.xlim(R*(-1.03),R*(1.03))
    plt.ylim(R*(-1.03),R*(1.03))
    plt.xlabel("x (units of R)")
    plt.ylabel("y (units of R)")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    
if generate_poincare:
    
    poincare_y, poincare_vy = poincare_section(x_traj, y_traj, velocities[0], dt)
    
    plt.scatter(poincare_y, poincare_vy, s=1, color='black')
    plt.xlabel('y (units of R)')
    plt.ylabel('dy/dt (units of R/s)')
    plt.show()
    
if calculate_lyapunov:
    
    times, distances = lyapunov_exponent(x_init, y_init, dt, total_time, velocities[0])
    log_dist = np.log(distances)

    slope, intercept = np.polyfit(times, log_dist, 1)
    print(f"Estimated Lyapunov Exponent: {slope}")
    plt.plot(times, log_dist)
    plt.plot(times, slope * times + intercept, 'r--')
    plt.xlabel('time')
    plt.ylabel('ln(distances)')
    plt.show()        
        

