'"Gravitational Collapse: simulating interactions between N bodies"'

import numpy as np
import foo

initial_positions = np.zeros((N,3))
initial_positions[0]=[1.496e11,0.,0.]     # earth
initial_positions[1]=[0.,0.,0.]           # sun



# Initially store the positions/velocities for each timestep but ultimately unecessary

all_positions = np.zeros((n,N,3))
all_positions[0] = initial_positions 

initial_velocities = np.zeros((N,3))
initial_velocities[0]=[2557.5,29668.52,0.]
initial_velocities[1]=[0.,0.,0.]


all_velocities = np.zeros((n,N,3))
all_velocities[0] = initial_velocities 


# Empty arrays to store values

kinetic_energy = np.zeros(n)
potential_energy = np.zeros(n)
total_energy = np.zeros(n)
all_energy_change = np.zeros(n)

# Functions to find force, velocity, position, kinetic/potential/total energy

def find_force(positions, mass, softening):
    force = np.zeros((N,3))
    for i in range(N):
        for j in range(N):
            if i!=j:
                d = positions[i]-positions[j]
                modulus_d2 = d[0]**2+d[1]**2+d[2]**2    # squared modulus d
                F = -G*mass[i]*mass[j]*d/(modulus_d2+softening**2)**1.5
                force[i] += F
    return force
   
def find_velocities(velocities, dt, force, mass):
    new_velocities = np.zeros((N,3))
    for i in range(N):
        new_velocities[i] = velocities[i] + force[i]*dt/mass[i]
    return new_velocities

def find_positions(positions, dt, new_velocites):
    new_positions = np.zeros((N,3))
    for i in range(N):
        new_positions[i] = positions[i] + new_velocities[i]*dt
    return new_positions

def calculate_kinetic_energy(mass,velocities):
    total_kinetic_energy = 0
    for i in range(N):
        modulus_v = np.sqrt(velocities[i,0]**2+velocities[i,1]**2+velocities[i,2]**2)
        k_e = 0.5*mass[i]*modulus_v**2
        total_kinetic_energy += k_e
    return total_kinetic_energy

def calculate_potential_energy(mass, positions, softening):
    total_potential_energy = 0
    for i in range(N):
        for j in range(N):
            if i!=j:
                d = positions[i]-positions[j]
                modulus_d2 = d[0]**2+d[1]**2+d[2]**2
                p_e = -G*mass[i]*mass[j]*0.5/np.sqrt(modulus_d2+softening**2)  # *0.5 because should only count pe from each pair of particles once
                total_potential_energy += p_e
    return total_potential_energy

def calculate_total_energy(total_potential_energy, total_kinetic_energy):
    total_energy = total_potential_energy + total_kinetic_energy
    return total_energy

current_positions = initial_positions
current_velocities = initial_velocities


for i in range(n):

    force = find_force(current_positions, m, eps)
    # print 'F = ', force
    print 'v = ', current_velocities
    print 'r = ',current_positions
    new_velocities = find_velocities(current_velocities, dt, force, m)
    current_velocities = new_velocities
    new_positions = find_positions(current_positions, dt, new_velocities)
    all_positions[i] = new_positions
    current_positions = new_positions
    
    
#    total_potential_energy = calculate_potential_energy(m, current_positions, eps)
#    potential_energy[i] = total_potential_energy
#    total_kinetic_energy = calculate_kinetic_energy(m, current_velocities)
#    kinetic_energy[i] = total_kinetic_energy 
#    energy = calculate_total_energy(total_potential_energy, total_kinetic_energy)
#    total_energy[i] = energy

# Stop program if change in total energy not conserved within 1%

#    all_energy_change[i] = abs(total_energy[0] - total_energy[-1])/total_energy[0]
#    if energy_change >= max_tolerance:
#        break

#print 'ke', kinetic_energy
#print 'pe', potential_energy
#print 'e', total_energy
#print all_energy_change

# Plot graphs of KE vs t, PE vs t, E vs t and sun-earth r vs t

#t=np.zeros(n)    # need to plot initial positions at t=0
#ni = np.arange(1, n+1)

#for i in range(n):              # better way to do time array??
#    t[i] = ni[i]*dt

#py.figure()
#py.plot(t, kinetic_energy)
#py.plot(t, potential_energy)
#py.plot(t, total_energy)
#py.xlabel('Time')               # units?
#py.ylabel('Energy (J)')
#py.title('Title')

#py.figure()
#py.plot(t, all_energy_change)
#py.xlabel('Time')               # units?
#py.ylabel('Label')
#py.title('Title')

#py.figure()
#py.plot(t, sun_earth_r)
#py.xlabel('Time')               # units?
#py.ylabel('Separation (km)')    # units?
#py.title('Title')

# plot radius of particle as function of time??

#py.show()

# Other things to think about:
# Softening?
# Initial velocities need to be set half a timestep before initial position
# Fij = -Fji?







