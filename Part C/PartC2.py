#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 16:19:32 2020

@author: amandaseger
"""
from fenics import *
from dolfin import *
import numpy as np
import random

##############################################################################

#class DirichletBoundary(SubDomain):
 #  def inside(self, x, on_boundary):
  #     return on_boundary
####################################################
   
# Create mesh and define function space
mesh = Mesh("circle_finest.xml")
# Construct the finite element space
V = VectorFunctionSpace(mesh, 'P', 1)

# Define parameters:
T = 1000
dt = 0.5
alpha = 0.4
beta = 2
gamma = 0.8 
delta1 = 1
delta2 = 1    
   
# Class representing the intial conditions
class InitialConditions(UserExpression): 
    def eval(self, values, x):
        r = random.uniform(0, 1)
        values[0] = (1/2)*(1-r)
        values[1] = (1/4)+ ((1/2)*r)
    def value_shape(self): 
        return (2,)

    
# Define initial condition
indata = InitialConditions(degree=2) 
u0 = Function(V)
u0 = interpolate(indata , V)

# Test and trial functions
u = TrialFunction(V)
v = TestFunction(V)
  

# Create bilinear and linear forms
a0 = u[0]*v[0]*dx + 1/2*dt*delta1*inner(grad(u[0]), grad(v[0]))*dx -1/2*dt*u[0]*v[0]*dx   
a1 = u[1]*v[1]*dx + 1/2*dt*delta2*inner(grad(u[1]), grad(v[1]))*dx + gamma*dt*1/2*u[1]*v[1]*dx
      
L0 = u0[0]*v[0]*dx -(1/2*dt*delta1*inner(grad(u0[0]), grad(v[0]))*dx) -\
    (dt*((u0[0]*u0[1])/(u0[0]+alpha)+u0[0]*u0[0])*v[0]*dx) + 1/2*dt*u0[0]*v[0]*dx  
    
L1 = u0[1]*v[1]*dx -(1/2*dt*delta2*inner(grad(u0[1]), grad(v[1]))*dx) -\
    (dt*(-(beta*u0[1]*u0[0])/(u0[0]+alpha))*v[1]*dx) - gamma*1/2*dt*u0[1]*v[1]*dx 
    
a = a0+a1
L = L0+L1

#Set upp boundary condition
bc = [] #Neumann

#solve(a==L, u, bc) #Could just solve directly

# Assemble matrix
#A = assemble(a)
#b = assemble(L)
#bc.apply(A)

# Set an output file
out_file = File("resultsPartC2finest/C2solution.pvd", "compressed")

# Set initial condition
u = Function(V)
u.assign(u0)

t = 0

out_file << (u,t)

u_initial = Function(V)
u_initial.assign(u0)

t_save = 0
num_samples = 20

#Initial population rate
Mprey= []
Mpred = []
# Define the integrals
M0 = u[0] * dx
M1 = u[1] * dx

time = []
# Time-stepping
while t <= (T):
# assign u0
    time.append(t)
    u0.assign(u)
    
    # Assemble vector and apply boundary conditions
    
    A = assemble(a)
    b = assemble(L)
    #bc.apply(b)
    
    # Solve linear system
    #solve(A, u.vector(), b, "bicgstab", "default")
    #solve(A, u.vector(), b, "gmres", "ilus")
    solve(A, u.vector(), b, "lu")
    #solve(a==L,u,bc)
    
    t_save += dt
    #Calculate population rate
    # compute the functional
    population_u = assemble(M0) 
    population_v = assemble(M1)

    
    Mprey.append(population_u)
    Mpred.append(population_v)

    if t==50:
        print('Saving at 50')
        plot(u)
        out_file << (u,t)

    if t==100:
        print('Saving at 100')
        plot(u)
        out_file << (u,t)
        
    if t==500:
        print('Saving at 500')
        plot(u)
        out_file << (u,t)
    if t == 1000:
        print('Saving at 1000')
        plot(u)
        out_file << (u, t)
        
    # Move to next interval and adjust boundary condition
    t += dt
np.savetxt('C2Mprey.txt',Mprey)
np.savetxt('C2Mpred.txt',Mpred)
np.savetxt('C2time.txt',time)
        
    

