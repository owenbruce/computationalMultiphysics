import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import sympy as sym

startTime = time.time()

FlexCode = """TITLE 'Design Project Owen Bruce'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u (threshold=1e-3)
  v (threshold=1e-3)
SELECT         { method controls }
	ngrid = 1
DEFINITIONS    { parameter definitions }
	MoonRad = 1737.4e3   !Constants
    Mm = 7.346e22	!Moon Mass
    Mr = 104535   !Rocket Mass
    Mp = 11e3	!Payload Mass
    gamma = 1.4
    Rg = 8.314
    p0 = 3*(10e-15)*100e3	!External Pressure
    G = 6.67e-11
    

	Astar = %s !Rocket Parameters
    Ae = %s
    Pt = %s
    Tt = %s
    Mfuel =%s
    Me = %s
    
    
	mdot = Astar*pt/sqrt(Tt)*sqrt(gamma/Rg)*((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)))   !Thrust Equations
    Te = Tt*(1+(gamma-1)/2*Me^2)^(-1)
    pe = pt*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1))
    ve = Me*sqrt(gamma*Rg*Te)
    
    Tfuel = Mfuel/mdot
    
    M = if t<tfuel then Mfuel + Mr + Mp - tfuel*mdot else Mr+Mp
    
	Ft = if t<tfuel then mdot*Ve+(pe-p0)*Ae else (pe-p0)*Ae   !Force Equations
    Fg = -G*M*mm/u^2
    Fnet = Ft + Fg
    a = Fnet/M
INITIAL VALUES
u = MoonRad
EQUATIONS        { PDE's, one for each variable }
	u: dt(u) = v
    v: dt(v) = a
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 3600 halt(u-MoonRad>3e3 or u<(MoonRad-1e3))    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for time = endtime
  history(u) at (0,0)
  history(Fg) at (0,0)
  history(Ft) at (0,0)
  history(if val(u,0,0)<MoonRad+2.9e3 then 0 else if t<tfuel then tintegral(mdot) else mdot*tfuel)	Export Format '#t#b#1' file = 'FuelConsumed.txt'
Summary
	report(if t<tfuel then tintegral(mdot) else mdot*tfuel) as 'Total mass of fuel consumed'
END
"""

FlexFileName = "D:\\EngPhys\\2CM4\\DesignProject\\designProjectScriptOwenBruce2CM4.pde"    #Change to whatever folder and file you want

outcome = 0 #Variable equals 0 if the next value is lower, variable equals 1 if the next value is higher, variable equals 2 if the rocket fails to reach the required height, variable equals 3 if the simulation crashes

#Initial Values for Parameters
Astar = 1
Ae = 10
Pt = 10*10**6
Tt = 5000
Mfuel = 780000

#Constants
gamma = 1.4

#Paramters
initDeltaFuel = 10
deltaFuel = 10
deltaFuelPrevious = 0
fuelArray = np.array([])
fuelConsumedNew = 0
fuelConsumedPrevious = 235364
accuracy = 1
maxAstar = 1
minAstar = 0.1
deltaAstar = 0.01
maxAe = 15
minAe = 1
deltaAe = 0.01
maxPt = 25*10**6
minPt = 2*10**6
deltaPt = 5*10**4
minTt = 2985
maxTt = 10000
deltaTt = 50
simulation = 0
line = 0

while initDeltaFuel > accuracy:
    
    line += 1
    
    while deltaFuel > accuracy:    #Continue to change the same amount while the fuel amount is changing 
        simulation += 1
        #First, we need to solve for Me based on Ae and Astar
        
        w = sym.symbols('w')
    
        f = sym.Eq(Ae/Astar,((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)*w**2/2)**((gamma+1)/(2*(gamma-1)))/w)
        a = sym.solve(f,w)
    
        Me = a[1]
        
        
        with open(FlexFileName, "w") as f:
            print(FlexCode%(Astar,Ae,Pt,Tt,Mfuel,Me), file=f)
        
        #completed = subprocess.run(["ls"], shell=True)
        completed =subprocess.Popen(["FlexPDE7n","-S", FlexFileName])   #Change to FlexPDE6s if that's what you are running
        time.sleep(0.2)
        print('Started')
        completed.wait()
        
        try:
            with open("D:\\EngPhys\\2CM4\\DesignProject\\designProjectScriptOwenBruce2CM4_output\\FuelConsumed.txt") as f:   #Change to whatever folder the output is in
                data=sp.loadtxt(f, skiprows=8)
            fuelArray=data[:,1]
            if(fuelArray[-1]==0):
                outcome = 2
                print("Simulation #",simulation," on line #",line," failed to reach the required height.")
            elif(fuelArray[-1] >= fuelConsumedPrevious):
                outcome = 1
                fuelConsumedNew = fuelArray[-1]
                print("Simulation #",simulation," on line #",line," used more fuel.")
            else:
                outcome = 0
                fuelConsumedNew = fuelArray[-1]
                print("Simulation #",simulation," on line #",line," used less fuel.")
                
        except:
            print("Failed")
            outcome = 3
        
        Mfuel = fuelConsumedNew
        deltaFuel = abs(fuelConsumedNew-fuelConsumedPrevious)
        
        if(outcome==2):
            Ae -= deltaAe
            Astar -= deltaAstar
            Pt -= deltaPt
            Tt -= deltaTt
            deltaAe *= 0.9
            deltaAstar *= 0.9
            deltaPt *= 0.9
            deltaTt *= 0.9
            
        elif(outcome==1):
            deltaAe *= -0.9
            deltaAstar *= -0.9
            deltaPt *= -0.9
            deltaTt *= -0.9
            
        elif(outcome==0):
            
            if(deltaFuel>=deltaFuelPrevious*0.9):
                deltaAe *= 1.1
                deltaAstar *= 1.1
                deltaPt *= 1.1
                deltaTt *= 1.1
                
            else:
                deltaAe *= 0.9
                deltaAstar *= 0.9
                deltaPt *= 0.9
                deltaTt *= 0.9
 
            
        if(Ae+deltaAe<maxAe and Ae+deltaAe>minAe):
            Ae+=deltaAe
            
        if(Astar+deltaAstar<maxAstar and Astar+deltaAstar>minAstar):
            Astar+=deltaAstar
            
        if(Pt+deltaPt<maxPt and Pt+deltaPt>minPt):
            Pt+=deltaPt
            
        if(Tt+deltaTt<maxTt and Tt+deltaTt>minTt):
            Tt+=deltaTt
            
        if(simulation == 1):
            initDeltaFuel = deltaFuel
            
        deltaFuelPrevious = deltaFuel
        fuelConsumedPrevious = fuelConsumedNew
            
        print("Ae =", Ae)
        print("Astar =", Astar)
        print("Pt =", Pt)
        print("Tt =", Tt)
        print("Mfuel =", Mfuel)
        print("Delta fuel =",deltaFuel)
        
        
        
    if(line==5):
        line -= 4
        
    if(line==1):
        deltaAe *= -2
        deltaAstar *= 2
        deltaPt *= 2
        deltaTt *= 2
    
    elif(line==2):
        deltaAe *= 2
        deltaAstar *= -2
        deltaPt *= 2
        deltaTt *= 2
    
    elif(line==3):
        deltaAe *= 2
        deltaAstar *= 2
        deltaPt *= -2
        deltaTt *= 2
        
    elif(line==4):
        deltaAe *= 2
        deltaAstar *= 2
        deltaPt *= 2
        deltaTt *= -2
        
    if(Ae+deltaAe<maxAe and Ae+deltaAe>minAe):
            Ae+=deltaAe
            
    if(Astar+deltaAstar<maxAstar and Astar+deltaAstar>minAstar):
            Astar+=deltaAstar
            
    if(Pt+deltaPt<maxPt and Pt+deltaPt>minPt):
            Pt+=deltaPt
            
    if(Tt+deltaTt<maxTt and Tt+deltaTt>minTt):
            Tt+=deltaTt
            
    simulation = 0
    deltaFuel = 10
        
    
print("Minimized Fuel Consumed: ", fuelConsumedNew)
            
print("Full simulation takes",time.time()-startTime,"s")
