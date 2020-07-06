import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import time

startTime = time.time()

FlexCode = """TITLE 'Static Beam Temp Change'
COORDINATES
cartesian3
VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

Select
quadratic

DEFINITIONS
mag=0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))

Lx=0.1/1000
Lya=.005/1000
Lyg=%s
Lz=.02/1000
E
nu
rho

c = 2*v/x^2

vtip = 1/2*c*Lx^2

G=E/(2*(1+nu))

DeltaTemp = 5

alpha

alphax = alpha
alphay = alpha
alphaz = alpha

!Stiffness matrix components (isotropic material)
C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C22 = C11
C33 = C11

C12 = E*nu/(1+nu)/(1-2*nu)
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
gyz=(dy(w)+dz(v))
gxz=(dz(u)+dx(w))

!mechanical strain
exm=ex - alphax*DeltaTemp
eym=ey - alphay*DeltaTemp
ezm=ez - alphaz*DeltaTemp

!!Stress via Hooke's law
!Axial Stress
sx = C11*exm+C12*eym+C13*ezm
sy = C21*exm+C22*eym+C23*ezm
sz = C31*exm+C32*eym+C33*ezm
!Shear stress
sxy=G*gxy
sxz=G*gxz
syz=G*gyz

EQUATIONS
!FNet = 0; !apply gravity in the direction you want a mode just to get it started
u:	dx(sx)+dy(sxy)+dz(sxz) = 0
v:	dx(sxy)+dy(sy)+dz(syz) = 0
w:	dx(sxz)+dy(syz)+dz(sz) = 0 

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES
surface 'bottom'
		load(u)=0
		load(v)=0 
		load(w)=0
surface 'top'
		load(u)=0
		load(v)=0 
		load(w)=0
        
    Region 'Gold'
    E = 79E9
    nu = 0.4
    alpha = 14.2E-6
    rho = 19300
    START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0 
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		load(u)=0
		load(v)=0 
		load(w)=0
		LINE TO (Lx,Lyg) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,Lyg) !x=0 surface
		value(u)=0
		value(v)=0 
		value(w)=0
		LINE TO CLOSE

  	REGION 'Aluminum'
    E = 68.9E9
    nu = 0.33
    alpha = 2.32e-5
    rho = 2700
    	START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0 
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		load(u)=0
		load(v)=0 
		load(w)=0
		LINE TO (Lx,-Lya) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,-Lya) !x=0 surface
		value(u)=0
		value(v)=0 
		value(w)=0
		LINE TO CLOSE
        
PLOTS
	!grid(x+u, y+v, z+w)
    elevation(vtip) from(0,0,0) to (Lx,0,0) PrintOnly Export Format '#t#b#1' file = '%s_thick.txt'
SUMMARY
	!report val(ex,Lx,0,0)
    !report val(ey,Lx,0,0)
    !report val(ez,Lx,0,0)
    !report val(c,Lx/2,-Lya,Lz/2)
end

"""

FlexFileName = "D:\\EngPhys\\2CM4\\Week6\\assignment6_4OwenBruce2CM4.pde"    #Change to whatever folder and file you want

import matplotlib.pyplot as plt
vtip = np.array([])
maxThick = 0
maxvtip = 0
thickmin = 10**-6
thickmax = (10**-4)
thickstep = (10**-6)
ThickRange = sp.arange(thickmin,thickmax,thickstep)

processes=[]

while thickstep>=10**-10:
    for Thick in ThickRange:
        with open(FlexFileName, "w") as f:
            print(FlexCode%(Thick,Thick), file=f)
        
        #completed = subprocess.run(["ls"], shell=True)
        completed =subprocess.Popen(["FlexPDE7n","-S", FlexFileName])   #Change to FlexPDE6s if that's what you are running
        time.sleep(0.2)
        processes.insert(0,completed)
        print('Started ', Thick)
    
    for p in processes:
        p.wait()
        
    #input("Press enter to continue once processes are done")  #Press enter to resume after processes are finished
    time.sleep(240) #This is approximately how long mine needs to sleep for at maximum, it may vary from computer to comptuer
    
    for Thick in ThickRange:
        try:
            with open("D:\\EngPhys\\2CM4\\Week6\\assignment6_4OwenBruce2CM4_output\\"+str(Thick)+"_thick.txt") as f:   #Change to whatever folder the output is in
                data=sp.loadtxt(f, skiprows=8)
            vtip=data[:,1]
            if vtip[-1] > maxvtip:
                maxvtip = vtip[-1]
                maxThick = Thick
            print('With a Gold Thikness of {Thick}m a tip displacement of {vtip}m'.format(Thick=Thick, vtip=vtip[-1]))
        except:
            print("missing: ", Thick)
    
    thickmin = maxThick - thickstep
    thickmax = maxThick + thickstep
    thickstep = thickstep/10
    ThickRange = sp.arange(thickmin, thickmax, thickstep)

print("")
print("The maximum vtip displacement is", maxvtip, "m with a gold thickness of:", maxThick)
print("")
    

print("Full simulation takes",time.time()-startTime,"s")
