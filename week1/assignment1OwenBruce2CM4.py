import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import time

startTime = time.time()

FlexCode = """TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
 xd(threshold=0.01)
 vx(threshold=0.01)
 yd(threshold=0.01)
 vy(threshold=0.01)
SELECT         { method controls }
ngrid=1
DEFINITIONS    { parameter definitions }
vi = 20
thetai = %s*pi/180
mi=1000
q=50
vfuel=2000
g=9.81

mfueli = 0.8*mi
mDry=mi-mfueli

tfuel = mfueli/q
mfuel = if(t<tfuel) then mfueli-q*t else 0

Ft = if(t<tfuel) then q*vfuel else 0
m = mDry + mfuel

r = sqrt(xd^2+yd^2)
vwindx = -25
vwindy = 10
vrelx = vx-vwindx
vrely = vy-vwindy
vrel = sqrt(vrelx^2+vrely^2)
v= sqrt(vx^2+vy^2)+1e-6

Area = 5
CD = 0.6
rhoAir = 2
Fdmag = 0.5*rhoAir*CD*Area*vrel^2

Fdx = Fdmag*(-vrelx/vrel)
Fdy = Fdmag*(-vrely/vrel)

Ftx = Ft*vx/v
Fty = Ft*vy/v
ay = (Fty+Fdy)/m-g
ax = (Ftx+Fdx)/m
INITIAL VALUES
xd=0
vx=vi*cos(thetai)
vy=vi*sin(thetai)

EQUATIONS        { PDE's, one for each variable }
  yd: dt(yd)=vy
  vy: dt(vy)=ay
   xd: dt(xd)=vx
  vx: dt(vx)=ax
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 200 halt(yd<0)   { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by 1 to endtime
  history(yd,vy,ay) at (0,0)
  history(yd) at (0,0) vs xd
  history(xd,yd,v) at (0,0) PrintOnly Export Format '#t#b#1#b#2#b#3' file = '%s_traj.txt'
  SUMMARY
  report val(xd,0,0)
END

"""

FlexFileName = "D:\\EngPhys\\2CM4\\Week1\\assignment1OwenBruce2CM4.pde"    #Change to whatever folder and file you want

import matplotlib.pyplot as plt
xfinal = np.array([])
tfinal = np.array([])
keimpact = np.array([])
maxkei = 0
maxx = 0
anglemin = 0
anglemax = 181
anglestep = 1
AngleRange = sp.arange(anglemin,anglemax,anglestep)

processes=[]

while anglestep>=0.01:
    for Angle in AngleRange:
        with open(FlexFileName, "w") as f:
            print(FlexCode%(round(Angle,2),round(Angle,2)), file=f)
        
        #completed = subprocess.run(["ls"], shell=True)
        completed =subprocess.Popen(["FlexPDE7n","-S", FlexFileName])   #Change to FlexPDE6s if that's what you are running
        time.sleep(0.2)
        processes.insert(0,completed)
        print('Started ', round(Angle,2))
    
    for p in processes:
        p.wait()
        
    #input("Press enter to continue once processes are done")  #Press enter to resume after processes are finished
    time.sleep(3) #This is approximately how long mine needs to sleep for at maximum, it may vary from computer to comptuer
    
    for Angle in AngleRange:
        try:
            with open("D:\\EngPhys\\2CM4\\Week1\\assignment1OwenBruce2CM4_output\\"+str(round(Angle,2))+"_traj.txt") as f:   #Change to whatever folder the output is in
                data=sp.loadtxt(f, skiprows=8)
            t=data[:,0]
            xd=data[:,1]
            yd=data[:,2]
            vy=data[:,3]
            plt.plot(xd,yd)
            xfinal = np.append(xfinal,[xd[-1]])
            tfinal = np.append(tfinal, [t[-1]])
            keimpact = np.append(keimpact,[vy[-1]**2*0.5*200])
            if keimpact[-1] > maxkei:
                maxkei = keimpact[-1]
                keangle = round(Angle,2)
            if xd[-1] > maxx:
                maxx = xd[-1]
                maxangle = round(Angle,2)
            print('{angle}\N{DEGREE SIGN} lands at xd = {xdfinal} after {time} seconds'.format(angle=round(Angle,2), xdfinal=xd[-1], time=t[-1]))
        except:
            print("missing: ", Angle)
    
    anglemin = maxangle - anglestep
    anglemax = maxangle + anglestep
    anglestep = anglestep/10
    AngleRange = sp.arange(anglemin, anglemax, anglestep)

print("")
print("The maximum x-displacement is:", maxx, "m with an angle of:", maxangle)
print("The maximum impact kinetic energy is:", maxkei, "N with an angle of:", keangle)
print("")
    
#df = pd.DataFrame({'Angle': AngleRange,
                   #'Landing Position': xfinal,
                   #'Flight Time': tfinal,
                   #'Impact Kinetic Energy': keimpact})    
    
#print(df)
plt.title('Trajectory for various launch angles')
plt.legend(AngleRange)
plt.show()

print("Full simulation takes",time.time()-startTime,"s")
