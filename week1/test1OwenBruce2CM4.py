import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import time

startTime = time.time()

FlexCode = """TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
r(threshold=1e-3) = vector(rx,ry)
v(threshold=1e-3) = vector(vx,vy)
 SELECT         { method controls }
ngrid = 1
quadratic
DEFINITIONS    { parameter definitions }
ag = vector(0, -9.81)

thetai=45*pi/180
vi = 30

mfueli = 500
mdry = 100

q = %s
tfuel = mfueli/q

mfuel = if(t<tfuel) then mfueli-q*t else 0
vfuel = 1200
Ftmag = if(t<tfuel) then q*vfuel else 0
Ft = Ftmag*v/magnitude(v)

m = mdry + mfuel

Fg = ag*m


rho = 1.2
area = 2
Cd = 0.8
vwind = vector(0,0)
vrel = v - vwind
Fdmag = 0.5*rho*Cd*area*(magnitude(vrel))^2
Fd = -Fdmag*vrel/magnitude(vrel)

Fnet =  Fg + Ft + Fd
a = Fnet/m

INITIAL VALUES
v = vi*vector(cos(thetai),sin(thetai))
rx = -3000
ry = 500
EQUATIONS        { PDE's, one for each variable }
r: dt(r)=v
v: dt(v) = a

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 100 halt(ry<0)    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for t = 0 by endtime/50 to endtime
history(ry) at (0,0) vs rx
history(v, a) at (0,0)
history(Fd, Ft) at (0,0)
history(rx, ry) at (0,0) PrintOnly Export Format '#t#b#1#b#2' file = '%s_traj.txt'
summary
report eval(rx, 0,0)
report eval(ry, 0,0)
END


"""

FlexFileName = "D:\\EngPhys\\2CM4\\Week1\\test1bOwenBruce2CM4.pde"

import matplotlib.pyplot as plt
xfinal = np.array([])
tfinal = np.array([])
maxx = -3000
flowmin = 1
flowmax = 61
flowstep = 1
FlowRange = sp.arange(flowmin,flowmax,flowstep)

processes=[]

for Flow in FlowRange:
    with open(FlexFileName, "w") as f:
        print(FlexCode%(Flow,Flow), file=f)
    
    #completed = subprocess.run(["ls"], shell=True)
    completed =subprocess.Popen(["FlexPDE7n","-S", FlexFileName])
    time.sleep(0.1)
    processes.insert(0,completed)
    print('Started flow =', round(Flow,2), "kg/s")

for p in processes:
    p.wait()
    
#input("Press enter to continue once processes are done")
time.sleep(2)

for Flow in FlowRange:
    try:
        with open("D:\\EngPhys\\2CM4\\Week1\\test1bOwenBruce2CM4_output\\"+str(Flow)+"_traj.txt") as f:
            data=sp.loadtxt(f, skiprows=8)
        t=data[:,0]
        xd=data[:,1]
        yd=data[:,2]
        plt.plot(xd,yd)
        xfinal = np.append(xfinal,[xd[-1]])
        tfinal = np.append(tfinal,[t[-1]])
        if (xd[-1]+3000) > maxx:
            maxx = xd[-1] + 3000
            maxFlow = Flow
        print('Flow rate of {Flow} kg/s lands at xd = {xdfinal}m after {time} seconds'.format(Flow=round(Flow,2), xdfinal=xd[-1], time=t[-1]))
    except:
        print("missing: ", Flow)


print("")
print("The maximum range is:", maxx, "m with a flow of:", maxFlow, "kg/s.")
print("")
    
df = pd.DataFrame({'Flow Rate': FlowRange,
                   'Landing Position': xfinal,
                   'Flight Time': tfinal,
                   'Range': 3000 + xfinal})    
    

print(df)
plt.title('Trajectory for various flow rates in kg/s')
plt.legend(FlowRange)
plt.show()

print("Total simulation time:", (time.time()-startTime))
