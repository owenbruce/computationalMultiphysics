{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=1e-6)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
baconw=0.02
baconl=0.08
ironw=0.008
ironl = 0.1
r = sqrt(x^2+y^2)
w = 0.01
q0 = 4E6
h = 20
tFluid = 22

k
cp
rho 
Tempi

qdot=-k*grad(Temp)
qdotv
avgIronTopTemp = integral(Temp,'Iron_Top')/integral(1,'Iron_Top')
avgBaconTemp = integral(Temp,'Bacon')/integral(1,'Bacon')
avgBaconMid = bintegral(Temp,'bacon_mid')/bintegral(1,'bacon_mid')

baconConvHeat = tintegral(bintegral(normal(qdot),'bacon_conv'))
totalBaconConvHeat = 4*baconConvHeat
ironConvHeat = tintegral(bintegral(normal(qdot),'iron_conv'))
totalIronConvHeat = 4*ironConvHeat
maserHeat = tintegral(integral(qdotv,'Bacon'))
totalMaserHeat = 4*maserHeat
INITIAL VALUES
Temp = Tempi
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) + qdotv = rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }

  	REGION 'Iron_Top'       { For each material region }
    
  	k = 80.4
  	cp = 450
  	rho = 7874
  	Tempi = 180
  	qdotv = if(AvgIronTopTemp<180) then 7.15E6 else 0
    
    START(0,baconw/2)
    Line to (baconl/2,baconw/2) load(Temp) = -h*(Temp-tFluid)
    LINE TO (ironl/2,baconw/2)
    line TO (ironl/2,baconw/2+ironw) load(Temp) = 0
    line TO (0,baconw/2+ironw)
    line TO CLOSE
    
  	REGION 'Bacon'
    
  	k = 0.4
  	cp = 4200
  	rho = 1000
  	Tempi = 4
    qdotv = q0*exp(-2*r^2/w^2)
    
  	START 'contact' (0,0)
    Line to (baconl/2,0) load(Temp) = -h*(Temp-tFluid)
    LINE TO (baconl/2,baconw/2) load(Temp) = 0
    LINE TO (0,baconw/2)
    LINE TO CLOSE
    
FEATURE
	start 'iron_conv' (baconl/2,baconw/2) line to (ironl/2,baconw/2) line to (Ironl/2, baconw/2+ironw)
    start 'bacon_conv' (baconl/2,0) line to (baconl/2,baconw/2)
    start 'bacon_mid' (0,0) line to (baconl/2,0)
    
TIME 0 TO 600 halt(avgBaconMid>63)  { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for time=0 by 5 to endtime
  contour(Temp) painted
  surface(Temp)
  vector(qdot) norm
  
SUMMARY
report avgBaconTemp
report avgBaconMid
report baconConvHeat
report totalBaconConvHeat
report ironConvHeat
report totalIronConvHeat
report maserHeat
report totalMaserHeat
END