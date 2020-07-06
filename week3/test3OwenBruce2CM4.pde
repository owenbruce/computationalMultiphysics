{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=1e-6)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
Rb = 0.04
Rm = 0.03

r = sqrt(x^2+y^2)
w = 0.03
q0 = 7E6
h = 15
tFluid = 20

k
gamma
cp
rho 
Tempi = -20

qdot=-k*grad(Temp)
qdotv = q0*gamma*exp(-r^2/w^2)

avgMeatTemp = integral(Temp,'MysteryMeat')/integral(1,'MysteryMeat')
avgBreadTemp = integral(Temp,'Breading')/integral(1,'Breading')
avgMeatTempCheck = tintegral(integral(qdotv,'MysteryMeat')+bintegral(normal(qdot),'FillingConv'))/integral(cp*rho,'MysteryMeat')+Tempi
INITIAL VALUES
Temp = Tempi
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) + qdotv = rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
    
  	REGION 'Breading'
    
  	k = 0.02
  	cp = 1100
  	rho = 200
    gamma = 0.03
    
    start(Rb,0) load(Temp) = -h*(Temp-tFluid) 
    arc(center=0,0) angle 0.1 load(Temp) = 0
    line to (0,0)
    line to close
    
    REGION 'MysteryMeat'       { For each material region }
    
  	k = 0.4
  	cp = 4200
  	rho = 1000
  	gamma = 0.8
    
    start(Rm,0)
    arc(center=0,0) angle 0.1
    line to (0,0)
    line to close
    
FEATURE
	start 'FillingConv' (Rm,0)arc(center=0,0) angle 0.1
    
TIME 0 TO 60!6000 halt(avgBreadTemp>100) !commented out part 1.c)
MONITORS         { show progress }
PLOTS            { save result displays }
  for time=endtime!351.825 by 0.0001 to endtime !commented out part 1.c)
  contour(Temp) painted
  surface(Temp)
  vector(qdot) norm
  
SUMMARY
report avgMeatTemp
report avgBreadTemp
report avgMeatTempCheck
report t
END