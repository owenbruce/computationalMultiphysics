TITLE 'AdvancedExample1_Regions'     { the problem identification }
COORDINATES cartesian3  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=1E-3)
SELECT         { method controls }
!ngrid = 1
DEFINITIONS    { parameter definitions }
!k = if(Temp<1300) then 100/(11.75+.0235*Temp) else 2.2
tfac=(Temp-273.15)/1000
k=115.8/(7.5408+17.692*tfac+3.6142*tfac^2)+7410.5*tfac^(-5/2)*exp(-16.35/tfac)

rho = 10600
cp = 350

qdotvol = 279e6
R = 12e-3
qdot = -k*grad(Temp)
INITIAL VALUES
Temp = 300

EQUATIONS        { PDE's, one for each variable }
dt(rho*cp*Temp) = qdotvol-div(qdot)
!0 = qdotvol-div(qdot)

! CONSTRAINTS    { Integral constraints }
EXTRUSION
	surface 'back' z=0
    surface 'front' z=R
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
  
    START'edge'(R,0)   load(Temp) = 3000*(Temp-300) arc(center=0,0) angle 90
    load(Temp) = 0 line to (0,0) line to (R,0)

TIME 0 TO 3    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for t =  endtime
  CONTOUR(Temp) painted on z=0
vector(qdot) norm on z=0
summary
report pi*R^2/4
report integral(Temp,1)/integral(1, 1)
report tintegral(integral(qdotvol,1)-surf_integral(normal(qdot),'edge'))/integral(cp*rho,1)+300
END
