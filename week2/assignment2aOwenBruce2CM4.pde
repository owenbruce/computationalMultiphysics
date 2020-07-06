{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=1e-6)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
k
cp
rho 
Tempi
qdotv = 0
INITIAL VALUES
Temp = Tempi
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) + qdotv = rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Iron'       { For each material region }
  k = 80.4
  cp = 450
  rho = 7874
  Tempi = 180
    START(0,0)
    LINE TO (0.1,0)
    line TO (0.1,0.008)
    line TO (0,0.008)
    line TO CLOSE
  REGION 'Bacon'
  k = 0.4
  cp = 4200
  rho = 1000
  Tempi = 4
  	START 'Contact' (0.01,0.008)
    Line to (0.09,0.008)
    LINE TO (0.09,0.028)
    LINE TO (0.01,0.028)
    LINE TO CLOSE
    
TIME 0 TO 10*60  { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by endtime/10 to endtime
  contour(Temp) painted
  surface(Temp)
  
SUMMARY

END