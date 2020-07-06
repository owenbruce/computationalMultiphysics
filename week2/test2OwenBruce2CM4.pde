{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=0.1)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
k
cp
rho = 4000
qdot = -k*grad(Temp)
INITIAL VALUES
Temp = 20
EQUATIONS        { PDE's, one for each variable }
  div(grad(k*Temp)) = rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Section_1'       { For each material region }
  k = 80
  cp = 400
    START(0,0) value(Temp) = 100
    LINE TO (4,0) load(Temp) = 0
    line TO (4,2)
    line TO (0,2)
    line TO CLOSE
  REGION 'Section_2'
  k = 600
  cp = 800
  	START (2,2)
    Line to (3,2)
    LINE TO (3,4)
    LINE TO (2,4)
    LINE TO CLOSE
REGION 'Section_3'
  k = 90
  cp = 200
  	START (3,3)
    Line to (5,3) value(Temp) = 0
    LINE TO (5,4) load(Temp) = 0
    LINE TO (3,4)
    LINE TO CLOSE
    
TIME 0 TO 120  { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=120
  contour(Temp) painted
  surface(Temp)
  vector(qdot) norm
SUMMARY
END