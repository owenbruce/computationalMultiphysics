{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=0.0001)
SELECT         { method controls }
!quadratic
DEFINITIONS    { parameter definitions }
k
cp
rho = 4000
qdot = -k*grad(Temp)
rateheat12 = bintegral(dot(qdot, vector(0,1)),'12line')
rateheat23 = bintegral(dot(qdot, vector(1,0)),'23line')
heatflux = line_integral(normal(-qdot),'Boundary_1') + line_integral(normal(-qdot),'Boundary_2') + line_integral(normal(-qdot),'Boundary_3')
INITIAL VALUES
Temp = 20
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) = 0
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Section_1'       { For each material region }
  k = 80!160
  cp = 400
    START 'Boundary_1' (0,0)  value(Temp) = 100
    LINE TO (4,0) load(Temp) = 0
    line TO (4,2)
    line TO (0,2)
    line TO CLOSE
    start '12line' (2,2) line to (3,2)
  REGION 'Section_2'
  k = 600!1200
  cp = 800
  	START 'Boundary_2' (2,2)
    Line to (3,2)
    LINE TO (3,4)
    LINE TO (2,4)
    LINE TO CLOSE
    start '23line' (3,3) line to (3,4)
REGION 'Section_3'
  k = 90!180
  cp = 200
  	START 'Boundary_3' (3,3)
    Line to (5,3) value(Temp) = 0
    LINE TO (5,4) load(Temp) = 0
    LINE TO (3,4)
    LINE TO CLOSE
    
TIME 0 TO 120 { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for t=120
  contour(Temp) painted
  surface(Temp)
  vector(qdot) norm
SUMMARY
report rateheat12
report rateheat23
report val(heatflux,0,0)
END