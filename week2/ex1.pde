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
INITIAL VALUES
Temp = 20
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) = 0!rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
  k = 401
  cp = 385
  rho = 8960
    START(0,0)   { Walk the domain boundary }
    load(Temp) = 0 LINE TO (2,0)  value(Temp) = 100 line TO (2,1) value(Temp) = if(x<1) then 200*x else 200-100*(x-1) line TO (0,1) value(Temp)=0 line TO CLOSE
  REGION 2
  k = 7000
  cp = 400
  rho = 2000
  	START(1,0) load(Temp)=0
    Line to (1,-1) value(Temp) = 400
    LINE TO (1.5,-1) load(Temp) = 0
    LINE TO (1.5,0)
    LINE TO CLOSE
    
TIME 0 TO 1  { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t = endtime
  contour(Temp)
  surface(Temp)
END