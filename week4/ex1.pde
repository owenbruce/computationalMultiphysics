TITLE 'New Problem'     
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  V
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
E = -grad(V)
ep0 = 8.85E-12
rho0 = 1
rho = 0
chi_e = 0
epsilon = ep0*(1+chi_e)
D = epsilon*E

r = sqrt((x-0.5)^2+(y-0.5)^2)
RadCharge = 0.2
Etheory = if(r<RadCharge) then rho0*r/(2*ep0) else RadCharge^2*rho0/(2*r*ep0)
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(D) = rho
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(1,0.5) Value(V) = 0 arc(center = 0.5,0.5) angle = 360
    region 2
    rho = rho0
    start(0.5+RadCharge,0.5) arc(center = 0.5,0.5) angle = 360
    region 3 
    chi_e = 10
    start(0.8,0.5) arc(center = 0.5,0.5) angle = 30
    line to (0.5+0.4*cos(Pi/6), 0.5+0.4*sin(Pi/6))
    arc(center = 0.5,0.5) angle = -30
    line to close
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(V) painted
  vector(E)
  vector(D)
  surface(V)
  elevation(V,dot(E,if(x>0.5) then vector(1,0) else vector(-1,0)), Etheory) from (0,0.5) to (1,0.5)
  elevation(abs(dot(E,if(x>0.5) then vector(1,0) else vector(-1,0))-Etheory)/Etheory) from (0,0.5) to (1,0.5)
END
