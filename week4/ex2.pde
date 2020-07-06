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
sigma = 2
epsilon = ep0*(1+chi_e)
D = epsilon*E
J = sigma*E
I = J
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(D) = rho
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
  chi_e = 8E-6
    start(0,0.5) point load(V) = 10E-12 load(V) = 0
    arc(center=0,0) angle=180 point value(V) = 7 load(V) = -10E-42
    arc(center = 0,0) angle=180
    line to close
   
    
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(V) painted
  surface(V)
  vector(E) norm
  vector(D)
  vector(J)
  vector(I)
END
