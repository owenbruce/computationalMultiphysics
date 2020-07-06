TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { choose your own names }
SELECT         { method controls }
modes = 3
DEFINITIONS    { parameter definitions }
L = 0.6285
m = 1.09/1e3
T = 105.3
cw=sqrt(T*L/m)
f = sqrt(lambda)/(2*pi)
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  cw^2*del2(u)=-lambda*u
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    Start(0) point value(u) = 0 line to (L) point value(u) = 0
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  !CONTOUR(u)
  elevation(u)
  !surface(u)
  SUMMARY
  report f
END
