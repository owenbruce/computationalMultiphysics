TITLE 'Test5'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { choose your own names }
SELECT         { method controls }
modes = 7
DEFINITIONS    { parameter definitions }
cw=1000
f = sqrt(lambda)/(2*pi)
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  cw^2*del2(u)=-lambda*u
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       
  REGION 1       
    START(0,0)   value(u) = 0
    LINE TO (1,1.6) 
	arc(center=0.5,2) angle=180 
	arc(center = -0.5,2) angle=180 
	line to (-1,1.6) 
	line TO CLOSE
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(u) painted
  !elevation(u)
  surface(u)
  SUMMARY
  report f
END
