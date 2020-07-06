TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  xd (threshold=0.1)
  yd (threshold=0.1)
  vx (threshold=0.1)
  vy (threshold=0.1){ choose your own names }
SELECT         { method controls }
ngrid=1
DEFINITIONS    { parameter definitions }
vi = 21
theta0 = 60*pi/180

ax = 0
ay = -9.81
INITIAL VALUES

xd = 0
yd = 0
vx = vi*cos(theta0)
vy = vi*sin(theta0)


EQUATIONS        { PDE's, one for each variable }
  xd: dt(xd)=vx
  yd: dt(yd)=vy
  vx: dt(vx) = ax
  vy: dt(vy) = ay
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 50 halt(yd<0)    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by 0.1 to endtime
  history(yd,xd) at (0,0)
  history(xd,yd) at (0,0) PrintOnly Export Format '#t#b#1#b#2' file = 'test.txt'
END
