TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
 yd(threshold=0.1)
 yv(threshold=0.1)
SELECT         { method controls }
ngrid=1
DEFINITIONS    { parameter definitions }
mi=2
q = 0.1
vfuel = 250
g=9.81

mfuel=0.5*mi

tfuel=mfuel/q

Ft = if(t<tfuel) then q*vfuel else 0

m = if(t<tfuel) then mi-q*t else mi-mfuel

ya=Ft/m-g
!INITIAL VALUES

EQUATIONS        { PDE's, one for each variable }
  yd: dt(yd)=yv
  yv: dt(yv)=ya
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 50 halt(yd<-0.1)   { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by 1 to endtime
  history(yd,yv,ya) at (0,0)
  history(yd,yv,ya) at (0,0) PrintOnly Export Format '#t#b#1#b#2' file = 'test.txt'
END
