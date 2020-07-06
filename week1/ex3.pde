TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
 xd(threshold=0.1)
 vx(threshold=0.1)
 yd(threshold=0.1)
 vy(threshold=0.1)
SELECT         { method controls }
ngrid=1
DEFINITIONS    { parameter definitions }
vi = 20
thetai = 60*pi/180
mi=800
q=40
vfuel=1000
g=9.81

mfueli = 0.8*mi
mDry=mi-mfueli

tfuel = mfueli/q
mfuel = if(t<tfuel) then mfueli-q*t else 0

Ft = if(t<tfuel) then q*vfuel else 0
m = mDry + mfuel


vwindx = 7
vrelx = vx-vwindx
vrely = vy-0
vrel = sqrt(vrelx^2+vrely^2)
v= sqrt(vx^2+vy^2)+1e-6

Area = 0.8
CD = 0.05
rhoAir = 1.2
Fdmag = 0.5*rhoAir*CD*Area*v^2

Fdx = Fdmag*(-vrelx/vrel)
Fdy = Fdmag*(-vrely/vrel)

Ftx = Ft*vx/v
Fty = Ft*vy/v
ay = (Fty+Fdy)/m-g
ax = (Ftx+Fdx)/m
INITIAL VALUES
vx=  vi*cos(thetai)
vy=vi*sin(thetai)

EQUATIONS        { PDE's, one for each variable }
  yd: dt(yd)=vy
  vy: dt(vy)=ay
   xd: dt(xd)=vx
  vx: dt(vx)=ax
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 200 halt(yd<-0.1)   { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by endtime/100 to endtime
  history(yd,vy,ay) at (0,0)
  history(yd) at (0,0) vs xd
  history(yd,vy,ay) at (0,0) PrintOnly Export Format '#t#b#1#b#2' file = 'test.txt'
END
