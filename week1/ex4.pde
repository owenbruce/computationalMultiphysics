TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
 r(threshold=0.1) = vector(rx,ry)
 v(threshold=0.1) = vector(vx,vy)
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

vhat = v/magnitude(v)
Ftvec = Ft*vhat

agrav = vector(0,-g)

a = Ftvec/m + agrav
INITIAL VALUES
v = vi*vector(cos(thetai),sin(thetai))

EQUATIONS        { PDE's, one for each variable }
  r: dt(r) = v
  v: dt(v) = a
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 200 halt(ry<-0.1)   { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for t=0 by endtime/100 to endtime
  history(ry) at (0,0) vs rx
  SUMMARY
  report val(normal(Ftvec),0,0)
END
