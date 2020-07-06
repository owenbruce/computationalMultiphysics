TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
r(threshold=1e-3) = vector(rx,ry)
v(threshold=1e-3) = vector(vx,vy)
 SELECT         { method controls }
ngrid = 1
quadratic
DEFINITIONS    { parameter definitions }
ag = vector(0, -9.81)

mi = 600
thetai=50*pi/180
vi = 10

mfueli = 0.8*mi
mdry = mi - mfueli

q = 20
tfuel = mfueli/q

mfuel = if(t<tfuel) then mfueli-q*t else 0
vfuel = 1000
Ftmag = if(t<tfuel) then q*vfuel else 0
Ft = Ftmag*v/magnitude(v)

m = mdry + mfuel

Fg = ag*m


rho = 1.2
area = 5
Cd = 0.4
vwind = vector(-20,0)
vrel = v - vwind
Fdmag = 0.5*rho*Cd*area*(magnitude(vrel))^2
Fd = -Fdmag*vrel/magnitude(vrel)

Fnet =  Fg + Ft + Fd
a = Fnet/m

INITIAL VALUES
v = vi*vector(cos(thetai),sin(thetai))
EQUATIONS        { PDE's, one for each variable }
r: dt(r)=v
v: dt(v) = a

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 50 halt(ry<0)    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for t = 0 by endtime/50 to endtime
history(ry) at (0,0) vs rx
history(vx) at (0,0)
summary
report eval(rx, 0,0)

report eval(ry, 0,0)
END
