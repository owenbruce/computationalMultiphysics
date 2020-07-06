TITLE 'Two-Path Conduction'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  V (threshold=1e-6)
  Temp (threshold=1e-6)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
sigma
E = -grad(V)
J = -sigma*grad(V)


sigWire = 1/(26.5e-9)
sigRes = 1/(5E-6)

ww = 1e-3 !width of the wire
Lwx = 2e-2 !length of the wire in x
Lwy = 2e-2 !length of the wire in y

Lr = 0.01
arl = 0.003*0.001
arr = 0.003*0.001

thickness = 1e-3
right = vector(1,0)
down = vector(0,-1)
left = vector(-1,0)
I1 = thickness * surf_integral(dot(J, right), 'Entrance')
I2 = thickness * surf_integral(dot(J, right), 'top path')
I3 = thickness * surf_integral(dot(J, down), 'mid top')
v2 = bintegral(V,'top path')/ww
v3 = bintegral(V,'mid top')/ww

Lw = 2e-2*2+4*(5e-3-ww)
aw = ww*thickness
RTheoryW = Lw/(aw*sigWire)
RTheoryR = Lr/(arr*sigRes)
RTheoryL = Lr/(arl*sigRes)
ITheoryR = 10/RTheoryR
ITheoryL = 10/RTheoryL

h = 20
tFluid = 20

BPG = 400
MPA = 660.3

MaxGrTemp1 = GLOBALMAX(Temp, 'res_1')
MaxGrTemp2 = GLOBALMAX(Temp, 'res_2')
MaxAlTemp = GLOBALMAX(Temp, 'wire')

k
cp
rho 
Tempi = 20

qdot=-k*grad(Temp)

pow
qdotv = dot(E,J)

INITIAL VALUES
Temp = Tempi
EQUATIONS        { PDE's, one for each variable }
  V: div(J)=0 { one possibility }
  Temp: -div(-k*grad(Temp)) + qdotv = rho * cp*dt(Temp)
CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'wire'       { For each material region }
  k = 171
  rho = 2700
  cp = 0.9*1000
  pow = V^2/RTheoryW
	sigma = sigWire
    START(0,0) load(Temp) = -h*(Temp-tFluid)     { Walk the domain boundary }
	value(V) = 0 LINE TO (0,ww)
	load(V) = 0 LINE TO (Lwx/2,ww) to (Lwx/2,Lwy) to (0,Lwy)
	value(V) = 10 line to (0,Lwy+ww)
	load(V) = 0 line to (Lwx+ww,Lwy+ww) to (Lwx+ww,0) to close
	start(Lwx/2+ww,ww) line to (Lwx,ww) to (Lwx,Lwy) to (Lwx/2+ww, Lwy) to close
    REGION 'res_1'
    k = 130
    rho = 2150
    cp = 0.7078*1000
    pow = V^2/RTheoryL
    sigma = sigRes
    start(Lwx/2+ww/2-0.0010,(Lwy-ww)/2+ww-0.005) load(Temp) = -h*(Temp-tFluid) 
    line to (Lwx/2+ww/2+0.0010,(Lwy-ww)/2+ww-0.005)
    line to (Lwx/2+ww/2+0.0010,(Lwy-ww)/2+ww+0.005)
    line to (Lwx/2+ww/2-0.0010,(Lwy-ww)/2+ww+0.005)
    line to close
    REGION 'res_2'
    k = 130
    rho = 2150
    cp = 0.7078*1000
    sigma = sigRes
    pow = V^2/RTheoryR
    start(Lwx+ww/2-0.0015/2,(Lwy-ww)/2+ww-0.005) load(Temp) = -h*(Temp-tFluid) 
    line to (Lwx+ww/2+0.0015/2,(Lwy-ww)/2+ww-0.005)
    line to (Lwx+ww/2+0.0015/2,(Lwy-ww)/2+ww+0.005)
    line to (Lwx+ww/2-0.0015/2,(Lwy-ww)/2+ww+0.005)
    line to close
feature
	start 'Entrance'(Lwx/4,Lwy) line to (Lwx/4,Lwy+ww)
	start 'top path' (Lwx*.75, Lwy) line to (Lwx*.75,Lwy+ww)
	start 'mid top' (Lwx/2, 0.9*Lwy) line to (Lwx/2+ww, 0.9*Lwy)
TIME 0 TO 1e-6 halt(MaxGrTemp1>BPG or MaxGrTemp2>BPG or MaxAlTemp>MPA)     { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for time =  endtime
  CONTOUR(V) painted
  contour(Temp) painted
	vector(J)
SUMMARY
report I1
report I2
report I3
report v2
report v3
report RTheoryR
report ITheoryR
report RTheoryL
report ITheoryL
report MaxGrTemp1
report MaxGrTemp2
report MaxAlTemp
report abs(I2-ITheoryR) as 'Right Current Error'
report abs(I3-ITheoryL) as 'Left Current Error'
report (I1-(I2+I3))/I1 as 'KCL check top'

END
