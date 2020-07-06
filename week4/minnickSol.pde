TITLE 'Two-Path Conduction'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  V (threshold=1e-3)
Temp (threshold = 1e-3)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
sigma
E = -grad(V)
J = sigma*E

cp = 720
rho = 2150
k =140
Tlimit = 800

sigWire = 1/(26.5e-9)
sigR = 1/(5e-6)

BPG = 400
MPA = 660.3

ww = 1e-3 !width of the wire
Lwx = 2e-2 !length of the wire in x
Lwy = 2e-2 !length of the wire in y

wR = 3e-3 !width of [right] resistor in x
LR = 1e-2 !length of resistor in y
wR3 = 3e-3 !width of left resistor

DeltaV = 10
thickness = 1e-3
right = vector(1,0)
down = vector(0,-1)
left = vector(-1,0)
I1 = thickness * surf_integral(dot(J, right), 'Entrance')
I2 = thickness * surf_integral(dot(J, right), 'top path')
I3 = thickness * surf_integral(dot(J, down), 'mid top')

R2 = LR/(sigR*thickness*wR) !right resistor
R3 = LR/(sigR*thickness*wR3) !middle resistor
I2theory = DeltaV/R2
I3theory = DeltaV/R3

MaxGrTemp1 = GLOBALMAX(Temp, 'Left resistor (R3)')
MaxGrTemp2 = GLOBALMAX(Temp, 'Right resistor (R2)')
MaxAlTemp = GLOBALMAX(Temp, 'wire')

qdotvol = dot(J,E)
qdot = -k*grad(Temp)
INITIAL VALUES
Temp = 20
EQUATIONS        { PDE's, one for each variable }
V:  div(J)=0
Temp: dt(rho*cp*Temp) = qdotvol - div(qdot)
CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'wire'       { For each material region }
	cp = 902
	rho = 2700
	k = 237
	Tlimit = 660
	sigma = sigWire
    START(0,0)   { Walk the domain boundary }
	load(Temp) = 20*(Temp-20)
	value(V) = 0 LINE TO (0,ww)
	load(V) = 0 LINE TO (Lwx/2,ww) to (Lwx/2,Lwy) to (0,Lwy)
	value(V) = DeltaV line to (0,Lwy+ww)
	load(V) = 0 line to (Lwx+ww,Lwy+ww) to (Lwx+ww,0) to close
	start(Lwx/2+ww,ww) line to (Lwx,ww) to (Lwx,Lwy) to (Lwx/2+ww, Lwy) to close
	REGION 'Left resistor (R3)'
		sigma = sigR
		start(Lwx/2-wR3/2+ww/2, Lwy/2+ww/2-LR/2) 	load(Temp) = 20*(Temp-20) line to (Lwx/2+wR3/2+ww/2, Lwy/2+ww/2-LR/2) to (Lwx/2+wR3/2+ww/2, Lwy/2+ww/2+LR/2) to (Lwx/2-wR3/2+ww/2, Lwy/2+ww/2+LR/2) to close
	REGION 'Right resistor (R2)'
		sigma = sigR
		start(Lwx-wR/2+ww/2, Lwy/2+ww/2-LR/2) 	load(Temp) = 20*(Temp-20) line to (Lwx+wR/2+ww/2, Lwy/2+ww/2-LR/2) to (Lwx+wR/2+ww/2, Lwy/2+ww/2+LR/2) to (Lwx-wR/2+ww/2, Lwy/2+ww/2+LR/2) to close

feature
	start 'Entrance'(Lwx/4,Lwy) line to (Lwx/4,Lwy+ww)
	start 'top path' (Lwx*.75, Lwy) line to (Lwx*.75,Lwy+ww)
	start 'mid top' (Lwx/2, 0.9*Lwy) line to (Lwx/2+ww, 0.9*Lwy)
TIME 0 TO 6e-3 halt(MaxGrTemp1>BPG or MaxGrTemp2>BPG or MaxAlTemp>MPA)  { if time dependent }
MONITORS         { show progress }
for t = 0 by 5e-5 to endtime
contour(Temp) painted
PLOTS            { save result displays }
for t= endtime
  CONTOUR(V) painted
	vector(J)
  contour(Temp) painted
SUMMARY
report I1
report I2
report I3
report (I1-(I2+I3))/I1 as 'KCL check top'
report MaxGrTemp1
report MaxGrTemp2
report MaxAlTemp
report I2theory
report I3theory
END
