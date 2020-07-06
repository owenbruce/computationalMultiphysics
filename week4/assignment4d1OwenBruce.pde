TITLE 'Two-Path Conduction'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  V
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
sigma
J = -sigma*grad(V)


sigWire = 1/(26.5e-9)
sigRes = 1/(5E-6)

ww = 1e-3 !width of the wire
Lwx = 2e-2 !length of the wire in x
Lwy = 2e-2 !length of the wire in y

Lr = 0.01
arr = 0.002*0.001
arl = 0.0015*0.001

thickness = 1e-3
right = vector(1,0)
down = vector(0,-1)
left = vector(-1,0)
I1 = thickness * surf_integral(dot(J, right), 'Entrance')
I2 = thickness * surf_integral(dot(J, right), 'top path')
I3 = thickness * surf_integral(dot(J, down), 'mid top')
v2 = bintegral(V,'top path')/ww
v3 = bintegral(V,'mid top')/ww

RTheoryR = Lr/(arr*sigRes)
RTheoryL = Lr/(arl*sigRes)
ITheoryR = 10/RTheoryR
ITheoryL = 10/RTheoryL
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(J)=0 { one possibility }
CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'wire'       { For each material region }
	sigma = sigWire
    START(0,0)   { Walk the domain boundary }
	value(V) = 0 LINE TO (0,ww)
	load(V) = 0 LINE TO (Lwx/2,ww) to (Lwx/2,Lwy) to (0,Lwy)
	value(V) = 10 line to (0,Lwy+ww)
	load(V) = 0 line to (Lwx+ww,Lwy+ww) to (Lwx+ww,0) to close
	start(Lwx/2+ww,ww) line to (Lwx,ww) to (Lwx,Lwy) to (Lwx/2+ww, Lwy) to close
    REGION 'res_1'
    sigma = sigRes
    start(Lwx/2+ww/2-0.0015/2,(Lwy-ww)/2+ww-0.005)
    line to (Lwx/2+ww/2+0.0015/2,(Lwy-ww)/2+ww-0.005)
    line to (Lwx/2+ww/2+0.0015/2,(Lwy-ww)/2+ww+0.005)
    line to (Lwx/2+ww/2-0.0015/2,(Lwy-ww)/2+ww+0.005)
    line to close
    REGION 'res_2'
    sigma = sigRes
    start(Lwx+ww/2-0.0010,(Lwy-ww)/2+ww-0.005)
    line to (Lwx+ww/2+0.0010,(Lwy-ww)/2+ww-0.005)
    line to (Lwx+ww/2+0.0010,(Lwy-ww)/2+ww+0.005)
    line to (Lwx+ww/2-0.0010,(Lwy-ww)/2+ww+0.005)
    line to close
feature
	start 'Entrance'(Lwx/4,Lwy) line to (Lwx/4,Lwy+ww)
	start 'top path' (Lwx*.75, Lwy) line to (Lwx*.75,Lwy+ww)
	start 'mid top' (Lwx/2, 0.9*Lwy) line to (Lwx/2+ww, 0.9*Lwy)
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(V) painted
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
report abs(I2-ITheoryR) as 'Right Current Error'
report abs(I3-ITheoryL) as 'Left Current Error'
report (I1-(I2+I3))/I1 as 'KCL check top'

END
