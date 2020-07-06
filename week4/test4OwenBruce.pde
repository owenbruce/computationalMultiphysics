TITLE 'Two-Path Conduction'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  V
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
sigma
J = -sigma*grad(V)


wR = 1e-2
wR2 = wR/2
wR3 = wR*2
LR1 = 10e-2
LR2 = 4e-2
LR3 = 6e-2

sigR1 = 2e-3
sigR2 = 4e-3
sigR3 = 1e-3
sigwire = 1e6

LW = 5e-2 !wire length at each junction
WW = .5e-2 !wire width

thickness = 1e-3
right = vector(1,0)
down = vector(0,-1)
left = vector(-1,0)
up = vector(0,1)

I0 = thickness*surf_integral(dot(J,up),'I0')
I1 = thickness*surf_integral(dot(J,right),'I1')
I2 = thickness*surf_integral(dot(J,right),'I2')

KCL_Check = I0-I1-I2

!Theoretical
aR = thickness*wR
aR2 = thickness*wR2
aR3 = thickness*wR3
aW = thickness*WW
Rwire = 4*LW/(aW*sigwire)
R1 = LR1/(aR*sigR1)
R2 = LR2/(aR*sigR2) !R2 = LR2/(aR2*sigR2) !for part d)
R3 = LR3/(aR*sigR3) !R3 = LR3/(aR3*sigR3) !for part d)

I1Theory = 10/R1
I23Theory = 10/(R2+R3)
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(J)=0 { one possibility }
CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'R1'
	sigma = sigR1
	start (0,-LW) line to (LR1,-LW) to (LR1,-LW+wR) to (0,-LW+WR) to close
REGION 'R2'
	sigma = sigR2
	start (0,0) line to (LR2,0) to (LR2,0+wR) to (0,0+WR) to close !start (0,0) line to (LR2,0) to (LR2,0+wR2) to (0,0+WR2) to close !for part d)
REGION 'R3'
	sigma = sigR3
	start (LR2,0) line to (LR1,0) to (LR1,0+wR) to (LR2,0+WR) to close !start (0,0) line to (LR2,0) to (LR2,0+wR3) to (0,0+WR3) to close !for part d)
REGION 'wire'
	sigma = sigwire
	start (-ww,-2*LW) value(V) = 10 line to (0,-2*LW) load(V) = 0 line to (0,wr) to (-ww,wr) to close
	start (LR1,-2*LW) value(V) = 0 line to (LR1+ww,-2*LW) load(V) = 0 line to (LR1+ww,wr) to (LR1,wr) to close
FEATURE
	start 'I0' (-ww-0.01, -2*Lw*0.9) line to (0+0.01,-2*Lw*0.9)
    start 'I1' (0+0.01,-LW-0.01) line to (0+0.01,-LW+wR+0.01)
    start 'I2' (0+0.01, 0-0.01) line to (0+0.01, 0+wR+0.01)

! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(V) painted
	vector(J)
SUMMARY
report I0
report I1
report I2
report KCL_Check
report I1Theory
report I23Theory
report abs(I1-I1Theory)
report abs(I2-I23Theory)
report R1
report R2
report R3
report Rwire
END
