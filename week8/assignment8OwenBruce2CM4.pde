TITLE 'EfieldConductor & Dielectric'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Voltage
rho_free
u
v
Temp
SELECT         { method controls }
!errlim = 1e-12
DEFINITIONS    { parameter definitions }
mag = 0.1*globalmax(magnitude(x,y))/globalmax(magnitude(U,V))
hbot = 1e-6
hmid = 1e-6
htop = 2e-6
L = 24e-6
Ltip = 4e-6


rho
sigma = 1/rho !1/(Ohm-m)

Efield = -grad(Voltage) !V/m = N/C
J = sigma*Efield !N/(C-Ohm-m) = V/(m^2-Ohm)

epr = 1
k
ep = 8.85e-12*epr

Dfield = ep*Efield

qdot = -k*grad(Temp) !heat flux density
qdotvol = dot(J, Efield) !Joule's law


Efieldx = xcomp(Efield)
Efieldy = ycomp(Efield)
Efieldz = 0

E 
nu
alpha
h = 20
Tfluid = 20
Vmax = 0.06165



C11 =E*(1-nu)/(1+nu)/(1-2*nu) C12 = E*nu/(1+nu)/(1-2*nu) C13 = C12 C14 = 0 C15 = 0 C16 = 0
C21 = C12 C22 = C11 C23 = C12 C24 = 0 C25 = 0 C26 = 0
C31 = C12 C32 = C12 C33 = C11 C34 = 0 C35 = 0 C36 = 0
C41 = 0 C42 = 0 C43 = 0 C44 = E/(2*(1+nu)) C45 = 0 C46 = 0
C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = C44 C56 = 0
C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = C44

!thermal expansion coefficients
alphax = alpha !for e.g.
alphay = alpha 
alphaz = 0
alphayz = 0 !0 unless monoclinic or triclinic
alphaxz = 0
alphaxy = 0

DeltaTemp =Temp - 20 !If starting (unstrained) temperature is 20

!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = 0!dz(w)
gyz = 0!dy(w) + dz(v)
gxz = 0!dx(w) + dz(u)
gxy = dx(v) + dy(u)

!Mechanical strain
exm  = ex  - alphax *DeltaTemp
eym  = ey  - alphay *DeltaTemp
ezm  = ez  - alphaz *DeltaTemp
gyzm = gyz - alphayz*DeltaTemp
gxzm = gxz - alphaxz*DeltaTemp
gxym = gxy - alphaxy*DeltaTemp

!Hookes Law
sx  = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy  = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz  = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym

c = 2*v/x^2

vtip = 1/2*c*L^2

 INITIAL VALUES !added to give a starting point for the simulation (even though not time-dependent) in attempt to make it start more smoothly
Voltage = Vmax/2
Temp = 20
rho_free = 0

EQUATIONS        { PDE's, one for each variable }
  !div(grad(Voltage))=0 { one possibility }
	!div(J+Dfield) = 0
Voltage: div(J) = 0
rho_free: div(Dfield) = rho_free
u: dx(sx) + dy(sxy) = 0
v: dx(sxy) + dy(sy)  = 0
Temp: div(qdot) = qdotvol


! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
	REGION 'Aluminum'
    E = 68.9E9
    nu = 0.33
    alpha = 2.32e-5
    k = 151
    rho = 32.5e-9
    START(0,hmid+hbot)
     load(Temp) = -h*(Temp-tFluid)
     Line to (L-Ltip,hmid+hbot)
     load(Temp) = 0
    Line to (L,hmid+hbot)
     load(Temp) = -h*(Temp-tFluid)
    Line to (L,hmid+hbot+htop)
    Line to (0,hmid+hbot+htop)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    value(Voltage) = Vmax
    Line to Close
    
    REGION 'Gold'
    E = 79e9
    nu = 0.4
    alpha = 14.2e-6
    k = 318
    rho = 22.14e-9
    START(0,0)
     load(Temp) = -h*(Temp-tFluid)
    Line to (L,0)
    Line to (L,hbot+hmid)
    load(Temp) = 0
    Line to (L-Ltip,hbot+hmid)
     load(Temp) = -h*(Temp-tFluid)
    Line to (L-Ltip,hbot)
    Line to (0,hbot)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    value(Voltage) = 0
    Line to Close
FEATURE
	'Right Side' Start(L,0) line to (L, hbot+hmid+htop)
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
 grid(x+mag*u, y+mag*v)
  CONTOUR(Voltage) painted
    vector(J) norm
    vector(Efield) norm
  Contour(Temp) painted
  vector(qdot) norm
  elevation(v) from (L,0) to (L,hbot+hmid+htop)
  elevation(Globalmax(Temp)) from (0,0) to (L,0)
summary
report Vmax
report Globalmax(Temp)
report val(v,L,0)
report val(v,L,hbot+hmid+htop)
END
