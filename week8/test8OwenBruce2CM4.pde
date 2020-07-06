TITLE 'EfieldConductor & Dielectric'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Voltage
rho_free
u
v
Temp
SELECT         { method controls }
errlim = 1e-9
DEFINITIONS    { parameter definitions }
mag = 0.1*globalmax(magnitude(x,y))/globalmax(magnitude(U,V))
hbot = 1e-6
hmid = 1e-6
htop = 3e-6
L = 30e-6
Ltip = 1e-6


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
Vmax = 1

d11 = 0	d12 = 0 	d13 = 0 	d14 = 0 	d15 = 0	d16 = 0 
d21 = 0	d22 = 0	d23 = 0 	d24 = 0	d25 = 0	d26 = 0
d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0

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
exm  = ex  - alphax *DeltaTemp - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym  = ey  - alphay *DeltaTemp - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm  = ez  - alphaz *DeltaTemp - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz - alphayz*DeltaTemp - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz - alphaxz*DeltaTemp - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy - alphaxy*DeltaTemp - (d16*Efieldx+d26*Efieldy+d36*Efieldz)


!Hookes Law
sx  = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy  = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz  = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym

c = 2*v/x^2

vtip = 1/2*c*L^2

chargetop = integral(rho_free,'Grey') + integral(rho_free,'Piezo Top')
chargebot = integral(rho_free,'Piezo Bot') + integral(rho_free,'Yellow')
chargemid = bintegral(-ycomp(Dfield), 'midline')

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
	REGION 'Grey'
    E = 70E9
    nu = 0.33
    alpha = 24e-6
    k = 100
    rho = 35e-6
    START(0,hmid+hbot)
     Line to (L-Ltip,hmid+hbot)
    Line to (L,hmid+hbot)
    Line to (L,hmid+hbot+htop)
    Line to (0,hmid+hbot+htop)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    value(Voltage) = Vmax
    Line to Close
    
    REGION 'Piezo Top'
    E = 0
    nu = 0
    alpha = 10e-6
    k = 6
    rho = 1
    epr = 700
    C11 = 108E9 C12 = 20E9 C13 = 21E9 C14 = 0 C15 = 0 C16 = 0
	C21 = C12 C22 = C11 C23 = 21E9 C24 = 0 C25 = 0 C26 = 0
	C31 = C12 C32 = C12 C33 = 90E9 C34 = 0 C35 = 0 C36 = 0
	C41 = 0 C42 = 0 C43 = 0 C44 = 50E9 C45 = 0 C46 = 0
	C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = 50E9 C56 = 0
	C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = 40E9
    
    d11 = 0	d12 = 0 	d13 = 0 	d14 = 0 	d15 = 0	d16 = 0 
	d21 = -20e-12	d22 = 20e-12	d23 = 0 	d24 = 0	d25 = 0	d26 = 0
	d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0
    
    START(0,hbot+hmid/2)
    line to (L-Ltip,hbot + hmid/2)
    line to (L-Ltip,hbot+hmid)
    line to (0, hbot+hmid)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    line to Close
    
    REGION 'Piezo Bot'
    E = 0
    nu = 0
    alpha = 10e-6
    k = 6
    rho = 1
    epr = 700
    C11 = 108E9 C12 = 20E9 C13 = 21E9 C14 = 0 C15 = 0 C16 = 0
	C21 = C12 C22 = C11 C23 = 21E9 C24 = 0 C25 = 0 C26 = 0
	C31 = C12 C32 = C12 C33 = 90E9 C34 = 0 C35 = 0 C36 = 0
	C41 = 0 C42 = 0 C43 = 0 C44 = 50E9 C45 = 0 C46 = 0
	C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = 50E9 C56 = 0
	C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = 40E9
    
    d11 = 0	d12 = 0 	d13 = 0 	d14 = 0 	d15 = 0	d16 = 0 
	d21 = -20e-12	d22 = 20e-12	d23 = 0 	d24 = 0	d25 = 0	d26 = 0
	d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0
    
    START(0,hbot)
    line to (L-Ltip,hbot)
    line to (L-Ltip,hbot+hmid/2)
    line to (0, hbot+hmid/2)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    line to Close
    	
    
    REGION 'Yellow'
    E = 70e9
    nu = 0.33
    alpha = 14e-6
    k = 300
    rho = 40e-9
    START(0,0)
    Line to (L,0)
    Line to (L,hbot+hmid)
    Line to (L-Ltip,hbot+hmid)
    Line to (L-Ltip,hbot)
    Line to (0,hbot)
    value(u) = 0
    value(v) = 0
    value(Temp) = 20
    value(Voltage) = 0
    Line to Close
FEATURE
	'Midline' Start(0,hbot+hmid/2) line to (L-Ltip, hbot+hmid/2)
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
 grid(x+u, y+v)
  CONTOUR(Voltage) painted
  contour(rho_free) painted
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
report chargetop
report chargemid
report chargebot
END
