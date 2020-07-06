{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New code'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
Voltage
rho_free
Temp
u
v
           { choose your own names }
!SELECT   
     { method controls }
 DEFINITIONS 
 hbot = 1e-6
 hmid = 1e-6
 Ltip = 1e-6
 htop = 3e-6
 L = 30e-6
mag=0.1*globalmax(magnitude(x,y))/globalmax(magnitude(U,V))
 K
 E
 nu
 Vmax =1
 alpha
 
 
 rhoe   !Ohm-m
sigmae = 1/rhoe !1/(Ohm-m)

hconv = 20
roomtemp = 20
 Efield = -grad(Voltage) !V/m = N/C
J = sigmae*Efield !N/(C-Ohm-m) = V/(m^2-Ohm)

epsilonr = 1! !basically zero, could put it to zero completely and it'll still work
epsilon = 8.85e-12*epsilonr

Dfield = epsilon*Efield


Efieldx = xcomp(Efield)
Efieldy = ycomp(Efield)
Efieldz = 0



DeltaTemp = Temp - 20

qdot = -k*grad(Temp) !heat flux density
qdotvol = dot(J, Efield) !Joule's law



!piezoelectric coupling coefficients
d11 =0		d12 =0	d13 = 0	d14 =0	d15 = 0				d16 = 0
d21 =0				d22 = 0		d23 = 0		d24 = 0				d25 =0 d26 =0
d31 = 0				d32 = 0		d33 = 0		d34 = 0				d35 = 0				d36 = 0


G=E/(2*(1+nu))

C11 =E*(1-nu)/(1+nu)/(1-2*nu) C12 = E*nu/(1+nu)/(1-2*nu) C13 = C12 C14 = 0 C15 = 0 C16 = 0
C21 = C12 C22 = C11 C23 = C12 C24 = 0 C25 = 0 C26 = 0
C31 = C12 C32 = C12 C33 = C11 C34 = 0 C35 = 0 C36 = 0
C41 = 0 C42 = 0 C43 = 0 C44 = E/(2*(1+nu)) C45 = 0 C46 = 0
C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = C44 C56 = 0
C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = C44

!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = 0!dz(w)
gyz = 0!dy(w) + dz(v)
gxz = 0!dx(w) + dz(u)
gxy = dx(v) + dy(u)

!Mechanical strain
exm  = ex  - alpha *DeltaTemp - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym  = ey  - alpha*DeltaTemp - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm  = ez  !- alphaz *DeltaTemp - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz !- alphayz*DeltaTemp - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz !- alphaxz*DeltaTemp - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy !- alphaxy*DeltaTemp - (d16*Efieldx+d26*Efieldy+d36*Efieldz)


!Hookes Law
sx  = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy  = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz  = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym



 
INITIAL VALUES
temp = 20

EQUATIONS        { PDE's, one for each variable }
Voltage: div(J) = 0
rho_free: div(Dfield) = rho_free
u: dx(sx) + dy(sxy) = 0
v: dx(sxy) + dy(sy)  = 0
Temp: div(qdot) = qdotvol
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Alu'   
  
  E = 70e9
  nu = 0.33
  Alpha = 24e-6
  k = 100
  rhoe = 35e-6
  
  
  d11 =0		d12 =0	d13 = 0	d14 =0	d15 = 0				d16 = 0
d21 =0				d22 =0		d23 = 0		d24 = 0				d25 =0 d26 =0
d31 = 0				d32 = 0		d33 = 0		d34 = 0				d35 = 0				d36 = 0
  
	start(0,hmid+hbot)
    
    value(voltage) = Vmax
    value(temp) = 20
    value(u) = 0
    value(v) = 0
    
    line to (0,hmid+hbot+htop)
    
    load(voltage) = 0
    load(u) = 0
    load(v) = 0
    load(Temp)=0
    
    line to (L,hmid+hbot+htop)
    
        load(voltage) = 0
    load(u) = 0
    load(v) = 0
    load(Temp)=0
    
    line to (L,hmid+hbot)
    
        load(voltage) = 0
    load(u) = 0
    load(v) = 0
    
    load(temp) = 0
    
    line to (L-Ltip,hmid+hbot)
     load(voltage) = 0
    load(u) = 0
    load(v) = 0
    
    load(Temp)=0
    
    
    line to close
    

    
    
REGION 'Gold'
E = 70e9
nu = 0.433
Alpha = 14e-6
k = 300
rhoe = 40e-9

d11 =0		d12 =0	d13 = 0	d14 =0	d15 = 0				d16 = 0
d21 =0				d22 = 0		d23 = 0		d24 = 0				d25 =0 d26 =0
d31 = 0				d32 = 0		d33 = 0		d34 = 0				d35 = 0				d36 = 0

start(0,0) 

load(u) = 0
    load(v) = 0
load(Temp)=0

line to (L, 0)

line to (L,hbot+hmid)
load(temp) = 0
line to (L-Ltip,hbot+hmid)

load(Temp)=0

line to (L-Ltip,Hbot)

line to (0,Hbot)
value(Voltage) = 0
value(temp) = 20
value(u) = 0
value(v) = 0
line to close 

Region 'peizo'
E = 0
nu = 0 
Alpha = 10e-6
k = 6
rhoe = 1
epsilonr = 700


d11 =0		d12 =0	d13 = 0	d14 =0	d15 = 0				d16 = 0
d21 =-20e-12				d22 = 20e-12		d23 = 0		d24 = 0				d25 =0 d26 =0
d31 = 0				d32 = 0		d33 = 0		d34 = 0				d35 = 0				d36 = 0


C11 =108e9 C12 = 20e9 C13 = 21e9 C14 = 0 C15 = 0 C16 = 0
C21 = 20e9 C22 = C11 C23 = 21e9 C24 = 0 C25 = 0 C26 = 0
C31 = 21e9 C32 = C12 C33 = 90e9 C34 = 0 C35 = 0 C36 = 0
C41 = 0 C42 = 0 C43 = 0 C44 = 50e9 C45 = 0 C46 = 0
C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = 50e9 C56 = 0
C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = 40e9







start(0,hbot)

line to (0,hmid+hbot)

line to (L-Ltip,hmid+hbot)

line to (L-Ltip, hbot)

line to close




! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  grid(x+u*mag, y+v*mag)
  contour(temp) painted
  contour(voltage) painted
  vector(J) norm
  vector(qdot) norm
  vector(E) norm
summary 
report val(v,L,htop+hmid+hbot) as 'tip disp'
report val(Vmax,0,0) as 'Vmax'

report integral(rho_free, 'alu') + integral(rho_free, 'peizo')
	report integral(rho_free, 'gold') + integral(rho_free, 'peizo')

END