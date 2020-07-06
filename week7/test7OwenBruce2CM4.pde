TITLE 'PiezoThermalVoigtNotation'     { the problem identification }
COORDINATES cartesian3  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
u
v
w
SELECT         { method controls }
quadratic
DEFINITIONS    { parameter definitions }
mag = 0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))
Lx = .03
Lz = .0015

hbot = 3e-4
htop = 2e-4

sapplied=0
DeltaTemp = 0
Vmax = 100

Efieldx = 0
Efieldy = if(y<0) then -Vmax/hbot else 0
Efieldz = 0

E = 0
nu = 0

!Stiffness Matrix components
C11 	C12 	C13 	C14 	C15 	C16 
C21 	C22 	C23 	C24 	C25 	C26 
C31 	C32 	C33 	C34 	C35 	C36 
C41 	C42 	C43 	C44 	C45  C46 
C51 	C52 	C53 	C54 	C55 	C56 
C61 	C62 	C63 	C64 	C65 	C66 

!thermal expansion coefficients
alphax = 0 !for e.g.
alphay = 0 
alphaz = 0
alphayz = 0 !0 unless monoclinic or triclinic
alphaxz = 0
alphaxy = 0

!piezoelectric coupling coefficients
d11	d12 	d13 	d14 	d15	d16 
d21	d22	d23 	d24	d25	d26 
d31	d32	d33 	d34	d35	d36 
!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
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

utip = 1/2*c*Lx^2


EQUATIONS        { PDE's, one for each variable }
u: dx(sx) + dy(sxy) + dz(sxz) = 0
v: dx(sxy) + dy(sy) + dz(syz) = 0
w: dx(sxz) + dy(syz) + dz(sz) = 0

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES       { The domain definition }

surface 'bottom'
	load(u) = 0
    load(w) = 0
    load(v) = 0
surface 'top'

    REGION 'PZT'       { For each material region }
    !Stiffness
    C11 = 109e9	C12 = 14e9		C13 = 21e9		C14 = 0		C15 = 0		C16 = 0
	C21 = C12		C22 = 49e9	C23 = 13e9		C24 = 0		C25 =0			C26 = 0
	C31 = C13		C32 = C23		C33 = 62e9	C34 = 0		C35 = 0		C36 = 0
	C41 = C14		C42 = C24		C43 = C34		C44 = 100e9	C45 = 0		C46 = 0
	C51 = C15		C52 = C25		C53 = C35		C54 = C45	C55 = 100e9	C56 = 0
	C61 = C16		C62 = C26		C63 = C36		C64 = C46	C65 = C56	C66 = 83e9
    !piezoelectric
    d11 = 0	d12 = 0	d13 = 0	d14 = 0			d15 = 2.3e-12			d16 = 0
    d21 = -2e-12		 		d22 = 2e-12				d23 = 0			d24 = 0	d25 = 0			d26 = 0
    d31 = 0				d32 = 0				d33 = 0			d34 = 0			d35 = 0	d36 = 1.4e-12
    START(0,0) !z=0
	value(u) = 0
	value(v) = 0
	value(w) = 0
    LINE TO (0,-hbot) !y=Ly
	load(u)= 0
	load(v) = 0
	load(w) = 0

	LINE TO (Lx,-hbot) !z=Lz

	LINE TO (Lx,0) !y=0

	LINE TO CLOSE
    
    REGION 'Chromium'       { For each material region }
    E = 279E9
    nu = 0.21
    !stiffness
	C11 =E*(1-nu)/(1+nu)/(1-2*nu) C12 = E*nu/(1+nu)/(1-2*nu) C13 = C12 C14 = 0 C15 = 0 C16 = 0
	C21 = C12 C22 = C11 C23 = C12 C24 = 0 C25 = 0 C26 = 0
	C31 = C12 C32 = C12 C33 = C11 C34 = 0 C35 = 0 C36 = 0
	C41 = 0 C42 = 0 C43 = 0 C44 = E/(2*(1+nu)) C45 = 0 C46 = 0
	C51 = 0 C52 = 0 C53 = 0 C54 = 0 C55 = C44 C56 = 0
	C61 = 0 C62 = 0 C63 = 0 C64 = 0 C65 = 0 C66 = C44

    !piezoelectric
    d11 = 0	d12 = 0	d13 = 0	d14 = 0	d15 = 0	d16 = 0
    d21 = 0	d22 = 0	d23 = 0	d24 = 0	d25 = 0	d26 = 0
    d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0
    START(0,0) !z=0
	value(u) = 0
	value(v) = 0
	value(w) = 0
    LINE TO (0,htop) !y=Ly
	load(u)= 0
	load(v) = 0
	load(w) = 0

	LINE TO (Lx,htop) !z=Lz

	LINE TO (Lx,0) !y=0

	LINE TO CLOSE

PLOTS            { save result displays }
 grid(x+mag*u, y+mag*v, z+mag*w)
 !CONTOUR(c) painted on x = 0
elevation(c) from (Lx/2,0,0) to (Lx,0,0)
elevation(utip) from (Lx/2,0,0) to (Lx,0,0)

SUMMARY
!report val(ex, Ly/2,Lz/2,0)
!report val(ey,Ly/2,Lz/2, 0)
!report val(ez,Ly/2,Lz/2, 0)
!report val(gyz,Ly/2,Lz/2, 0)
!report val(gxz,Ly/2,Lz/2, 0)
!report val(gxy,Ly/2,Lz/2, 0)
report val(utip, Lx, 0, Lz/2)
report val(c, Lx, 0, Lz/2)
END
