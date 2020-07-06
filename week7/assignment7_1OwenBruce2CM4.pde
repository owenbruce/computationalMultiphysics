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
Ly = .03
Lz = .004

hbot = 4e-4
htop = 6e-4

sapplied=0
DeltaTemp = 0
Vmax = 200

Efieldx = if(x<0) then Vmax/hbot else -Vmax/htop
Efieldy = 0
Efieldz = 0

!Stiffness Matrix components
C11 =82e9	C12 = 35e9		C13 = 34e9	C14 = 0		C15 = -2e9	C16 = 0
C21 = C12	C22 = 63e9		C23 = 35e9	C24 = 0		C25 = -8e9	C26 = 0
C31 = C13	C32 = C23		C33 = 58e9	C34 = 0		C35 = -3e9	C36 = 0
C41 = C14	C42 = C24		C43 = C34	C44 = 21e9	C45 = 0		C46 = -5e9
C51 = C15	C52 = C25		C53 = C35	C54 = C45	C55 = 28e9	C56 = 0
C61 = C16	C62 = C26		C63 = C36	C64 = C46	C65 = C56	C66 = 29e9

!thermal expansion coefficients
alphax = 0 !for e.g.
alphay = 0 
alphaz = 0
alphayz = 0 !0 unless monoclinic or triclinic
alphaxz = 0
alphaxy = 0

!piezoelectric coupling coefficients
d11 = -2.3e-12	 d12 = 2.3e-12	d13 = 0	d14 = -0.67e-12	d15 = 0				d16 = 0
d21 = 0			 d22 = 0			d23 = 0	d24 = 0				d25 = 0.67e-12	d26 = 4.6e-12
d31 = 0			 d32 = 0			d33 = 0	d34 = 0				d35 = 0				d36 = 0
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

c = 2*u/y^2

vtip = 1/2*c*Ly^2

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

  REGION 1       { For each material region }
    START(0,0) !z=0
	value(u) = 0
	value(v) = 0
	value(w) = 0
    LINE TO (-hbot,0) !y=Ly
	load(u)= 0
	load(v) = 0
	load(w) = 0

	LINE TO (-hbot,Ly) !z=Lz

	LINE TO (0,Ly) !y=0

	LINE TO CLOSE
    
    REGION 2       { For each material region }
    START(0,0) !z=0
	value(u) = 0
	value(v) = 0
	value(w) = 0
    LINE TO (htop,0) !y=Ly
	load(u)= 0
	load(v) = 0
	load(w) = 0

	LINE TO (htop,Ly) !z=Lz

	LINE TO (0,Ly) !y=0

	LINE TO CLOSE

PLOTS            { save result displays }
 grid(x+mag*u, y+mag*v, z+mag*w)
 CONTOUR(vtip) painted on x = 0
SUMMARY
!report val(ex, Ly/2,Lz/2,0)
!report val(ey,Ly/2,Lz/2, 0)
!report val(ez,Ly/2,Lz/2, 0)
!report val(gyz,Ly/2,Lz/2, 0)
!report val(gxz,Ly/2,Lz/2, 0)
!report val(gxy,Ly/2,Lz/2, 0)
report val(vtip, 0, Ly/2, Lz/2)
END
