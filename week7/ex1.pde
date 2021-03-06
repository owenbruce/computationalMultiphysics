TITLE 'PiezoThermalVoigtNotation'    
COORDINATES cartesian3  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u
  v
  w
! SELECT         { method controls }
DEFINITIONS
mag = 0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))
Lx = 0.05
Ly = 0.01
Lz = 0.001

sapplied = 0
DeltaTemp = 0

Efieldx = 0
Efieldy = 0
Efieldz = -1e5

!Stiffness Matrix components
C11 = 82e9	C12 = 35e9	C13 = 34e9	C14 = 0	C15 = -2e9	C16 = 0
C21 = C12	C22 = 63e9	C23 = 35e9	C24 = 0	C25 = -8e9	C26 = 0
C31 = C13	C32 = C23	C33 = 58e9	C34 = 0	C35 = -3e9	C36 = 0
C41 = C14	C42 = C24	C43 = C34	C44 = 21e9	C45 = 0	C46 = -5e9
C51 = C15	C52 = C25	C53 = C35	C54 = C45	C55 = 28e9	C56 = 0
C61 = C16	C62 = C26	C63 = C36	C64 = C46	C65 = C56	C66 = 29e9

!Thermal Expansion Coefficients
alphax = 0
alphay = 0
alphaz = 0
alphayz = 0
alphaxz = 0
alphaxy = 0

!piezoelectric coupling coefficients
d11 = 0	d12 = 0	d13 = 0	d14 = 0	d15 = 584e-12	d16 = 0
d21 = 0	d22 = 0	d23 = 0	d24 = d15	d25 = 0	d26 = 0
d31 = -171e-12	d32 = d31	d33 = 374e-12	d34 = 0	d35 = 0	d36 = 0

!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
gxy = dx(v) + dy(u)

!Mechanical Strain
exm = ex - alphax*DeltaTemp - (d11*Efieldx + d21*Efieldy + d31*Efieldz)
eym = ey - alphay*DeltaTemp - (d12*Efieldx + d22*Efieldy + d32*Efieldz)
ezm = ez - alphaz*DeltaTemp - (d13*Efieldx + d23*Efieldy + d33*Efieldz)
gyzm = gyz - alphayz*DeltaTemp - (d14*Efieldx + d24*Efieldy + d34*Efieldz)
gxzm = gxz - alphaxz*DeltaTemp - (d15*Efieldx + d25*Efieldy + d35*Efieldz)
gxym = gxy - alphaxy*DeltaTemp - (d16*Efieldx + d26*Efieldy + d36*Efieldz)

!Hookes Law
sx = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym

! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  u: dx(sx) + dy(sxy) + dz(sxz) = 0
  v: dx(sxy) + dy(sy) + dz(syz) = 0
  w: dx(sxz) + dy(syz) + dz(sz) = 0
  
EXTRUSION
surface 'bottom' z = 0
surface 'top' z=Lz
BOUNDARIES       { The domain definition }
	surface 'bottom'
    load(u) = 0
    load(v) = 0
    load(w) = 0
    surface 'top'
    load(u) = 0
    
    Region 1
    	Start(0,0)
        Line to (Lx,0)
        Line to (Lx,Ly)
        Line to (0,Ly)
        !value(u) = 0
        Line to close
        
PLOTS            { save result displays }
	grid(x+mag*u, y+mag*v, z+mag*w)
    contour(ez) painted on x = Lx/2
SUMMARY
report val(ex, Lx/2, Ly/2, Lz/2)
report val(ey, Lx/2, Ly/2, Lz/2)
report val(ez, Lx/2, Ly/2, Lz/2)
report val(gyz, Lx/2, Ly/2, Lz/2)
report val(gxz, Lx/2, Ly/2, Lz/2)
report val(gxy, Lx/2, Ly/2, Lz/2)
END
