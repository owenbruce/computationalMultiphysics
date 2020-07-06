TITLE 'Plate'
SELECT
errlim = .01
modes = 3

COORDINATES
cartesian3

VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

DEFINITIONS
Lz = 1e-2
R=51e-2

E = 68.9e9
nu = 0.33
rho = 2700

ftheory = if(Mode = 1) then 91.6256 else if(mode = 2) then 366.5022 else 824.63

G=E/(2*(1+nu))


C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C22 = C11
C33 = C11

C12 = E*nu/(1+nu)/(1-2*nu)
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
gyz=(dy(w)+dz(v))
gxz=(dz(u)+dx(w))

!!Stress via Hooke's law
!Axial Stress
sx = C11*ex+C12*ey+C13*ez
sy = C21*ex+C22*ey+C23*ez
sz = C31*ex+C32*ey+C33*ez
!Shear stress
sxy=G*gxy
sxz=G*gxz
syz=G*gyz

  { scaling factor for displacement plots } 
   Mt =0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W)) 
EQUATIONS
!FNet = 0; !apply gravity in the direction you want a mode just to get it started
u:	dx(sx)+dy(sxy)+dz(sxz) = -rho*lambda*u
v:	dx(sxy)+dy(sy)+dz(syz) =-rho*lambda*v
w:	dx(sxz)+dy(syz)+dz(sz) =-rho*lambda*w 

CONSTRAINTS
   integral(U)=0               { eliminate translations } 

integral(V)=0 

!   integral(W)=0 

   integral(dx(V)-dy(U)) = 0   { eliminate rotations } 

 !  integral(dy(W) - dz(V)) = 0 
!Note: if you don't eliminate these rotations you can find odd-torsional modes too; this makes it look less like the string but if you really pin supported it these modes would exist!
  ! integral(dz(U) - dx(W))  = 0 

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,R) 
	value(w) = 0
	arc(center=0,0) angle=360
! TIME 0 TO 1    { if time dependent }
PLOTS            { save result displays }
  grid(x+Mt*U,y+Mt*V,z+Mt*W) as "Shape" 
  contour(w) painted on z=Lz/2
SUMMARY
!report f !in modal analysis, the reports are automatically done for each mode
END
