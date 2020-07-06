TITLE 'Assignment5'
SELECT
errlim=0.01
modes = 3

COORDINATES
cartesian3

VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

DEFINITIONS
mag=7e6

Lx=1
Ly=0.04
Lz=0.2
E=68.9E9
nu=0.33
rho=2.7/1000*100^2

G=E/(2*(1+nu))

!omega = 3455 + .01*stage!2614.82+.01*stage !1547.6315 +.0001*stage !2309.2+.01*stage

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


mt = 0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))
EQUATIONS
!FNet = 0; !apply gravity in the direction you want a mode just to get it started
u:	dx(sx)+dy(sxy)+dz(sxz) = -rho*lambda*u
v:	dx(sxy)+dy(sy)+dz(syz) =-rho*lambda*v
w:	dx(sxz)+dy(syz)+dz(sz) =-rho*lambda*w 

CONSTRAINTS
 integral(U) = 0
 !integral(V) = 0
 integral(W) = 0
 integral(dx(V) - dy(U)) = 0
 integral(dy(W) - dz(V)) = 0
 integral(dz(U) - dx(W)) = 0
 
 
EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES
surface 'bottom'
		load(u)=0
		load(v)=0
		load(w)=0
surface 'top'
		load(u)=0 
		load(v)=0
		load(w)=0

  	REGION 1
    	START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0 
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		load(u)=0
		value(v)=0 
		value(w)=0
		LINE TO (Lx,Ly) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,Ly) !x=0 surface
		load(u)=0
		value(v)=0 
		value(w)=0
		LINE TO CLOSE

PLOTS
	!grid(x+u*mag, y+v*mag, z+w*mag)
	grid(x+mt*u, y+mt*v, z+mt*w)
 	!contour(sxz) on surface z=Lz
	!elevation(sx,sy,sz,syz,sxz,sxy) from (0,0,0) to (0,0,Lz)
	!elevation(u) from (Lx/2,Ly/2,0) to(Lx/2,Ly/2,Lz)
	elevation(v) from (Lx/2,Ly/2,0) to(Lx/2,Ly/2,Lz)
	!history(u) at (Lx/2,Ly/2,Lz/2) !middle for the x-displacement since end will be fixed
	!history(v) at (Lx/2,Ly/2,Lz) !tip should have the most y-displacement 
summary
report lambda
	
{	summary
		report val(u,Lx,Ly,Lz)
		report val(v,Lx,Ly,Lz)
		report val(w,Lx,Ly,Lz)
		report val(w,Lx/2,Ly/2,Lz)}
end
