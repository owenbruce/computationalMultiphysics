TITLE 'Static Beam Temp Change'
COORDINATES
cartesian3
VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

Select
quadratic

DEFINITIONS
mag=0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))

Lx=0.1
Ly=.005
Lz=.02
E
nu
rho

G=E/(2*(1+nu))

DeltaTemp = 5

alpha

alphax = alpha
alphay = alpha
alphaz = alpha

!Stiffness matrix components (isotropic material)
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

!mechanical strain
exm=ex - alphax*DeltaTemp
eym=ey - alphay*DeltaTemp
ezm=ez - alphaz*DeltaTemp

!!Stress via Hooke's law
!Axial Stress
sx = C11*exm+C12*eym+C13*ezm
sy = C21*exm+C22*eym+C23*ezm
sz = C31*exm+C32*eym+C33*ezm
!Shear stress
sxy=G*gxy
sxz=G*gxz
syz=G*gyz

EQUATIONS
!FNet = 0; !apply gravity in the direction you want a mode just to get it started
u:	dx(sx)+dy(sxy)+dz(sxz) = 0
v:	dx(sxy)+dy(sy)+dz(syz) = 0
w:	dx(sxz)+dy(syz)+dz(sz) = 0 

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

  	REGION 'Aluminum'
    E = 68.9E9
    nu = 0.33
    alpha = 2.32e-5
    rho = 2700
    	START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0 
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		load(u)=0
		load(v)=0 
		load(w)=0
		LINE TO (Lx,Ly) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,Ly) !x=0 surface
		value(u)=0
		value(v)=0 
		value(w)=0
		LINE TO CLOSE
        
PLOTS
	grid(x+u*mag, y+v*mag, z+w*mag)
SUMMARY
	report val(ex,Lx,0,0)
    report val(ey,Lx,0,0)
    report val(ez,Lx,0,0)
end
