

{  2D_MOVEPOINT.PDE

 

 This example is a variation  of 2D_STRETCH_XY.PDE demonstrating the use of

 moving and non-moving point declarations.

 

 A point defined by name as a MOVABLE POINT and used in the definition of the domain will move

 with the mesh.

 

 Any point declared explicitly or not used in the domain definition will remain fixed.

 

}

TITLE "stretching brick"

 

SELECT

 regrid=off    

 

VARIABLES

 u

 xm = move(x)

 ym = move(y)

 

DEFINITIONS

 Hl = 1/2

 gwid = 0.15

 u0= exp(-(x^2+y^2)/gwid^2)

 lmove = Hl + t

 ms = gwid^2/u0

 vx = dt(xm)

 vy = dt(ym)

 P = movable point(Hl,Hl)

 Q = movable point(0.1,0)

 R = point(-0.2,-0.2)

 

INITIAL VALUES

 u= u0

 dt(xm) = x/Hl

 dt(ym) = y/Hl

 

EQUATIONS

 U:  dt(u)=0

 Xm:  div(grad(vx))=0

 Ym:  div(grad(vy))=0

 

BOUNDARIES

REGION 1

  mesh_spacing = ms

  START(-Hl,-Hl)

  value(u) = 0 nobc(xm) value(ym)=-lmove

  Line to (Hl,-Hl)

  value(u)=0 value(xm)=lmove nobc(ym)

  line to P

  value(u)=0 nobc(xm) value(ym) = lmove

  line to (-Hl,Hl)

  value(u)=0 value(xm)=-lmove nobc(ym)

  line to close

 

NODE POINT Q

 

TIME 0 TO 0.5 by 0.01! 10

 

MONITORS

for cycle=1

  grid(x,y) zoom(-Hl-1/2,-Hl-1/2, 2*Hl+1,2*Hl+1)

  grid(x,y) zoom(-0.6,0.4, 0.2,0.2)

  contour(vx) zoom(-0.6,0.4, 0.2,0.2)

  contour(vy) zoom(-0.6,0.4, 0.2,0.2)

  contour(u)

  elevation(u,u0) from(-10*Hl,0) to (10*Hl,0) range (0,1)

  elevation(u,u0) from(0,-10*Hl) to (0,10*Hl) range (0,1)

 

PLOTS

for time=0.1 by 0.1 to endtime

  grid(x,y) zoom(-Hl-1/2,-Hl-1/2, 2*Hl+1,2*Hl+1)

      report(distance(P,(0.2,0)))

  contour(u)

  contour(u-u0) as "True Total Error"

  contour(space_error()) as "Estimated Step Error" painted

  elevation(u,u0)   from(-10*Hl,0) to (10*Hl,0) range (0,1)

  elevation(vx) from(-10*Hl,0) to (10*Hl,0) range (0,1)

  elevation(u,u0)   from(0,-10*Hl) to (0,10*Hl) range (0,1)

  elevation(vy) from(0,-10*Hl) to (0,10*Hl) range (0,1)

 

History(u) at P,Q, (0.2,0) as "Points a(P) and b(Q) move with the mesh, c(0.2,0) is fixed in space"

History(u,u0) at R,(-0.2,-0.2) as "both points are fixed in space"

History(distance(P,R)) ! both are movable points, so the distance changes as the mesh moves

              at P,R   ! name the points to get markers on the domain

History(time_error())

 

END
 
