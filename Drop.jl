# Basado en:
# Cleve Moler, Experiments with MATLAB,
# http://www.mathworks.com/moler/exm/chapters/water.pdf
using GLVisualize, GLAbstraction, Colors, Reactive, GeometryTypes, GLWindow

#Parameters

n = 10;                  # grid size
g = 9.8;                 # gravitational constant
dt = 0.01;               # hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           # plot interval
dropstep = 200;          # drop interval
npasos = 400

# Gotas:
function droplet(height::Real, width::Int64)
    # DROPLET  2D Gaussian
    # D = droplet(height,width)
    x = linspace(-1,1, width)
    D = [height*exp(-5*(i^2+j^2)) for i in x, j in x]
    return x, D
end

# Inicializamos arreglos:
x, D = droplet(1.5,5);
H = ones(n+2,n+2);   U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
#ndrop = ceil(rand*ndrops);
nstep = 0;

#First drop
w = size(D,1);
i = Int(ceil(rand()*(n-w)))+(1:w);
j = Int(ceil(rand()*(n-w)))+(1:w);
H[i,j] = H[i,j] + D;

window = glscreen()
timesignal = loop(1:npasos)

function surf(nstep, H,U,V)
  # Reflective boundary conditions
 H[:,1] = H[:,2];      U[:,1] = U[:,2];       V[:,1] = -V[:,2];
 H[:,n+2] = H[:,n+1];  U[:,n+2] = U[:,n+1];   V[:,n+2] = -V[:,n+1];
 H[1,:] = H[2,:];      U[1,:] = -U[2,:];      V[1,:] = V[2,:];
 H[n+2,:] = H[n+1,:];  U[n+2,:] = -U[n+1,:];  V[n+2,:] = V[n+1,:];

  # First half step

  # x direction
  i = 1:n+1;
  j = 1:n;

  # height
  Hx[i,j] = (H[i+1,j+1]+H[i,j+1])/2 - dt/(2*dx)*(U[i+1,j+1]-U[i,j+1]);

  # x momentum
  Ux[i,j] = (U[i+1,j+1]+U[i,j+1])/2 -
           dt/(2*dx)*((U[i+1,j+1].^2./H[i+1,j+1] + g/2*H[i+1,j+1].^2) -
                      (U[i,j+1].^2./H[i,j+1] + g/2*H[i,j+1].^2));

  # y momentum
  Vx[i,j] = (V[i+1,j+1]+V[i,j+1])/2 -
           dt/(2*dx)*((U[i+1,j+1].*V[i+1,j+1]./H[i+1,j+1]) -
                      (U[i,j+1].*V[i,j+1]./H[i,j+1]));

  # y direction
  i = 1:n;
  j = 1:n+1;

  # height
  Hy[i,j] = (H[i+1,j+1]+H[i+1,j])/2 - dt/(2*dy)*(V[i+1,j+1]-V[i+1,j]);

  # x momentum
  Uy[i,j] = (U[i+1,j+1]+U[i+1,j])/2 -
           dt/(2*dy)*((V[i+1,j+1].*U[i+1,j+1]./H[i+1,j+1]) -
                      (V[i+1,j].*U[i+1,j]./H[i+1,j]));
  # y momentum
  Vy[i,j] = (V[i+1,j+1]+V[i+1,j])/2 -
           dt/(2*dy)*((V[i+1,j+1].^2./H[i+1,j+1] + g/2*H[i+1,j+1].^2) -
                      (V[i+1,j].^2./H[i+1,j] + g/2*H[i+1,j].^2));

  # Second half step
  i = 2:n+1;
  j = 2:n+1;

  # height
  H[i,j] = H[i,j] - (dt/dx)*(Ux[i,j-1]-Ux[i-1,j-1]) -
                   (dt/dy)*(Vy[i-1,j]-Vy[i-1,j-1]);
  # x momentum
  U[i,j] = U[i,j] - (dt/dx)*((Ux[i,j-1].^2./Hx[i,j-1] + g/2*Hx[i,j-1].^2) -
                   (Ux[i-1,j-1].^2./Hx[i-1,j-1] + g/2*Hx[i-1,j-1].^2))
                 - (dt/dy)*((Vy[i-1,j].*Uy[i-1,j]./Hy[i-1,j]) -
                   (Vy[i-1,j-1].*Uy[i-1,j-1]./Hy[i-1,j-1]));
  # y momentum
  V[i,j] = V[i,j] - (dt/dx)*((Ux[i,j-1].*Vx[i,j-1]./Hx[i,j-1]) -
                   (Ux[i-1,j-1].*Vx[i-1,j-1]./Hx[i-1,j-1]))
                 - (dt/dy)*((Vy[i-1,j].^2./Hy[i-1,j] + g/2*Hy[i-1,j].^2) -
                   (Vy[i-1,j-1].^2./Hy[i-1,j-1] + g/2*Hy[i-1,j-1].^2));
  return H
end

t = timesignal

bb = Signal(AABB{Float32}(Vec3f0(0), Vec3f0(1)))

_view(visualize(const_lift(surf, t, H,U,V), :surface, boundingbox=bb))

renderloop(window)
