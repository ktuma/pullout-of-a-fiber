from dolfin import *
from mshr import *
import numpy as np

import sys
if len(sys.argv)==1:
   info("Input parameter is missing. Run the script using 'python3 axisym_neoHooke.py muf n [kapparatio] [b]'.")
   sys.exit()

muf = Constant(sys.argv[1])
mum = Constant(1.0)
n = Constant(sys.argv[2])

if len(sys.argv) > 3:
    kapparatio = Constant(sys.argv[3])
else:
    kapparatio = Constant(100.0)

if len(sys.argv) > 4:
    b = Constant(sys.argv[4])
else:
    b = Constant(1.0)

kappaf = kapparatio*muf
kappam = kapparatio*mum

PETScOptions.set('mat_mumps_icntl_14', 1000) # work array, multiple of estimate to allocate
PETScOptions.set('mat_mumps_icntl_24', 1)  # detect null pivots
PETScOptions.set('mat_mumps_cntl_1', 1.0)  # pivoting threshold, this solves to machine precision 
 
comm = MPI.comm_world
rank = MPI.rank(comm)
set_log_level(LogLevel.INFO if rank==0 else LogLevel.INFO)
parameters["std_out_all_processes"] = False
parameters["form_compiler"]["quadrature_degree"] = 8
parameters["form_compiler"]["cpp_optimize"] = True
ffc_opts = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

H = 0.5
R = 1.0
ratio = 10.0
Rf = R/ratio
Hext = 0.2*H

mesh = Mesh()
XDMFFile("meshes/data/mesh2.xdmf").read(mesh)

class Axis(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0) and on_boundary
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], R) and on_boundary
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0) and on_boundary
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[1] - (H + Hext)) < 1e-9 and on_boundary

facets = MeshFunction("size_t", mesh, 1)
facets.set_all(0)
Axis().mark(facets, 1)
Right().mark(facets, 2)
Bottom().mark(facets, 3)
Top().mark(facets, 4)
ds = Measure("ds", subdomain_data=facets)

class FiberDomain(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (0.0, Rf)) and between(x[1], (0.5*H, H + Hext))

class MatrixDomain(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (Rf, R)) or between(x[1], (0, 0.5*H))
axisbndry = Axis()
rightbndry = Right()
bottombndry = Bottom()
topbndry = Top()
fiberdomain = FiberDomain()
matrixdomain = MatrixDomain()

boundaries = MeshFunction("size_t", mesh, 1)
boundaries.set_all(0)
axisbndry.mark(boundaries, 1)
rightbndry.mark(boundaries, 2)
bottombndry.mark(boundaries, 3)
topbndry.mark(boundaries, 4)

boundaries_file = XDMFFile(mesh.mpi_comm(),"boundaries.xdmf")
boundaries_file.write(boundaries)

domains = MeshFunction("size_t", mesh, 2)
domains.set_all(1)
fiberdomain.mark(domains, 1)
matrixdomain.mark(domains, 2) 

domain_file = XDMFFile(mesh.mpi_comm(),"domains.xdmf")
domain_file.write(domains) 

dx = Measure("dx", subdomain_data = domains) 

resultsdir = f'results_case2_kapparatio_{float(kapparatio)}_n_{float(n)}_muf_{float(muf)}'
ufile = XDMFFile(resultsdir+'/u.xdmf')
ufile.parameters["flush_output"] = True
TRZfile = XDMFFile(resultsdir+'/TRZ.xdmf')
TRZfile.parameters["flush_output"] = True
TfRZfile = XDMFFile(resultsdir+'/TfRZ.xdmf')
TfRZfile.parameters["flush_output"] = True
TmRZfile = XDMFFile(resultsdir+'/TmRZ.xdmf')
TmRZfile.parameters["flush_output"] = True
normefile = XDMFFile(resultsdir+'/norme.xdmf')
normefile.parameters["flush_output"] = True

V = VectorFunctionSpace(mesh, 'CG', 2)
u = Function(V)
u_ = TestFunction(V)

V1 = FunctionSpace(mesh, 'CG', 2)
tfrz = Function(V1)
tmrz = Function(V1)
norme = Function(V1)

u0T = Expression("disp", disp=0.0, degree=1)
bc_axis = DirichletBC(V.sub(0), Constant(0), facets, 1)
bc_right = DirichletBC(V, Constant((0,0)), facets, 2)
bc_bottom = DirichletBC(V, Constant((0,0)), facets, 3)
bc_pull = DirichletBC(V.sub(1), u0T, facets, 4)
bcs = [bc_axis, bc_right, bc_pull]

x = SpatialCoordinate(mesh)

def axi_grad(v):
    return as_tensor([[v[0].dx(0), 0, v[0].dx(1)],
                          [0, v[0]/x[0], 0],
                          [v[1].dx(0), 0, v[1].dx(1)]])

I = Identity(3)
Ff = variable(I + axi_grad(u))
B = Ff*Ff.T
C = Ff.T*Ff
E = 0.5*(C - I)
normE = sqrt(tr(E*E))
I1 = tr(B)
I2 = 0.5*(tr(B)*tr(B)-tr(B*B))
I3 = det(B)
I1bar = I1*(I3**(-1.0/3.0))

Wf = kappaf*(I3**(0.5) - 1.0)**2.0 + 0.5*muf*(1.0 + I1bar - 3.0)
Pf = diff(Wf, Ff)

Wm = kappam*(I3**(0.5) - 1.0)**2.0 + 0.5*mum/b*(1.0 + b/n*(I1bar - 3.0))**n
Pm = diff(Wm, Ff)

PfZZ = Pf[2,2]
Tf = I3**(-0.5)*Pf*Ff.T
Tm = I3**(-0.5)*Pm*Ff.T
TfRZ = Tf[0,2]
TfZZ = Tf[2,2]
TmRZ = Tm[0,2]

Eq = inner(Pf, axi_grad(u_))*x[0]*dx(1) + inner(Pm, axi_grad(u_))*x[0]*dx(2)

J = derivative(Eq, u)

info("Solving problem of size: {0:d}".format(V.dim()))

class NSNonlinearProblem(NonlinearProblem):
    def __init__(self, Eq, w, bcs, J, *args, **kwargs):
        NonlinearProblem.__init__(self)
        self.PEq = Eq
        self.PJ = J
        self.bcs = bcs
        form_compiler_parameters = {"optimize": True}
        if "form_compiler_parameters" in kwargs:
            self.form_compiler_parameters = kwargs["form_compiler_parameters"]

        assembler = SystemAssembler(J, Eq, bcs)
        self.w = w
        self.assembler = assembler
                
    def F(self, b, x):
        self.assembler.assemble(b, x)
            
    def J(self, A, x):
        self.assembler.assemble(A)

problem=NSNonlinearProblem(Eq, u, bcs, J, form_compiler_parameters=ffc_opts)
solver=PETScSNESSolver()
prm = solver.parameters

prm['absolute_tolerance'] = 5e-9
prm['relative_tolerance'] = 1e-20
prm['maximum_iterations'] = 80
prm['linear_solver'] = 'mumps'
prm['report'] = True
prm['error_on_nonconvergence'] = False
prm['method']='newtonls'
prm['line_search']='bt' #[basic, bt, cp, l2, nleqerr]

def element_project(mesh, quad=2):
   """ make element-wise projection of function f to DG0 space"""
   
   e=FiniteElement("DG", mesh.ufl_cell(), 0)
   DG0=FunctionSpace(mesh, e)
   u = TrialFunction(DG0)
   v = TestFunction(DG0)
   A = u*v*dx(1, metadata={"quadrature_degree": quad}) + u*v*dx(2, metadata={"quadrature_degree": quad})
   b = TfRZ*v*dx(1, metadata={"quadrature_degree": quad}) + TmRZ*v*dx(2, metadata={"quadrature_degree": quad})
   solver = LocalSolver(A, b)
   solver.factorize()
   Pf=Function(DG0)
   solver.solve_local_rhs(Pf)
   return(Pf)

e = FiniteElement("DG", mesh.ufl_cell(), 0)
DG0 = FunctionSpace(mesh, e)
trz = Function(DG0)

def save_sol(i):
   u.rename("u", "displacement")
   ufile.write(u, i)
   
   trz.assign(element_project(mesh))
   trz.rename("TRZ", "stress")
   TRZfile.write(trz, i)
   
   tfrz.assign(project(TfRZ, V1))
   tfrz.rename("TfRZ", "stress")
   TfRZfile.write(tfrz, i)
   
   tmrz.assign(project(TmRZ, V1))
   tmrz.rename("TmRZ", "stress")
   TmRZfile.write(tmrz, i)

   norme.assign(project(normE, V1))
   norme.rename("normE", "normE")
   normefile.write(norme, i)


def parallel_eval(f, x):
   mesh = f.function_space().mesh()
   bb = mesh.bounding_box_tree()

   p=Point(x)
   value=np.array([0.0])
   ic=0
   cf=bb.compute_first_entity_collision(p)
   inside= cf < mesh.num_cells()
   if inside :
      f.eval_cell(value,x,Cell(mesh,cf))
      ic=1
       
   comm=MPI.comm_world
   v= MPI.sum(comm, value) / MPI.sum(comm, ic)
   return(v)


info("muf = {}".format(float(muf)))
info("n = {}".format(float(n)))
info("kapparatio = {}".format(float(kapparatio)))
info("b = {}".format(float(b)))

# Compute
for i in np.arange(0.000,0.12501,0.001): 
   i = round(i, 3)
   info("u0T = {}".format(i))
   u0T.disp = i
   solver.solve(problem, u.vector())
   Ftop = assemble(2.0*pi*x[0]*PfZZ*ds(4))
   info("u0T, Force = {{{0:f}, {1:f}}},".format(i, float(Ftop)))   
   if 1000*i % 25 == 0: save_sol(i)
