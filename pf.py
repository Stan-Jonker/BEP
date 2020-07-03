import sympy as sp
import numpy as np
# import sympy.vector as spv
import time
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from numpy.linalg import solve, norm
#from scipy.linalg import norm
from scipy.sparse import diags
import math

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 6),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)


def create_powermismatch(G_bus, B_bus, loads, gens, slack):


    # Put nodes in order of load, gen, slack

    # Voorbeeld Y_bus, componenten G en B zijn nodig dus die apart neerzetten. 

    # Variabelen opslaan in een dictionary zodat het makkelijker is om bij de entries te komen als de nodes niet 0, 1, 2, etc. heten
    VA = []
    VM = []
    P = []
    Q = []

    nr_of_nodes = len(loads) + len(gens) + len(slack)
    x_init = list(range(2*nr_of_nodes))

    Node_index = {}

    # Maak alle variablen aan. Gebruik index voor als de nodes niet als een index gebruikt kunnen worden. 
    index = 0
    name = 0
    for l in [loads, gens, slack]:
        bustype = ['l', 'g', 's'][name]
        for i in l:
            Node_index[i] = index

            VA.append(sp.Symbol("VA{}_{}".format(bustype,i)))
            VM.append(sp.Symbol("VM{}_{}".format(bustype,i)))
            P.append(sp.Symbol("P{}_{}".format(bustype,i)))
            Q.append(sp.Symbol("Q{}_{}".format(bustype,i)))

            # x_init[index] = VA[i]
            # x_init[index + nr_of_nodes] = VM[i]

            index +=1
        name += 1

    # Categorizeer alles in load, gen of slack. Dit later gebruiken.
    VAl = VA[0:len(loads)]
    VAg = VA[len(loads):(len(loads) + len(gens))]
    VAs = [VA[-1]]
    VMl = VM[0:len(loads)]
    VMg = VM[len(loads):(len(loads) + len(gens))]
    VMs = [VM[-1]]
    Pl = P[0:len(loads)]
    Pg = P[len(loads):(len(loads) + len(gens))]
    Ps = [P[-1]]
    Ql = Q[0:len(loads)]
    Qg = Q[len(loads):(len(loads) + len(gens))]
    Qs = [Q[-1]]

    # Maak vector x met eerst de delta_i en daarna de |V_i|
    # x = sp.Transpose(sp.Matrix(x_init))

    # Creeer componenten voor de power mismatch function
    PQ = list(range(2*nr_of_nodes))
    PQ_x = list(range(2*nr_of_nodes))

    # Maak jacobiaan specifiek. Rijen eerst P dan Q, kolommen eerst d/d delta, dan d/d |V|
    Jacobian = sp.eye(2*nr_of_nodes)

    
    for index_i in range(nr_of_nodes):
        
        Pi_x = 0
        Qi_x = 0

        for index_j in range(nr_of_nodes):
            # Loop door alle elementen heen
            if (G_bus[index_i, index_j] != 0) and (B_bus[index_i, index_j] != 0):
                delta_ij = VA[index_i] - VA[index_j]
                Pij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.cos(delta_ij) + B_bus[index_i,index_j]*sp.sin(delta_ij))
                Qij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.sin(delta_ij) - B_bus[index_i,index_j]*sp.cos(delta_ij))

                # Creeer de Pi_x door te sommeren.
                Pi_x = Pi_x + Pij_x_part*VM[index_j]
                Qi_x = Qi_x + Qij_x_part*VM[index_j]

                
                if index_i is not index_j:
                    # maak de dPi/dδj, i =/= j in de jacobiaan
                    Jacobian[index_i, index_j] = Qij_x_part*VM[index_j]

                    # maak de dQi/dδj, i =/= j in de jacobiaan
                    Jacobian[nr_of_nodes + index_i, index_j] = -Pij_x_part*VM[index_j]

                    # maak de dPi/d|Vj| i =/= j in de jacobiaan
                    Jacobian[index_i, nr_of_nodes + index_j] = Pij_x_part

                    # maak de dQi/d|Vj| i =/= j in de jacobiaan 
                    Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qij_x_part

        # Maak jacobaan dPi/dδi
        Jacobian[index_i, index_i] = -Qi_x - B_bus[index_i, index_i]*VM[index_i]**2
        
        # Maak jacobaan dQi/dδi
        Jacobian[nr_of_nodes + index_i, index_i] = Pi_x - G_bus[index_i, index_i]*VM[index_i]**2

        # dPi/d|Vi|
        Jacobian[index_i, nr_of_nodes + index_i] = Pi_x/VM[index_i] - G_bus[index_i, index_i]*VM[index_i]

        # dQi/d|Vi|
        Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qi_x/VM[index_i] - B_bus[index_i, index_i]*VM[index_i]


        PQ[index_i] = P[index_i]
        PQ[index_i + nr_of_nodes] = Q[index_i]

        PQ_x[index_i] = Pi_x
        PQ_x[index_i + nr_of_nodes] = Qi_x


    # Maak de power mismatch function
    F  = sp.Matrix(PQ) - sp.Matrix(PQ_x)
    
    known = VMg + VMs + VAs + Pl + Pg + Ql
    unknown = VMl + VAl + VAg + Ps + Qg + Qs

    #Fx = sp.lambdify([PQ, x_init], F, "numpy")
    Fx = sp.lambdify([unknown, known], F, 'numpy')

    # Maak de jacobiaan
    # Jacobian = F.jacobian(x_init)
    #J = sp.SparseMatrix(Jacobian)
    J = F.jacobian(unknown)
    #Jx = sp.lambdify([x_init], J, "numpy")
    Jx = sp.lambdify([unknown, known], J, 'numpy')

    #print(Jacobian)

    return F, Jacobian, Fx, Jx, VM, VA, P, Q, Node_index

def create_reduced_powermismatch_v1(G_bus, B_bus, loads, gens, slack):


    # Put nodes in order of load, gen, slack

    # Voorbeeld Y_bus, componenten G en B zijn nodig dus die apart neerzetten. 

    # Variabelen opslaan in een dictionary zodat het makkelijker is om bij de entries te komen als de nodes niet 0, 1, 2, etc. heten
    VA = []
    VM = []
    P = []
    Q = []

    nr_of_nodes = len(loads) + len(gens) + len(slack)
    #x_init = list(range(2*nr_of_nodes))

    Node_index = {}

    # Maak alle variablen aan. Gebruik index voor als de nodes niet als een index gebruikt kunnen worden. 
    index = 0
    name = 0
    for l in [loads, gens, slack]:
        bustype = ['l', 'g', 's'][name]
        for i in l:
            Node_index[i] = index

            VA.append(sp.Symbol("VA{}_{}".format(bustype,i)))
            VM.append(sp.Symbol("VM{}_{}".format(bustype,i)))
            P.append(sp.Symbol("P{}_{}".format(bustype,i)))
            Q.append(sp.Symbol("Q{}_{}".format(bustype,i)))

            # x_init[index] = VA[i]
            # x_init[index + nr_of_nodes] = VM[i]

            index +=1
        name += 1

    # Categorizeer alles in load, gen of slack. Dit later gebruiken.
    VAl = VA[0:len(loads)]
    VAg = VA[len(loads):(len(loads) + len(gens))]
    VAs = [VA[-1]]
    VMl = VM[0:len(loads)]
    VMg = VM[len(loads):(len(loads) + len(gens))]
    VMs = [VM[-1]]
    Pl = P[0:len(loads)]
    Pg = P[len(loads):(len(loads) + len(gens))]
    Ps = [P[-1]]
    Ql = Q[0:len(loads)]
    Qg = Q[len(loads):(len(loads) + len(gens))]
    Qs = [Q[-1]]

    # Maak vector x met eerst de delta_i en daarna de |V_i|
    # x = sp.Transpose(sp.Matrix(x_init))

    # Creeer componenten voor de power mismatch function
    PQ = list(range(2*nr_of_nodes))
    PQ_x = list(range(2*nr_of_nodes))

    # Maak jacobiaan specifiek. Rijen eerst P dan Q, kolommen eerst d/d delta, dan d/d |V|
    Jacobian = sp.eye(2*nr_of_nodes)

    
    for index_i in range(nr_of_nodes):
        
        Pi_x = 0
        Qi_x = 0

        for index_j in range(nr_of_nodes):
            # Loop door alle elementen heen
            if (G_bus[index_i, index_j] != 0) and (B_bus[index_i, index_j] != 0):
                delta_ij = VA[index_i] - VA[index_j]
                Pij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.cos(delta_ij) + B_bus[index_i,index_j]*sp.sin(delta_ij))
                Qij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.sin(delta_ij) - B_bus[index_i,index_j]*sp.cos(delta_ij))

                # Creeer de Pi_x door te sommeren.
                Pi_x = Pi_x + Pij_x_part*VM[index_j]
                Qi_x = Qi_x + Qij_x_part*VM[index_j]

                
                if index_i is not index_j:
                    # maak de dPi/dδj, i =/= j in de jacobiaan
                    Jacobian[index_i, index_j] = Qij_x_part*VM[index_j]

                    # maak de dQi/dδj, i =/= j in de jacobiaan
                    Jacobian[nr_of_nodes + index_i, index_j] = -Pij_x_part*VM[index_j]

                    # maak de dPi/d|Vj| i =/= j in de jacobiaan
                    Jacobian[index_i, nr_of_nodes + index_j] = Pij_x_part

                    # maak de dQi/d|Vj| i =/= j in de jacobiaan 
                    Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qij_x_part

        # Maak jacobaan dPi/dδi
        Jacobian[index_i, index_i] = -Qi_x - B_bus[index_i, index_i]*VM[index_i]**2
        
        # Maak jacobaan dQi/dδi
        Jacobian[nr_of_nodes + index_i, index_i] = Pi_x - G_bus[index_i, index_i]*VM[index_i]**2

        # dPi/d|Vi|
        Jacobian[index_i, nr_of_nodes + index_i] = Pi_x/VM[index_i] - G_bus[index_i, index_i]*VM[index_i]

        # dQi/d|Vi|
        Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qi_x/VM[index_i] - B_bus[index_i, index_i]*VM[index_i]


        PQ[index_i] = P[index_i]
        PQ[index_i + nr_of_nodes] = Q[index_i]

        PQ_x[index_i] = Pi_x
        PQ_x[index_i + nr_of_nodes] = Qi_x


    known = VMg + VMs + VAs + Pl + Pg + Ql
    unknown = VAl + VAg+ VMl + Ps + Qg + Qs
    argument_unknown = VAl + VAg + VMl

    # Maak de power mismatch function
    F  = sp.Matrix(PQ) - sp.Matrix(PQ_x)
    Fx = sp.lambdify([unknown, known], F, 'numpy')

    # Gebruik alleen de delen die nuttig zijn

    powermismatch = sp.Matrix(F[0:len(loads)] + F[len(loads):(len(loads)+len(gens))] + F[nr_of_nodes:(nr_of_nodes + len(loads))])
    #print(powermismatch)
    J_pmm = powermismatch.jacobian(argument_unknown)
    #print(J_pmm)

    powermismatch_x = sp.lambdify([argument_unknown, known], powermismatch, 'numpy')
    J_pmm_x = sp.lambdify([argument_unknown, known], J_pmm, 'numpy')

    # Maak functies voor de slack en gen P en Q evalueren
    
    Ps_x = sp.lambdify([argument_unknown, known], PQ_x[len(loads)+len(gens)], 'numpy')
    Qg_x = sp.lambdify([argument_unknown, known], PQ_x[(nr_of_nodes+len(loads)):(nr_of_nodes+len(loads)+len(gens))], 'numpy')
    Qs_x = sp.lambdify([argument_unknown, known], PQ_x[nr_of_nodes+len(loads)+len(gens)], 'numpy')

    PQ_f = sp.lambdify([argument_unknown, known], PQ_x, 'numpy')

    #print(Jacobian)

    return powermismatch, powermismatch_x, J_pmm, J_pmm_x, Ps_x, Qg_x, Qs_x, Fx, PQ_f, VM, VA, P, Q, Node_index

def create_reduced_powermismatch(G_bus, B_bus, loads, gens, slack):

    # Put nodes in order of load, gen, slack

    # Voorbeeld Y_bus, componenten G en B zijn nodig dus die apart neerzetten. 

    # Variabelen opslaan in een dictionary zodat het makkelijker is om bij de entries te komen als de nodes niet 0, 1, 2, etc. heten
    VA = []
    VM = []
    P = []
    Q = []

    nr_of_nodes = len(loads) + len(gens) + len(slack)
    x_init = list(range(2*nr_of_nodes))

    Node_index = {}

    # Maak alle variablen aan. Gebruik index voor als de nodes niet als een index gebruikt kunnen worden. 
    index = 0
    name = 0
    for l in [loads, gens, slack]:
        bustype = ['l', 'g', 's'][name]
        for i in l:
            Node_index[i] = index

            VA.append(sp.Symbol("VA{}_{}".format(bustype,i)))
            VM.append(sp.Symbol("VM{}_{}".format(bustype,i)))
            if bustype == 'l' or bustype == 'g':
                P.append(sp.Symbol("P{}_{}".format(bustype,i)))
            if bustype == 'l':
                Q.append(sp.Symbol("Q{}_{}".format(bustype,i)))

            # x_init[index] = VA[i]
            # x_init[index + nr_of_nodes] = VM[i]

            index +=1
        name += 1

    # Categorizeer alles in load, gen of slack. Dit later gebruiken.
    VAl = VA[0:len(loads)]
    VAg = VA[len(loads):(len(loads) + len(gens))]
    VAs = [VA[-1]]
    VMl = VM[0:len(loads)]
    VMg = VM[len(loads):(len(loads) + len(gens))]
    VMs = [VM[-1]]
    # Alleen deze drie zijn nodig 
    Pl = P[0:len(loads)]
    Pg = P[len(loads):(len(loads) + len(gens))]
    Ql = Q[0:len(loads)]
    


    # Maak vector x met eerst de delta_i en daarna de |V_i|
    # x = sp.Transpose(sp.Matrix(x_init))

    # Creeer componenten voor de power mismatch function
    PQ = list(range(2*nr_of_nodes))
    PQ_x = list(range(2*nr_of_nodes))

    # Maak jacobiaan specifiek. Rijen eerst P dan Q, kolommen eerst d/d delta, dan d/d |V|
    Jacobian = sp.eye(2*nr_of_nodes)

    
    for index_i in range(nr_of_nodes):
        
        Pi_x = 0
        Qi_x = 0

        for index_j in range(nr_of_nodes):
            # Loop door alle elementen heen
            if (G_bus[index_i, index_j] != 0) and (B_bus[index_i, index_j] != 0):
                delta_ij = VA[index_i] - VA[index_j]
                Pij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.cos(delta_ij) + B_bus[index_i,index_j]*sp.sin(delta_ij))
                Qij_x_part = VM[index_i]*(G_bus[index_i,index_j]*sp.sin(delta_ij) - B_bus[index_i,index_j]*sp.cos(delta_ij))

                # Creeer de Pi_x door te sommeren.
                Pi_x = Pi_x + Pij_x_part*VM[index_j]
                Qi_x = Qi_x + Qij_x_part*VM[index_j]

                
                if index_i is not index_j:
                    # maak de dPi/dδj, i =/= j in de jacobiaan
                    Jacobian[index_i, index_j] = Qij_x_part*VM[index_j]

                    # maak de dQi/dδj, i =/= j in de jacobiaan
                    Jacobian[nr_of_nodes + index_i, index_j] = -Pij_x_part*VM[index_j]

                    # maak de dPi/d|Vj| i =/= j in de jacobiaan
                    Jacobian[index_i, nr_of_nodes + index_j] = Pij_x_part

                    # maak de dQi/d|Vj| i =/= j in de jacobiaan 
                    Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qij_x_part

        # Maak jacobaan dPi/dδi
        Jacobian[index_i, index_i] = -Qi_x - B_bus[index_i, index_i]*VM[index_i]**2
        
        # Maak jacobaan dQi/dδi
        Jacobian[nr_of_nodes + index_i, index_i] = Pi_x - G_bus[index_i, index_i]*VM[index_i]**2

        # dPi/d|Vi|
        Jacobian[index_i, nr_of_nodes + index_i] = Pi_x/VM[index_i] - G_bus[index_i, index_i]*VM[index_i]

        # dQi/d|Vi|
        Jacobian[nr_of_nodes + index_i, nr_of_nodes + index_j] = Qi_x/VM[index_i] - B_bus[index_i, index_i]*VM[index_i]


        PQ[index_i] = P[index_i]
        PQ[index_i + nr_of_nodes] = Q[index_i]

        PQ_x[index_i] = Pi_x
        PQ_x[index_i + nr_of_nodes] = Qi_x


    # Maak de power mismatch function
    F  = sp.Matrix(PQ) - sp.Matrix(PQ_x)
    
    known = VMg + VMs + VAs + Pl + Pg + Ql
    unknown = VMl + VAl + VAg #+ Ps + Qg + Qs

    #Fx = sp.lambdify([PQ, x_init], F, "numpy")
    Fx = sp.lambdify([unknown, known], F, 'numpy')

    # Maak de jacobiaan
    # Jacobian = F.jacobian(x_init)
    #J = sp.SparseMatrix(Jacobian)
    J = F.jacobian(unknown)
    #Jx = sp.lambdify([x_init], J, "numpy")
    Jx = sp.lambdify([unknown, known], J, 'numpy')

    #print(Jacobian)

    return F, Jacobian, Fx, Jx, VM, VA, P, Q, Node_index

def test_speed():
    max_nodes = 20
    min_nodes = 1
    n = min_nodes

    times = []
    while n <= max_nodes:
        nodes = list(range(n))
        G_bus = sp.ones(n)
        B_bus = sp.ones(n)

        t_start = time.time()

        F, Jaco, Fxx, Jxx, VM, VA, P, Q, Node_index = create_powermismatch(G_bus, B_bus, nodes)

        times.append(time.time() - t_start)

        n+=1

    plt.plot(range(min_nodes, max_nodes+1), times, 'o')
    plt.show()

def reshapeFx(unknown, known, Fx, Jx):
    # Reshape voor gebruik van fsolve
    res = Fx(unknown, known)
    return np.reshape(res, res.shape[0])

def evalJx(unknown, known, Fx, Jx):
    # Evalueer met goede functie entries voor fsolve
    res = Jx(unknown, known)
    # print(res)
    return res

def apply_fsolve(loads, gens, slack, G_bus, B_bus, known, x0):

    #

    F, J, Fx, Jx, VM, VA, P, Q, Node_index = create_powermismatch(G_bus, B_bus, loads, gens, slack)
    
    res = fsolve(func = reshapeFx, x0 = x0, args = (known, Fx, Jx), fprime = evalJx, xtol = 1e-10)

    return res, -1, 1

def apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, known, x0):

    powermismatch, powermismatch_x, J_pmm, J_pmm_x, Ps_x, Qg_x, Qs_x, Fx, PQ_f, VM, VA, P, Q, Node_index = create_reduced_powermismatch_v1(G_bus, B_bus, loads, gens, slack)
    
    # res bestaat uit VAl + VAg + VMl
    res, infodict, _, _ = fsolve(func = reshapeFx, x0 = x0, args = (known, powermismatch_x, J_pmm_x), fprime = evalJx, full_output = 1)
    #Ps = Ps_x(res, known)
    #Qg = Qg_x(res, known)
    #Qs = Qs_x(res, known)
    #now_known = np.concatenate((res, np.array([Ps]), Qg, np.array([Qs])))
    #print("DEZE")
    #print(Fx(now_known, known))
    return res, infodict['nfev'], norm(powermismatch_x(res, known))

def apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, known, x0):
    powermismatch, powermismatch_x, J_pmm, J_pmm_x, Ps_x, Qg_x, Qs_x, Fx, PQ_f, VM, VA, P, Q, Node_index = create_reduced_powermismatch_v1(G_bus, B_bus, loads, gens, slack)
    max_iter = 10000
    tol = 1e-7
    x = x0
    curr_iter = 1

    while norm(powermismatch_x(x, known)) > tol and curr_iter < max_iter:
        A = J_pmm_x(x, known)
        B = reshapeFx(x, known, powermismatch_x, J_pmm_x)
        s = solve(A, -B)
        x += s

        curr_iter += 1
    #print('NR performed {} iterations'.format(curr_iter))
    curr_norm = norm(powermismatch_x(x, known))
    if curr_iter >= max_iter:
        
        print('Current norm was {}'.format(curr_norm))
    return x, curr_iter, curr_norm

def apply_line_search_to_reduced_system(method, loads, gens, slack, G_bus, B_bus, known, x0):
    # Setup
    powermismatch, powermismatch_x, J_pmm, J_pmm_x, Ps_x, Qg_x, Qs_x, Fx, PQ_f, VM, VA, P, Q, Node_index = create_reduced_powermismatch_v1(G_bus, B_bus, loads, gens, slack)
    # Definieer alle iteraties (horen eigenlijk een functieargument te zijn)
    max_iter = 10000
    max_lambda_iter = 1000
    tol = 1e-7
    x = x0
    curr_iter = 1
    # Constantes voor line search
    # vermenigvuldig hier lambda mee
    lambda_constant = .5
    # deze constante voor Armijo rule
    alpha = 1e-5

    funceval = 0

    curr_norm = norm(powermismatch_x(x, known))
    funceval += 1
    while curr_norm > tol and curr_iter < max_iter:
        curr_iter += 1
        # Solve de Newton step
        A = J_pmm_x(x, known)
        B = reshapeFx(x, known, powermismatch_x, J_pmm_x)
        funceval += 2
        s = solve(A, -B)
        # Gebruik hier de methode om lambda telkens met een vaste waarde te verminderen totdat we iets vinden
        # Indien dat niet mogelijk is, doe een newton stap %%%%
        if method == 'n':
            l = 1
            lambda_iter = 0
            while lambda_iter < max_lambda_iter:
                funceval += 1
                if norm(powermismatch_x(x+l*s, known)) <= math.sqrt(1-2*alpha*l)*curr_norm:
                    x = x + l*s
                    break
                l = lambda_constant*l
                lambda_iter += 1
            if lambda_iter == max_lambda_iter:
                x = x + s
                print('A newton step was done')
            
        # Gebruik hier een methode zoals in het boek beschreven staat, waar we zoeken naar het minimum van een parabool
        if method == 'p':
            # Initiële waarden voor lambda
            l = 1
            l_prev = 1
            lambda_iter = 0
            
            # Maak waarden voor parabool
            g_zero = (norm(powermismatch_x(x, known))**2)/2
            funceval += 1
            g_prime_0 = -2*g_zero

            while lambda_iter < max_lambda_iter:
                # Check voor Armijo regel
                funceval += 1
                if norm(powermismatch_x(x+l*s, known)) <= math.sqrt(1-2*alpha*l)*curr_norm:
                    x = x + l*s
                    # print(lambda_iter)
                    break
                    
                # Reken top uit
                g_lambda = (norm(powermismatch_x(x+l*s, known))**2)/2
                funceval += 1
                top = -g_prime_0*l/(2*((g_lambda-g_zero)/l - g_prime_0))
                # Safeguarding
                if top < 0.5*l and top > 0.1*l:
                    l_prev = l
                    l = top
                else:
                    l_prev = l
                    funceval += 1
                    if norm(powermismatch_x(x + 0.1*l*s, known)) < norm(powermismatch_x(x + 0.5*l*s, known)):
                        l = 0.1*l
                    else:
                        l = 0.5*l
                lambda_iter += 1
            # Als het niet goed genoeg is, doe een Newton stap
            if lambda_iter >= max_lambda_iter - 1:
                x = x + s
                print('A newton step was done') 

        funceval += 1
        temp_norm = norm(powermismatch_x(x, known))
        # Als het algoritme geen goede hoeveelheid progress kan maken, stop
        #if abs(temp_norm - curr_norm)<0.00001*tol:
        #    print('No significant progress was made, the current norm is {}'.format(temp_norm))
        #    curr_iter -= 1
        #    break
        #else:
        #    curr_norm = temp_norm
        curr_norm = temp_norm

    #print('line search performed {} iterations'.format(curr_iter))
    return x, curr_iter, curr_norm


#known = VMg + VMs + VAs + Pl + Pg + Ql
#unknown = VAl + VAg+ VMl + Ps + Qg + Qs
#argument_unknown = VAl + VAg + VMl
#print("___ NIET PER UNIT____")
loads = [1, 2]
gens = []
slack = [0]
#G_bus = np.array([[2*0.991746, -0.991746, -0.991746], [-0.991746, 0.991746, 0], [-0.991746, 0, 0.991746]])
#B_bus = np.array([[2*0.128216, -0.128216, -0.128216], [-0.128216, 0.128216, 0], [-0.128216, 0, 0.128216]])

#print('----------LINE-SEARCH-----------')
#print(apply_line_search_to_reduced_system('p', loads, gens, slack, G_bus, B_bus, [50000, 0, -1e6, -1e6, -0.25*1e6, -0.25*1e6], [0,0,5e4,5e4]))
#print(apply_line_search_to_reduced_system('n', loads, gens, slack, G_bus, B_bus, [50000, 0, -1e6, -1e6, -0.25*1e6, -0.25*1e6], [0,0,5e4,5e4]))
#print('----------NR-----------')
#
#for i in range(40):
#    x, iters_1[i], h = apply_line_search_to_reduced_system('p', loads, gens, slack, G_bus, B_bus, [50000, 0, -1e6*i/20, -1e6, -0.25*1e6*i/20, -0.25*1e6], [0,0,5e4,5e4])
#print(apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, [50000, 0, -1e6, -1e6, -0.25*1e6, -0.25*1e6], [0,0,5e4,5e4]))

#res, Ps, Qg, Qs, known = apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, [50000, 0, -1e6, -1e6, -0.25*1e6, -0.25*1e6], [0,0,5e4,5e4])
#print(res)
#print(Ps)
#print(Qg)
#print(Qs)

print("_____PER UNIT______")

#known = VMg + VMs + VAs + Pl + Pg + Ql
#unknown = VAl + VAg+ VMl + Ps + Qg + Qs
#argument_unknown = VAl + VAg + VMl

k = 1e8
print(k)
G_bus = np.array([[(1+k)*0.991746/4e-4, k*-0.991746/4e-4, -0.991746/4e-4], [k*-0.991746/4e-4, k*0.991746/4e-4, 0], [-0.991746/4e-4, 0, 0.991746/4e-4]])
B_bus = np.array([[(1+k)*0.128216/4e-4, k*-0.128216/4e-4, -0.128216/4e-4], [k*-0.128216/4e-4, k*0.128216/4e-4, 0], [-0.128216/4e-4, 0, 0.128216/4e-4]])


print('----------FSOLVE----------------')
print(apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, [1, 0, -1, -1, -0.1, -0.1], [0,0,1,1]))
print('----------LINE-SEARCH-----------')
print(apply_line_search_to_reduced_system('p', loads, gens, slack, G_bus, B_bus, [1, 0, -1, -1, -0.1, -0.1], [0,0,1,1]))
print(apply_line_search_to_reduced_system('n', loads, gens, slack, G_bus, B_bus, [1, 0, -1, -1, -0.1, -0.1], [0,0,1,1]))
print('----------NR-----------')
print(apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, [1, 0, -1, -1, -0.25, -0.25], [0,0,1,1]))
#k = 10
#factor = 40
#iters_1 = list(range(k))
#iters_2 = list(range(k))
#iters_3 = list(range(k))
#x = [0,0,1,1]
#x1 = x
#x2 = x
#x3 = x
#for i in range(k):
#    x1, h, iters_1[i] = apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, [1, 0, -1*i*factor, -1*i*factor, -0.25*i*factor, -0.25*i*factor], x1)
#    x2, h, iters_2[i] = apply_line_search_to_reduced_system('p', loads, gens, slack, G_bus, B_bus, [1, 0, -1*i*factor, -1*i*factor, -0.25*i*factor, -0.25*i*factor], x2)
#    x3, h, iters_3[i] = apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, [1, 0, -1*i*factor, -1*i*factor, -0.25*i*factor, -0.25*i*factor], x3)



#plt.plot(factor*np.arange(0, k), iters_1, 'o')
#plt.plot(factor*np.arange(0, k), iters_2, 'o')
#plt.plot(factor*np.arange(0, k), iters_3, 'o')
#plt.xlabel('Power demand bus 1')
#plt.ylabel('Residual norm')
#plt.yscale('log')
#plt.show()

#n = 10+1

#loads = np.arange(1, n)
#gens = []
#slack = [0]

#A = diags([-1, 2, -1], [-1, 0, 1], shape=(n, n))
#G_bus = (0.991746/4e-4*A).toarray()
#B_bus = (0.128216/4e-4*A).toarray()
#known = np.array([1.0, 0.0])
#known = np.append(known, -1.0*np.ones((1, n-1)))
#known = np.append(known, -0.1*np.ones((1, n-1)))
#known = np.reshape(known, known.shape[0])
#known = known.tolist()
#x0 = np.append(np.zeros((1, n-1)), np.ones((1, n-1)))
#x0 = np.reshape(x0, x0.shape[0])
#x0 = x0.tolist()

#x, iters, norm = apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, known, x0)
#print(x)
#print(norm)

#plt.plot(np.arange(3, 20+1, 1), norms)
#plt.show()
#res, Ps, Qg, Qs, known = apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, [1, 0, -1, -1, -0.25, -0.25], [0,0,1,1])
#print(res)
#print(Ps)
#print(Qg)
#print(Qs)








# Hier fsolve gebruiken
#t_start = time.time()
#print(fsolve(func = reshapeFx, x0 = [1.,1.,1.,-1.,-1.,1.], args = ([1.,1.,1.,1.,1.,1.], Fx, Jx), fprime = evalJx))
#print(Fx([1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]))
#F, Jacobian, Fx, Jx, powermismatch, powermismatch_x, VM, VA, P, Q, Node_index = create_reduced_powermismatch_v1(G_bus, B_bus, loads, gens, slack)
#print(powermismatch)

#print(time.time()-t_start)

#print(F)
#print(Jaco)
#print(Fxx)
#print(Jxx)
#print(VM)
#print(VA)
#print(P)
#print(Q)

#print(x)
#print(VA)
#print(Y_bus)

#loads = [0, 1]
#gens = [2]
#slack = [3]
#G_bus = np.array([[3, 0, -1, -2], [0, 3, -1, -2], [-1, -1, 2, 0], [-2, -2, 0, 4]])
#B_bus = np.array([[3, 0, -1, -2], [0, 3, -2, -1], [-1, -2, 3, 0], [-2, -1, 0, 3]])

#print('fsolve')
#apply_fsolve_to_reduced_system(loads, gens, slack, G_bus, B_bus, [230, 230, 0, 2, 3, 4, 0.5, 0.5], [0, 0, 0, 230, 230])
#print('NR')
#print(apply_NR_to_reduced_system(loads, gens, slack, G_bus, B_bus, [230, 230, 0, 2, 3, 4, 0.5, 0.5], [0, 0, 0, 230, 230]))
#print('Line seach (normal)')
#print(apply_line_search_to_reduced_system('n', loads, gens, slack, G_bus, B_bus, [230, 230, 0, 2, 3, 4, 0.5, 0.5], [0, 0, 0, 230, 230]))
#print('Line seach (parabola)')
#print(apply_line_search_to_reduced_system('p', loads, gens, slack, G_bus, B_bus, [230, 230, 0, 2, 3, 4, 0.5, 0.5], [0, 0, 0, 230, 230]))