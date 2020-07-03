import pandapower as pp
import pandapower.networks
import pandapower.topology
import pandapower.plotting
import pandapower.converter
import pandapower.estimation
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

nr_of_steps = 20

contin = True


iterations = np.empty([nr_of_steps, 1])

while contin:
    



    mogelijke_vars = ["q1", "q2", "pload1", "pload2", "phi1", "phi2"]
    # modussen
    
    q1 = 0.25
    q2 = 0.25
    p_load_1 = 1
    p_load_2 = 1
    q_load_1 = .25
    q_load_2 = .25
    s_load_1 = 1 #p_load_1**2 + q_load_1**2
    s_load_2 = 1 #p_load_2**2 + q_load_2**2

    while 1:
        var = str(input("Which variable? Choose from q1/q2/pload1/pload2/phi1/phi2: "))
        if var in mogelijke_vars:
            break
        else:
            print("Input something else")
    while 1:
        try:   
            start_var = float(input("Start value for variable: "))
        except:
            print("Invalid input")
        break
    while 1:
        try:   
            end_var = float(input("End value for variable: "))
        except:
            print("Invalid input")
        break
    
    varlist = np.linspace(start = start_var, stop = end_var, num = nr_of_steps, endpoint = True)
    
    while 1:
        try:
            modus = int(input("Welke modus?\n 1 - 2 loads los aan een generator\n 2 - een generator met 2 nodes aan een tussennode\n 3 - loads en generator in een lijn \n input (1/2/3):"))
            if modus in [1, 2, 3]:
                break
            else:
                print("input another number")
        except:
            print("Input invalid")
    i = 0

    
    try:
        max_iteration = int(input("Max aantal iteraties (default naar 10)"))
    except:
        print("Dat was geen goede input, algoritme neemt default")
        max_iteration = 'auto'

    if not modus == 3:
        connect23input = str(input("Connect the two loads? (y/n): "))
    
    if modus == 2:
        nr_of_nodes = 4
    else:
        nr_of_nodes = 3
    
    vm_pu_res = np.empty([nr_of_nodes, nr_of_steps])
    va_degree_res = np.empty([nr_of_nodes, nr_of_steps])
    p_mw_res = np.empty([nr_of_nodes, nr_of_steps])
    q_mvar_res = np.empty([nr_of_nodes, nr_of_steps])


    while i < nr_of_steps:
        
        if var == "q1":
            q_load_1 = varlist[i]
            varlabel = "$Q_1$"
        elif var == "q2":
            q_load_2 = varlist[i]
            varlabel = "$Q_2$"
        elif var == "pload1":
            p_load_1 = varlist[i]
            varlabel = "$P_1$"
        elif var == "pload2":
            p_load_2 = varlist[i]
            varlabel = "$P_2$"
        elif var == "phi1":
            p_load_1 = s_load_1*math.cos(varlist[i])
            q_load_1 = s_load_1*math.sin(varlist[i])
            varlabel = "power phase angle $\phi$ bus 1"
        elif var == "phi2":
            p_load_2 = s_load_2*math.cos(varlist[i])
            q_load_2 = s_load_2*math.sin(varlist[i])
            varlabel = "power phase angle $\phi$ bus 2"

        
        p_gen_1 = p_load_1 + p_load_2

        net = pp.create_empty_network(name = "net", f_hz = 50, sn_mva = 1, add_stdtypes=True)
        if modus == 1:
            
            
            bus_1 = pp.create_bus(net, vn_kv = 50, name = "Bus 1", geodata = (0, 0), type = "b", zone = None, in_service = True)
            bus_2 = pp.create_bus(net, vn_kv = 50, name = "Bus 2", geodata = (1, 1), type = "b", zone = None, in_service = True)
            bus_3 = pp.create_bus(net, vn_kv = 50, name = "Bus 3", geodata = (1, -1), type = "b", zone = None, in_service = True)

            line_12 = pp.create_line(net = net, from_bus = bus_1, to_bus = bus_2, length_km = 1, std_type = "NAYY 4x50 SE")
            line_13 = pp.create_line(net = net, from_bus = bus_1, to_bus = bus_3, length_km = 1, std_type = "NAYY 4x50 SE")

            
            if connect23input == "y":
                # Voor als de lijn tussen 2 en 3 wel meedoet
                line_23 = pp.create_line(net = net, from_bus = bus_2, to_bus = bus_3, length_km = 1, std_type = "NAYY 4x50 SE")
            


        elif modus == 2:
            bus_0 = pp.create_bus(net, vn_kv = 50, name = "Bus 0", geodata = (-1, 0), type = "b", zone = None, in_service = True)
            bus_1 = pp.create_bus(net, vn_kv = 50, name = "Bus 1", geodata = (0, 0), type = "b", zone = None, in_service = True)
            bus_2 = pp.create_bus(net, vn_kv = 50, name = "Bus 2", geodata = (math.sqrt(2)/2, math.sqrt(2)/2), type = "b", zone = None, in_service = True)
            bus_3 = pp.create_bus(net, vn_kv = 50, name = "Bus 3", geodata = (math.sqrt(2)/2, -math.sqrt(2)/2), type = "b", zone = None, in_service = True)


            line_01 = pp.create_line(net = net, from_bus = bus_0, to_bus = bus_1, length_km = 1, std_type = "NAYY 4x50 SE")
            line_12 = pp.create_line(net = net, from_bus = bus_1, to_bus = bus_2, length_km = 1, std_type = "NAYY 4x50 SE")
            line_13 = pp.create_line(net = net, from_bus = bus_1, to_bus = bus_3, length_km = 1, std_type = "NAYY 4x50 SE")
            
            
            if connect23input == "y":
                # Voor als de lijn tussen 2 en 3 wel meedoet
                line_23 = pp.create_line(net = net, from_bus = bus_2, to_bus = bus_3, length_km = 1, std_type = "NAYY 4x50 SE")

           
        elif modus == 3:
            bus_1 = pp.create_bus(net, vn_kv = 50, name = "Bus 1", geodata = (0, 0), type = "b", zone = None, in_service = True)
            bus_2 = pp.create_bus(net, vn_kv = 50, name = "Bus 2", geodata = (1, 0), type = "b", zone = None, in_service = True)
            bus_3 = pp.create_bus(net, vn_kv = 50, name = "Bus 3", geodata = (2, 0), type = "b", zone = None, in_service = True)

            line_12 = pp.create_line(net = net, from_bus = bus_1, to_bus = bus_2, length_km = 1, std_type = "NAYY 4x50 SE")
            line_23 = pp.create_line(net = net, from_bus = bus_2, to_bus = bus_3, length_km = 1, std_type = "NAYY 4x50 SE")




        if type(q1) is float or type(q1) is int:
            load_2 = pp.create_load(net = net, bus = bus_2, p_mw = p_load_1, q_mvar = q1*p_load_1, name = "load 2")
        elif q1:    
            load_2 = pp.create_load(net = net, bus = bus_2, p_mw = p_load_1, q_mvar = q_load_1, name = "load 2")
        else:
            load_2 = pp.create_load(net = net, bus = bus_2, p_mw = p_load_1, name = "load 2")
        
        
        if type(q2) is float or type(q2) is int:
            load_3 = pp.create_load(net = net, bus = bus_3, p_mw = p_load_2, q_mvar = q2*p_load_2, name = "load 3")
        elif q2:    
            load_3 = pp.create_load(net = net, bus = bus_3, p_mw = p_load_2, q_mvar = q_load_2, name = "load 3")
        else:
            load_3 = pp.create_load(net = net, bus = bus_3, p_mw = p_load_2, name = "load 3")
        #load_3 = pp.create_load(net = net, bus = bus_3, p_mw = p_tot - p_load_1, q_mvar = q_load_2, name = "load 3")

        if modus == 2:
            gen_1 = pp.create_gen(net = net, bus = bus_0, p_mw = p_gen_1, vm_pu = 1, slack = True)
        else:
            gen_1 = pp.create_gen(net = net, bus = bus_1, p_mw = p_gen_1, vm_pu = 1, slack = True)

        # Run the network
        pp.runpp(net, max_iteration = max_iteration, tolerance_mva = 1e-10)

        iterations[i] = net["_ppc"]["iterations"]

        vm_pu_res[:, i] = net["res_bus"].loc[:, "vm_pu"]
        va_degree_res[:, i] = net["res_bus"].loc[:, "va_degree"]
        p_mw_res[:, i] = net["res_bus"].loc[:, "p_mw"]
        q_mvar_res[:, i] = net["res_bus"].loc[:, "q_mvar"]
        i += 1

    
    if nr_of_nodes == 4:
        extra_index = 1
    else:
        extra_index = 0

    plt.subplot(2, 2, 1)
    plt.plot(varlist, vm_pu_res[0, :], label = "Generator")
    if extra_index == 1:
        plt.plot(varlist, vm_pu_res[1, :], label = "Tussennode")
    plt.plot(varlist, vm_pu_res[1+extra_index, :], label = "Load 1")
    plt.plot(varlist, vm_pu_res[2+extra_index, :], label = "Load 2")
    plt.xlabel(varlabel)
    plt.ylabel('$|V_i|$ (p.u.)')
    plt.legend()

    txt = "variable %s, q1 %s, q load 1 %s (alleen als q1 True), q2 %s, q load 2 %s, p load 1 %s, p load 2 %s" % (var, q1, q_load_1, q2, q_load_2, p_load_1, p_load_2)
    plt.title(txt)

    plt.subplot(2, 2, 2)
    plt.plot(varlist, 2*math.pi/360*va_degree_res[0, :], label = "Generator")
    if extra_index == 1:
        plt.plot(varlist, 2*math.pi/360*va_degree_res[1, :], label = "Tussennode")
    plt.plot(varlist, 2*math.pi/360*va_degree_res[1+extra_index, :], label = "Load 1")
    plt.plot(varlist, 2*math.pi/360*va_degree_res[2+extra_index, :], label = "Load 2")
    plt.xlabel(varlabel)
    plt.ylabel('$\delta_i$ (rad)')
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(varlist, -p_mw_res[0, :], label = "Generator")
    if extra_index == 1:
        plt.plot(varlist, -p_mw_res[1, :], label = "Tussennode")
    plt.plot(varlist, -p_mw_res[1+extra_index, :], label = "Load 1")
    plt.plot(varlist, -p_mw_res[2+extra_index, :], label = "Load 2")
    plt.xlabel(varlabel)
    plt.ylabel('$P_i$ (MW)')
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(varlist, -q_mvar_res[0, :], label = "Generator")
    if extra_index == 1:
        plt.plot(varlist, -q_mvar_res[1, :], label = "Tussennode")
    plt.plot(varlist, -q_mvar_res[1+extra_index, :], label = "Load 1")
    plt.plot(varlist, -q_mvar_res[2+extra_index, :], label = "Load 2")
    plt.xlabel(varlabel)
    plt.ylabel('$Q_i$ (MVAr)')
    plt.legend()


    plt.show()

    plt.plot(varlist, iterations, 'o', label = "Number of iterations", color = 'black')
    plt.xlabel(var)
    plt.ylabel('Iterations')
    plt.show()

    if input("Print network? (y/n)") == 'y':
        pandapower.plotting.simple_plot(net)

    ans = input("Continue (y/n)?: ")
    if ans != "y":
        contin = False
#print(net)
#print("------ Results bus ------")
#print(net["res_bus"])
#print(net["res_bus"].at[1, "va_degree"])
#print("------ Results line ------")
#print(net["res_line"])
#print("------ Results load ------")
#print(net["res_load"])
#print("------ Results generator ------")
#print(net["res_gen"])
print(net['res_bus'])
#