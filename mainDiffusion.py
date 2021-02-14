from functionsDiffusion import *

# fixed global parameters from Doong et. al. 2017
x_s       = 20              # synthase production rate(au/min)
x_a       = 0.1
x_g       = .1
lambda_a  = 2.3             # hill coeff ahl2
K_a       = 40              # M and M constant for ahl2
D         = 0.03            # AHL diffusion coeff (#mm2/min)
D_a       = 0.3
w         = 0.75  # dx, dy
rho_n     = 3
rc        = 6 * 0.0001
Dc        = 0.0001
rho       = 5 * 0.001
lambda_n  = 2.0
K_n       = 80

def model_small(t, U_flat, shape):
    U_grid = U_flat.reshape(shape)

    # 0 LuxI
    # 1 Arabinose
    # 2 Nutrients
    # 3 Sender
    # 4 C6
    # 5 Receiver
    # 6 GFP
    
    N = hill(U_grid[2], K_n, lambda_n)
        
    LuxI_ficks      = ficks(U_grid[0], w)
    arabinose_ficks = ficks(U_grid[1], w)
    n_ficks         = ficks(U_grid[2], w)
    S_ficks         = ficks(U_grid[3], w)
    c6_ficks        = ficks(U_grid[4], w)
    R_ficks         = ficks(U_grid[5], w)

    S         = Dc * S_ficks + rc * N * U_grid[3]
    R         = Dc * R_ficks + rc * N * U_grid[5]
    LuxI      = x_s * N * hill(U_grid[1], lambda_a, K_a) * U_grid[3]
    c6        = D_a * c6_ficks + (x_a * U_grid[0]) - rho * U_grid[4]
    arabinose = D * arabinose_ficks
    n         = D * n_ficks - rho_n * N * (U_grid[3] + U_grid[5])
    gfp       = x_g * N * hill(U_grid[4], lambda_a, K_a) * U_grid[5]

    return(np.concatenate((LuxI.flatten(),
                           arabinose.flatten(),
                           n.flatten(),
                           S.flatten(),
                           c6.flatten(),
                           R.flatten(),
                           gfp.flatten()) ) )

def run_cross_setup():
    # This is Luca's orginal function simulating a single setup of the form:
    #             R
    #             R
    #             R
    #       R R R S R R R
    #             R
    #             R
    #             R
    
    # setup 7 concentrations on a grid 30 x 30
    U = np.zeros([7,30,30])
    
    shape = U.shape
    size = U.size

    U[2] = 100

    U[3][14, 14] = 0.5
    U[3][14, 15] = 0.5
    U[3][15, 14] = 0.5
    U[3][15, 15] = 0.5

    U[1][14, 14] = 5
    U[1][14, 15] = 5
    U[1][15, 14] = 5
    U[1][15, 15] = 5

    i = 0
    for j in np.arange(1, 4):
        U[5][11 - i, 14] = 0.5
        U[5][11 - i, 15] = 0.5
        U[5][12 - i, 14] = 0.5
        U[5][12 - i, 15] = 0.5
        i = i + 3
    i = 0
    for j in np.arange(1, 4):
        U[5][17 + i, 14] = 0.5
        U[5][17 + i, 15] = 0.5
        U[5][18 + i, 14] = 0.5
        U[5][18 + i, 15] = 0.5
        i = i + 3

    i = 0
    for j in np.arange(1, 4):
        U[5][14, 11 - i] = 0.5
        U[5][14, 12 - i] = 0.5
        U[5][15, 11 - i] = 0.5
        U[5][15, 12 - i] = 0.5
        i = i + 3

    i = 0
    for j in np.arange(1, 4):
        U[5][14, 17 + i] = 0.5
        U[5][14, 18 + i] = 0.5
        U[5][15, 17 + i] = 0.5
        U[5][15, 18 + i] = 0.5
        i = i + 3

    t_final = 1000  # mins
    dt = .1
    t_points = int(t_final / dt)

    t = np.arange(0, t_final, dt)
    U_init = U.flatten()
    start_time = time.time()

    sim_ivp = solve_ivp(model_small, [0, t_final], U_init,
                        t_eval=t, args=(shape,))

    sim_ivp = sim_ivp.y.reshape(7, 30, 30, t_points)

    if os.path.isdir(os.getcwd() + "/output_cross") is False:
        print("ciao")
        os.mkdir(os.getcwd() + "/output_cross")

    with PdfPages("output_cross/simulation-cross-" + time.strftime("%H%M%S") + ".pdf") as pdf:

        for i in np.arange(0, 901, 120):
            t = int(i / 60)
            f1 = multi_plots(sim_ivp[:, :, :, i], str(t) + " hours")
            pdf.savefig()
            # plt.show()
            plt.close()

    end_time = time.time()
    print(end_time - start_time)


def run_simple_setup():
    # This simulates a single sender and receiver pair
    # Also extracts time information by summing over colony locations

    # setup 7 concentrations on a grid 30 x 30
    U = np.zeros([7,30,30])
    
    shape = U.shape
    size = U.size

    # 0 LuxI
    # 1 Arabinose
    # 2 Nutrients
    # 3 Sender
    # 4 C6
    # 5 Receiver
    # 6 GFP

    # Nutrients
    U[2] = 100

    # Senders
    U[3][14, 14] = 0.5
    U[3][14, 15] = 0.5
    U[3][15, 14] = 0.5
    U[3][15, 15] = 0.5

    # Arabinose
    U[1][14, 14] = 5
    U[1][14, 15] = 5
    U[1][15, 14] = 5
    U[1][15, 15] = 5

    # Receivers
    print("Receivers at:")
    i = 0
    for j in np.arange(1, 2):
        U[5][11 - i, 14] = 0.5
        U[5][11 - i, 15] = 0.5
        U[5][12 - i, 14] = 0.5
        U[5][12 - i, 15] = 0.5
        print( "(",11 - i, 14,")",  "(",11 - i, 15,")", "(",12 - i, 14,")",  "(",12 - i, 15,")" )
        i = i + 3

    t_final = 960  # mins
    dt = .1
    t_points = int(t_final / dt)

    t = np.arange(0, t_final, dt)
    U_init = U.flatten()
    start_time = time.time()

    sim_ivp = solve_ivp(model_small, [0, t_final], U_init,
                        t_eval=t, args=(shape,))

    sim_ivp = sim_ivp.y.reshape(7, 30, 30, t_points)

    if os.path.isdir(os.getcwd() + "/output_simple") is False:
        os.mkdir(os.getcwd() + "/output_simple")


    tp = np.arange(0, 18, 2)    
        
    with PdfPages("output_simple/simulation-simple.pdf") as pdf:

        for i in np.arange(0,tp.size):        
            f1 = multi_plots(sim_ivp[:, :, :, tp[i]*60], str(tp[i]) + " hours")
            pdf.savefig()
            # plt.show()
            plt.close()
    
    with PdfPages("output_simple/timecourse-simple.pdf") as pdf:

        # calculate statistics over time
      
        # coordinates to sum over
        coords = []
        for xi in range(10,14):
            for yi in range(13,17):
                coords.append( [xi, yi] )

        print("Outputs at:")
        print(coords)
        
        # GFP
        x_gfp_t = np.zeros(tp.size)
         
        for i in np.arange(0,tp.size):
            #x_gfp_t[i] = sim_ivp[6, 11, 14, tp[i]*60]
            for jc in coords:
                x_gfp_t[i] += sim_ivp[6, jc[0], jc[1], tp[i]*60]
    
        plt.plot( tp, x_gfp_t )
        plt.xlabel("time (hours)")
        plt.ylabel("GFP")
        pdf.savefig()
        plt.close()

        # R (amount of receiver strain)
        x_r_t = np.zeros(tp.size)
       
        for i in np.arange(0,tp.size):
            for jc in coords:
                #print( "\t", sim_ivp[5, jc[0], jc[1], tp[i]*60]) 
                x_r_t[i] += sim_ivp[5, jc[0], jc[1], tp[i]*60]
    
        plt.plot( tp, x_r_t )
        plt.xlabel("time (hours)")
        plt.ylabel("Receiver")
        pdf.savefig()
        plt.close()
    
    end_time = time.time()
    print(end_time - start_time)


def run_simple(tp, ara, setup, outputdir):
    # this version of the simple setup takes in parameters so that it can be run repeatedly
    # it also contains three different sender arrangements all surrounding a single receiver

    # setup 7 concentrations on a grid 30 x 30
    U = np.zeros([7,30,30])
    
    shape = U.shape
    size = U.size

    # Nutrients
    U[2] = 50

    # Arabinose
    U[1] = ara
    
    # Senders
    if setup == 1:
        U[3][14, 14] = 0.5
        U[3][14, 15] = 0.5
        U[3][15, 14] = 0.5
        U[3][15, 15] = 0.5
        
    if setup == 2:
        U[3][14, 13] = 0.5
        U[3][14, 14] = 0.5
        U[3][15, 13] = 0.5
        U[3][15, 14] = 0.5

        U[3][14, 16] = 0.5
        U[3][14, 17] = 0.5
        U[3][15, 16] = 0.5
        U[3][15, 17] = 0.5

    if setup == 3:
        U[3][14, 11] = 0.5
        U[3][14, 12] = 0.5
        U[3][15, 11] = 0.5
        U[3][15, 12] = 0.5

        U[3][14, 14] = 0.5
        U[3][14, 15] = 0.5
        U[3][15, 14] = 0.5
        U[3][15, 15] = 0.5
        
        U[3][14, 17] = 0.5
        U[3][14, 18] = 0.5
        U[3][15, 17] = 0.5
        U[3][15, 18] = 0.5

    # Receivers
    print("Receivers at:")
    i = 0
    for j in np.arange(1, 2):
        U[5][11 - i, 14] = 0.5
        U[5][11 - i, 15] = 0.5
        U[5][12 - i, 14] = 0.5
        U[5][12 - i, 15] = 0.5
        print( "(",11 - i, 14,")",  "(",11 - i, 15,")", "(",12 - i, 14,")",  "(",12 - i, 15,")" )
        i = i + 3

    t_final = tp[len(tp)-1]*60  # mins
    dt = .1
    t_points = int(t_final / dt)

    t = np.arange(0, t_final, dt)
    U_init = U.flatten()
    start_time = time.time()

    #sim_ivp = solve_ivp(model_small, [0, t_final], U_init, method='LSODA',
    #                    t_eval=t, args=(shape,))

    sim_ivp = solve_ivp(model_small, [0, t_final], U_init,
                        t_eval=t, args=(shape,))
    
    sim_ivp = sim_ivp.y.reshape(7, 30, 30, t_points)

    with PdfPages(outputdir+"/simulation_"+str(setup)+".pdf") as pdf:    

        for i in np.arange(0,tp.size):        
            f1 = multi_plots(sim_ivp[:, :, :, tp[i]*60], str(tp[i]) + " hours")
            pdf.savefig()
            plt.close()
           
    # coordinates to sum over
    coords = []
    for xi in range(10,14):
        for yi in range(13,17):
            coords.append( [xi, yi] )

    print("Outputs at:")
    print(coords)
        
    # GFP
    x_gfp_t = np.zeros(tp.size)
         
    for i in np.arange(0,tp.size):
        for jc in coords:
            x_gfp_t[i] += sim_ivp[6, jc[0], jc[1], tp[i]*60]
    
    return x_gfp_t

def make_dist_plots():
    
    
    if os.path.isdir(os.getcwd() + "/output_dist") is False:
        os.mkdir(os.getcwd() + "/output_dist")

    tp = np.arange(0, 18, 2) 
    x1 = run_simple(tp, 5, 1, "output_dist")
    x2 = run_simple(tp, 5, 2, "output_dist")
    x3 = run_simple(tp, 5, 3, "output_dist")

    with PdfPages("output_dist/timecourse.pdf") as pdf:
        plt.plot( tp, x1 )
        plt.plot( tp, x2 )
        plt.plot( tp, x3 )
        plt.xlabel("time (hours)")
        plt.ylabel("GFP")
        pdf.savefig()
        plt.close()
    

def main():
    run_cross_setup()

    run_simple_setup()

    make_dist_plots()

    

    

main()
