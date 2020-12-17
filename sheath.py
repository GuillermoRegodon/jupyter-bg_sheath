import functions as f
import numpy as np


def getNrk(beta, Ip):
    '''
    These values of Nrk have been found after checking the convergence
    of the model and the duration of the calculations for many Ip values.'''
    Nrk = 1000;
    if beta < 0.1:
        if Ip<=5:
            Nrk = 1000
        elif Ip<=14:
            Nrk = 4000
        elif Ip<=25:
            Nrk = 8000
        #else:
        #    Nrk = 12000

    else:
        if Ip <=60:
            Nrk = 250
        elif Ip<100:
            Nrk = 1000
        else:
            Nrk = 2000
    
    return Nrk



def cylindrical(bg):
    xl = bg.g.xl
    kappa = bg.g.kappa      # must be 2 or 3 for these expresions
    sol_Euler = bg.g.sol_Euler
    beta = bg.g.beta
    alpha0 = bg.g.alpha0
    gamma = bg.g.gamma
    Ip = bg.g.Ip
    x0a = bg.g.x0a
    xMAX = bg.g.xMAX
    
    numprec = bg.config.numprec
    graph = bg.config.graph
    save = bg.config.save
    warn = bg.config.warn
    time_calculation_reduction = bg.config.time_calculation_reduction
    xinf = bg.config.xinf
  
    __xMAX = xMAX + 20; # local of xMAX
    frac_tail = 0.5;
  
    xd = 0;
    xu = 0;
    xd_prev = 0;
    xu_prev = 0;
  
    __xp = bg.g.xp; #local of xp

# to reduce calculation time if such precission is unneeded
    Nrk = getNrk(beta, Ip)/bg.config.time_calculation_reduction

    if bg.config.warn:
        print("Nrk = {0}; Ip = {1}; beta = {2}".format(Nrk, Ip, beta))

# for beta not equals 0, the singularity must be found
    if beta != 0:
        bg.g.sol_Euler = +1
      
        #Looks for the singularity
        xl = f.getxl(bg)
        xd = 0
        xu = xl
        
        while True:

            xtest = 0.5*(xu+xd)     # Calculate xtest for the iteration

            [x0, y0, ydot0, xl] = f.getInitFloat(xtest, bg)
            if not(bg.config.graph):
                __xp = x0

            rkincrement = min([1, xl])/Nrk     # key line for calculation time
            
# bg_runge_kutta stops if non validity of the solution, valid only for this problem
            [y, x] = f.runge_kutta(f.poisson_cyl, [y0, ydot0], [x0, xinf+x0], rkincrement, bg)

            if y[-1][0] > 0:
                xd_prev = xd
                xu_prev = None
                xd = xtest
                if xd_prev == None:
                    print("663192")
                
            else:
                xd_prev = None
                xu_prev = xu
                xu = xtest
                if xu_prev == None:
                    print("663193")
          
            if bg.config.warn:
                print("bg_sheath looking for singularity loop, xu - xd = {0}; xtest = {1}".format(xu-xd, xtest))
            
            if (xtest<0.001):
                print("Ip too low for sheath to form, Ip = {0}".format(Ip))
                xu = xd
                return None
            
            if abs(xu-xd) < max([bg.config.numprec*0.001, 1e-12]):
                break
                
        x0a = xtest;
        
        if bg.config.graph:
            print("Singularity at x0a = ", x0a)
        
        if bg.config.graph or bg.config.save:
            __xp = max(__xMAX,4*x0a)
        else:
            __xp = xp

# if xp is to the right, to the plasma, more values have to be calculated
        if len(y) > 1:

            if xd_prev == None:
                xu = xu_prev
            elif xu_prev == None:
                xd = xd_prev
            
# We have to calculate where the solutions u and d are too distint
            print("Calculate solutions diverging down")
            [x0u, y0u, ydot0u, xl] = f.getInitFloat(xu, bg)
            [y_d, x_d] = f.runge_kutta(f.poisson_cyl, [y0u, ydot0u], [x0u, xinf+x0u], rkincrement, bg)
        
            print("Calculate solutions diverging up")
            [x0d, y0d, ydot0d, xl] = f.getInitFloat(xd, bg)
            [y_u, x_u] = f.runge_kutta(f.poisson_cyl, [y0d, ydot0d], [x0d, xinf+x0d], rkincrement, bg)
            
            l_y_dif = min(len(y_u),len(y_d))
            y_dif = np.subtract(y_u[1:l_y_dif], y_d[1:l_y_dif])
            for i in range(l_y_dif):
                if y_dif[i][0] > bg.config.numprec*1000:     # Here, if precission is too low, convergence fails
                    l_y_dif = i # We take the value out of the loop
                    break
            i = l_y_dif

# We take the point were we trust the solution and calculate again from there
            while True:

                # Position and potential are trusted
                x0w = x[i] 
                y0w = y[i][0]

                # Electric field numeric search
                ydotd = y[i][1]-0.01
                ydotu = y[i][1]+0.01

                x = x[1:(i-1)]
                y = y[1:(i-1)]

                while True:
                    ydotwin = 0.5*(ydotu+ydotd)
                    [ywin, xwin] = f.runge_kutta(f.poisson_cyl, [y0w, ydotwin], [x0w, xinf+x0w], rkincrement, bg)

                    if ywin[-1][0] < 0:
                        ydotd = ydotwin
                    else:
                        ydotu = ydotwin
                
                    if abs(ydotu-ydotd) < bg.config.numprec:
                        break

                if bg.config.warn:
                    print('quasineutral search 0, xwin = {0} , i = {1}'.format(xwin[0], i))

#Calculate the index where the solutions u and d are too distint
                [y_u, x_u] = f.runge_kutta(f.poisson_cyl, [y0w, ydotu], [x0w, xinf+x0w], rkincrement, bg)
                [y_d, x_d] = f.runge_kutta(f.poisson_cyl, [y0w, ydotd], [x0w, xinf+x0w], rkincrement, bg)
        
                l_y_dif = min(len(y_u),len(y_d))
                y_dif = np.subtract(y_u[1:l_y_dif], y_d[1:l_y_dif])
                for j in range(1, len(y_dif)):
                    if y_dif[j][0] > bg.config.numprec*1000:     # Here, if precission is too low, convergence fails
                        l_y_dif = j
                        break

                i = i + l_y_dif

                x = x + xwin     # Concatenate
                y = y + ywin     # Concatenate

                if x[-1] > __xp:
                    break
            
        try:
            i = x.index(__xp)
        except ValueError:
            pass
        else:
            x = x[0:i]
            y = y[0:i]

        N = [f.NnormEuler(x[i], y[i][0], bg) for i in range(len(x))]
            
# Calculates towards the axis
        _xp = 0
      
        if (__xp < x0+rkincrement): # Calculates to the left if necessary
            bg.g.sol_Euler = -1
        
            rkincrement = min(1, xl)/Nrk     # this line is key to control de duration of the calculations
        
            if rkincrement < x0-__xp:     # If there is something to calculate. Note redinition of rkincrement
                x_up = __xp-rkincrement     # Set x_up below __xp
                if x_up < 0:                     # If negative,
                    x_up = bg.config.numprec     # Set a minimum nonzero value

                [y2, x2] = f.runge_kutta(f.poisson_cyl, [y0, ydot0], [x0, x_up], -rkincrement, bg)
                # This one always converges, no need to check anything

                N2 = [f.NnormEuler(x2[i], y2[i], bg) for i in range(len(x2))]

            else:     # __xp is to the left but very close, no need to calculate
                y2 = []
                x2 = []

            if length(x2)!=0:
                
                x_array = [x2(i) for i in reversed(range(len(x2)))] + x
                
                yz = [y2(i) for i in reversed(range(len(y2)))] + y
                y_array = [yz[i][0] for i in range(len(yz))]
                z_array = [yz[i][1] for i in range(len(yz))]

                N_array = [N2(i) for i in reversed(range(len(N2)))] + N

            else:
                x_array = x
                y_array = [y[i][0] for i in range(len(y))]
                z_array = [y[i][1] for i in range(len(y))]
                N_array = N
      
        else:
            x_array = x
            y_array = [y[i][0] for i in range(len(y))]
            z_array = [y[i][1] for i in range(len(y))]
            N_array = N

# To calculate the sheath for beta equals 0, we follow the ABR squeme
    else: #beta==0
      
        if bg.config.graph:
            __xp = 0
        else:
            __xp = xp
   
        beta_0 = 4
        [x0, y0, ydot0, xl] = f.getInitFloat(1, bg)     # beta == 0 in bg

# if Ip is high, the calculation is long, we try to accelerate without losing precission
        rkincrement = 0.005*(max(200, x0))/Nrk;
      
        [y2, x2] = f.runge_kutta(f.poisson_cyl, [y0, ydot0], [x0, __xp+rkincrement], -rkincrement, bg);
      
        N2 = [f.NnormEuler(x2[i], y2[i][0], bg) for i in range(len(x2))]
      
        x_array = [x2[i] for i in reversed(range(len(x2))) if x2[i] > 0]
        y_array = [y2[i][0] for i in reversed(range(len(x2))) if x2[i] > 0]
        z_array = [y2[i][1] for i in reversed(range(len(x2))) if x2[i] > 0]
        N_array = [N2[i] for i in reversed(range(len(x2))) if x2[i] > 0]
        
    return [x_array, y_array, z_array, N_array]
    













