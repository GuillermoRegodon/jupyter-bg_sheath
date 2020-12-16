import numpy as np
from scipy import optimize



def defineI0tuple(bg=[]):
    '''
    defineI0tuple(bg): Defines the values of I01, and I02, the constants
    before the effusion negative particles currents collected by the
    probe.'''
    I01 = bg.g.I01
    I02 = bg.g.I02
    A = bg.p.A
    B = bg.p.B
    m_e = bg.c.m_e
    alpha0 = bg.g.alpha0
    gamma = bg.g.gamma
    
    I02 = np.sqrt(A/(4*np.pi*m_e))
    I01 = I02*alpha0*np.sqrt(m_e/(B*gamma))
    
    bg.g.I01, bg.g.I02 = I01, I02

    
    
def fN(y, bg=[]):
    '''
    fN(y): Auxiliary function for getydoCI, the electric field according
    to the solution to the sheath potential according to the quasineutral
    assumption which is only valid in the presheath, far from the probe.
    
    It returns the negative charge density in a region with potential y.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    
    return np.e**(-y) + alpha0*np.e**(-gamma*y)



def fNp(y, bg=[]):
    '''
    fNp(y): Auxiliary function for getydoCI, the electric field according
    to the solution to the sheath potential according to the quasineutral
    assumption which is only valid in the presheath, far from the probe.
    
    It returns the derivative of the negative charge density in a region
    with potential y with respect to the potential.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    
    return -np.e**(-y) -gamma*alpha0*np.e**(-gamma*y)



def getxl(bg=[]):
    '''
    getxl(bg): When beta is not zero, returns the theoretical maximum value
    of the singularity x.
    
    When beta equals zero, returns a big singularity x for which the potential
    is small, y0.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    xp = bg.g.xp
    xl = bg.g.xl
    kappa = bg.g.kappa
    
    if beta != 0:
        if kappa==2:
            xl = ((1.5)**(1.5))*Ip/(np.sqrt(beta)*(1+alpha0))
        elif kappa==3:
            xl = np.sqrt(8/(3*beta))*Ip/(1+alpha0)
          
        y0asin = 1.1*getylimit(bg)
        x0 = getxCI(y0asin, bg)
        
        y0 = y0asin
        ydot0 = getydotCI(x0, y0, bg)
        
    else: #beta == 0:
        y0 = Ip*0.0003/(1+alpha0)
        N0 = np.e**(-y0) + alpha0*np.e**(-gamma*y0)
        N0 = fN(y0, bg)
        x0 = Ip/(np.sqrt(y0)*N0)
        ydot0 = 1/(x0*((np.e**(-y0) + alpha0*gamma*np.e**(-gamma*y0))/N0 - 0.5/y0))
        ydot0 = 1/(x0*((-fNp(y0, bg))/N0 - 0.5/y0))
        xl = x0
    
    bg.g.xl = xl
    
    return xl

        

def getInitFloat(x0a, bg=[]):
    '''
    getInitFloat(bg, x0a): When beta is not zero, returns the singularity x,
    y and ydot, as well as it calculates xl and saves to globals, maybe not necessary.
    
    When beta equals zero, returns a big singularity x.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    xp = bg.g.xp
    xl = bg.g.xl
    kappa = bg.g.kappa
    
    if beta != 0:
        if kappa==2:
            x0 = x0a
            if x0a == None:
                print("Warn")
            y0 = 3*(beta*Ip/(x0*(1+alpha0)))**(2/3) - 2*beta
            xl = ((1.5)**(1.5))*Ip/(np.sqrt(beta)*(1+alpha0))
            ydot0 = -2*(y0 + 2*beta)/(3*x0)
        elif kappa==3:
            x0 = x0a;
            y0 = np.sqrt(6*beta)*(Ip/(x0*(1+alpha0))) - 1.5*beta
            xl = np.sqrt(8/(3*beta))*Ip/(1+alpha0)
            ydot0 = -np.sqrt(6*beta)*Ip/((x0**2)*(1+alpha0))
        N0 = NnormEuler(x0, y0, bg)

    else: # beta == 0:
        y0 = Ip*0.0003
        #N0 = e^(-y0) + alpha0*e^(-gamma*y0)
        N0 = fN(y0, bg)
        x0 = Ip/(np.sqrt(y0)*N0)
        ydot0 = 1/(x0*((np.e**(-y0) + alpha0*gamma*np.e**(-gamma*y0))/N0 - 0.5/y0))
        xl = x0
        if x0>40:
            print('x0 very high, it may take time, x0 =', x0)
    
    bg.g.xl = xl
    
    return [x0, y0, ydot0, xl]



def getxCI(y, bg=[]):
    '''
    getxCI(bg, y): Returns the x for a given y according to the quasineutral
    solution of the model.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    kappa = bg.g.kappa

    c = fN(y, bg)
    
    if (kappa == 2):
        x = Ip/(c*np.sqrt(y-2*beta*(c/(1+alpha0)-1)))
    elif (kappa == 3):
        x = Ip/(c*np.sqrt(y-1.5*beta*(c^2/(1+alpha0)**2 - 1)))
    else:
        x = None
    
    return x
    


def getydotCI(x0, y0, bg=[]):
    '''
    getydotCI(bg, y): Returns the ydot, the electric fiels for a given x and y
    according to the quasineutral solution of the model.'''
    Ip = bg.g.Ip
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    kappa = bg.g.kappa
  
    if kappa == 2:
        a = fNp(y0, bg)/fN(y0, bg)
        b = (1-2*beta*fNp(y0, bg)/(1+alpha0))/(y0 - 2*beta*(fN(y0, bg)/(1+alpha0)-1))
        ydot0_1 = -x0*(a + 0.5*b)
        ydot0 = 1/ydot0_1

    elif kappa == 3:
        a = fNp(y0, bg)/fN(y0, bg)
        bnum = (1 - 1.5*beta*(2*fN(y0, bg)/(1+alpha0)**2)*fNp(y0, bg))
        bden = (y0 - 1.5*beta*((fN(y0, bg))**2/(1+alpha0)**2 - 1))
        ydot0_1 = -x0*(a + 0.5*bnum/bden)
        ydot0 = 1/ydot0_1
    
    return ydot0



def getylimit(bg=[]):
    '''
    getylimit(bg=[]): Returns the maximum value of the potential according
    to the quasineutral solution, before it returns back to the plasma.
    The quasineutral solution becomes invalid for lower values of y. This
    limit is useful in plots and in estimations of quasineutral validity'''
    kappa = bg.g.kappa
    
#    if (kappa == 2):
#        retval = optimize.root(limitk2, 0, args=bg)
#    elif (kappa == 3):
#        retval = optimize.root(limitk2, 0, args=bg)

    retval = optimize.root(limitkappa, 0, args=bg)

    return retval.x



def limitkappa(y, bg=[]):
    '''
    limitkappa(y, bg=[]): derivative of the position with respect to the
    potential in the quasineutral solution, to find the turning point.'''
    gamma = bg.g.gamma
    alpha0 = bg.g.alpha0
    beta = bg.g.beta
    kappa = bg.g.kappa
    
    if kappa == 2:
        if ((beta!=0) or (y!=0)):
            f_y = (np.e**y)/(1 + alpha0*np.e**(-y*(gamma-1))) - 3*beta/((1+alpha0)*(y+2*beta))
        else:
            f_y = 0

    elif kappa == 3:
        if ((beta!=0) or(y!=0)):
            f_y = (np.e**(2*y)*(1+alpha0)**2)/(( 1 + alpha0*np.e**(-y*(gamma-1)) )**2) - 3*beta/(y + 1.5*beta);
        else:
            f_y = 0

    return f_y



def NnormEuler(x, y, bg=[]):
    '''
    NnormEuler(x, y, bg=[]): Solves momentum conservation equation in the
    model. It is a cubic or a quartic, depending on model parameter kappa.'''
    Ip = bg.g.Ip
    beta = bg.g.beta
    alpha0 = bg.g.alpha0
    kappa = bg.g.kappa
    sol_Euler = bg.g.sol_Euler
    __config__cylindrical = bg.config.cylindrical
  
    a = (kappa/(kappa-1))*beta/(1+alpha0)**(kappa-1)
    b = -y - (kappa/(kappa-1))*beta
    c = 0
    if __config__cylindrical == 1:
        d = Ip**2/x**2
    else:
        d = Ip**2

    if (beta != 0):

        if kappa==2:
            # trigonometric solution for cubic
            p = - b**2/(3*a**2)
            q = (2*b**3 + 27*a**2*d)/(27*a**3)

            m=3*q*np.sqrt(-3/p)/(2*p)
            if abs(m)>1:
                m = 1*np.sign(m)

            rm = np.real((1/3)*np.arccos(m))
            im = np.imag((1/3)*np.arccos(m))

            if sol_Euler == -1:
                n = rm- 2*np.pi/3 #Smaller positive solution
            elif sol_Euler == +1:
                n = rm #Higher positive solution

            t0 = 2*np.sqrt(-p/3)*np.cos(n)

            result = t0 - b/(3*a)

        elif kappa == 3:
            if sol_Euler == -1:
                result = np.real(np.sqrt((-b - np.sqrt(b**2-4*a*d))/(2*a))) #Smaller positive solution
            elif sol_Euler == +1:
                result = np.real(np.sqrt((-b + np.sqrt(b**2-4*a*d))/(2*a))) #Higher positive solution

    else: #return the ABR quasineutral solution
        if y<0:
            print("Warn: ", y, x)
        result = Ip/(x*np.sqrt(y))
        
    return result



def poisson_cyl(y, x, bg=[]):
    '''
    poisson_cyl(y, x, bg=[]): Returns the derivatives of the potential and 
    the electric field according to Poisson equation as stated in the model'''
    beta = bg.g.beta
    gama = bg.g.gamma
    alpha0 = bg.g.alpha0
    Ip = bg.g.Ip
    
    ydot = [0, 0]
    ydot[0] = y[1]
  
    C = NnormEuler(x, y[0], bg) - fN(y[0], bg)
    ydot[1] = -y[1]/x + C
    
    return ydot
    
    
    
def runge_kutta(f, y_0, tlim, h, bg=[]):
    '''
    runge_kutta(f, y_0, tlim, h, bg=[]): Runge kutta with a special exit
    condition for this model. It exits if potential crosses the y = 0 axis
    or if the electric field becames positive.'''
    y = []
    y.append(y_0)
    
    if (tlim[1]-tlim[0])*h > 0:
        t = []
        i = 0
        while True:
            t.append(tlim[0] + h*i)
            if (t[-1]-tlim[1])*(tlim[0]-tlim[1]) <= 0:
                break
            i += 1
    else:
        raise RuntimeError from exc

    for i in range(1, len(t)):
        k1 = f(y[i-1], t[i-1], bg)
        
        a = [h/2*k for k in k1]
        b = [y[i-1][j] + a[j] for j in range(len(a))]
        k2 = f(b, t[i-1] + h/2, bg)
        
        a = [h/2*k for k in k2]
        b = [y[i-1][j] + a[j] for j in range(len(a))]
        k3 = f(b, t[i-1] + h/2, bg)
        
        a = [h*k for k in k3]
        b = [y[i-1][j] + a[j] for j in range(len(a))]
        k4 = f(b, t[i-1] + h, bg)

        a = [h/6*k for k in k1]
        b = [h/3*k for k in k2]
        c = [h/3*k for k in k3]
        d = [h/6*k for k in k4]
        e = [y[i-1][j] + a[j] + b[j] + c[j] + d[j] for j in range(len(a))]
        y.append(e)
        #warning('vuelta en bucle dentro de RK') %optional loop warning
        
        if h > 0:
            # exit condition
            if (i<len(t) and (t[i] <= 0)) or (y[-1][0]*y[-1][1] > 0):
                break
                
        else:
            beta = bg.g.beta
            Ip = bg.g.Ip
            alpha0 = bg.g.alpha0
            kappa = bg.g.kappa
            if beta!=0:
                N = (2*Ip**2*((1+alpha0)**(kappa-1))/(kappa*beta*t[i]**2))**(1/(kappa+1));
                if (y[-1][0] < (Ip/(t[i]*N))**2 + (kappa/(kappa-1))*beta*((N/(1+alpha0))**(kappa-1)-1)):
                    y = y[1:]
                    break

    t = t[0:len(y)]
    
    return [y, t]

