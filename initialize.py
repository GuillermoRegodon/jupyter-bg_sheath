import numpy as np
import functions as f
import copy

def none_initializer():
    for i in range(100):
        yield None
        
        
        
class Constants:
    kb, k, m_e, q_e, uma, *_ = none_initializer()

constants = Constants()

constants.kb = 1.380648e-23     #Boltzmann constant, SI
constants.k = 8.61732739e-5     #k_B/e, SI
constants.m_e = 5.4858e-4       #electron mass, umas
constants.q_e = 1.602176e-19    #electron absolute charge, SI
constants.uma = 1.66054e-27     #SI
  



plasma_options = ["Argon and electrons", "Neon and electrons", "Helium and electrons", "Oxygen O+ and O-"]
        
        
        
def get_int(s):
    while True:
        try:
            x = int(input(s))
            break
        except ValueError:
            if x == "m_e":
                x = constants.m_e
            else:
                print("No valid number.  Try again...")
    return x

class Plasma:
    A, B, *_ = none_initializer()
    
    def __init__(self, kind):
        if kind in plasma_options:
            if kind == "Argon and electrons":
                self.A = 39.948        #Argon cation mass, umas
                self.B = constants.m_e           #defined above, electron mass, umas
            elif kind == "Neon and electrons":
                self.A = 20.18         #Neon cation mass, umas
                self.B = constants.m_e           #defined above, electron mass, umas
            elif kind == "Helium and electrons":
                self.A = 4.0026        #Helium cation mass, umas
                self.B = constants.m_e           #defined above, electron mass, umas
            elif kind == "Oxygen O+ and O-":
                self.A = 15.9989       #Oxygen single charged cation mass, umas
                self.B = 15.9999       #Oxygen single charged anion mass, umas
            else:
                self.A = get_int("Enter positive ion mass, in umas")
                self.B = get_int("Enter negative ion mass, in umas")


plasma = Plasma("Argon and electrons")
        
        

class Globals:
    beta, gamma, alpha0, kappa, xp, Ip, I01, I02, yp, yf, sol_Euler, xl, x0a, xMAX, *_ = none_initializer()

glob = Globals()

glob.beta = 0.1            #positive ion temperature to electron temperature ratio
glob.gamma = 10.0          #electron temperature to negative ion temperature ratio
glob.alpha0 = 0.5          #negative ion density to electron density ratio
glob.kappa = 2.0           #thermodynamic adiabatic index
#glob.vmas0 = 0.5           #Bohm velocity
glob.xp = 1.0              #probe radius to Debye length ratio

glob.Ip = 5.0              #adimensional ion current
glob.I01 = None            #negative ion current constant
glob.I02 = None            #electron current constant
glob.yp = None             #adimensional probe potential
glob.yf = None             #adimensional floating potential

glob.sol_Euler = 1         #index to select Euler's equation branch for the solution

glob.xl = None             #maximum x for singularity, according to theory
glob.x0a = None            #singularity x for bounded quasineutral solution
glob.xMAX = 50.0           #Runge-Kutta integration limit if singularity x is small



class Configuration:
    cylindrical, numprec, graph, save, warn, time_calculation_reduction, *_ = none_initializer()

configuration = Configuration()

configuration.cylindrical = 1
configuration.numprec = 1e-9
configuration.graph = 1
configuration.save = 0
configuration.warn = 0
configuration.time_calculation_reduction = 5     #reduces compute time, 3 is ok
configuration.xinf = 20
configuration.save_approx = False
#Ap = np.pi*6e-3*2e-4 ??



class BetaGammaSheath:
    c = copy.deepcopy(constants)
    p = copy.deepcopy(plasma)
    g = copy.deepcopy(glob)
    config = copy.deepcopy(configuration)
    
    
    
    
    
