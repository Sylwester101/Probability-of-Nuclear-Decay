import numpy as np
import pandas as pd 
import random
from mpmath import mp
import math
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt


#https://www.davidhbailey.com//dhbpapers/dhb-tanh-sinh.pdf
def numericalIntegrate(func,min,max,d):
    #change in interval
    multiple= (max-min)/2
    offset = (max+min)/2
    #int_-1^1 f(multiple*x+offset)*multiple dx
    
    total=0
    
   
    placeholder=0
    for j in range(0,20):
        h=2**(-j)
        value_inside = 1/np.pi * np.log(2 * 10**d * np.min([1, multiple]))

        if value_inside < 1:
            raise ValueError(f"Invalid value for arccosh: {value_inside} < 1")
        N = math.floor(1/h * np.arccosh(value_inside))

        for k in range(-N-1,N+1):
            
            if j==0 or (k%2==1):
                tk=k*h
                
                xk=np.tanh(np.pi/2*np.sinh(tk))
                
                wk=(1/2)*np.pi*np.cosh(tk)/(np.cosh(np.pi/2*np.sinh(tk))**2)
                
            
                if(k>0):
                    total=total+func(multiple*xk+offset)*wk
                else:
                    total=total+func(multiple*xk+offset)*wk
        if abs(h*total-placeholder)<10**(-d/2)*placeholder:
            break     
        
        
        placeholder=h*total
        
    return h*multiple*total

def neutronwidthfunctionunintegrated(x,an,xn):
    return np.exp(2*np.sqrt(an*x))*(xn-x)
def fissionunintegrated(x,af):
    return np.exp(2*np.sqrt(af*x))
def gammaemission(T,A):
    return 0.624*10.**(-9.)*A**(1.6)*T**(5.)

def neutronCDF(k,T,xn):
    return (T-np.exp(-k/T)*(T+k))/(T-np.exp(-xn/T)*(T+xn))


def betterneutroned(k): 
    func = lambda k: np.exp(2*np.sqrt(an*xn-an*k))*k
    return np.exp(2*np.sqrt(an*xn-an*k))*k/numericalIntegrate(func,0,xn,10)


#https://notebook.community/tommyogden/quantum-python-lectures/11_Monte-Carlo-Maxwell-Boltzmann-Distributions
def cdfinversion(func,max,stepsize):
    x=np.arange(0,max+stepsize,stepsize)
    
    cdf=func(x)
    inv_cdf=interp(cdf,x)
    
    return inv_cdf

#The monte carlo function proper
def montecarlofission(Z,weight,m):
    frequencies={key:0 for key in actinide_data[Z]}
    #the for loop, loops through all the samples that the user specifies, can be further optimized through parallel processing since the processes are independent of each other, however some variables I think need to be even more local than they are at present
    for i in range(m):
        A=weight
        energy=excitedenergy
        #The while loop serves to cycle through the various processes that an individual sample can undergo
        while (1>0):
            
            
            Bf,Bn=actinide_data[Z][A]
            xf=energy-Bf
            xn=energy-Bn
            an=A/8.5
            af=an   
            T=np.sqrt(energy/an)
            
            #Ensures that fission is possible, if not the loop is stopped
            if(xf<0):  
                break
            def rho0pi(x): return 1/(np.pi*2*np.exp(2.*np.sqrt(an*x)))
            def fu(x): return fissionunintegrated(x,af)
            def nu(x): return neutronwidthfunctionunintegrated(x,an,xn)
            integratedfission=numericalIntegrate(fu,0,xf,6)
            #Ensures that neutron emission is possible, if not the neutron decay width is set to 0
            if (xn<0):
                integratedneutron=0
            else:
                integratedneutron=numericalIntegrate(nu,0,xn,6)

            #calculation of the probabilities, neutron and fission
            probfission= rho0pi(energy)*integratedfission/(rho0pi(energy)*integratedfission+(A**(2/3)/5.286)*rho0pi(energy)*integratedneutron+gammaemission(T,A))
            probneutron= (A**(2/3)/5.286)*rho0pi(energy)*integratedneutron/((rho0pi(energy)*integratedfission+(A**(2/3)/5.286)*rho0pi(energy)*integratedneutron+gammaemission(T,A)))
            


            #random number generator to determine which process the atom undergoes 
            fatenumber=random.random()
            #fission
            if (fatenumber<probfission):
                frequencies[A]=frequencies[A]+1
                break
            #neutron emission
            elif(fatenumber<probneutron+probfission):
                rand_num=1
                while rand_num>=neutronCDF(xn,T,xn):
                    rand_num=random.random()
                A=A-1
                def nemission(x): return neutronCDF(x,T,xn)
                
                energy=energy-Bn-cdfinversion(nemission,xn,0.1)(rand_num)
                continue
            #gamma radiation
            else:
                energy=energy-2
                
            
    #returns the map of how many particles underwent fission at a given atomic weight        
    return frequencies


#Map of binding energies for Uranium
actinide_data = {
    92: {240:(6.59,5.64), 239:(7.05,4.99),238:(6.27,5.96), 237:(6.49,5.31),236:(6.13,6.36),235:(6.14,5.48),234:(6.16,6.49),233:(6.23,5.97),232:(5.95,6.84),231:(5.84,10)}

}

excitedenergy=60

samples=1
frequencies=montecarlofission(92,240,samples)
probabilities={key:frequencies[key]/samples for key in frequencies}
print(probabilities)
'''rand_nums = np.random.random(1000000)
T=np.sqrt(60/(240/8.5))
xn=60-5.64
def nem(x): return neutronCDF(x,T,xn)
speeds = cdfinversion(nem,60-5.64,0.1)(1)
print(np.average(speeds))'''







def probdistribution(x,Z,A):
    Bf,Bn=actinide_data[Z][A]
    
    an=A/8.5
    af=an   
    Bf=5.71
    Bn=6.2
    
    results =np.zeros_like(x)
    for i, energy in enumerate(x):
        integratedfission=0
        integratedneutron=0
        xf=energy-Bf
        xn=energy-Bn
        T=np.sqrt(energy/an)
        def rho0pi(x): return 1/(np.pi*2*np.exp(2.*np.sqrt(an*x)))
        def fu(x): return fissionunintegrated(x,af)
        def nu(x): return neutronwidthfunctionunintegrated(x,an,xn)
        if xf>=0:
            integratedfission=numericalIntegrate(fu,0,xf,6)
        if xn>=0:
            integratedneutron=numericalIntegrate(nu,0,xn,6)
        results[i]=rho0pi(energy)*integratedfission/(rho0pi(energy)*integratedfission+(A**(2/3)/5.286)*rho0pi(energy)*integratedneutron+gammaemission(T,A))
    return results



x = np.linspace(0, 20, 500)
plt.plot(x, probdistribution(x,92,238))

plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Function with Zero Below Threshold')

plt.grid(True)
plt.show()