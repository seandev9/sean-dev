#Python code for Chemical concentration - and approximation using Euler method
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

#starting conditions
CK = 1.0 #coef
CQ = 2.5 #Rate 
CV = 100.0 #Volume 
CC0 = 40.0 # initial concentration

#need to pick run time to capture the exponential decay (look at a few graphs)
runtime = 40 #seconds

#calculate baseline concentration and then Euler approx

#baseline concentration calculation
numpts = 125 # number of points to approximate the curve
temp = 0.00 #temporary variable to check exp calculation (not totally necessary)
CT = np.zeros(numpts) #matrix of inverse concentrations
time = np.zeros(numpts)
#print (CT,time)
f = open("ChemOutput.txt","w+")
f.write(" Volume (L)=" + str(CV) + "\n")
f.write(" Flow Rate (L/s)=" + str(CQ) + "\n")
f.write(" Initial concentration (g/L)=" + str(CC0) + "\n")
print ("Time    Concentration")
#Calculate Concentration vs time using Bernouli Eqn solution
for m in range (0,numpts):
    time[m]= (m)* runtime/numpts
    expvar1=(CQ*time[m])/CV
    if expvar1 > 700: expvar1=700
    temp = math.exp(expvar1)
    CT[m] = (CK*CV)/CQ+(1/CC0+(CK*CV)/CQ)*(temp)
    output1=("{:.2f}".format(time[m]) + "      " + "{:.8f}".format(1/CT[m])+"\n")
    print(output1)
    f.write(output1)
    
#plot of baseline Chemical reaction rate
plt.style.use('dark_background')    
plt.figure(figsize=(12,6)) # 10 is width, 6 is height
plt.plot(time,1/CT,'white', label='Concentration ')  # blue dots
plt.plot(time,1/CT, 'white') #overplot black lines
plt.title('Concentration vs. Time')  
plt.xlabel('Time (sec)')
plt.ylabel('Concentration (g/L)')
plt.text(0.76*max(time), 0.77*max(1/CT)," Volume (L)=" + str(CV),color='red')
plt.text(0.76*max(time), 0.81*max(1/CT)," Flow Rate (L/s)=" + str(CQ),color='red')
plt.text(0.76*max(time), 0.85*max(1/CT), " Initial concentration (g/L)=" + str(CC0),color='red')
plt.xlim(0, max(time))
plt.ylim(0, max(1/CT))
plt.grid(True)

plt.show()

#Compute Eulers approximation of curve
numit = 9 #number of interations of step sizes
hstep = runtime/(2*numit+1) #smallest step size in seconds
X0 = 0 #initial x value
Y0 = 1/CT[X0] #initial y value
hstepit = np.zeros(numit) #hstep iterative values
DIM = int(numpts/hstep) #number of data points = dim of arry
DCT = np.zeros((numit,DIM)) #derivative of CT
ACT = np.zeros((numit,DIM)) #approximate CT value
DTime = np.zeros((numit,DIM)) #time scale for euler function
Dtemp = np.zeros(DIM)
f.write('Euler Section' + "\n")
f.write('Iterations = ' + str(numit) + "\n")
f.write("Volume (L)=" + str(CV) + " Flow Rate (L/s)=" + str(CQ) + "Concentration (g/L)=" + str(CC0))
print('Euler section')
print("Volume (L)=" + str(CV) + " Flow Rate (L/s)=" + str(CQ) + "Concentration (g/L)=" + str(CC0))
print('Iterations = ' + str(numit))
for i in range(0, numit):
    hstepit[i]= hstep*(numit-i)  #run through step sizes
    print('H Step size =',"{:.3f}".format(hstepit[i]))
    print('Time  Concentration  Euler Approx')
    f.write('Time   Concentration   Euler Approx' + "\n")
    for n in range (X0, DIM):
        expvar =  (CQ*n*hstepit[i])/CV
        if expvar > 700: expvar=700
        Dtemp[n] = math.exp(expvar) #temporary variable EXP function
        DCT[i][n]=(CQ/(CC0*CV) + (CK))*Dtemp[n] #compute derivative or slope of CT
        DTime[i][n] = n*hstepit[i]
        if n == 0 :
            ACT[i][n] = 1/Y0 #initial value 
        else:
            ACT[i][n] = ACT[i][n-1] + (DCT[i][n-1]*hstepit[i]) #Euler approximation = prior value times slope and delta
        if (DTime[i][n] < runtime) :
            output2=str("{:.2f}".format(DTime[i][n]) + '   ' + "{:.6f}".format(1/CT[int(DTime[i][n]/(runtime/numpts))]) + '        ' + "{:.6f}".format(1/ACT[i][n]) + "\n")
            print(output2)
            f.write(output2)

#plot of derviative vs function
plt.style.use('dark_background')            
plt.figure(figsize=(12,6)) # 10 is width, 6 is height
plt.plot(time,1/CT,'white', linewidth= 3.0,linestyle = ':')  
plt.title('Concentration and Euler Approx vs. Time')  
plt.xlabel('Time (sec)')
plt.ylabel('Concentration (g/L)')
plt.text(0.72*max(time), 0.80*max(1/CT)," Volume (L)=" + str(CV),color='red')
plt.text(0.72*max(time), 0.84*max(1/CT)," Flow Rate (L/s)=" + str(CQ),color='red')
plt.text(0.72*max(time), 0.88*max(1/CT), " Initial concentration (g/L)=" + str(CC0),color='red')
plt.xlim(0, max(time))
plt.ylim(0, max(1/CT))
plt.grid(True)

colors=cm.rainbow(np.linspace(0,1,numit))
for j in range(0,numit):
    #plt.plot(DTime[j],1/ACT[j],'*')#color=colors[j])
    plt.plot(DTime[j],1/ACT[j],color=colors[j],linewidth=(0.8+(2+numit)/(j+2)),alpha=(0.30 + 0.69*(j/numit))) 

plt.show()

f.close   
