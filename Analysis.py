import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size
from scipy.ndimage.measurements import standard_deviation
from scipy.sparse.extract import find
import tifffile as tf
from scipy.ndimage import rotate
from scipy.special import erf
from scipy.optimize import curve_fit

#Global variables ----------------------------
pixlength = 500*10**(-6)/215    #meters per pixel
channelStart = 920              #pixels, lower is further out in channel
plt.style.use("seaborn")        #Stylish
# --------------------------------------------

def generateYvalData(filename = "110ul_4xBF_startCanal.tif"):
    #  Crop and rotate image data array
    rotPic = rotate(tf.imread(filename),1.6)
    pictif = rotPic[130:315,30:800]          

    #  generates array to put relative luminance-values
    Yvals = np.zeros((len(pictif[:,1]),len(pictif[1,:]))) 
    #  Creates X, Y arrays for looping
    Ys,Xs = np.arange(len(pictif[:,1])), np.arange(len(pictif[1,:]))

    for j in Ys:
        for i in Xs:
            Yvals[j,i] = linearizeLuminocity(pictif[j,i]) #generates array of values for relative intensity
    return Yvals

def findDs(interval = (50,-50), startx = 220, filename = "110ul_4xBF_startCanal.tif" ):
    Yvals = generateYvalData(filename)[interval[0]:interval[1],startx:] # generates intervals of interest
    Dlist = []
    Yaxis = np.arange(len(Yvals[:,1]))*pixlength
    Xaxis = np.arange(len(Yvals[1,:]))
    for i in Xaxis:
        #  finds lowest and highest intensity that minimizes squared distance to data
        #  in the interval of interest, 30 outer points where concentration should be unaffectd by diffusion
        lowInt = findAv(Yvals[:30,i])
        highInt = findAv(Yvals[-30:,i]) 
        #  calculates C(x = X,y)/C_max according to Beer-Lamberts' law.
        cRelative = np.log10(highInt/Yvals[:,i])/np.log10(highInt/lowInt) 
        #  curve fits data with parameters D*t and x0 where x0 is the location of the interface, using the
        #  1D infinite boundary conditions with step function initial condition. pars are the fitted parameters as tuple.
        #pars, covs = curve_fit(conc,Yaxis, cRelative) 
        #  curve fits data with parameters D and x0 instead, using average time at position x to estimate the diffusiontime.
        pars, covs = curve_fit(conc,Yaxis,cRelative)
        dT,X0 = pars[0], pars[1]
        Dlist.append(dT/avgtime(i))
        # alternatively using D solo.
        #Dlist.append(dT) 
        """
        plt.plot(Yaxis,conc(Yaxis,dT,X0), label = "fit x = "+str(i+startx))
        plt.plot(Yaxis,cRelative, label = "actual data x ="+str(i+startx))
        plt.text(10**(-6)*120,0.3,f"$D_{str(i+1)} = {round(dT/avgtime((770-220)*pixlength),11)}\pm {round(np.sqrt(covs[0,0]),14)} $")
        plt.legend()
        plt.show()"""
    fig, axs = plt.subplots(1,2)
    axs[0].hist(Dlist[200:300], bins = 20 ,label = "D coeffs.")
    axs[1].hist(R_list(Dlist[200:300]), bins = 20, label = "Radii")    
    #axs[0].plot(Xaxis[200:300],np.asarray(Dlist[200:300]), label = "D coeffs.")
    #axs[1].plot(Xaxis[200:300],R_list(Dlist[200:300]), label = "Radii")
    fig.text(0.05, 0.5, 'numer of particles', va='center', rotation='vertical', size = 14)
    fig.text(0.65, 0.03, 'particle radius [$nm$]', ha='center', size = 14)
    fig.text(0.25, 0.03, 'Diffusion coefficient [$m^2/s$]', ha='center' , size = 14)
    print(sum(Dlist[200:300])/len(Dlist[200:300]), np.std(Dlist[200:300]))
    print(sum(R_list(Dlist[200:300]))/len(R_list(Dlist[200:300])), np.std(R_list(Dlist[200:300])))
    plt.show()
    
    return Yvals

def R_list(D_list): #calculates R from Stokes-Einstain eq and calculated Ds
    radii = []
    T = 273 + 20 #Kelvin
    K_b = 1.38 * 10 **(-23) #J/K
    nu = 1.0016 * 10**(-3) #Pa*s
    for elem in D_list:
        try:
            radii.append(K_b*T*10**9/(elem*6*np.pi*nu)) #radii in nm
        except:
            ZeroDivisionError  
    return radii


def findAv(xarray):
    def f(x,a):
        return a +0*x
    av,varAv = curve_fit(f,np.linspace(0,1,len(xarray)),xarray,p0 = 0.45)
    return av

def linearizeLuminocity(RGBarr: np.ndarray):
    GreyValWeight = [0.2126, 0.7152,0.0722] #from wikipedia
    RGBlin = np.array([RGBarr[0]/255, RGBarr[1]/255, RGBarr[2]/255])
    for i,Cs in enumerate(RGBlin):
        if Cs < 0.04045:
            RGBlin[i] = Cs/12.92
        else:
            RGBlin[i] = ((Cs+0.055)/(1.055))**2.4
    Y = np.dot(GreyValWeight,RGBlin)
    return Y

def avgtime(x):
    pixlength = 500*10**(-6)/215 #meters per pixel
    flowrate = 2*110/60 *10**(-9) # m^3 per sec
    crossArea = 75*10**(-6)*500*10**(-6) #m^2
    xtrue = (channelStart-x)*pixlength
    speed = flowrate/crossArea
    return xtrue/speed

def conc(ydata,dt,x0):
        return 1/2 - 1/2*erf((ydata-x0)/(2*np.sqrt(dt)))

def conc2(xdata,ydata,d,x0):
        return 1/2 - 1/2*erf((ydata-x0)/(2*np.sqrt(d*avgtime(xdata))))

def curvefitSci(xdata,ydata):
    return curve_fit(conc,xdata, ydata)


findDs()



"""for i in range(X):
    S[:,i] = np.dot((picTif[:,i]),[0.1140, 0.5870,0.2989])
    #print(len(picTif[:,i]))
    #S[i,:] = np.dot(picTif[i,:],[0.1140, 0.5870,0.2989])
plt.imshow(S,cmap="Greys")
def ru(x,y):
    return S[x,y]

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(Xs,Ys,ru(Xs,Ys), cmap="Blues")
"""


