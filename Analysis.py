import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.extract import find
import trackpy as pt
import tifffile as tf
from scipy.ndimage import rotate
from scipy.special import erf
from scipy.optimize import curve_fit

#Global variables ----------------------------
pixlength = 500*10**(-6)/210    #meters per pixel
channelStart = 920              #pixels, lower is further out in channel
plt.style.use("seaborn")        #Stylish
# --------------------------------------------

def generateYvalData(filename = "fluidics\\110ul_4xBF_startCanal.tif"):
    #  Crop and rotate image data array
    if filename == "fluidics\\110ul_4xBF_startCanal.tif":
        rotPic = rotate(tf.imread(filename),1.6)
        pictif = rotPic[130:315,20:-10]#slicing
    else:
        pictif = rotate(tf.imread(filename)[670:-100,0:-300],1.18)[28:-40,5:-15]
    #  generates array to put relative luminance-values
    Yvals = np.zeros((len(pictif[:,1]),len(pictif[1,:]))) 
    #  Creates X, Y arrays for looping
    Ys,Xs = np.arange(len(pictif[:,1])), np.arange(len(pictif[1,:]))
    for j in Ys:
        for i in Xs:
            Yvals[j,i] = linearizeLuminocity(pictif[j,i]) #generates array of values for relative intensity
    return Yvals

def findDs(interval = (50,-50),pltD=False,endx=-1, startx = 220, filename = "fluidics\\110ul_4xBF_startCanal.tif", showplots = [] ):
    Yvals = generateYvalData(filename)[interval[0]:interval[1],startx:endx] # generates intervals of interest
    Dlist = []
    Yaxis = np.arange(len(Yvals[:,0]))*pixlength
    Xaxis = np.arange(len(Yvals[0,:]))
    sentinel = 0
    if len(showplots):
        fig, axarr = plt.subplots(len(showplots)//2,len(showplots)//2)
        axarr = axarr.flatten()
        axind = len(showplots)-1
        fig.tight_layout()
        labels = ["(a)","(b)","(c)","(d)"]
        lowInt = findAv(Yvals[10,10:-10])
        highInt = findAv(Yvals[-10,10:-10]) 
    for i in Xaxis:
        #  finds lowest and highest intensity that minimizes squared distance to data
        #  in the interval of interest, 30 outer points where concentration should be unaffectd by diffusion
        lowInt = findAv(Yvals[0:20,i])
        highInt = findAv(Yvals[-20:,i]) 
        #  calculates C(x = X,y)/C_max according to Beer-Lamberts' law.
        cRelative = np.log10(highInt/Yvals[:,i])/np.log10(highInt/lowInt) 
        #  curve fits data with parameters D*t and x0 where x0 is the location of the interface, using the
        #  1D infinite boundary conditions with step function initial condition. pars are the fitted parameters as tuple.
        #pars, covs = curve_fit(conc,Yaxis, cRelative) 
        #  curve fits data with parameters D and x0 instead, using average time at position x to estimate the diffusiontime.
        pars, covs = curve_fit(conc,Yaxis,cRelative)
        Dt,X0 = pars[0], pars[1]
        D = Dt/avgtime(i)
        #if D > 10**(-9): #filter out outliers
        Dlist.append(D)
        sentinel+=1
        if i+startx in showplots:
            X = channelStart-i-startx
            Y0 = Yaxis - X0
            axarr[axind].plot(-(Y0)*10**6,conc(Yaxis,Dt,X0), label = "Fit")
            #axarr[axind].plot([],[])
            axarr[axind].plot(-(Y0)*10**6,cRelative, label = "Data")
            axarr[axind].text(2,0.01,labels[axind]+": $D "+f"= {round(D*10**12,0)}" + " \mathrm{\mu m^2 /s} $",fontsize=12) #\pm {round(np.sqrt(covs[0,0])*10**(12),0)}
            axarr[axind].legend()
            #axarr[axind].set_title(f"X = {X}")
            axarr[axind].set_xlabel("Y coordinate [$\mathrm{\mu m}]$")
            axarr[axind].set_ylabel("Relative concentration")
            axind-=1
    plt.show()
    if pltD:
        plt.plot(Xaxis[:sentinel],Dlist)
        fig, axs = plt.subplots(2,1)
        rs =np.asarray((R_list(Dlist))) 
        ds = np.asarray(Dlist)*10**12
        axs[0].hist(ds, label = "D coeffs.", bins = 20)
        axs[0].set_xlabel("Diffusion coefficient [$\mathrm{\mu m^2/s}$]")
        axs[0].set_ylabel("Number of measurements")
        axs[0].set_xlim([1000,3800])
        axs[0].axvline(ds.mean(),color = "black", label="$\langle D \\rangle$")
        axs[0].legend()
        print(f"avg D is {ds.mean()} +- {ds.std(ddof=1)}, avg r is {rs.mean()} +- {rs.std(ddof=1)}")
        axs[1].hist(rs, label = "Stokes radii",bins = 18)
        axs[1].set_ylabel("Number of measurements")
        axs[1].set_xlabel("Calculated Stokes' radius [$\mathrm{nm}$]")
        axs[1].set_xlim([0.05,0.165])
        axs[1].axvline(rs.mean(),color = "black", label="$\langle r \\rangle$")
        axs[1].legend()
        plt.show()
    axarr[0].set_xlabel("")  
    axarr[1].set_xlabel("") 
    axarr[3].set_ylabel("")
    axarr[1].set_ylabel("")

    plt.show()
    return Dlist

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
    pixlength = 500*10**(-5)/215 #meters per pixel
    flowrate = 2*110/60 *10**(-9) # m^3 per sec
    crossArea = 75*10**(-6)*500*10**(-6) #m^2
    xtrue = (channelStart-x)*pixlength #x-axis is flipped.
    speed = flowrate/crossArea
    return xtrue/speed

def conc(ydata,dt,x0):
        return 1/2 - 1/2*erf((ydata-x0)/(2*np.sqrt(dt)))

def conc2(xdata,ydata,d,x0):
        return 1/2 - 1/2*erf((ydata-x0)/(2*np.sqrt(d*avgtime(xdata))))

def curvefitSci(xdata,ydata):
    return curve_fit(conc,xdata, ydata)

def showimage(filename):
    rotPic = rotate(tf.imread(filename),1.6)[130:315,20:-10]
    plt.imshow(rotPic)
    plt.show()
#showimage("fluidics\\110ul_4xBF_startCanal.tif")    
findDs(interval=(60,-60),showplots=[160,325,490,555],startx=160, pltD=True,endx=900)# filename="fluidics\\Vik og Johansen.tif")
#findDs(interval=(60, -60),startx=160, pltD=False,endx=650)# filename="fluidics\\Vik og Johansen.tif")


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


