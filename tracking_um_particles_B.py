import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from math import log, pi
import numpy as np
# Globals
t = 100 / 7.5 #s
pixel_size = 2.5 * 10**(-7) #m/pixel
T = 273 + 20 #Kelvin
K_b = 1.38 * 10 **(-23) #J/K
nu = 1.0016 * 10**(-3) #Pa*s
# Globals

def parse_file(file):
    tree = ET.parse(file)
    root = tree.getroot()
    return tree, root

def remove_noise(file): #removes particles tracked in 10 or less frames
    tree, root = parse_file(file)
    for particle in root.findall('particle'):
        if int(particle.attrib['nSpots'] ) <= 10:
            root.remove(particle)
    tree.write(file)


treeA1, rootA1 = parse_file('ResultofDF40x_7_Tracks_B1_3-3 – Kopi.xml')

treeA2, rootA2 = parse_file('ResultofDF40x_7_Tracks_B2_3-3 – Kopi.xml')

treeA3, rootA3 = parse_file('ResultofDF40x_7_Tracks_B3_3-3 – Kopi.xml')

def move(i,para):
    return float(i[-1].attrib[para]) -float(i[0].attrib[para])

def make_D_list(particle): #calculates D from displacement
    D_coeff_list = []
    for i in particle:
        xmove = move(i,"x") 
        ymove = move(i,"y") 
        tmove = move(i,"t") 
        displace2 = (xmove**2+ymove**2)*pixel_size**2
        D = displace2/(4*tmove*t)
        if displace2 > 16*pixel_size**2:
            D_coeff_list.append(D)
    return D_coeff_list

def make_R_list(D_list): #calculates R from Stokes-Einstain eq and calculated Ds
    radii = []
    for elem in D_list:
        try:
            radii.append(K_b*T*10**6/(elem*6*pi*nu)) #radii in um
        except:
            ZeroDivisionError  
    return radii

def avg_D_const(D_list):
    return sum(D_list)/len(D_list)

def avg_R(R_list): #lagde vist to funksjoner for å regne gjennomsnitt:)
    return sum(R_list)/len(R_list)

def SD(datalist):
    calc_list = []
    avg = avg_D_const(datalist)

    for data in datalist:
        calc_list.append((data-avg)**2)
    return np.sqrt(sum(calc_list)/(len(datalist)-1))

def dataDict(filenames : tuple, hshshs : int) -> dict:
    fig, ax = plt.subplots(len(filenames), sharex= False, constrained_layout = False )
    figD, axD = plt.subplots(len(filenames), sharex = True)
    # create loopable structure
    fileDict = dict({})
    for file in filenames:
        tree, root = parse_file(file)
        fileDict[file[-5:-3]] = dict({"Dlist" :  make_D_list(root),
                                      "Rlist" :  make_R_list(root),
                                      "avgD"  :  avg_D_const(make_D_list(root)),
                                      "avgR"  :  avg_R(make_R_list(root)),
                                      "particles"  :  root,
                                      "tree"  :  tree})
    return fileDict     

def plotPaths(tree):
    plt.xlim(0,1360)
    plt.ylim(-1024,0)
    for j, particle in enumerate(tree.getroot()):
        data = np.zeros((len(particle),2))
        for i, detection in enumerate(particle):
            data[i] = np.asarray([detection.attrib["x"],detection.attrib["y"]])
        plt.plot(data[:,0],-data[:,1])
    plt.show()

#plotPaths(treeA1)
#plotPaths(treeA2)
#plotPaths(treeA3)

fig, ax = plt.subplots(3,1, sharex= False, constrained_layout = False )
figD , axD = plt.subplots(3,1 , sharex= True)

D_coeff_list_A1 = make_D_list(rootA1)
R_list_A1 = make_R_list(D_coeff_list_A1)
avg_D_A1 = avg_D_const(D_coeff_list_A1)
avg_R_A1 = avg_R(R_list_A1)
print(f'Average D in B1 is {avg_D_A1} +- {SD(D_coeff_list_A1)}, {len(D_coeff_list_A1)}')
print(f'Average radius in B1 is {avg_R_A1} +- {SD(R_list_A1)}')

D_coeff_list_A2 = make_D_list(rootA2)
R_list_A2 = make_R_list(D_coeff_list_A2)
avg_D_A2 = avg_D_const(D_coeff_list_A2)
avg_R_A2 = avg_R(R_list_A2)
print(f'Average D in B2 is {avg_D_A2} +- {SD(D_coeff_list_A2)}, {len(D_coeff_list_A2)}')
print(f'Average radius in B2 is {avg_R_A2} +- {SD(R_list_A2)}')

D_coeff_list_A3 = make_D_list(rootA3)
R_list_A3 = make_R_list(D_coeff_list_A3)
avg_D_A3 = avg_D_const(D_coeff_list_A3)
avg_R_A3 = avg_R(R_list_A3)
print(f'Average D in B3 is {avg_D_A3} +- {SD(D_coeff_list_A3)}, {len(D_coeff_list_A3)}')
print(f'Average radius in B3 is {avg_R_A3} +- {SD(R_list_A3)}') 

totlistD = D_coeff_list_A1 + D_coeff_list_A2 + D_coeff_list_A3
totlistR = R_list_A1 + R_list_A2 + R_list_A3

print(f"Total avg D in A: {avg_D_const(totlistD)} +- {SD(totlistD)} ")
print(f"total avg R in A: {avg_R(totlistR)} +- {SD(totlistR)}")
print(len(D_coeff_list_A1), len(D_coeff_list_A2), len(D_coeff_list_A3))

#print(np.percentile(totlistR, [25, 50, 75], interpolation= 'midpoint'))
#print(np.percentile(totlistD, [25, 50, 75], interpolation= 'midpoint'))

print(min(totlistD), max(totlistD))
print(min(totlistR), max(totlistR))

_, bins_r = np.histogram(R_list_A1, bins = 40)
logbins_r = np.logspace(np.log10(bins_r[0]), np.log10(bins_r[-1]), len(bins_r))
_, bins_D = np.histogram(D_coeff_list_A1, bins = 40)
logbins_D = np.logspace(np.log10(bins_D[0]), np.log10(bins_D[-1]), len(bins_D))

ax[0].hist(R_list_A1, bins = logbins_r)
ax[0].set_xscale('log')
ax[1].hist(R_list_A2, bins = logbins_r)
ax[1].set_xscale('log')
ax[2].hist(R_list_A3, bins = logbins_r)
ax[2].set_xscale('log')

axD[0].hist(D_coeff_list_A1, bins = logbins_D)
axD[0].set_xscale('log')
axD[1].hist(D_coeff_list_A2, bins = logbins_D)
axD[1].set_xscale('log')
axD[2].hist(D_coeff_list_A3, bins = logbins_D)
axD[2].set_xscale('log')

ax[0].annotate('B1', (0.5, 1.05), xycoords = 'axes fraction')
ax[1].annotate('B2', (0.5, 1.05), xycoords = 'axes fraction')
ax[2].annotate('B3', (0.5, 1.05), xycoords = 'axes fraction')
fig.text(0.05, 0.5, 'numer of tracked particles', va='center', rotation='vertical')
fig.text(0.5, 0.03, 'particle radius [$\mu m$]', ha='center')

axD[0].annotate('B1', (0.5, 1.05), xycoords = 'axes fraction')
axD[1].annotate('B2', (0.5, 1.05), xycoords = 'axes fraction')
axD[2].annotate('B3', (0.5, 1.05), xycoords = 'axes fraction')
figD.text(0.05, 0.5, 'numer of tracked particles', va='center', rotation='vertical')
figD.text(0.5, 0.03, 'Diffusion coefficient [$m^2/s$]', ha='center')

for ax in ax.flat:
    ax.label_outer()

for axD in axD.flat:
    ax.label_outer()   
#plt.show()

fig.savefig('histR_B.jpg')
figD.savefig('histD_B.jpg')



