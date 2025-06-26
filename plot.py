import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
plt.rcParams.update({"text.usetex": True,'font.size' : 14,})
timestep=.05
nyqTime=2*timestep
We=.05
rho=We
sigma=2*2**.5/3
diam=2*np.pi/3
vis=.219 #from measurements of rmsVel and Re_lambda=58
tsig=(rho*diam**3/sigma)**.5
td=diam**(2/3.)
cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0.8, 0.2, 3)] 

def energyVsFreq(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_ylabel('$E$', rotation=0)
  #ax.set_ylabel('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}\omega T /(2\\pi\\epsilon)$')#, rotation=0, labelpad=15)
  ax.set_xlabel('$\\omega T /2\\pi$', rotation=0)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([1,50])
  legs = ['rectangular', 'Hann', 'tanh']
  direcs = ['we_05nTime100', 'we_05windowing', 'we_05tanhWindow']
  for leg, direc in zip(legs, direcs):
    fname = f'../'+direc+'/ESpec_run_break_002.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    E = np.array([np.sum(df[:,i]) for i in range(df.shape[1])])
    #E = [E[i]*i for i in range(len(E))]
    if 'Hann' in leg :
      print('maxEHann',np.max(E))
      E*=2
    E[0]=0
    print('maxE',np.max(E))
    ax.plot(E,label=leg)
  ax.legend()
  fname = 'plots/energyVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyFreqVsWav(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_ylabel('$\\frac{kL}{2\\pi}$', rotation=0,size=22,labelpad=15)
  ax.set_xlabel('$\\omega T /2\\pi$', rotation=0)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([1,50])
  ax.set_ylim([1,85])
  df=np.zeros((86,51))
  for i in range(200):
    fname = f'../we_05windowing/ESpec_run_break_{i:03}.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df += np.loadtxt(f)
        print(f"found: {fname}")
    except FileNotFoundError:
      continue
  tot=np.sum(df)
  E = [[df[i, j] * i * j/tot for j in range(df.shape[1])] for i in range(df.shape[0])]
  #mappable = ax.pcolormesh(np.log10(E))#, vmin=-8, vmax=0)
  lels=range(-11,1,1)
  print('lels',lels)
  mappable = ax.contourf( np.log10(E), levels=lels)
  cBar = plt.colorbar(mappable, ax=ax)
  #cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0, labelpad=15)
  cBar.set_label('$\\log(\\mathbf{\\hat u\\cdot\\hat u^*}kL\\omega T /(4\\pi^2E)$')#, rotation=0, labelpad=15)
  fname = 'plots/energyFreqVsWav.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsWav(): 
  import matplotlib.colors as mcolors
  from matplotlib.ticker import ScalarFormatter
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_ylabel('$\\frac{kL}{2\\pi}$', rotation=0,size=22,labelpad=15)
  ax.set_xlabel('$\\omega T /2\\pi$', rotation=0)
  #ax.set_ylim([1,int(timeLen/2+1)])
  #spec=np.zeros((int(timeLen/2+1),85))
  df=np.zeros((86,51))
  count=0
  for i in range(200):
    fname = f'../we_05windowing/FSpec_run_break_{i:03}.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df += np.loadtxt(f)
        count+=1
        print(f"found: {fname}")
    except FileNotFoundError:
      continue
  for i in range(df.shape[0]):
    for j in range(df.shape[1]): 
      x = df[i,j]*i*j/count 
      #df[i,j] = np.log10(-x) 
      #df[i,j] = np.sign(x) * np.log1p(np.abs(x)) 
  maxV=np.max(np.abs(df))
  #mappable = ax.pcolormesh(df, cmap='bwr_r', vmin=-maxV, vmax=maxV)# levels=levels,
  #lels = [-1e5,-1e4,-1e2,-1e0,-1e-2,1e-2,1e0,1e2,1e4,1e5]
  lels = [-1e5,-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,1e-1,1e0,1e1,1e2,1e3,1e4,1e5]
  colors = plt.cm.bwr_r(np.linspace(0,1,len(lels)-1))  # or define your own list of colors
  cmap = mcolors.ListedColormap(colors)
  norm = mcolors.BoundaryNorm(lels, ncolors=cmap.N)
  mappable = ax.contourf(df, levels=lels, cmap=cmap, norm=norm)
  #mappable = ax.contourf(df, levels=[-3e4,-3e2,-3e0,-3e-2,-3e-4,3e-4,3e-2,3e0,3e2,3e4])#, cmap='bwr_r')
  cBar = plt.colorbar(mappable, ax=ax)
  #formatter = ScalarFormatter(useMathText=True)
  #formatter.set_powerlimits((-1, 1))
  #cBar.ax.yaxis.set_major_formatter(formatter)
  cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}kL\\omega T /(4\\pi^2\\epsilon)$')#, rotation=0, labelpad=15)
  cBar.set_ticks(lels)
  cBar.set_ticklabels(lels)
  ax.set_xlim([1,50])
  ax.set_ylim([1,85])
  ax.set_xscale('log')
  ax.set_yscale('log')
  fname = 'plots/surPowFreqVsWav.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#energyVsFreq()
#energyFreqVsWav()
surPowFreqVsWav()
