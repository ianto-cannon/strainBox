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
  from scipy.interpolate import griddata
  for case in 'we_02  we_05  we_05T100  we_05nTime100  we_05tanhWindow  we_05window  we_05windowing  we_10'.split():
    print(case)
    fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$\\frac{kd}{2\\pi}$', rotation=0,size=22)#,labelpad=15)
    ax.set_ylabel('$\\frac{\\omega t_\\sigma }{2\\pi}$', rotation=0,size=22)
    df=np.zeros((51,86))
    count=0
    allData=[]
    for i in range(600):
      fname = '../'+case+f'/FSpec_run_break_{i:03}.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          #print('loadin ',fname)
          data = np.loadtxt(f).T
          try: df += data
          except ValueError: continue
          allData.append(data)
          count+=1
      except FileNotFoundError: continue
    print('count',count)
    if count<1:continue
    data = np.stack(allData, axis=0)
    l = np.nanpercentile(data, 75, axis=0)
    wavNumb = np.array([i/3 for i in range(86)])
    freq = np.array([i*td/timestep/100 for i in range(51)])
    for i in range(51):
      for j in range(86): 
        x = -df[i,j]*freq[i]*wavNumb[j]/count 
        #x = l[i,j]*freq[i]*wavNumb[j]
        df[i,j]=x
    maxV=np.max(np.abs(df))
    print('maxV',maxV)
    lels = np.array([-1e6,-1e5,-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])*4e-4
    colors = plt.cm.RdBu(np.linspace(0.05,.95,len(lels)-1))  # or define your own list of colors
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(lels, ncolors=cmap.N)
    X, Y = np.meshgrid(wavNumb,freq)
    if False:
      x=wavNumb
      y=freq
      Z=df
      # Interpolate to finer grid
      xi = np.linspace(x.min(), x.max(), 200)
      yi = np.linspace(y.min(), y.max(), 200)
      XI, YI = np.meshgrid(xi, yi)
      ZI = griddata((X.ravel(), Y.ravel()), Z.ravel(), (XI, YI), method='cubic')
    else: mappable = ax.contourf(X,Y,df, levels=lels, cmap=cmap, norm=norm)
    #mappable = ax.contourf(XI,YI,ZI, levels=lels, cmap=cmap, norm=norm)
    cBar = plt.colorbar(mappable, ax=ax)
    cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}kd\\omega t_\\sigma /(4\\pi^2\\epsilon)$')#, rotation=0, labelpad=15)
    cBar.set_ticks(lels)
    cBar.set_ticklabels(lels)
    cBar.ax.tick_params(which='both', direction='in', top=True, left=True)
    ax.set_xlim([wavNumb[1],wavNumb[-1]])
    ax.set_ylim([freq[1],freq[-1]])
    ax.set_xscale('log')
    ax.set_yscale('log')
    fname = 'plots/surPowFreqVsWav_'+case+'.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#energyVsFreq()
#energyFreqVsWav()
surPowFreqVsWav()
