import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
plt.rcParams.update({"text.usetex": True,'font.size' : 14,})
timestep=.05
nu=1 #kinematic viscosity
diss = 1
sigma=2*2**.5/3
diam=2*np.pi/3
vis=.219 #from measurements of rmsVel and Re_lambda=58
td=diam**(2/3.)
cmap = plt.get_cmap('plasma')
freq = np.array([i*td/timestep/100 for i in range(51)])
wavNumb = np.array([i/3 for i in range(86)])
X, Y = np.meshgrid(wavNumb,freq)
wavLabel='$\\frac{kd}{2\\pi}$'
freqLabel='$\\frac{\\omega t_d}{2\\pi}$'
cases = {
    'we_02'          :{'We':.02, 'col':cmap(0)},
    'we_05PSpec'     :{'We':.05, 'col':cmap(.25)},
    'we_08'          :{'We':.08, 'col':cmap(.5)},
    'we_10'          :{'We':.10, 'col':cmap(.75)}
       }
for dirName, case in cases.items():
  case['rho']=case['We']
  case['WeJfm21'] = case['We']*3.63/.1
  case['tsig'] = ( case['rho'] * diam**3 / sigma ) **.5
  case['Oh'] = nu * case['rho']**.5  / ( diam * sigma) ** .5
  case['dHinze'] = .725 * sigma**(3/5) * case['rho']**(-3/5) * diss**(-2/5)

def energyVsFreq(): 
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70))
  #fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xscale('log')
    a.set_xlim([1,100])
    a.set_xlabel('$\\frac{kL}{2\\pi},\\,\\, \\frac{\\omega T_L}{2\\pi}$', size=20)
  ax[0].set_ylabel("$\\frac{2\\pi E(k)}{L u'^2},\\,\\, \\frac{2\\pi E(\\omega)}{T_L u'^2}$", size=20)
  ax[0].set_yscale('log')
  #ax[0].set_ylim([1e-8,1])
  ax[1].set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}\omega t_d}{2\\pi\\epsilon}$', rotation=0, labelpad=15, size=22)
  #ax[1].set_ylabel('$\\Pi$', rotation=0)
  secax = ax[1].secondary_xaxis(0, transform=ax[1].transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(labelbottom=False, labeltop=False, bottom=True, top=True, direction='inout')
  uRms = 1.25
  L=2*np.pi
  TL=L/uRms
  wav = np.array([i for i in range(86)])
  omega = np.array([i*TL/timestep/100 for i in range(51)])
  ax[0].plot(wav, wav**(-5/3), c='grey')
  for dirName, case in cases.items():
    print(dirName)
    count=0
    allSpec=[]
    E=np.zeros(51)
    Ewav=np.zeros(86)
    for i in range(1000):
      fname = '../'+dirName+f'/ESpec_run_break_{i:03}.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          df = np.loadtxt(f)
      except FileNotFoundError:
        continue
      count+=1
      E += (np.array([np.sum(df[:,i]) for i in range(df.shape[1])]) -E)/count
      Ewav += (np.array([np.sum(df[i,:]) for i in range(df.shape[0])]) -Ewav)/count
      
      fname = '../'+dirName+f'/FSpec_run_break_{i:03}.txt'
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
      F = np.array([np.sum(df[:,i]*omega[i]) for i in range(51)])
      #df[:,0]=0
      #F = np.array([np.sum(df[:,:i]) for i in range(51)])
      F/=case['We']
      F/=256**2
      allSpec.append(F)
    
    print('count',count)
    specData = np.stack(allSpec, axis=0)
    specData[:,0]=np.nan
    specData[:,-1]=np.nan
    l = np.nanpercentile(specData, 25, axis=0)
    u = np.nanpercentile(specData, 75, axis=0)
    m = np.nanmean(specData, axis=0)
    #diss = 2*E*omega**2
    #D = np.array([np.sum(diss[:i]) for i in range(51)])
    #ax[1].plot([td/case['tsig'], td/case['tsig']], [-8E4,8E4], c=case['col'], alpha=.3)
    ax[1].fill_between(omega, l, u, color=case['col'], alpha=0.3, edgecolor='none') 
    ax[1].plot(omega, m, c=case['col'], alpha=0.85)
    #ax[1].plot(omega, D, c=case['col'], alpha=0.85)
    
    invUMeanSq = 1.5 / sum(Ewav) 
    #invUMeanSq = 1 / 256**3
    Ewav *= invUMeanSq
    E *= invUMeanSq
    E[0]=np.nan
    E[-1]=np.nan
    Ewav[0]=np.nan
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[0].plot(omega, E, ls='dashed', c=case['col'])
    ax[0].plot(wav, Ewav, label=lbl, c=case['col'])
    
  ax[0].legend()
  fname = 'plots/energyVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowVsFreq(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  secax = ax.secondary_xaxis(0, transform=ax.transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(labelbottom=False, labeltop=False, bottom=True, top=True, direction='inout')
  ax.set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}\omega t_d}{2\\pi\\epsilon}$', rotation=0, labelpad=15, size=22)
  ax.set_xlabel(freqLabel, rotation=0, size=22)
  ax.set_xscale('log')
  ax.set_ylim([-8e4,8e4])
  #freq = np.array([i*td/nyqTime/50 for i in range(1,51)])
  ax.set_xlim([freq[1],freq[-1]])
  #dFreq=100*timestep/td
  #ax.plot([dFreq,dFreq],[-2E5,2E5],c='k',alpha=.3)
  for dirName, stat in cases.items():
    #if 'PSpec' in dirName: continue
    print(dirName)
    count=0
    allSpec=[]
    E=np.zeros(51)
    for i in range(1000):
      fname = '../'+dirName+f'/FSpec_run_break_{i:03}.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          df = np.loadtxt(f)
      except FileNotFoundError:
        continue
      count+=1
      F = np.array([np.sum(df[:,i]*freq[i]) for i in range(51)])
      F/=stat['We']
      if 'we_05windowing'in dirName: F*=4
      E += (F -E)/count
      allSpec.append(F)
    specData = np.stack(allSpec, axis=0)
    l = np.nanpercentile(specData, 25, axis=0)
    u = np.nanpercentile(specData, 75, axis=0)
    print('count',count)
    ax.plot([td/stat['tsig'],td/stat['tsig']],[-8E4,8E4],c=stat['col'],alpha=.3)
    ax.fill_between(freq,l,u,color=stat['col'],alpha=0.3,edgecolor='none') 
    ax.plot(freq,E,label=dirName,c=stat['col'],alpha=0.85)
  ax.legend()
  fname = 'plots/surPowVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyWindowsVsFreq(): 
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
  fname = 'plots/energyVsFreq_.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyFreqVsWav(): 
  for case in 'we_02 we_05windowing  we_10 we_05PSpec'.split():
    if 'PSpec' in case: 
      specTyp='/P'
    else: 
      continue#specTyp='/E'
    fig, ax = plt.subplots(1)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$\\frac{kL}{2\\pi}$', rotation=0,size=22)
    ax.set_ylabel('$\\frac{\\omega T}{2\\pi}$', rotation=0,size=22,labelpad=15)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1,50])
    ax.set_xlim([1,85])
    df=np.zeros((51,86))
    for i in range(1000):
      fname = '../'+case+specTyp+f'Spec_run_break_{i:03}.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          df += np.loadtxt(f).T
          print(f"found: {fname}")
      except FileNotFoundError:
        continue
    tot=np.sum(df)
    #E = [[df[i, j] * i * j/tot for j in range(df.shape[1])] for i in range(df.shape[0])]
    E = df/tot 
    #mappable = ax.pcolormesh(np.log10(E))#, vmin=-8, vmax=0)
    #lels=range(-11,1,1)
    #mappable = ax.contourf( np.log10(E))#, levels=lels)
    mappable = ax.contourf( E, levels=10)
    cBar = plt.colorbar(mappable, ax=ax)
    #cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0, labelpad=15)
    cBar.set_label('$\\log(\\mathbf{\\hat u\\cdot\\hat u^*}kL\\omega T /(4\\pi^2E)$')#, rotation=0, labelpad=15)
    fname = 'plots/energyFreqVsWav_'+case+'.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def phaseFreqVsWav(): 
  for dirName, stat in cases.items():
    if 'we_05PSpec' not in dirName: continue
    fig, ax = plt.subplots(1)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel(wavLabel, rotation=0,size=22)
    ax.set_ylabel(freqLabel, rotation=0,size=22)#,labelpad=15)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([freq[1],freq[-1]])
    ax.set_xlim([wavNumb[1],wavNumb[-1]])
    df=np.zeros((51,86))
    count=0
    for i in range(1000):
      fname = '../'+dirName+f'/PSpec_run_break_{i:03}.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          count+=1
          df += np.loadtxt(f).T
          #print(f"found: {fname}")
      except FileNotFoundError:
        continue
    print('count',count)
    tot=np.sum(df)
    E = df/count/np.pi/(256/3)**2/sigma#/tot 
    print('tot',np.sum(df))
    ax.plot([wavNumb[1],wavNumb[-1]],[td/stat['tsig'],td/stat['tsig']],c=stat['col'],alpha=.3)
    mappable = ax.contourf(X, Y, E, levels=10, cmap='Greys')
    cBar = plt.colorbar(mappable, ax=ax)#, orientation='horizontal')
    cBar.set_label('$\\alpha\\mathbf{k\\cdot k}\\hat c^2 /(\\pi d^2\\sigma)$')
    fname = 'plots/phaseFreqVsWav.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsWav(): 
  import matplotlib.colors as mcolors
  from scipy.interpolate import griddata
  #for case in 'we_02  we_05  we_05T100  we_05nTime100  we_05tanhWindow  we_05window  we_05windowing  we_10'.split():
  for dirName, stat in cases.items():
    if 'PSpec' in dirName: continue
    print(dirName)
    fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel(wavLabel, rotation=0,size=22)#,labelpad=15)
    ax.set_ylabel(freqLabel, rotation=0,size=22)
    df=np.zeros((51,86))
    count=0
    allData=[]
    for i in range(600):
      fname = '../'+dirName+f'/FSpec_run_break_{i:03}.txt'
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
    for i in range(51):
      for j in range(86): 
        x = -df[i,j]*freq[i]*wavNumb[j]/count 
        #x = l[i,j]*freq[i]*wavNumb[j]
        df[i,j]=x
    maxV=np.max(np.abs(df))
    print('maxV',maxV)
    lels = np.array([-1e6,-1e5,-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])*1e-2
    colors = plt.cm.RdBu(np.linspace(0.05,.95,len(lels)-1))  # or define your own list of colors
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(lels, ncolors=cmap.N)
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
    ax.plot([wavNumb[1],wavNumb[-1]],[td/stat['tsig'],td/stat['tsig']],c=stat['col'],alpha=.3)
    ax.set_xlim([wavNumb[1],wavNumb[-1]])
    ax.set_ylim([freq[1],freq[-1]])
    ax.set_xscale('log')
    ax.set_yscale('log')
    fname = 'plots/surPowFreqVsWav_'+dirName+'.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

energyVsFreq()
#surPowVsFreq()
#energyWindowsVsFreq()
#energyFreqVsWav()
#phaseFreqVsWav()
#surPowFreqVsWav()
