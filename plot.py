import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
plt.rcParams.update({"text.usetex": True,'font.size' : 14,})
timestep=.05
diss = 1
sigma=2*2**.5/3
L=2*np.pi
diam=L/3#*.3
uRms = 1.25
TL=L/uRms
#ud=(diam/L)**(1/3.)*uRms
ud=(diss*diam/L)**(1/3.)
Kmax=2*np.pi*256/3/L
lKol=4/Kmax
print(f"lKol {lKol/L:.3g}")
nu=(diss*lKol**4)**(1/3)
print(f"nu {nu:.3g}")
#td=diam**(2/3.)*L**(1/3.)/uRms
td=diam/ud
print(f"td {td/TL:.3g}")
tvis=diam**2/nu
print(f"tvis {tvis/TL:.3g}")
tKol=(nu/diss)**.5
print(f"tKol {tKol/TL:.3g}")
lTay=(15*nu*uRms**2/diss)**.5
print(f"lTay {lTay/L:.3g}")
ReTay=uRms*lTay/nu
print(f"ReTay {ReTay:.3g}")
ReDrop=(diss*diam)**(1/3)*diam/nu
print(f"ReDrop {ReDrop:.3g}")
print(f"normDiss {diss*.42*L/uRms**3:.3g}")
cmap = plt.get_cmap('plasma')
freq = np.array([i*td/timestep/100 for i in range(51)])
wavNumb = np.array([i/3 for i in range(86)])
X, Y = np.meshgrid(wavNumb,freq)
wavLabel='$\\frac{\kappa d}{2\\pi}$'
freqLabel='$\\frac{\\omega t_d}{2\\pi}$'
cases = {
    #'we_02'          :{'rho':.2, 'col':cmap(0)},
    #'we_05windowing' :{'rho':.5, 'col':'k'},
    #'we_08'          :{'rho':.8, 'col':cmap(.5)},
    #'we_10'          :{'rho':1, 'col':cmap(.75)},
    'we_02PSpec'     :{'rho':.2, 'col':cmap(0)},
    'we_05PSpec'     :{'rho':.5, 'col':'k'},
    'we_08PSpec'     :{'rho':.8, 'col':cmap(.5)},
    'we_10PSpec'     :{'rho':1, 'col':cmap(.75)},
    'we_inf'         :{'rho':1e8, 'col':'k'},
       }
for dirName, case in cases.items():
  case['We'] = case['rho'] * (diss*diam)**(2/3) * diam / sigma #gives We from VelaSciAdv
  case['tsig'] = ( case['rho'] * diam**3 / sigma ) **.5
  case['Oh'] = nu * case['rho']**.5  / ( diam * sigma) ** .5
  case['dHinze'] = .725 * sigma**(3/5) * case['rho']**(-3/5) * diss**(-2/5)
  print(f"{dirName} We {case['We']:.5g} Oh {case['Oh']:.3g} dHinze {case['dHinze']/2/np.pi:.3g} tsig {case['tsig']/TL:.3g}")

def energyVsFreq(): 
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70))
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.5, hspace=0.05)
  import cmasher as cmr
  cmap = cmr.lavender
  #colors = [cmap(i) for i in [.8,.5,0]]#np.linspace(0.8, 0.2, 3)] 
  blu = cmap(.8)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xscale('log')
    a.set_xlim([1,100])
    a.text(.4, -.17, '$\\frac{kL}{2\\pi},$', size=20, transform=a.transAxes)
    a.text(.53, -.17, '$\\frac{\\omega T}{2\\pi}$', size=20, transform=a.transAxes, c=blu)
  ax[0].text(-0.36, 0.55, r"$\frac{2\pi E(\kappa)}{L u'^2},$", fontsize=20, transform=ax[0].transAxes)
  ax[0].text(-0.36, 0.35, r"$\frac{2\pi E(\omega)}{T u'^2}$", fontsize=20, transform=ax[0].transAxes, c=blu)
  ax[0].set_yscale('log')
  ax[0].set_ylim([1e-8,1])
  ax[1].text(-0.36, 0.55, r"$\frac{\Phi(\kappa)\kappa L}{2\pi u'^2},$", fontsize=20, transform=ax[1].transAxes)
  ax[1].text(-0.36, 0.35, r"$\frac{\Phi(\omega)\omega T}{2\pi u'^2}$", fontsize=20, transform=ax[1].transAxes, c=blu)
  ylim=[-.8,.6]
  ax[1].set_ylim(ylim)
  ax[1].set_yticks([-.8,-.4,0,.4])
  secax = ax[1].secondary_xaxis(0, transform=ax[1].transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(labelbottom=False, labeltop=False, bottom=True, top=True, direction='inout')
  insax = False#ax[1].inset_axes((.5,.05,.45,.35))
  if insax:
    insax.xaxis.set_label_position("top")
    insax.xaxis.tick_top()
    insax.set_xlabel('$\\mathrm{We}$', rotation=0)
    insax.text(-.4, .65, r"$\frac{l_\Phi}{d},$", fontsize=20, transform=insax.transAxes)
    insax.text(-.4, .25, r"$\frac{t_\Phi}{t_d}$", fontsize=20, transform=insax.transAxes)
    insax.set_xlim([0,4])
    insax.set_ylim([.6,1.4])
    insax.set_xticks([0,1,2,3,4])
    insax.set_yticks([.6,1,1.4])
    insax.tick_params(which='both', direction='in', top=True, right=True, bottom=True)
  wav = np.array([1,100])
  ax[0].plot([3,3], [1e-8,1], c='k', lw=0.5)
  ax[1].plot([3,3], ylim, c='k', lw=0.5)
  ax[0].plot([TL/td,TL/td], [1e-8,1], c=blu, lw=0.5)
  ax[1].plot([TL/td,TL/td], ylim, c=blu, lw=0.5)
  ax[0].plot(wav, wav**(-5/3), ls='dotted', c='k')
  wav = np.array([i for i in range(86)])
  omega = np.array([i*TL/timestep/100 for i in range(51)])
  for dirName, case in cases.items():
    if '05' not in dirName and 'inf' not in dirName: continue
    print(dirName)
    count=0
    allSpec=[]
    allWav=[]
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
      F = np.array([np.sum(-df[:,i]*omega[i]) for i in range(51)])
      F/=case['rho']
      F/=256**2
      allSpec.append(F)
      F = np.array([np.sum(-df[i,:]*wav[i]) for i in range(86)])
      F/=case['rho']
      F/=256**2
      allWav.append(F)
    print('count',count)
    specData = np.stack(allSpec, axis=0)
    specData[:,0]=np.nan
    specData[:,-1]=np.nan
    m = np.nanmean(specData, axis=0)
    ax[1].plot(omega, m, c=blu)
    sigInd = int(TL/case['tsig'])
    ax[1].plot(sigInd, m[sigInd], 'o', c=blu)#case['col'])#, mfc='None')
    for i in range(len(m)-1):
      if m[i]*m[i+1]<0:
        flipTime = 1 / ( omega[i] - m[i] * (omega[i+1] - omega[i]) / (m[i+1] - m[i]) )
        if insax: insax.plot(case['We'], TL*flipTime/td, 'o', c=case['col'], mfc='None', clip_on=False)
    specData = np.stack(allWav, axis=0)
    specData[:,0]=np.nan
    specData[:,-1]=np.nan
    m = np.nanmean(specData, axis=0)
    ax[1].plot(wav, m, c=case['col'])
    hinzInd = int(2*np.pi/case['dHinze'])
    if hinzInd<len(m): ax[1].plot(hinzInd, m[hinzInd], 'o', c=case['col'])
    for i in range(len(m)-1):
      if m[i]*m[i+1]<0:
        flipLen = 1 / ( wav[i] - m[i] * (wav[i+1] - wav[i]) / (m[i+1] - m[i]) )
        if insax: insax.plot(case['We'], flipLen*3, 'o', c=case['col'], clip_on=False)
    invUMeanSq = 1.5 / sum(Ewav) 
    Ewav *= invUMeanSq
    E *= invUMeanSq
    E[0]=np.nan
    E[-1]=np.nan
    Ewav[0]=np.nan
    lbl = rf"$\mathrm{{We}}={case['We']:.1f}$"
    if 'inf' in dirName: liSty='dashed'
    else: liSty='solid'
    ax[0].plot(omega, E, c=blu, ls=liSty)
    ax[0].plot(sigInd, E[sigInd], 'o', c=blu)#, mfc='None')
    ax[0].plot(wav, Ewav, label=lbl, c=case['col'], ls=liSty)
    if hinzInd<len(Ewav): ax[0].plot(hinzInd, Ewav[hinzInd], 'o', c=case['col'])
  fname = '../plots/energyVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def interfaceVsFreq():
  fig, ax = plt.subplots(1, 3, figsize=(11, 3), constrained_layout=True)
  insax = [ax[0].inset_axes((.5, .05, .45, .35)), 
           ax[1].inset_axes((.5, .05, .45, .35)),
           ax[2].inset_axes((.07, .05, .45, .35))]
  ax[0].set(yscale='log')
  ax[2].set(yscale='log')
  ax[0].set_ylim(1e-8,1e4)
  ax[1].set_ylim(-65e-4,4e-3)
  #ax[2].set_ylim(1e-3,2e-1)
  ax[0].set_ylabel(r"$\tilde E(\kappa),\,\,\,\,\tilde E(\omega)$")
  ax[1].set_ylabel(r"$\tilde F(\kappa)\kappa d/2\pi,\,\,\,\, \tilde F(\omega)\omega t_d/2\pi$")
  #ax[1].set_ylabel(r"$\frac{1}{2}F(\kappa)\kappa d\pi^{-1}u'^{-2}, \frac{1}{2}F(\omega)\omega t_d\pi^{-1}u'^{-2}$")
  #ax[1].set_ylabel(r"$F(\kappa)\kappa d/(2\pi u'^2), F(\omega)\omega t_d/(\pi u'^2)$")
  #ax[2].set_ylabel(r"$2A(\kappa)d^{-3}, 2A(\omega)d^{-2}t_d^{-1}$")
  ax[2].set_ylabel(r"$\tilde A(\kappa),\,\,\,\, \tilde A(\omega)$")
  ax[0].plot(wavNumb, 1e3*wavNumb**(-5/3), ls='solid', c='k', lw=.7)
  ax[0].text(.6, .8, r'$10^3\left(\frac{\kappa d}{2\pi}\right)^{-\frac{5}{3}}$', transform=ax[0].transAxes, rotation=-8)
  insax[0].set_visible(False)
  insax[1].text(-.5, .75, r"$\frac{l_\Phi}{d},$", fontsize=20, transform=insax[1].transAxes)
  insax[1].text(-.5, .15, r"$\frac{t_\Phi}{t_d}$",fontsize=20, transform=insax[1].transAxes)
  insax[2].text(1.3, .5,  r"$\frac{A}{\pi d^2}$", fontsize=20, transform=insax[2].transAxes)
  insax[2].yaxis.set_label_position("right")
  insax[2].yaxis.tick_right()
  insax[1].set_ylim(.5,1)
  #insax[2].set_ylim(1,2.3)
  #insax[2].set_yticks([1,1.5,2,2.5])
  secax = ax[1].secondary_xaxis(0, transform=ax[1].transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(labelbottom=False, labeltop=False, bottom=True, top=True, direction='inout')
  for a, ins in zip(ax,insax):
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xscale('log')
    a.set_xlim(wavNumb[1], wavNumb[-1])
    a.text(.37, -.14, r'$\frac{\kappa d}{2\pi},$', size=20, transform=a.transAxes)
    a.text(.56, -.14, r'$\frac{\omega t_d}{2\pi}$', size=20, transform=a.transAxes)
    ins.set(xlim=(0, 4), xticks=[0,1,2,3,4])#ylim=(.5, 1.5), 
    ins.xaxis.set_label_position("top")
    ins.xaxis.tick_top()
    ins.tick_params(which='both', direction='in', left=True, top=True, right=True, bottom=True)
    ins.set_xlabel(r'$\mathrm{We}$')
  def load(name):
    try: return np.loadtxt(f'../plottedData/{name}.txt')
    except FileNotFoundError: return None
  def crossings(x, y):
    i = np.where(y[:-1]*y[1:] < 0)[0]
    return 1/(x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i]))
  for dirName, case in reversed(cases.items()):
    print(dirName)
    c = 'grey' if 'inf' in dirName else plt.cm.PiYG((case['We']-.72)/10)
    sig = int(TL/case['tsig'])
    hinz = int(2*np.pi/case['dHinze'])
    for i, (tag, func) in enumerate([
      ("ESpecdf", lambda df: (
        df.sum(1)*2/ud**2*2*np.pi/td, 
        df.sum(0)*2/ud**2*2*np.pi/diam)),
      ("FSpecdf", lambda df: (
        -(df/case['rho'] * freq[:,None]).sum(1)*2/ud**2,
        -(df/case['rho'] * wavNumb[None,:]).sum(0)*2/ud**2*td/diam
      )),
      ("PSpecdf", lambda df: ( 
        df.sum(1)/(np.pi*(diam/2/np.pi)**2)*2*np.pi/td, 
        df.sum(0)/(np.pi*(diam/2/np.pi)**2)*2*np.pi/diam ) )
    ]):
      df = load(dirName + tag)
      if df is None: continue
      y, yk = func(df)
      y[[0,-1]] = yk[[0,-1]] = np.nan
      a = ax[i]
      ins = insax[i]
      a.plot(freq, y, c=c, ls='dashed')
      a.plot(freq[sig], y[sig], 'o', c=c, mfc='none')
      a.plot(wavNumb, yk, c=c)
      if hinz < len(yk): a.plot(wavNumb[hinz], yk[hinz], 'o', c=c)
      if False and tag == "ESpecdf":
        We = case['We']
        if case['We']>5: We=5
        ins.plot(We, np.nansum( yk[1:] * ( wavNumb[1:] - wavNumb[:-1] ) ),'o', c=c, clip_on=False)
      if tag == "FSpecdf":
        for z in crossings(freq, y): ins.plot(case['We'], z, 'o', c=c, mfc='none', clip_on=False)
        for z in crossings(wavNumb, yk): ins.plot(case['We'], z, 'o', c=c, clip_on=False)
      if tag == "PSpecdf" and 'inf' not in dirName:
        #ins.plot(case['We'], np.nansum(yk),'o', c=c, clip_on=False)
        ins.plot(case['We'], np.nansum( yk[1:] * ( wavNumb[1:] - wavNumb[:-1] ) ),'o', c=c, clip_on=False)
  fig.savefig('../plots/interfaceVsFreq.pdf', bbox_inches='tight', pad_inches=0.01, transparent=True)

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
  for dirName, case in cases.items():
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
      F/=case['rho']
      if 'we_05windowing'in dirName: F*=4
      E += (F -E)/count
      allSpec.append(F)
    specData = np.stack(allSpec, axis=0)
    l = np.nanpercentile(specData, 25, axis=0)
    u = np.nanpercentile(specData, 75, axis=0)
    print('count',count)
    ax.plot([td/case['tsig'],td/case['tsig']],[-8E4,8E4],c=case['col'],alpha=.3)
    ax.fill_between(freq,l,u,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.plot(freq,E,label=dirName,c=case['col'],alpha=0.85)
  ax.legend()
  fname = '../plots/surPowVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyWindowsVsFreq(): 
  fig, ax = plt.subplots()
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
  fname = '../plots/energyVsFreq_.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyFreqVsWav(var): 
  import matplotlib.colors as mcolors
  import matplotlib.ticker as mticker
  for dirName, case in cases.items():
    if 'inf' in dirName and 'ESpec' not in var: continue
    print(dirName)
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.set_box_aspect(1)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$\\frac{\kappa d}{2\\pi}$', rotation=0,size=22,labelpad=-5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([wavNumb[1],wavNumb[-1]])
    ax.set_ylim([freq[1],freq[-1]])
    savname = '../plottedData/'+dirName+var+"df.txt"
    try:
      with open(savname, encoding = 'utf-8') as f:
        print('load ', savname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print('not found ', savname)
      df=np.zeros((51,86))
      count=0
      for i in range(1000):
        fname = '../'+dirName+'/'+var+f'_run_break_{i:03}.txt'
        try:
          with open(fname, encoding = 'utf-8') as f:
            count+=1
            df += (np.loadtxt(f).T - df) / count
        except FileNotFoundError: continue
      print(fname,'count',count)
      if count<1: continue
      if 'we_05windowing'in dirName: df*=4
      df /= 256**3
      print('save ', savname)
      np.savetxt(savname, df)
    if 'ESpec' in var: 
      df /= ud**2*diam*td/2**3/np.pi**2
      lels=[10**e for e in range(-7, 5, 1)]
      colors = plt.cm.PiYG(np.linspace(.5,0,len(lels)-1))
    if 'FSpec' in var: 
      df *= -4*np.pi**2/ud**2/diam/case['rho']
      print('rho',case['rho'])
      lels = [-1*10**e for e in range(-1, -9, -1)] + [10**e for e in range(-8, -0)]
      print('lels',lels)
      colors = plt.cm.PiYG(np.linspace(1,0,len(lels)-1))
    if 'PSpec' in var: 
      df *= 4*np.pi/diam**3/td
      lels=[10**e for e in range(-13, -1, 1)]
      colors = plt.cm.PiYG(np.linspace(.5,0,len(lels)-1))
    print('df',np.min(df),np.max(df))
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(lels, ncolors=cmap.N)
    mappable = ax.contourf(X,Y,df, levels=lels, cmap=cmap, norm=norm)
    if 'inf' in dirName or '10' in dirName and 'ESpec' not in var: 
      ax.set_ylabel('$\\frac{\\omega t_d}{2\\pi}$', rotation=0,size=22)#,labelpad=15)
    else: ax.set_yticklabels([])
    if '_02' in dirName: 
      cBar = plt.colorbar(mappable, ax=ax)
      cBar.ax.tick_params(which='both', direction='in', top=True, left=True)
      cBar.ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext())
      if 'ESpec' in var: cBar.set_label("$8\pi^2E(\kappa,\omega)/(u_d^2 d t_d)$")
      if 'FSpec' in var: 
        cBar.set_label("$4\pi^2F(\kappa,\omega)/(u_d^2 d)$")
        lels=np.array(lels)
        cBar.set_ticks(lels[[1,3,5,7,8,10,12,14]])
      if 'PSpec' in var: 
        cBar.set_label("$4\pi A(\kappa,\omega)/(d^3 t_d)$")
    ax.plot([wavNumb[1],wavNumb[-1]],[td/case['tsig'],td/case['tsig']],c='k',lw=.7, ls='dotted')
    ax.plot([diam/case['dHinze'],diam/case['dHinze']],[freq[1],freq[-1]],c='k',lw=.7, ls='dotted')
    grad = ud**.5*td/diam
    ax.plot( [wavNumb[1], wavNumb[-1]], [wavNumb[1]*grad, wavNumb[-1]*grad], c='k', lw=.7)
    if 'inf' not in dirName: ax.plot( wavNumb, ( wavNumb**3 * (8*6/128)**2 / case['rho'] ) **.5 * td/diam, c='k', lw=.7, ls='dashed')
    if 'inf' in dirName: lbl = rf"$\mathrm{{We}}=\infty$"
    else: lbl = rf"$\mathrm{{We}}={case['We']:.2f}$"
    ax.text(.01, .97, lbl, transform=ax.transAxes, ha='left', va='top')
    fname = '../plots/'+var+'LocalFreqVsWav_'+dirName+'.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', pad_inches=0.01, transparent=True)
  return

def phaseFreqVsWav(): 
  for dirName, case in cases.items():
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
    ax.plot([wavNumb[1],wavNumb[-1]],[td/case['tsig'],td/case['tsig']],c=case['col'],alpha=.3)
    mappable = ax.contourf(X, Y, E, levels=10, cmap='Greys')
    cBar = plt.colorbar(mappable, ax=ax)#, orientation='horizontal')
    cBar.set_label('$\\alpha\\mathbf{k\\cdot k}\\hat c^2 /(\\pi d^2\\sigma)$')
    fname = '../plots/phaseFreqVsWav.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsWav(): 
  import matplotlib.colors as mcolors
  from scipy.interpolate import griddata
  #for case in 'we_02  we_05  we_05T100  we_05nTime100  we_05tanhWindow  we_05window  we_05windowing  we_10'.split():
  for dirName, case in cases.items():
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
    colors = plt.cm.PiYG(np.linspace(0.05,.95,len(lels)-1))  # or define your own list of colors
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
    cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}\kappa d\\omega t_\\sigma /(4\\pi^2\\epsilon)$')#, rotation=0, labelpad=15)
    cBar.set_ticks(lels)
    cBar.set_ticklabels(lels)
    cBar.ax.tick_params(which='both', direction='in', top=True, left=True)
    ax.plot([wavNumb[1],wavNumb[-1]],[td/case['tsig'],td/case['tsig']],c=case['col'],alpha=.3)
    ax.set_xlim([wavNumb[1],wavNumb[-1]])
    ax.set_ylim([freq[1],freq[-1]])
    ax.set_xscale('log')
    ax.set_yscale('log')
    fname = '../plots/surPowFreqVsWav_'+dirName+'PG.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#energyVsFreq()
interfaceVsFreq()
#surPowVsFreq()
#energyWindowsVsFreq()
energyFreqVsWav('ESpec')
energyFreqVsWav('FSpec')
energyFreqVsWav('PSpec')
#phaseFreqVsWav()
#surPowFreqVsWav()
