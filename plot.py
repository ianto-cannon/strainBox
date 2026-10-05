import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcdefaults()
#plt.rcParams.update({"text.usetex": True,'font.size' : 14,})
plt.rcParams['pdf.fonttype'] = 42
timestep=.05
nyqTime=2*timestep
nu=1 #kinematic viscosity
diss = 1
sigma=2*2**.5/3
diam=2*np.pi/3
rad=85
vis=.219 #from measurements of rmsVel and Re_lambda=58
td=diam**(2/3.)
#cmap = plt.get_cmap('plasma')
import cmasher as cmr
cmap = cmr.lavender
colors = [cmap(i) for i in [.8,.5,0]]#np.linspace(0.8, 0.2, 3)] 
listy=['solid','dashed','dotted']
cases = {
    'we_02'   :{'output':'',        'We':.02, 'col':cmap(0)},
    'we_05'   :{'output':'/output', 'We':.05, 'col':'k'},#cmap(.25)},
    'we_08'   :{'output':'/output', 'We':.08, 'col':cmap(.5)},
    'we_10'   :{'output':'/output', 'We':.10, 'col':cmap(.75)}
       }
for dirName, case in cases.items():
  case['rho'] = case['We']
  case['WeJfm21'] = case['We']*3.63/.1
  case['tsig'] = ( case['rho'] * diam**3 / sigma ) **.5
  case['Oh'] = nu * case['rho']**.5  / ( diam * sigma) ** .5
  case['dHinze'] = .725 * sigma**(3/5) * case['rho']**(-3/5) * diss**(-2/5)

def areaAndForcingVsTime(): 
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$t$', rotation=0)
  #ax.set_ylabel('$A$', rotation=0)
  fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_024/FSpec.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    df = np.loadtxt(f)
  print('shapedf',np.shape(df))
  #ax.plot(df[1:,0]-df[:-1,0]) to get timestep
  #energy injected to drop at each timestep
  surP = np.array([np.sum(df[i,1:]) for i in range(1,len(df[:,0]))])
  surP*=2 #forgot that half of kx is missing in real FFT
  ax.plot(df[1:,0], surP, label='$\\int \\mathbf{\\hat f_\\sigma\\cdot\\hat u}d^3k$')
  fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_024/statDrops.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    df = np.loadtxt(f)
  #ax.plot(df[:,0],df[:,2])#*1e-4*2*np.sqrt(2)/3)
  ax.plot(df[2:,0],(df[2:,2]-df[1:-1,2])*np.sqrt(2)/3,label='$\\sigma dA$')#/timestep)
  ax.legend()
  fname = 'clusterPlots/areaAndForcingVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyVsWaveNumber(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$kL/2\\pi$', rotation=0)
  ax.set_ylabel('$\\mathbf{\\hat u\\cdot\\hat u^*}$', rotation=0, labelpad=15)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([1e0,1e2])
  #ax.set_ylim([1e-2,1e2])
  spec=np.zeros(85)
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/ESpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        print('loadin ',fname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    print('times',times)
    if times<=100:continue
    dfre = np.fft.rfft(df[50:-50,2:], axis=0) / len(df[50:-50,0])
    count+=1
    spec += (dfre[0,:].real-spec)/count
  wavNumb = np.array([i+1 for i in range(len(spec))])
  ax.plot(wavNumb, spec,'.',c='b')#,clip_on=False)
  ax.plot(wavNumb,-spec,'.',c='r')#,clip_on=False)
  fname = 'clusterPlots/ESpec.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowVsWaveNumber(): 
  fig, ax = plt.subplots(nrows=2,ncols=1, figsize=(6.4/1.5, 4.8))
  for i in range(2):
    ax[i].tick_params(which='both', direction='in', top=True, right=True)
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
  #ax.set_ylabel('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0)#, labelpad=10)
  ax[0].set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[1].set_ylabel('$\\frac{-\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[0].tick_params(labelbottom=False)
  ax[1].set_xlabel('$kd/2\\pi$', rotation=0)
  ax[0].set_ylim([1e-1,2e2])
  ax[1].set_ylim([2e2,1e-1])
  spec=np.zeros(85)
  allSpec=[]
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        #print('loadin ',fname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    #print('times',times)
    if times<=100:continue
    dfre = np.fft.rfft(df[50:-50,2:], axis=0) / len(df[50:-50,0])
    count+=1
    spec += (dfre[0,:].real-spec)/count
    allSpec.append(dfre[0,:].real)
  print('count',count)
  wavNumb = np.array([(i+1)/3 for i in range(len(spec))])
  specData = np.stack(allSpec, axis=0)
  l = np.nanpercentile(specData, 25, axis=0)
  m = np.nanpercentile(specData, 50, axis=0)
  u = np.nanpercentile(specData, 75, axis=0)
  #ax.plot(wavNumb, spec,'.',c='b')#,clip_on=False)
  #ax.plot(wavNumb,-spec,'.',c='r')#,clip_on=False)
  mSpec=-spec
  spec[spec < 0] = np.nan
  mSpec[mSpec < 0] = np.nan
  #ax.plot(wavNumb, spec,c='b')#,clip_on=False)
  #ax.plot(wavNumb,mSpec,c='r')#,clip_on=False)
  ax[0].plot(wavNumb,m,c='b')#,clip_on=False)
  ax[1].plot(wavNumb,-m,c='r')#,clip_on=False)
  #print('l',l)
  #print('u',u)
  #ax.plot(wavNumb,l)#,color='b',alpha=0.3,edgecolor='none') 
  #ax.plot(wavNumb,u)#,color='b',alpha=0.3,edgecolor='none') 
  ax[0].fill_between(wavNumb,l,u,color='b',alpha=0.3,edgecolor='none') 
  ax[1].fill_between(wavNumb,-l,-u,color='r',alpha=0.3,edgecolor='none') 
  #ax.plot(wavNumb,mSpec,c='r')#,clip_on=False)
  for i in range(2): ax[i].set_xlim([wavNumb[0],wavNumb[-1]])
  plt.tight_layout()
  fname = 'clusterPlots/FSpecLinLog.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowWavVsWaveNumber(): 
  fig, ax = plt.subplots()
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xscale('log')
  ax.set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}kd}{2\\pi\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax.set_xlabel('$kd/2\\pi$', rotation=0)
  #ax.spines['bottom'].set_position('zero')
  secax = ax.secondary_xaxis(0, transform=ax.transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(
      labelbottom=False,  # Remove bottom labels
      labeltop=False,     # Remove top labels (secondary axes often default to top)
      bottom=True,        # Show bottom ticks
      top=True,           # Show top ticks
      direction='inout',  # Extend both in and out
      #length=6            # Length of ticks (adjust as needed)
  )
  ax.set_ylim([-40,80])
  spec=np.zeros(85)
  allSpec=[]
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    if times<=100:continue
    dfre = np.fft.rfft(df[50:-50,2:], axis=0) / len(df[50:-50,0])
    count+=1
    spec += (dfre[0,:].real-spec)/count
    allSpec.append(dfre[0,:].real)
  print('count',count)
  wavNumb = np.array([(i+1)/3 for i in range(len(spec))])
  specData = np.stack(allSpec, axis=0)
  l = np.nanpercentile(specData, 25, axis=0)*wavNumb
  m = np.nanpercentile(specData, 50, axis=0)*wavNumb
  u = np.nanpercentile(specData, 75, axis=0)*wavNumb
  mSpec=-spec
  spec[spec < 0] = np.nan
  mSpec[mSpec < 0] = np.nan
  ax.plot(wavNumb,m,c='b')#,clip_on=False)
  ax.fill_between(wavNumb,l,u,color='b',alpha=0.3,edgecolor='none') 
  ax.set_xlim([wavNumb[0],wavNumb[-1]])
  fname = 'clusterPlots/FKSpec.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowVsFreg(): 
  fig, ax = plt.subplots(nrows=2,ncols=1, figsize=(6.4/1.5, 4.8))
  for i in range(2):
    ax[i].tick_params(which='both', direction='in', top=True, right=True)
    ax[i].set_yscale('log')
  #ax.set_ylabel('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0)#, labelpad=10)
  ax[0].set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[1].set_ylabel('$\\frac{-\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[0].tick_params(labelbottom=False)
  ax[1].set_xlabel('$\\omega t_\\sigma /2\\pi$', rotation=0)
  #ax.set_xscale('log')
  #ax.set_xlim([1e0,1e2])
  ax[0].set_ylim([1e-2,1e2])
  ax[1].set_ylim([1e2,1e-2])
  cut=50
  timeLen=40
  spec=np.zeros(int(timeLen/2+1))
  count=0
  allSpec=[]
  for i in range(1000):
    #fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        #print('loadin ',fname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    if times<=timeLen+cut+10:continue
    df=df[-timeLen-cut:-cut,2:]
    #Han window
    #for j in range(timeLen):
    #  df[j,:]=df[j,:]*np.sin( np.pi*j / (timeLen-1) )**2
    dfre = np.fft.rfft(df, axis=0) / timeLen
    count+=1
    spec += (dfre[:,0].real-spec)/count
    allSpec.append(dfre[:,0].real)
  print('count',count)
  specData = np.stack(allSpec, axis=0)
  l = np.nanpercentile(specData, 45, axis=0)
  m = np.nanpercentile(specData, 50, axis=0)
  u = np.nanpercentile(specData, 55, axis=0)
  freq = np.array([i*tsig/nyqTime/int(timeLen/2+1) for i in range(int(timeLen/2+1))])
  #ax.plot(freq, spec,'.',c='b')#,clip_on=False)
  #ax.plot(freq,-spec,'.',c='r')#,clip_on=False)
  ax[0].plot(freq,spec, c='b')#,clip_on=False)
  ax[1].plot(freq,-spec,c='r')#,clip_on=False)
  ax[0].fill_between(freq,0,u,color='b',alpha=0.3,edgecolor='none') 
  ax[1].fill_between(freq,-l,0,color='r',alpha=0.3,edgecolor='none') 
  for i in range(2): ax[i].set_xlim([freq[0],freq[-1]])
  plt.tight_layout()
  fname = 'clusterPlots/surPowVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsFreg(): 
  fig, ax = plt.subplots()
  ax.tick_params(which='both', direction='in', top=True, right=True)
  secax = ax.secondary_xaxis(0, transform=ax.transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(
      labelbottom=False,  # Remove bottom labels
      labeltop=False,     # Remove top labels (secondary axes often default to top)
      bottom=True,        # Show bottom ticks
      top=True,           # Show top ticks
      direction='inout',  # Extend both in and out
      #length=6            # Length of ticks (adjust as needed)
  )
  ax.set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}\\omega t_\\sigma}{2\\pi\\epsilon}$', rotation=0, size=20, labelpad=20)
  ax.set_xlabel('$\\omega t_\\sigma /2\\pi$', rotation=0)
  ax.set_ylim([-1.5,2])
  ax.set_xscale('log')
  cut=50
  timeLen=40
  count=0
  allSpec=[]
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    if times<=timeLen+cut+10:continue
    df=df[-timeLen-cut:-cut,2:]
    dfre = np.fft.rfft(df, axis=0) / timeLen
    count+=1
    allSpec.append(dfre[1:,0].real)
  print('count',count)
  specData = np.stack(allSpec, axis=0)
  freq = np.array([(i+1)*tsig/nyqTime/int(timeLen/2+1) for i in range(int(timeLen/2))])
  l = np.nanpercentile(specData, 45, axis=0)*freq
  m = np.nanpercentile(specData, 50, axis=0)*freq
  u = np.nanpercentile(specData, 55, axis=0)*freq
  ax.plot(freq,m, c='b')
  ax.fill_between(freq,l,u,color='b',alpha=0.3,edgecolor='none') 
  #ax.set_xlim([freq[0],6])
  ax.set_xlim([0.4,6])
  ax.set_xticks([0.4,0.5,0.7,1,2,3,4,5])
  ax.set_xticklabels(['$0.4$','$0.5$','$0.7$','$1$','$2$','$3$','$4$','$5$'])
  fname = 'clusterPlots/surPowFreqVsFreq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsWaveNumber(): 
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$kd/2\\pi$', rotation=0)
  #ax.set_ylabel('$\\omega t_\\sigma /2\\pi$', rotation=0)
  ax.set_ylabel('$\\omega t_d /2\\pi$', rotation=0)
  ax.set_xscale('log')
  #ax.set_yscale('log')
  cut=50
  timeLen=40
  #ax.set_xlim([1,85])
  #ax.set_ylim([1,int(timeLen/2+1)])
  spec=np.zeros((int(timeLen/2+1),85))
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    if times<=timeLen+cut+10:continue
    df=df[-timeLen-cut:-cut,2:]
    dfre = np.fft.rfft(df, axis=0) / timeLen
    count+=1
    spec += (dfre.real-spec)/count
  print('count',count)
  #dfre = np.fft.rfft(df[:,1:], axis=0) / len(df[:,0])
  #freqs = np.rfft.fftfreq(len(df[:,0]))
  #wavNumb = np.array([i/6 for i in range(len(df[0,:]))])
  #wavNumb = range(1,100)
  #ax.pcolormesh(wavNumb, freqs, dfre, shading='auto')
  #mappable = ax.pcolormesh(np.log(np.abs(dfre)**2))
  wavNumb = np.array([(i+1)/6 for i in range(85)])
  freq = np.array([i*tsig/nyqTime/int(timeLen/2+1) for i in range(int(timeLen/2+1))])
  X, Y = np.meshgrid(wavNumb, freq)
  mappable = ax.pcolormesh(X,Y,dfre.real, cmap='bwr_r', vmin=-4, vmax=4)
  cBar = plt.colorbar(mappable, ax=ax)
  cBar.set_label('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0, labelpad=15)
  fname = 'clusterPlots/surPowFreqVsWaveNumberNyq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def MoIAlignStrainVsTime(): 
  fig, ax = plt.subplots(1, 3, figsize=(6.4*100/70,  4.8*50/70), sharey=True)
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    a.set_xlim([-3,0])
  ax[0].set_ylabel('$|\\mathbf{\\hat s_i\\cdot \hat l_3}|$', rotation=0, labelpad=10)
  ax[0].set_ylim([0,1])
  cut=30
  radStr=f'{rad:03}'
  time = np.array([(i-199)*timestep/td for i in range(200)])
  j=0
  for dirName, case in cases.items():
    print(dirName)
    allStrain = []
    allS1 = []
    allS3 = []
    allS2 = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/strain/SphereAvR'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<cut+4:continue
      strain = strain[cut:,:]
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/MoInDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        MoI = np.loadtxt(f)
      MoI = MoI[cut:,:]
      for i in range(len(MoI[:,0])):
        if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
      inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      inProdS2 = [ abs(sum(MoI[i,4:7]*strain[i,7:10])) for i in range(len(MoI[:,0])) ]
      inProdS1 = [ abs(sum(MoI[i,4:7]*strain[i,4:7])) for i in range(len(MoI[:,0])) ]
      pad_rows = np.full(200+cut-strainLen, np.nan)
      inProdS1 = np.concatenate((pad_rows, inProdS1))
      inProdS2 = np.concatenate((pad_rows, inProdS2))
      inProdS3 = np.concatenate((pad_rows, inProdS3))
      s=np.stack([inProdS1,inProdS2,inProdS3], axis=-1)
      allStrain.append(s)#, axis=1))
      allS1.append(inProdS1)#, axis=1))
      allS2.append(inProdS2)#, axis=1))
      allS3.append(inProdS3)#, axis=1))
      count+=1
    print('count',count)
    if count<5: continue
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    s1Data = np.stack(allS1, axis=0)
    s2Data = np.stack(allS2, axis=0)
    s3Data = np.stack(allS3, axis=0)
    ms1 = np.nanmean(s1Data, axis=0)
    ms2 = np.nanmean(s2Data, axis=0)
    ms3 = np.nanmean(s3Data, axis=0)
    ls1 = np.nanpercentile(s1Data, 25, axis=0)
    us1 = np.nanpercentile(s1Data, 75, axis=0)
    ls2 = np.nanpercentile(s2Data, 25, axis=0)
    us2 = np.nanpercentile(s2Data, 75, axis=0)
    ls3 = np.nanpercentile(s3Data, 25, axis=0)
    us3 = np.nanpercentile(s3Data, 75, axis=0)
    ax[j].fill_between(time,ls1,us1,color=colors[0],alpha=0.3,edgecolor='none') 
    ax[j].fill_between(time,ls2,us2,color=colors[1],alpha=0.3,edgecolor='none') 
    ax[j].fill_between(time,ls3,us3,color=colors[2],alpha=0.3,edgecolor='none') 
    ax[j].plot(time, ms1, color=colors[0])# label=labels[0],
    ax[j].plot(time, ms2, color=colors[1])# label=labels[1],
    ax[j].plot(time, ms3, color=colors[2])# label=labels[2],
    for i in range(-1):
      ax[j].fill_between(time,p10[:,i],p90[:,i],color=colors[i-1],alpha=0.3,edgecolor='none') 
      ax[j].plot(time, p50[:,i], color=colors[i-1])
    j+=1
  fname = 'clusterPlots/MoIAlignStrainsVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def dropAxDotExtensionVsT(): 
  fig, ax = plt.subplots(1, 3, figsize=(6.4*100/70,  4.8*50/70), sharey=True)
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    a.set_xlim([-2.5,0])
  ax[0].set_ylabel('$|\\mathbf{\hat l^{(\\alpha)} \\cdot \\hat s^{(3)} }|$')#, rotation=0, labelpad=10)
  ax[0].set_ylim([0,1])
  ax[0].set_yticks([0,0.25,0.5,0.75,1])
  radStr=f'{rad:03}'
  for c, (dirName, case) in enumerate(cases.items()):
    print(dirName)
    for a in range(2):
      if a==0:continue
      allS3 = []
      allS2 = []
      allS1 = []
      count=0
      if a==0:
        time = np.array([i*timestep/td for i in range(200)])
        cut=50
      else:
        time = np.array([(i-199)*timestep/td for i in range(200)])
        cut=30
      for i in range(1000):
        fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/strain/SphereAvR'+radStr+'.txt'
        try:
          with open(fname, encoding = 'utf-8') as f:
            strain = np.loadtxt(f)
        except FileNotFoundError:
          continue
        strainLen=np.shape(strain)[0]
        if strainLen>198:continue
        if strainLen<cut+4:continue
        fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/statDrops.txt'
        with open(fname, encoding = 'utf-8') as f:
          defo = np.loadtxt(f)
        for t in range(len(defo)):
          #Exclude data when drop spans domain
          if defo[t,3] < 0: strain[t,:] = np.nan
        strain = strain[cut:,:]
        fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/MoInDrops.txt'
        with open(fname, encoding = 'utf-8') as f:
          MoI = np.loadtxt(f)
        MoI = MoI[cut:,:]
        for i in range(len(MoI[:,0])):
          if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
        #inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
        #inProdS2 = [ abs(sum(MoI[i,4:7]*strain[i,7:10])) for i in range(len(MoI[:,0])) ]
        #inProdS1 = [ abs(sum(MoI[i,4:7]*strain[i,4:7])) for i in range(len(MoI[:,0])) ]
        inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
        inProdS2 = [ abs(sum(MoI[i,7:10]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
        inProdS1 = [ abs(sum(MoI[i,10:13]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
        pad_rows = np.full(200+cut-strainLen, np.nan)
        if a==0:
          inProdS3 = np.concatenate((inProdS3, pad_rows))
          inProdS2 = np.concatenate((inProdS2, pad_rows))
          inProdS1 = np.concatenate((inProdS1, pad_rows))
        else:
          inProdS3 = np.concatenate((pad_rows, inProdS3))
          inProdS2 = np.concatenate((pad_rows, inProdS2))
          inProdS1 = np.concatenate((pad_rows, inProdS1))
        allS3.append(inProdS3)
        allS2.append(inProdS2)
        allS1.append(inProdS1)
        count+=1
      print('count',count)
      if count<5: continue
      for s, allS in enumerate([allS1, allS2, allS3]):
        data = np.stack(allS, axis=0)
        l = np.nanpercentile(data, 25, axis=0)
        m = np.nanpercentile(data, 50, axis=0)
        u = np.nanpercentile(data, 75, axis=0)
        #m = np.nanmean(data, axis=0)
        #std = np.nanstd(data, axis=0)
        #ax[c-1].fill_between(time,m-.25*std,m+.25*std,color=case['col'],alpha=0.3,edgecolor='none') 
        ax[c-1].fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none') 
        ax[c-1].plot(time, m, ls=listy[2-s], c=case['col'])
        diff=0
        pk=0
        for t in range(len(time)):
          if time[t]<-2.5:continue
          if (m[t]-.5)**2>diff: 
            diff=(m[t]-.5)**2
            pk=t
        print(dirName,c-1,'pk',pk,time[pk],m[pk])
        #ax[c-1].plot(time[pk], m[pk], 'o', c=case['col'], mfc='None')
  fname = 'clusterPlots/dropAxDotExtensionVsT.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def dropAxDotExtWe05VsT(): 
  fig, ax = plt.subplots(1, 1, figsize=(6.4*100/70,  4.8*50/70), sharey=True)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_xlim([-2.5,0])
  ax.set_ylabel('$|\\mathbf{\hat l^{(\\alpha)} \\cdot \\hat s^{(3)} }|$')#, rotation=0, labelpad=10)
  ax.set_ylim([0,1])
  ax.set_yticks([0,0.25,0.5,0.75,1])
  radStr=f'{rad:03}'
  for c, (dirName, case) in enumerate(cases.items()):
    if '05' not in dirName: continue
    print(dirName)
    allS3 = []
    allS2 = []
    allS1 = []
    count=0
    time = np.array([(i-199)*timestep/td for i in range(200)])
    cut=50
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/strain/SphereAvR'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<cut+4:continue
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/statDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        defo = np.loadtxt(f)
      for t in range(len(defo)):
        #Exclude data when drop spans domain
        if defo[t,3] < 0: strain[t,:] = np.nan
      strain = strain[cut:,:]
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/MoInDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        MoI = np.loadtxt(f)
      MoI = MoI[cut:,:]
      for i in range(len(MoI[:,0])):
        if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
      inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      inProdS2 = [ abs(sum(MoI[i,7:10]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      inProdS1 = [ abs(sum(MoI[i,10:13]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      pad_rows = np.full(200+cut-strainLen, np.nan)
      inProdS3 = np.concatenate((pad_rows, inProdS3))
      inProdS2 = np.concatenate((pad_rows, inProdS2))
      inProdS1 = np.concatenate((pad_rows, inProdS1))
      allS3.append(inProdS3)
      allS2.append(inProdS2)
      allS1.append(inProdS1)
      count+=1
    print('count',count)
    if count<5: continue
    for s, allS in enumerate([allS1, allS2, allS3]):
      data = np.stack(allS, axis=0)
      l = np.nanpercentile(data, 25, axis=0)
      m = np.nanpercentile(data, 50, axis=0)
      u = np.nanpercentile(data, 75, axis=0)
      ax.fill_between(time,l,u,color=colors[s],alpha=0.3,edgecolor='none') 
      ax.plot(time, m, ls=listy[2-s], c=colors[s])
      diff=0
      pk=0
      for t in range(len(time)):
        if time[t]<-2.5:continue
        if (m[t]-.5)**2>diff: 
          diff=(m[t]-.5)**2
          pk=t
      print(f"peak {m[pk]:.3g} angle {180*np.acos(m[pk])/np.pi:.3g} time {time[pk]:.3g}")
  fname = 'clusterPlots/dropAxDotExtWe05VsT.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def sanBernado(): 
  fig, ax = plt.subplots(nrows=1,ncols=2)
  h = 2.5
  img = mpimg.imread('SanBernardo_AlonsoCano2.jpg')
  for i in range(2):
    ax[i].tick_params(which='both', direction='in', top=True, right=True)
    ax[i].set_xlabel('$x$', rotation=0)
    ax[i].set_ylabel('$y$', rotation=0)
    ax[i].imshow(img, extent=[0, h*564/800, 0, h])
  ax[0].set_axis_off()
  y=np.array([137,250,390])
  y=(800-y)/800*h
  x=np.array([57,227,397])
  x=x*h/800
  ax[1].plot(x, y,'o',mfc='None')
  C = np.polyfit(x, y, 2)
  xAr=np.linspace(x[0],x[-1],100)
  ax[1].plot(xAr, C[2] + C[1]*xAr + C[0]*xAr**2,'--')
  print('C',C)
  print('v',np.sqrt(-9.81/2/C[0]))
  fname = 'sanBernado.jpg'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='jpg', dpi=600)
  return

def plotStrainVsTime(dirName,case):
  import os
  cut=30
  timeLen=100
  #for matrix in 'ReStress  strain'.split():
  for matrix in 'strain'.split():
    #for region in 'BoxR  FarBoxR  FarModesR  FarSphereR  Modes  SphereR  Drop  Ellipse'.split():
    for region in 'BoxR  SphereR  SphereAvR Drop  Ellipse'.split():
      #if 'ReStress' in matrix and 'Modes' in region:continue
      #for rad in [21,43,64,85,107]:
      for rad in [64,85,107]:
        radStr=f'{rad:03}'
        if 'Drop' in region or 'Ellipse' in region: 
          if rad==85: radStr=''
          if rad!=85: continue
        allStrain = []
        count=0
        lenCount=0
        meanStrainLen=0
        for i in range(1000):
          fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/'+matrix+'/'+region+radStr+'.txt'
          try:
            with open(fname, encoding = 'utf-8') as f:
              #print('loadin ',fname)
              strain = np.loadtxt(f)
          except FileNotFoundError:
            #print(f"File not found: {fname}, skipping.")
            continue
          strainLen=np.shape(strain)[0]
          lenCount+=1
          meanStrainLen+=(strainLen-meanStrainLen)/lenCount
          if strainLen>198:continue
          if strainLen<timeLen+cut:continue
          strain = strain[-timeLen:-1, :4]
          strain=strain*(256/3.)**(2/3.)
          #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
          #times=np.shape(strain)[0]
          #if times < timeLen:
          #pad_rows = np.full((200-times, 4), np.nan)
          #strain = np.vstack((pad_rows, strain))
          #else:
          #  strain = strain[-timeLen:]  # trim to last `timeLen` rows
          count+=1
          allStrain.append(strain)
        print(dirName+matrix+region+radStr,'count',count,'meanStrainLen',meanStrainLen)
        if count==0: return
        strainData = np.stack(allStrain, axis=0)
        fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
        labels = ['$S_1d^{2/3}\\epsilon^{-1/3}$', r'$S_2$', r'$S_3$']
        if 'strain' in matrix: 
          ax.set_ylabel('$\\frac{S^2}{\\epsilon}$', rotation=0, size=20)
        if 'ReStress' in matrix: ax.set_ylabel('$\\frac{Rd^{2/3}}{\\epsilon^{1/3}}$', rotation=0, size=20)
        #ax.set_ylim([-.4,.4])
        ax.set_xlim([-3,0])
        #stDev = np.sqrt(avgStrainSq - avgStrain**2)
        p10 = np.nanpercentile(strainData, 25, axis=0)
        p50 = np.nanpercentile(strainData, 50, axis=0)
        p90 = np.nanpercentile(strainData, 75, axis=0)
        m = np.nanmean(strainData, axis=0)
        std = np.nanstd(strainData, axis=0)
        #time = np.array([(i-199)*timestep/td for i in range(200)])
        time = np.array([(i+2-timeLen)*timestep/td for i in range(timeLen-1)])
        #for i in range(1, 4):
        for i in range(1, 2):
          #ax.plot(time, strainData[1,:,i], color=colors[i-1], alpha=.3)
          ax.plot(time, m[:,i], label=labels[i-1], color=colors[i-1])
          ax.fill_between(
            time,p10[:,i],p90[:,i],
            #time,m[:,i]-std[:,i],m[:,i]+std[:,i],
            #avgStrain[:, i] - stDev[:, i],
            #avgStrain[:, i] + stDev[:, i],
            color=colors[i-1],
            alpha=0.3,edgecolor='none') 
        #ax.legend()
        fname = 'clusterPlots/dissPlots/' + dirName + '/' + matrix + region + radStr + '.pdf'
        directory = os.path.dirname(fname)
        os.makedirs(directory, exist_ok=True)
        print('savin ',fname)
        fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
        plt.close(fig)

def strainVsTime(): 
   for dirName, case in cases.items():
     if '05' in dirName: continue
     print(dirName)
     plotStrainVsTime(dirName,case)
   return

def strainAvVsTime(): 
  fig, ax = plt.subplots(1, 5, figsize=(6.4*100/70,  4.8*50/70), sharey=True)
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for i, a in enumerate(ax):
    a.tick_params(which='both', direction='in', top=True, right=True)
    if i==0:continue
    if i>3:continue
    a.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    a.set_xlim([-2.5,0])
    a.set_xticks([-2,-1,0])
  ax[0].set_ylabel('$S^{(\\alpha)}\\sqrt{\\frac{\\nu}{\\epsilon}}$', labelpad=-5)#size=20)
  ax[0].set_ylim([-.017,.017])
  ax[0].set_xlim([0,2.5])
  ax[0].set_xticks([0,1,2])
  ax[0].set_xlabel('$t/t_d$', rotation=0)
  ax[-1].set_xlabel('$\\mathrm{We}$', rotation=0)
  ax[-1].set_xlim([1.5,4])
  ax[-1].set_xticks([2,3,4])
  m1=np.zeros(3)
  m3=np.zeros(3)
  Wes=[]
  #shapes = ['SphereAv','FarSphereAv']
  for j, (dirName, case) in enumerate(cases.items()):
    print(dirName)
    radStr=f'{rad:03}'
    if '02' in dirName:
      shape = 'FarSphereAv'
      a = 0
      col = 'grey'
      we = 4
      time = np.array([i*timestep/td-3.5 for i in range(200)])
      cut=0
    else:
      shape = 'SphereAv'
      a = j
      col = case['col']
      we = case['WeJfm21']
      time = np.array([(i-199)*timestep/td for i in range(200)])
      cut=20
    allStrain = []
    count=0
    for run in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{run:03}/strain/'+shape+'R'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError: continue
      strainLen=np.shape(strain)[0]
      if '02' not in dirName and strainLen>198: continue
      if '02' not in dirName and strainLen<cut+4: continue
      strain = strain[cut:, :4]
      #I forgot to normalise the FFT in fortran
      if '05' in dirName: strain/=256**4
      pad_rows = np.full([200+cut-strainLen,4], np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      if count==1: print(fname)
      allStrain.append(strain)
      if '02' not in dirName: continue
      if run != 6: continue
      for i in range(1, 4):
        ax[a].plot(time, strain[:,i], ls=listy[3-i], c=col, lw=1.5, alpha=0.3)
    print(f"count {count}")
    if count==0: 
      print("skipping")
      continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[a].fill_between( time, l[:,1], u[:,1], color=col, alpha=0.3, edgecolor='none', label=lbl) 
    for i in range(1, 4):
      if i>1:ax[a].fill_between( time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none') 
      ax[a].plot(time, m[:,i], ls=listy[3-i], c=col, lw=1.5)
    if '02' in dirName: 
      print('mean m1', np.nanmean(m[:,1]))
      print('mean m2', np.nanmean(m[:,2]))
      print('mean m3', np.nanmean(m[:,3]))
      continue
    for t in range(len(time)):
      if time[t]<-2:continue
      if m[t,1]<m1[a-1]: 
        m1[a-1]=m[t,1]
        time1=time[t]
      if m[t,3]>m3[a-1]: 
        m3[a-1]=m[t,3]
        time3=time[t]
    Wes.append(we)
    ax[-1].plot(we, m1[a-1], 'o', c=case['col'], mfc='None')
    ax[-1].plot(we, m3[a-1], 'o', c=case['col'])
    ax[a].plot(time1, m1[a-1], 'o', c=case['col'], mfc='None')
    ax[a].plot(time3, m3[a-1], 'o', c=case['col'])
  we=np.linspace(1,4)
  print('Wes',Wes)
  print('m1',m1)
  print('m3',m3)
  print('Ca',m3*Wes/58/diam**(2/3))
  coeffs = np.polyfit(np.log(Wes), np.log(-m1), 1)  # degree 1 = straight line
  print('coeffs',-np.exp(coeffs[1]),coeffs[0])
  ax[-1].plot(we, -np.exp(coeffs[1]) * we**coeffs[0], c='grey', zorder=-1)
  coeffs = np.polyfit(np.log(Wes), np.log(m3), 1)  # degree 1 = straight line
  print('coeffs',np.exp(coeffs[1]),coeffs[0])
  ax[-1].plot(we, np.exp(coeffs[1]) * we**coeffs[0], c='grey', zorder=-1)
  insax = ax[-1].inset_axes((-.02,.15,.9,.7))
  insax.set_xlim([-3,3])
  insax.set_ylim([-3,3])
  insax.set_aspect('equal', adjustable='box')
  insax.set_axis_off()
  xC=0
  yC=0
  l=.2
  h=2
  w=2
  x=[]
  y=[]
  for phi in np.linspace(0,2*np.pi,100):
    x.append(xC+w*np.cos(phi))
    y.append(yC+h*np.sin(phi))
  for a in np.linspace(0.1,2,8):
    x=[]
    y=[]
    for phi in np.linspace(0,3,500):
      #yP = (2-a)/phi
      yP = a/phi
      if yP**2 + phi**2 > (h/2 + w/2)**2:continue
      x.append(phi)
      y.append(yP)
    for quad in range(4): insax.plot( np.array(x)*(-1)**quad, np.array(y)*(-1)**(quad//2),lw=.5,color='k')
  for s in [0]:
    insax.plot(xC+w+l,yC+s,'>',c='k')
    insax.plot(xC-w-l,yC+s,'<',c='k')
    insax.plot(xC+s,yC-h+l,'^',c='k')
    insax.plot(xC+s,yC+h-l,'v',c='k')
  insax.text(xC+w+2*l,yC,'$\\mathbf{\\hat s^{(3)}}$',color='k', ha='left', va='center')
  insax.text(xC,yC+w+2*l,'$\\mathbf{\\hat s^{(1)}}$',color='k', ha='center', va='center')
  fname = 'clusterPlots/strainAvVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def strainWe05VsTime(): 
  fig, ax = plt.subplots(3, 2, figsize=(6.4*100/70,  2*4.8*50/70), sharey='row', sharex='col')
  fig.subplots_adjust(left=0.1, right=0.97, top=0.97, bottom=0.1, wspace=0.1, hspace=0.07)
  for j, (dirName, case) in enumerate(cases.items()):
    if '05' not in dirName: continue
    print(dirName)
    #for region in 'BoxR  FarBoxR  FarModesR  FarSphereR  Modes  SphereR  Drop  Ellipse'.split():
    for shp in range(3):
     shpLbl = 'box sphere wave'.split()[shp]
     for far in range(2):
      if far: 
        cut=20
        time = np.array([(i-199)*timestep/td for i in range(200)])
      else: 
        cut=0
        time = np.array([i*timestep/td-3 for i in range(200)])
        ax[shp,far].set_ylabel('$\\sqrt{\\frac{\\nu}{\\epsilon}}S_{\\mathrm{'+shpLbl+'}}^{(\\alpha)}$', labelpad=-5)#, rotation=0, size=20)
      shape = 'FarBoxR  FarSphereR FarModesR BoxR SphereAvR Modes'.split()[3*far+shp]
      for r, rad in enumerate([64,85,107]):
        radStr=f'{rad:03}'
        if 'Drop' in shape or 'Ellipse' in shape: 
          if rad==85: radStr=''
          if rad!=85: continue
        ax[shp,far].tick_params(which='both', direction='in', top=True, right=True)
        col = colors[r]
        allStrain = []
        count=0
        for run in range(1000):
          fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{run:03}/strain/'+shape+radStr+'.txt'
          try:
            with open(fname, encoding = 'utf-8') as f:
              strain = np.loadtxt(f)
          except FileNotFoundError: continue
          strainLen=np.shape(strain)[0]
          if strainLen>198: continue
          if strainLen<cut+4: continue
          strain = strain[cut:, :4]
          #I forgot to normalise the FFT in fortran
          if '05' in dirName and 'SphereAv' in shape: strain/=256**4
          pad_rows = np.full([200+cut-strainLen,4], np.nan)
          strain = np.concatenate((pad_rows,strain))
          count+=1
          if count==1: print(fname)
          allStrain.append(strain)
          if '02' not in dirName: continue
          if run != 6: continue
          print('run6',strain[:,3])
          for i in range(1, 4):
            ax[shp,far].plot(time, strain[:,i], ls=listy[3-i], c=col, lw=1.5, alpha=0.3, label=lbl)
        print(f"count {count}")
        if count==0: 
          print(fname,"skipping")
          continue
        strainData = np.stack(allStrain, axis=0)
        l = np.nanpercentile(strainData, 25, axis=0)
        u = np.nanpercentile(strainData, 75, axis=0)
        m = np.nanmean(strainData, axis=0)
        lbl=f"$l={1.5+r/2:.2}d$"
        #print('ax',shp,far,ax)
        for i in range(1, 4):
          if i==1: 
            ax[shp,far].plot(time, m[:,i]*0, c='r')
            ax[shp,far].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none', label=lbl) 
          else: ax[shp,far].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none') 
          ax[shp,far].plot(time, m[:,i], ls=listy[3-i], c=col, lw=1.5)
  if False:
    import matplotlib.patches as patches
    from operator import sub
    figW, figH = ax[0,0].get_figure().get_size_inches()
    _, _, w, h = ax[0,0].get_position().bounds
    disp_ratio = (figH * h) / (figW * w)
    data_ratio = sub(*ax[0,0].get_ylim()) / sub(*ax[0,0].get_xlim())
    rat =disp_ratio / data_ratio
    rect = patches.Rectangle((-2.5-.35,-.033-.35/rat), .7, .7/rat, edgecolor='grey', facecolor='none', lw=2)
    ax[0,0].add_patch(rect)
    circle = patches.Ellipse((-2.5,-.033), .7, .7/rat, linewidth=2, edgecolor='None', facecolor='grey', alpha=0.7)
    ax[1,0].add_patch(circle)
    ax[0,0].plot(-2.5,-.033,'+',c='k')
    ax[1,0].plot(-2.5,-.033,'+',c='k')
  ax[2,1].legend(loc='upper left', ncol=3, columnspacing=.7, handletextpad=0.1, frameon=False)
  ax[0,0].set_ylim([-.02,.02])
  ax[1,0].set_ylim([-.03,.03])
  ax[2,0].set_ylim([-.08,.08])
  ax[2,0].set_xlabel('$t/t_d$', rotation=0)
  ax[2,0].set_xlim([0,3])
  ax[2,1].set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax[2,1].set_xlim([-3,0])
  fname = 'clusterPlots/strainWe05VsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, transparent=True)
  return

def QRVsTime(): 
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70))
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.2, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    a.set_xlim([-3,0])
  ax[0].set_ylabel('$\\frac{Q\\nu}{\\epsilon}$', rotation=0, size=20)
  ax[1].set_ylabel('$\\frac{R\\nu^{3/2}}{\\epsilon^{3/2}}$', rotation=0, size=20)
  time = np.array([(i-199)*timestep/td for i in range(200)])
  cut=50
  for dirName, case in cases.items():
    print(dirName)
    radStr=f'{rad:03}'
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/dVeldx/SphereAvR'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError: continue
      strainLen=np.shape(strain)[0]
      if strainLen>198: continue
      if strainLen<cut+4: continue
      strain = strain[cut:, :3]
      #I forgot to normalise the FFT in fortran
      if '05' in dirName: 
        strain[:,1]/=256**(4*2)
        strain[:,2]/=256**(4*3)
      pad_rows = np.full([200+cut-strainLen, 3], np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      if count==1: print(fname)
      allStrain.append(strain)
    print(f"count {count}")
    if count==0: 
      print("skipping")
      continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[0].fill_between(time, l[:,1], u[:,1], color=case['col'], alpha=0.3, edgecolor='none', label=lbl) 
    for i in range(1, 3):
      if i>1:ax[i-1].fill_between( time, l[:,i], u[:,i], color=case['col'], alpha=0.3, edgecolor='none') 
      ax[i-1].plot(time, m[:,i], color=case['col'])
  ax[0].legend()
  fname = 'clusterPlots/QRVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def vortVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\frac{\\omega^2\\nu}{\\epsilon}$', rotation=0, size=20)
  ax.set_xlim([-3,0])
  time = np.array([(i-199)*timestep/td for i in range(200)])
  cut=50
  for dirName, case in cases.items():
    print(dirName)
    radStr=f'{rad:03}'
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/vortic/SphereAvR'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError: continue
      strainLen=np.shape(strain)[0]
      if strainLen>198: continue
      if strainLen<cut+4: continue
      strain = strain[cut:,1]**2 + strain[cut:,2]**2 + strain[cut:,3]**2
      #I forgot to normalise the FFT in fortran
      if '05' in dirName: 
        strain/=256**(4*2)
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      allStrain.append(strain)
      if count==1: print(fname)
    print(f"count {count}")
    if count==0: 
      print("skipping")
      continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax.fill_between(time, l, u, color=case['col'], alpha=0.3, edgecolor='none', label=lbl) 
    ax.plot(time, m, color=case['col'])
  ax.legend()
  fname = 'clusterPlots/vortVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$E$', rotation=0)#, size=20)
  ax.set_xlim([-3,0])
  radStr=f'{rad:03}'
  time = np.array([(i-199)*timestep/td for i in range(200)])
  cut=50
  for dirName, case in cases.items():
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/ReStress/SphereAvR'+radStr+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError: continue
      strainLen=np.shape(strain)[0]
      if strainLen>198: continue
      if strainLen<cut+4: continue
      strain = 0.5*case['rho']**0*( strain[cut:,1] + strain[cut:,2] + strain[cut:,3] )
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      if count==1: print(fname)
      allStrain.append(strain)
    print(f"count {count}")
  for i in range(2):
    ax[i].tick_params(which='both', direction='in', top=True, right=True)
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
  #ax.set_ylabel('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0)#, labelpad=10)
  ax[0].set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[1].set_ylabel('$\\frac{-\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}}{\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax[0].tick_params(labelbottom=False)
  ax[1].set_xlabel('$kd/2\\pi$', rotation=0)
  ax[0].set_ylim([1e-1,2e2])
  ax[1].set_ylim([2e2,1e-1])
  spec=np.zeros(85)
  allSpec=[]
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        #print('loadin ',fname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    #print('times',times)
    if times<=100:continue
    dfre = np.fft.rfft(df[50:-50,2:], axis=0) / len(df[50:-50,0])
    count+=1
    spec += (dfre[0,:].real-spec)/count
    allSpec.append(dfre[0,:].real)
  print('count',count)
  wavNumb = np.array([(i+1)/3 for i in range(len(spec))])
  specData = np.stack(allSpec, axis=0)
  l = np.nanpercentile(specData, 25, axis=0)
  m = np.nanpercentile(specData, 50, axis=0)
  u = np.nanpercentile(specData, 75, axis=0)
  #ax.plot(wavNumb, spec,'.',c='b')#,clip_on=False)
  #ax.plot(wavNumb,-spec,'.',c='r')#,clip_on=False)
  mSpec=-spec
  spec[spec < 0] = np.nan
  mSpec[mSpec < 0] = np.nan
  #ax.plot(wavNumb, spec,c='b')#,clip_on=False)
  #ax.plot(wavNumb,mSpec,c='r')#,clip_on=False)
  ax[0].plot(wavNumb,m,c='b')#,clip_on=False)
  ax[1].plot(wavNumb,-m,c='r')#,clip_on=False)
  #print('l',l)
  #print('u',u)
  #ax.plot(wavNumb,l)#,color='b',alpha=0.3,edgecolor='none') 
  #ax.plot(wavNumb,u)#,color='b',alpha=0.3,edgecolor='none') 
  ax[0].fill_between(wavNumb,l,u,color='b',alpha=0.3,edgecolor='none') 
  ax[1].fill_between(wavNumb,-l,-u,color='r',alpha=0.3,edgecolor='none') 
  #ax.plot(wavNumb,mSpec,c='r')#,clip_on=False)
  for i in range(2): ax[i].set_xlim([wavNumb[0],wavNumb[-1]])
  plt.tight_layout()
  fname = 'clusterPlots/FSpecLinLog.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowWavVsWaveNumber(): 
  fig, ax = plt.subplots()
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xscale('log')
  ax.set_ylabel('$\\frac{\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}kd}{2\\pi\\epsilon}$', rotation=0, size=20, labelpad=10)
  ax.set_xlabel('$kd/2\\pi$', rotation=0)
  #ax.spines['bottom'].set_position('zero')
  secax = ax.secondary_xaxis(0, transform=ax.transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(
      labelbottom=False,  # Remove bottom labels
      labeltop=False,     # Remove top labels (secondary axes often default to top)
      bottom=True,        # Show bottom ticks
      top=True,           # Show top ticks
      direction='inout',  # Extend both in and out
      #length=6            # Length of ticks (adjust as needed)
  )
  ax.set_ylim([-40,80])
  spec=np.zeros(85)
  allSpec=[]
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    allStrain.append(strain)
    print('count',count)
    if count==0: continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    print('m',m[0,:])
    std = np.nanstd(strainData, axis=0)
    for i in range(1, 4):
      if insax[0]: col = case['col']
      else: col = colors[3-i]
      if i==1:
        lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
        ax[0].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none', label=lbl) 
      else: ax[0].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none') 
      ax[0].plot(time, m[:,i], c=col, ls=listy[i-1])
    Av=0
    count=0
    for t in range(len(time)):
      if time[t]<1:continue
      if time[t]>1.5:continue
      count+=1
      Av += (m[t,1]/m[t,3]-Av)/count
    if insax[0]: insax[0].plot(case['WeJfm21'], Av, 'o', c=case['col'], clip_on=False)
    Wes.append(case['WeJfm21'])
    Avs.append(Av)
    #Wes.append(np.log(case['WeJfm21']))
    #Avs.append(np.log(Av))
  coeffs = np.polyfit(Wes, Avs, 1)  # degree 1 = straight line
  print('coeffs',coeffs)
  we=np.linspace(0,4)
  #insax[0].plot(we, coeffs[1] + we*coeffs[0], c='grey', zorder=-1)
  #print('coeffs',np.exp(coeffs[1]),coeffs[0])
  #we=np.linspace(.7,4)
  #insax[0].plot(we, np.exp(coeffs[1]) * we**coeffs[0], c='grey', zorder=-1)
  cut=30
  time = np.array([(t-199)*timestep/td for t in range(200)])
  for dirName, case in reversed(cases.items()):
    if not insax[0] and '05' not in dirName: continue
    print(dirName)
    allStrain = []
    count=0
    for run in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{run:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<cut+4:continue
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{run:03}/statDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        defo = np.loadtxt(f)
      for t in range(len(defo)):
        #Exclude data when drop spans domain
        if defo[t,3] < 0: strain[t,:] = np.nan
      strain = strain[cut:-1, :4]
      defo = defo[cut:-1, :4]
      MoI1 = strain[:,2] + strain[:,3] - strain[:,1]
      MoI2 = strain[:,1] + strain[:,3] - strain[:,2]
      MoI3 = strain[:,1] + strain[:,2] - strain[:,3]
      strain[:,1] = MoI1 * 2.5 / defo[:,1]
      strain[:,2] = MoI2 * 2.5 / defo[:,1]
      strain[:,3] = MoI3 * 2.5 / defo[:,1]
      strain[:,1:] = 6/256 * strain[:,1:] ** .5 /.918
      #true drop radius is R=0.918L/6
      #strain=(strain*15/8/np.pi/(256/6)**5)**-.5 / 1.234
      pad_rows = np.full([201+cut-strainLen, 4], np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count==0: continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)
    print('m',m[-1,:])
    print('m',m[-1,1]/m[-1,3])
    for i in range(1, 4):
      if insax[0]: col = case['col']
      else: col = colors[3-i]
      if i==1:
        lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
        ax[1].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none', label=lbl) 
      else: ax[1].fill_between(time, l[:,i], u[:,i], color=col, alpha=0.3, edgecolor='none') 
      ax[1].plot(time, m[:,i], c=col, ls=listy[i-1])
  if len(ax)>2: #cylinder
    ax[2].set_xlim([0,1])
    theta=np.linspace(0,np.pi)
    x=.5+.5*np.cos(theta)
    y=.05*np.sin(theta)
    ax[2].plot(x,-y, c='k', clip_on=False)
    ax[2].plot(x,y, c='lightgrey', clip_on=False)
    ax[2].plot(x,np.pi+y, c='k', clip_on=False)
    ax[2].plot(x,np.pi-y, c='k', clip_on=False)
    y=np.linspace(0,np.pi)
    x=.1-.1*np.cos(2*y)
    ax[2].plot(x,y, c='k', clip_on=False)
    ax[2].plot(1-x,y, c='k', clip_on=False)
    ax[2].set_axis_off()
  elif insax[1]: #ellipse
    sz=.3
    insax[1] = ax[1].inset_axes((.35-sz/2,.78-sz/2,sz,sz))
    insax[1].set_aspect('equal', adjustable='box')
    insax[1].set_xlim([-.5,.5])
    insax[1].set_ylim([-.5,.5])
    insax[1].set_axis_off()
    h=.4
    w=h*np.pi
    l=.05
    insax[1].plot([-w,w],[0,0], c='k', clip_on=False)
    insax[1].plot([-w,-w],[-l,l], c='k', clip_on=False)
    insax[1].plot([w,w],[-l,l], c='k', clip_on=False)
    insax[1].plot([0,0],[-h,h], c=colors[0], ls='dotted', clip_on=False)
    insax[1].plot([-l,l],[h,h], c=colors[0], clip_on=False)
    insax[1].plot([-l,l],[-h,-h], c=colors[0], clip_on=False)
    insax[1].text(.8*w,0,'$L^{(3)}$',c='k', ha='right', va='bottom')
    insax[1].text(0.05*w,.8*h,'$L^{(1)}$',c=colors[0], ha='left', va='top')
    x=[]
    y=[]
    for phi in np.linspace(0,2*np.pi,100):
      x.append(w*np.cos(phi))
      y.append(h*np.sin(phi))
    insax[1].plot(x, y, lw=.5, c='k', clip_on=False)
  fname = 'clusterPlots/axesLenVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def axesLenVsStrain(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_ylabel('$\\frac{L_3}{d}$', rotation=0, size=20, labelpad=10)
  ax.set_xlabel('$\\frac{S_3^2d^{2/3}}{\\epsilon^{1/3}}$', rotation=0, size=20, labelpad=20)
  #ax.set_xlim([0,.4])
  #ax.set_ylim([1,2])
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for run in range(80):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{run:03}/strain/SphereR085.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strain=strain[:,3] * (256/3.)**(2/3.)
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{run:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          MoI = np.loadtxt(f)
      except FileNotFoundError:
        continue
      for i in range(len(MoI[:,0])):
        if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
      MoI =(MoI[:,1]*15/8/np.pi/(256/6)**5)**-.5 / 1.234
      count+=1
      for t in range(1,len(strain)-2):
        if t==len(strain)-3:sym='x'
        else:sym='.'
        #ax.plot(strain[t], MoI[t]-MoI[t-1], sym, color=case['col'])
        ax.plot(strain[t], MoI[t], sym, color=case['col'])
      #ax.plot(strain, MoI, color=case['col'])
  fname = 'clusterPlots/axesLen_VsStrain.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def aspectRatioVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\frac{2L_3}{L_1+L_2}$', rotation=0, size=20, labelpad=20)
  ax.set_ylim([1.4,1.95])
  ax.set_xlim([-3,0])
  time = np.array([(i-199)*timestep/td for i in range(200)])
  cut=30
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      #fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/MoInDrops.txt'
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        #print(f"File not found: {fname}, skipping.")
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<4+cut:continue
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/statDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        defo = np.loadtxt(f)
      for t in range(len(defo)):
        #Exclude data when drop spans domain
        if defo[t,3] < 0: strain[t,:] = np.nan
      strain = 2*strain[cut:-1,3]/(strain[cut:-1,2]+strain[cut:-1,1])
      #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
      pad_rows = np.full(201+cut-strainLen, np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5: 
      print('count',count,fname)
      continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    #ax.plot([-3,0], [np.pi/4,np.pi/4], color='k') #above this aspect ratio, the streching is positive
    #ax.plot(time, strainData[31,:,i], color=colors[3-i])#, alpha=.3)
    ax.fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.plot(time, m, color=case['col'])
  fname = 'clusterPlots/aspectRatioVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def dissVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\frac{\\rho S^2 d^3}{\\sigma}$', rotation=0, size=20, labelpad=20)
  #ax.set_ylabel('$\\frac{\\nu S^2 }{\\epsilon}$', rotation=0, size=20, labelpad=20)
  ax.set_xlim([-2,0])
  cut=20
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/strain/BoxR085.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        #print(f"File not found: {fname}, skipping.")
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<case['timeLen']: continue
      strain = strain[cut-case['timeLen']:-1,:4]
      strain = strain[:,1]**2 + strain[:,2]**2 + strain[:,3]**2
      strain=strain*(256/3.)**3*case['rho']/sigma
      #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<2: continue
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)/5
    time = np.array([(t+2+cut-case['timeLen'])*timestep/td for t in range(case['timeLen']-cut-1)])
    #ax.plot([-3,0], [np.pi/4,np.pi/4], color='k') #above this aspect ratio, the streching is positive
    #ax.plot(time, strainData[31,:,i], color=colors[3-i])#, alpha=.3)
    ax.plot(time, m, color=case['col'], label=f"$\\mathrm{We}={case['WeJfm21']}$")
    #ax.fill_between(time,p10,p90,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.fill_between(time, m-std, m+std, color=case['col'],alpha=0.3,edgecolor='none') 
  ax.legend()
  fname = 'clusterPlots/dissVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def areaVsTimePad(): 
  print('areaVsTimePad')
  normByWeber = False
  saveTxtArea = False
  insax = [None]*2
  dTimes=np.array([1,100,165,179])
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70), gridspec_kw={'width_ratios': [1, 1.5]}, sharey=True)
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
  if normByWeber: 
    ax[0].set_ylabel('$\\left(\\frac{A}{\\pi d^2}-1\\right)\\frac{1}{\\mathrm{We}}$', size=20)
    ax[0].set_ylim([0,.3])
  else: 
    ax[0].set_ylabel('$\\frac{A}{\\pi d^2}$', rotation=0, size=20, labelpad=15)
    if insax: ax[0].set_ylim([1,2])
    else: ax[0].set_ylim([1,1.8])
  ax[0].set_xlim([0,2])
  ax[0].set_xlabel('$t/t_d$', rotation=0)
  ax[1].set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax[1].set_xlim([-3,0])
  #ax[1].set_xticks([-2.5,-2,-1.5,-1,-.5,0])
  if insax:
    insax[0] = ax[0].inset_axes((.06,.6,.66,.35))
    insax[1] = ax[1].inset_axes((.04,.6,.44,.35))
    for a in insax:
      a.set_xlabel('$\\mathrm{We}$', rotation=0)
      a.set_xlim([0,4])
      a.set_xticks([0,1,2,3,4])
      a.set_ylim([1,2])
      a.yaxis.set_label_position("right")
      a.yaxis.tick_right()
      a.tick_params(which='both', direction='in', top=True, right=True, left=True, bottom=True)
    insax[0].text(1.25, .5, r"$\frac{\langle A\rangle}{\pi d^2}$", fontsize=20, transform=insax[0].transAxes, ha='center', va='center')
    insax[1].text(1.25, .5, r"$\frac{A_b}{\pi d^2}$", fontsize=20, transform=insax[1].transAxes, ha='center', va='center')
  cut=50
  time = np.array([i*timestep/td for i in range(200)])
  Wes=[]
  AAvs=[]
  errs=[]
  for dirName, case in reversed(cases.items()):
    #if not insax and '05' not in dirName:continue
    allArea = []
    headers = []
    allStrain = []
    count=0
    for i in range(1000):
      #fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      #fname='/home/ianto.cannon/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/statDrops.txt'
      fname='/home/ianto.cannon/fsmBackup/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/statDrops.txt'
      if i==0: print(fname)
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if saveTxtArea:
        area=strain[:,:3]
        pad_area = np.full([200-strainLen, 3], 0.0)
        area = np.concatenate((area,pad_area))
        allArea.append(area)
        headers.extend([f'run{i:03}time', f'run{i:03}volu', f'run{i:03}area' ])
      if strainLen<cut+4:continue
      #strain = strain[:-cut,2]/strain[:-cut,1]**(2/3)/(4*9*np.pi)**(1/3)
      strain = strain[:,2]/strain[:,1]**(2/3)#/(4*9*np.pi)**(1/3)
      #strain /= strain[1]*.971
      strain /= strain[1]*.96 #as the area of the sphere is underestimated
      strain = strain[1:1-cut]
      #strain = strain[1:1-cut,2]/4/np.pi/(256./6*.92)**(2)
      #IC 2025 Sep 21: why does the high surface tension case have smaller initial drop?
      #if 'we_02' in dirName: strain/=0.96596215 
      if normByWeber: strain = ( strain - 1 ) / case['WeJfm21']
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((strain,pad_rows))
      if '05' in dirName and i==146:
        ax[0].plot(time, strain[:], c=case['col'], alpha=0.25, ls='dashed')
        ax[0].plot(time[dTimes], strain[dTimes], 'o', c=case['col'], ls='None', mfc='None')
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5:continue
    strainData = np.stack(allStrain, axis=0)
    if saveTxtArea:
      areaData = np.concatenate(allArea, axis=1)
      header_str = " ".join(f"{h:>15}" for h in headers)
      fname='areaVsTime/'+dirName+'.txt'
      print('save',fname)
      np.savetxt(fname, areaData, fmt='%15.6e', header=header_str, comments='')
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[0].fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none', label=lbl) 
    ax[0].plot(time, m, c=case['col'])
    count=0
    AAv=0
    lAv=0
    uAv=0
    for t in range(len(time)):
      if time[t]<1:continue
      if time[t]>1.5:continue
      count+=1
      AAv += (m[t]-AAv)/count
      lAv += (l[t]-lAv)/count
      uAv += (u[t]-uAv)/count
    insax.plot(case['WeJfm21'], AAv, 'o', c=case['col'], clip_on=False, zorder=4)
    print('AAv',case['WeJfm21'], AAv)
    if insax: insax[0].errorbar(case['WeJfm21'], AAv, yerr=[[AAv-lAv], [uAv-AAv]], fmt='o', c=case['col'], mfc='None', capsize=2, clip_on=False)

    Wes.append(case['WeJfm21'])
    #case['AAv']=AAv
    AAvs.append(AAv)
    errs.append(uAv-lAv)
  we=np.linspace(0,4)
  coeffs = np.polyfit(Wes, AAvs, 1)  # degree 1 = straight line
  print('coeffs',coeffs)
  fit_fn = np.poly1d(coeffs)
  if insax: insax[0].plot(we, fit_fn(we), c='grey', zorder=-1)
  Wes=np.array(Wes)
  AAvs=np.array(AAvs)
  errs=np.array(errs)
  num = np.sum(Wes * (AAvs - 1))# / errs**2)
  den = np.sum(Wes**2)# / errs**2)
  a = num / den
  print('a',a)
  insax[0].plot(we, a*we+1, c='grey', zorder=-1)
  Ab=[]
  Wes=[]
  time = np.array([(i-199)*timestep/td for i in range(200)])
  for dirName, case in reversed(cases.items()):
    #if not insax and '05' not in dirName:continue
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      #fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      fname='/home/ianto.cannon/fsmBackup/drops/boxStrain/'+dirName+'velGAv/output'+f'/run_break_{i:03}/statDrops.txt'
      if i==0: print(fname)
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<cut+4:continue
      strain = strain[:,2]/strain[:,1]**(2/3)/(4*9*np.pi)**(1/3)
      strain /= strain[1]*.96 #set the area of the sphere to 1
      strain = strain[cut:]
      #strain = strain[cut:,2]/4/np.pi/(256./6*.92)**(2)
      #if 'we_02' in dirName: strain/=0.96596215 
      if normByWeber: strain = ( strain - 1 ) / case['WeJfm21']
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((pad_rows,strain))
      if i==146 and '05' in dirName: 
        ax[1].plot(time, strain, c=case['col'], alpha=0.25, ls='dashed')
        print('sl',strainLen)
        print('dt',dTimes)
        ax[1].plot(time[dTimes-strainLen], strain[dTimes-strainLen], 'o', c=case['col'], ls='None', mfc='None')
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5:continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[1].fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none', label=lbl) 
    ax[1].plot(time, m[:], color=case['col'])
    insax.plot(case['WeJfm21'], m[-1], 'o', c=case['col'], clip_on=False, zorder=4)
    Ab.append(m[-1])#-case['AAv'])
    print('Ab',case['WeJfm21'], Ab[-1])
    if insax: insax[1].errorbar(case['WeJfm21'], Ab[-1], yerr=[[m[-1]-l[-1]], [u[-1]-m[-1]]], fmt='o', c=case['col'], mfc='None', capsize=2, clip_on=False, zorder=4)
    Wes.append(case['WeJfm21'])
  coeffs = np.polyfit(Wes, Ab, 1)  # degree 1 = straight line
  print('coeffs',coeffs)
  fit_fn = np.poly1d(coeffs)
  if insax: insax[1].plot(we, fit_fn(we), c='grey', zorder=-1)
  Wes=np.array(Wes)
  Ab=np.array(Ab)
  num = np.sum(Wes * (Ab - 1))# / errs**2)
  den = np.sum(Wes**2)# / errs**2)
  a = num / den
  print('a',a)
  insax[1].plot(we, a*we+1, c='grey', zorder=-1)
  ax[0].legend()
  fname = 'clusterPlots/areaVsTimePad.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf')
  return

def survivalVsTime(): 
  fig, ax = plt.subplots(1, figsize=(6.4*100/70,  4.8*50/70))
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$t/t_d$', rotation=0)
  ax.set_xlim([0,6])
  ax.set_ylabel('$\\frac{N}{N_0}$', rotation=0, size=20, labelpad=15)
  ax.set_ylim([0,1])
  insax = ax.inset_axes((.63,.63,.32,.32))
  #ax.fill_between([0,1],[0,0],[1,1],color='lightgrey',edgecolor='none', zorder=-2) 
  #ax.fill_between([4,6],[0,0],[1,1],color='lightgrey',edgecolor='none', zorder=-2) 
  #ax.plot(np.linspace(1,4), np.sin(np.pi*np.linspace(0,1))**2, ls='dashed', c='grey',zorder=-1)#, alpha=0.5) 
  ax.fill_between(np.linspace(0,6), [0]*50, [1]*50, zorder=-2, color=(.93,.93,.93), edgecolor='None')#, alpha=0.7) 
  ax.fill_between(np.linspace(1,4), [0]*50, np.sin(np.pi*np.linspace(0,1))**2, zorder=-2, color='white', edgecolor='None') 
  insax.set_xlabel('$\\mathrm{We}$', rotation=0)
  insax.set_ylabel('$\\kappa t_d$', rotation=0, labelpad=15)
  insax.set_ylim([0,2])
  insax.set_xlim([0,4])
  insax.set_xticks([0,1,2,3,4])
  insax.tick_params(which='both', direction='in', top=True, right=True, left=True, bottom=True)
  #insax.set_yscale('log')
  we=np.linspace(0,4)
  #insax.plot(1/we,14.8 * np.exp( -2*3.9/we ),c='grey')
  insax.plot(we,14.8 * np.exp( -2*3.9/we ),c='grey')
  rates = [0,0.24,1,1.65]
  time = np.array([i*timestep/td for i in range(200)])
  for i, (dirName, case) in enumerate(cases.items()):
    allStrain = []
    count=0
    for r in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{r:03}/statDrops.txt'
      if r==0: print(fname)
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strain = strain[1:,2]**0
      strainLen=np.shape(strain)[0]
      pad_rows = np.full(200-strainLen, 0)
      strain = np.concatenate((strain,pad_rows))
      count+=1
      allStrain.append(strain)
    print('count',count,case['WeJfm21'])
    strainData = np.stack(allStrain, axis=0)
    m = np.sum(strainData, axis=0)
    m /= m[0]
    mask = time < 6
    ax.plot(time[mask], np.exp( rates[i] * (1-time[mask]) ), c=case['col'], zorder=-2) 
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax.plot(time[mask], m[mask], 'o', c=case['col'], label=lbl, clip_on=False, mfc='None', zorder=-1) 
    #mask = (m > 0) & (time < 6)
    #if case['WeJfm21']>1: 
    insax.plot(case['WeJfm21'], rates[i], 'o', c=case['col'], clip_on=False)
  fname = 'clusterPlots/survivalVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def massVsTime(): 
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70))
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xlabel('$t/t_d$', rotation=0)
    a.set_xlim([0,6])
  ax[0].set_ylabel('$N$', rotation=0, labelpad=5)
  ax[0].set_ylim([0,200])
  ax[1].set_ylabel('$\\frac{d}{L}$', rotation=0, size=20)#, labelpad=15)
  ax[1].set_ylim([.97,1])
  ax[0].fill_between(np.linspace(0,6), [0]*50, [1]*50, zorder=-2, color=(.93,.93,.93), edgecolor='None')#, alpha=0.7) 
  ax[0].fill_between(np.linspace(1,4), [0]*50, np.sin(np.pi*np.linspace(0,1))**2, zorder=-2, color='white', edgecolor='None') 
  insax = [None]*2
  insax[0] = ax[0].inset_axes((.06,.04,.35,.25))
  insax[1] = ax[0].inset_axes((.47, .04,.35,.25))
  insax[0].set_xlabel('$\\mathrm{We}$', rotation=0)
  insax[1].set_xlabel('$\\mathrm{We}$', rotation=0)
  insax[1].text(1.2, .5, '$\\kappa t_d$', ha='left', va='center', transform=insax[1].transAxes)
  insax[1].yaxis.set_label_position("right")
  insax[1].yaxis.tick_right()
  insax[0].xaxis.tick_top()
  insax[1].xaxis.tick_top()
  insax[0].xaxis.set_label_position("top")
  insax[1].xaxis.set_label_position("top")
  insax[0].tick_params(which='both', direction='in', top=True, right=True, left=True, bottom=True)
  insax[1].tick_params(which='both', direction='in', top=True, right=True, left=True, bottom=True)
  insax[1].set_yticks([0,5,10,15])
  insax[0].set_yticks([0,1,2])
  insax[0].yaxis.set_minor_locator(MultipleLocator(1))
  insax[1].set_xticks([0,10,20,30])
  insax[0].set_xticks([0,1,2,3,4])
  we=np.linspace(0,4)
  rate = 14.8 * np.exp( -2*3.9/we)
  mask = rate<2
  insax[0].plot(we[mask], rate[mask], c='grey', clip_on=False, zorder=4)
  insax[1].plot(we, rate, c='grey',clip_on=False, zorder=4)
  we=np.linspace(4,30)
  insax[1].plot(we,14.8 * np.exp( -2*3.9/we ),c='grey',clip_on=False, zorder=4)
  #N=256**3 Re=58 measurements from velamartin22memoryless:
  wetom1=np.array([0.68661,0.54859,0.45733,0.39161,0.34307,0.30281,0.27486,0.22637,0.21239,0.18122,0.15987,0.13517])
  kappa=[0.096773,0.193114,0.453956,0.612631,1.184390,1.517420,1.857449,2.747264,3.171361,2.763910,2.949794,4.568208]
  insax[0].plot(1/wetom1[:-5], kappa[:-5], 'v', c='c', clip_on=False, mfc='None')
  insax[1].plot(1/wetom1, kappa, 'v', c='c', ms=3, mfc='None')
  insax[1].plot([0,30], [14.8,14.8], c='k', zorder=4, ls='dashed')
  rates = [0,0.24,1,1.65]
  insax[1].plot([0,4], [2,2], c='k',clip_on=False, zorder=4, lw=.8)
  insax[1].plot([4,4], [0,2], c='k',clip_on=False, zorder=4, lw=.8)
  insax[1].plot([0,-5], [2,15], c='k',clip_on=False, zorder=4, lw=.8, ls='dotted')
  insax[1].plot([0,-5], [0,0], c='k',clip_on=False, zorder=4, lw=.8, ls='dotted')
  insax[0].set_ylim([0,2])
  insax[0].set_xlim([0,4])
  insax[1].set_xlim([0,30])
  insax[1].set_ylim([0,15])
  time = np.array([i*timestep/td for i in range(200)])
  for i, (dirName, case) in enumerate(cases.items()):
    if '05' not in dirName: continue
    allStrain = []
    allMass = []
    count=0
    for r in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{r:03}/statDrops.txt'
      if r==0: print(fname)
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      mass = strain[1:,1]**(1/3)#/256
      mass /= mass[0]
      strain = strain[1:,2]**0
      strainLen=np.shape(strain)[0]
      pad_rows = np.full(200-strainLen, 0)
      strain = np.concatenate((strain,pad_rows))
      pad_rows = np.full(200-strainLen, np.nan)
      mass = np.concatenate((mass,pad_rows))
      count+=1
      allStrain.append(strain)
      allMass.append(mass)
    print('count',count,case['WeJfm21'])
    strainData = np.stack(allStrain, axis=0)
    m = np.sum(strainData, axis=0)
    mask = time < 6
    print('m',time[mask], m[mask]) 
    ax[0].plot(time[mask], m[mask], 'o', c='grey', clip_on=False, zorder=-2, mfc='None') 
    ax[0].plot(time[mask], 200*np.exp( rates[i] * (1-time[mask]) ), c='k', zorder=-1) 
    massData = np.stack(allMass, axis=0)
    l = np.nanpercentile(massData, 25, axis=0)
    u = np.nanpercentile(massData, 75, axis=0)
    m = np.nanmean(massData, axis=0)
    ax[1].fill_between(time[mask], l[mask], u[mask], color='k', alpha=0.3, edgecolor='none') 
    ax[1].plot(time[mask], m[mask], c='k', clip_on=False, zorder=-1) 
    insax[0].plot(case['WeJfm21'], rates[i], 'o', c='k', clip_on=False, zorder=4, mfc='None')
    insax[1].plot(case['WeJfm21'], rates[i], 'o', c='k', clip_on=False, zorder=4, mfc='None', ms=3)
  fname = 'clusterPlots/massVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def strainDiagram():
  fig, ax = plt.subplots(figsize=(4.8, 4.8))
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlim([-3,3])
  ax.set_ylim([-3,3])
  ax.set_axis_off()
  xC=0
  yC=0
  h=.4
  w=1.4
  l=.2
  if True:
    ax.plot([xC-w,xC+w],[yC,yC], color=colors[2])
    ax.plot([xC,xC],[yC-h,yC+h], color=colors[0])
    #ax.plot([xC+w,xC+w+l],[yC,yC], color=colors[2])
    #ax.plot([xC-w,xC-w-l],[yC,yC], color=colors[2])
    #ax.plot([xC,xC],[yC+h,yC+h+l], color=colors[0])
    #ax.plot([xC,xC],[yC-h,yC-h-l], color=colors[0])
    ax.plot(xC,yC,'o',color=colors[1],ms=3)
    ax.text(xC+w-.1,yC+.15,'$L^{(3)}\\mathbf{\\hat l^{(3)}}$',color=colors[2], ha='center', va='center')
    ax.text(xC,yC+h+.1,'$L^{(1)}\\mathbf{\\hat l^{(1)}}$',color=colors[0], ha='center', va='bottom')
    x=[]
    y=[]
    for phi in np.linspace(0,2*np.pi,100):
      x.append(xC+w*np.cos(phi))
      y.append(yC+h*np.sin(phi))
    ax.plot(x,y,lw=.5,color='k')
  if True:
    x=[]
    y=[]
    h=2
    w=2
    for phi in np.linspace(0,2*np.pi,100):
      x.append(xC+w*np.cos(phi))
      y.append(yC+h*np.sin(phi))
    ax.plot(x,y,lw=.5,color='k',ls='--')
    for s in [-.3,0,.3]:
      ax.plot([xC+w-l,xC+w+l],[yC+s,yC+s], color=colors[2])
      ax.plot([xC-w-l,xC-w+l],[yC+s,yC+s], color=colors[2])
      ax.plot([xC+s,xC+s],[yC+h-l,yC+h+l], color=colors[0])
      ax.plot([xC+s,xC+s],[yC-h-l,yC-h+l], color=colors[0])
      ax.plot(xC+w+l,yC+s,'>',color=colors[2])
      ax.plot(xC-w-l,yC+s,'<',color=colors[2])
      ax.plot(xC+s,yC-h+l,'^',color=colors[0])
      ax.plot(xC+s,yC+h-l,'v',color=colors[0])
    ax.text(xC+w+2*l,yC,'$S^{(3)}\\mathbf{\\hat s^{(3)}}$',color=colors[2], ha='left', va='center')
    ax.text(xC,yC+w+2*l,'$S^{(1)}\\mathbf{\\hat s^{(1)}}$',color=colors[0], ha='center', va='center')
  fname = 'clusterPlots/strainDiagram.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def aspecRatioVsDiss(): 
  fig = []
  ax = []
  backLen=50
  for back in range(backLen):
    figTemp, axTemp = plt.subplots(1)
    axTemp.tick_params(which='both', direction='in', top=True, right=True)
    axTemp.set_ylabel('$\\frac{2L_3}{L_1+L_2}$', rotation=0, size=20, labelpad=20)
    axTemp.set_xlabel('$\\frac{S^2d^{4/3}}{\\epsilon^{2/3}}$', rotation=0, size=20, labelpad=20)
    #axTemp.set_xlabel('$|\\mathbf{\\hat s_3\\cdot \hat l_3}|$', rotation=0, labelpad=10)
    #axTemp.set_xlim([0,1])
    axTemp.set_xlim([0,.4])
    axTemp.set_ylim([1,2])
    fig.append(figTemp)
    ax.append(axTemp)
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/strain/SphereR085.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      strain = strain[:,1]**2 + strain[:,2]**2 + strain[:,3]**2
      strain=strain*(256/3.)**(4/3.)
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          MoI = np.loadtxt(f)
      except FileNotFoundError:
        continue
      for i in range(len(MoI[:,0])):
        if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
      #inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      MoI = 2*MoI[:,3]/(MoI[:,2]+MoI[:,1])
      count+=1
      if strainLen<backLen+30: continue
      if strainLen==200: continue
      for back in range(backLen):
        if back==0:sym='x'
        else:sym='.'
        try: 
          #ax[back].plot(inProdS3[-2-back], MoI[-2-back], sym, color=case['col'])
          #ax[back].plot(inProdS3[-5-back:-1-back], MoI[-5-back:-1-back], color=case['col'], alpha=.5)
          ax[back].plot(strain[-2-back], MoI[-2-back], sym, color=case['col'])
          ax[back].plot(strain[-5-back:-1-back], MoI[-5-back:-1-back], color=case['col'], alpha=.5)
        except: IndexError
      allStrain.append(strain)
  for back in range(backLen):
    #fname = 'clusterPlots/back/alignmentVsDissBack'+str(back)+'.pdf'
    fname = f'clusterPlots/aspectRatioVsDiss/back{back:03}.pdf'
    print('savin ',fname)
    fig[back].savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def correlation(x, y):
  #n = len(x)
  #nfft = 2 ** (int(np.ceil(np.log2(2 * n - 1))))
  fx = np.fft.fft(x, n=4*len(x))
  fy = np.fft.fft(y, n=4*len(x))
  cross_corr = np.fft.ifft( fx * np.conj(fy) ).real
  return cross_corr#[:n]#/(n-np.arange(n))

def MoIDotStrainVsDelay(): 
  fig, ax = plt.subplots(1, 3, figsize=(6.4*100/70,  4.8*50/70), gridspec_kw={'width_ratios': [1,1,2]})
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.32, hspace=0.05)
  for i, a in enumerate(ax):
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_ylim([0,1])
    #a.set_yticks([0,.25,.5,.75,1])
    a.set_yticks([0,.5,1])
    a.yaxis.set_minor_locator(MultipleLocator(0.1))
    a.set_xlabel('$\\tau/t_d$')
  secax = a.secondary_yaxis(0, transform=ax[2].transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False, bottom=True, top=True, direction='inout')
  secax.set_yticks([0,.5,1])
  secax.yaxis.set_minor_locator(MultipleLocator(0.1))
  ax[0].set_xlim([0,2])
  ax[1].set_xlim([0,2])
  ax[2].set_xlim([-2,2])
  ax[0].set_ylabel('$\\langle|\\mathbf{\\hat s^{(3)}}(t+\\tau)\\cdot \\mathbf{\\hat s^{(3)}}(t)|\\rangle$',labelpad=3)
  ax[1].set_ylabel('$\\langle|\\mathbf{\\hat l^{(3)}}(t+\\tau)\\cdot \\mathbf{\\hat l^{(3)}}(t)|\\rangle$',labelpad=3)
  ax[2].set_ylabel('$\\langle|\\mathbf{\\hat l^{(3)}}(t+\\tau)\\cdot \\mathbf{\\hat s^{(3)}}(t)|\\rangle$',labelpad=3)
  insax = False#[None]*4
  if insax:
    insax[0] = ax[0].inset_axes((.06,.03,.6,.25))
    insax[1] = ax[1].inset_axes((.06,.03,.6,.25))
    insax[2] = ax[2].inset_axes((.03,.03,.3,.25))
    insax[3] = ax[2].inset_axes((.03,.03,.3,.25))
    #insax[3] = ax[2].inset_axes((.53,.03,.3,.25))
    insax[0].set_axis_off()
    insax[2].set_axis_off()
    for i, a in enumerate(insax):
      a.xaxis.set_label_position("top")
      a.xaxis.tick_top()
      a.set_xlabel('$\\mathrm{We}$', rotation=0)
      a.set_xlim([0,4])
      a.set_xticks([0,1,2,3,4])
      a.yaxis.tick_right()
      if i==3: 
        a.set_ylabel('$\\frac{\\tau_p}{t_\\sigma}$', rotation=0, size=20)
        a.set_ylim([-.5,.5])
        #a.set_yticks([-.3,0,.3])
        a.yaxis.set_label_coords(x=1.42, y=0.25)
        #a.yaxis.set_minor_locator(MultipleLocator(0.1))
      else:
        a.set_ylabel('$\\frac{\\tau_w}{t_\\sigma}$', rotation=0, size=20)
        a.set_ylim([0,2])
        a.set_yticks([0,1,2])
        a.yaxis.set_label_coords(x=1.3, y=0.25)
        a.yaxis.set_minor_locator(MultipleLocator(0.5))
      a.tick_params(which='both', direction='in', top=True, right=True, left=True, bottom=True)
  cutSta=5
  cutEnd=5
  delayMax=200-cutSta-cutEnd-1
  time = np.array([t*timestep for t in range(-delayMax, delayMax+1)])
  for dirName, case in (cases.items()):
    for far in ['Far','']:
      if not insax and '05' not in dirName: continue
      print(dirName)
      #datname = 'clusterPlots/MoIDotStrainMedianVsDelay_'+dirName+'.txt'
      datname = 'clusterPlots/MoIDotStrainVsDelay_'+dirName+far+'.txt'
      print('load',datname)
      try:
        with open(datname, encoding = 'utf-8') as f:
          m = np.loadtxt(f)
      except FileNotFoundError:
        print('not found')
        m = np.full([2*delayMax+1, 3], np.nan)
        listCorr0 = [[] for _ in range(2*delayMax+1)]
        listCorr1 = [[] for _ in range(2*delayMax+1)]
        listCorr2 = [[] for _ in range(2*delayMax+1)]
        count=0
        for i in range(1000):
          fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/strain/'+far+'SphereR085.txt'
          try:
            with open(fname, encoding = 'utf-8') as f:
              strain = np.loadtxt(f)
          except FileNotFoundError:
            continue
          strainLen = np.shape(strain)[0] - cutSta - cutEnd
          if strainLen < 4: continue
          strain = strain[cutSta:-cutEnd,:]
          fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/MoInDrops.txt'
          with open(fname, encoding = 'utf-8') as f:
            MoI = np.loadtxt(f)
          MoI = MoI[cutSta:-cutEnd,:]
          fname='/home/ianto.cannon/drops/boxStrain/'+dirName+f'velGAv/output/run_break_{i:03}/statDrops.txt'
          with open(fname, encoding = 'utf-8') as f:
            defo = np.loadtxt(f)
          defo = defo[cutSta:-cutEnd,3]
          for t in range(len(defo)):
            #Exclude data when drop spans domain
            if defo[t] < 0: 
              MoI[t,10:13] = np.nan
              strain[t,4:7] = np.nan
          count+=1
          for lag in range(1-strainLen,strainLen):
            t_start = max(0, -lag)
            t_end = min(strainLen, strainLen-lag)
            for t in range(t_start, t_end):
              align = abs( sum(strain[t,4:7] * MoI[t+lag,10:13] ) )
              if align!=align: continue
              listCorr2[lag+delayMax].append(align)
              align = abs( sum(strain[t,4:7] * strain[t+lag,4:7] ) )
              listCorr0[lag+delayMax].append(align)
              align = abs( sum(MoI[t,10:13] * MoI[t+lag,10:13] ) )
              listCorr1[lag+delayMax].append(align)
        print('count',count)
        count=0
        for i, lst0 in enumerate(listCorr0):
          if len(lst0) < 10*delayMax: continue    
          m[i,0] = np.mean(lst0)
          m[i,1] = np.mean(listCorr1[i])
          m[i,2] = np.mean(listCorr2[i])
          #m[i,0] = np.nanpercentile(lst0, 50)
          #m[i,1] = np.nanpercentile(listCorr1[i],50)
          #m[i,2] = np.nanpercentile(listCorr2[i],50)
          count+=1
        if count==0:continue
        np.savetxt(datname, m, fmt='%.18e', header=f'Correlation data for {dirName}\nColumns: corr0, corr1, corr2')
        print(f'Saved {dirName} to {datname}')
      #ax[1].plot(time/case['tsig'], m[:,1], c=case['col'])
      #ax[2].plot(time/td, m[:,2], c=case['col'])
      for a in range(3):
        #ax[a].plot(time/td, m[:,a], '.', c=case['col'], ls='None')
        if 'Far' in far: 
          col = colors[0]
          if a>0: continue
        else: col = case['col']
        ax[a].plot(time/td, m[:,a], c=col)
        tPeak = np.nanargmax(m[:,a])
        #if a!=1:continue
        tCorr = 0
        for t in range(delayMax,0,-1):
          if m[t,a]<0.5: break
          if m[t,a]!=m[t,a]: break
          #print('t',t,m[t,a])
          tCorr += (m[t,a] - 0.5)*timestep/(m[tPeak,a]-0.5)
        for t in range(delayMax+1,len(time)):
          if m[t,a]<0.5: break
          if m[t,a]!=m[t,a]: break
          #print('t',t,m[t,a])
          tCorr += (m[t,a] - 0.5)*timestep/(m[tPeak,a]-0.5)
        print('a',a,'tCorr',tCorr/td)
        if insax:insax[a].plot(case['WeJfm21'], tCorr/case['tsig'], 'o', c=case['col'], clip_on=False)
      #insax[3].errorbar(case['WeJfm21'], time[tPeak]/case['tsig'], yerr=(time[tPeak+1]-time[tPeak-1])/case['tsig'], fmt='o', c=case['col'], mfc='None', capsize=2, clip_on=False)
      if insax: insax[3].plot(case['WeJfm21'], time[tPeak]/case['tsig'], 'o', c=case['col'], clip_on=False)
      print('delay', time[tPeak]/td)
      #insax[3].errorbar(case['WeJfm21'], time[tPeak]/case['tsig'], yerr=(time[tPeak+1]-time[tPeak-1])/case['tsig'], fmt='o', c=case['col'], mfc='None', capsize=2, clip_on=False)
  fname = 'clusterPlots/MoIDotStrainVsDelay.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#surPowVsWaveNumber()
#surPowWavVsWaveNumber()
#surPowVsFreg()
#surPowFreqVsFreg()
#energyVsWaveNumber()
#surPowFreqVsWaveNumber()
areaVsTimePad()
#survivalVsTime()
#massVsTime()
#areaAndForcingVsTime()
#MoIAlignStrainVsTime()
#dropAxDotExtensionVsT()
#dropAxDotExtWe05VsT()
#sanBernado()
#strainVsTime()
#strainAvVsTime()
#strainWe05VsTime()
#QRVsTime()
#vortVsTime()
#energyVsTime()
#axesLenVsTime()
#axesLenVsStrain()
#aspectRatioVsTime()
#dissVsTime()
#strainDiagram()
#aspecRatioVsDiss()
#MoIDotStrainVsDelay()
