import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
plt.rcParams.update({"text.usetex": True,'font.size' : 14,})
plt.rcParams['pdf.fonttype'] = 42
timestep=.05
nyqTime=2*timestep
nu=1 #kinematic viscosity
We=.05
rho=We
sigma=2*2**.5/3
diam=2*np.pi/3
vis=.219 #from measurements of rmsVel and Re_lambda=58
td=diam**(2/3.)
cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0.8, 0.2, 3)] 
cases = {
    'we_02'   :{'output':'',        'We':.02, 'col':cmap(.2), 'timeLen':120},
    'we_05'   :{'output':'/output', 'We':.05, 'col':cmap(.4), 'timeLen':120},
    'we_08'   :{'output':'/output', 'We':.08, 'col':cmap(.6), 'timeLen':100},
    'we_10'   :{'output':'/output', 'We':.10, 'col':cmap(.8), 'timeLen':90}
       }
for dirName, case in cases.items():
  case['rho'] = case['We']
  case['WeJfm21'] = case['We']*3.63/.1
  case['tsig'] = ( case['rho'] * diam**3 / sigma ) **.5
  case['Oh'] = nu * case['rho']**.5  / ( diam * sigma) ** .5
  print(dirName, 'Oh', case['Oh'])

def areaAndForcingVsTime(): 
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$t$', rotation=0)
  #ax.set_ylabel('$A$', rotation=0)
  #fname = 'FSpec.txt'
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
  #fname = 'statDrops.txt'
  fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_024/statDrops.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    df = np.loadtxt(f)
  #ax.plot(df[:,0],df[:,2])#*1e-4*2*np.sqrt(2)/3)
  ax.plot(df[2:,0],(df[2:,2]-df[1:-1,2])*np.sqrt(2)/3,label='$\\sigma dA$')#/timestep)
  ax.legend()
  fname = 'plots/areaAndForcingVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyVsWaveNumber(): 
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
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
  fname = 'plots/ESpec.pdf'
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
  fname = 'plots/FSpecLinLog.pdf'
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
  fname = 'plots/FKSpec.pdf'
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
  fname = 'plots/surPowVsFreq.pdf'
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
  fname = 'plots/surPowFreqVsFreq.pdf'
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
  fname = 'plots/surPowFreqVsWaveNumberNyq.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def MoIAlignStrainVsTime(): 
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  labels = ['|$S_1d^{2/3}\\epsilon^{-1/3}|$', r'$S_2$', r'$S_3$']
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$|\\mathbf{\\hat s_i\\cdot \hat l_3}|$', rotation=0, labelpad=10)
  ax.set_ylim([0,1])
  ax.set_xlim([-3,0])
  cut=30
  timeLen=100
  allStrain = []
  allS1 = []
  allS3 = []
  allS2 = []
  count=0
  for i in range(1000):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/strainSphereR085.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        strain = np.loadtxt(f)
    except FileNotFoundError:
      continue
    strainLen=np.shape(strain)[0]
    if strainLen>198:continue
    if strainLen<timeLen+cut:continue
    strain = strain[-timeLen:,:]
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/MoInDrops.txt'
    with open(fname, encoding = 'utf-8') as f:
      MoI = np.loadtxt(f)
    MoI = MoI[-timeLen:,:]
    for i in range(len(MoI[:,0])):
      if abs(sum(MoI[i,4:7]**2)) > 2: MoI[i,4] = np.nan 
    inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
    inProdS2 = [ abs(sum(MoI[i,4:7]*strain[i,7:10])) for i in range(len(MoI[:,0])) ]
    inProdS1 = [ abs(sum(MoI[i,4:7]*strain[i,4:7])) for i in range(len(MoI[:,0])) ]
    s=np.stack([inProdS1,inProdS2,inProdS3], axis=-1)
    allStrain.append(s)#, axis=1))
    allS1.append(inProdS1)#, axis=1))
    allS2.append(inProdS2)#, axis=1))
    allS3.append(inProdS3)#, axis=1))
    count+=1
  print('count',count)
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
  time = np.array([(i+1-timeLen)*timestep/td for i in range(timeLen)])
  ax.plot(time, ms1, label=labels[0], color=colors[0])
  ax.plot(time, ms2, label=labels[1], color=colors[1])
  ax.plot(time, ms3, label=labels[2], color=colors[2])
  ax.fill_between(time,ls1,us1,color=colors[0],alpha=0.3,edgecolor='none') 
  ax.fill_between(time,ls2,us2,color=colors[1],alpha=0.3,edgecolor='none') 
  ax.fill_between(time,ls3,us3,color=colors[2],alpha=0.3,edgecolor='none') 
  for i in range(-1):
    ax.plot(time, p50[:,i], label=labels[i-1], color=colors[i-1])
    ax.fill_between(
      time,p10[:,i],p90[:,i],
      #time,m[:,i]-std[:,i],m[:,i]+std[:,i],
      color=colors[i-1],
      alpha=0.3,edgecolor='none') 
  fname = 'plots/MoIAlignStrainsVsTime.pdf'
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
    for region in 'BoxR  SphereR  Drop  Ellipse'.split():
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
        fname = 'dissPlots/' + dirName + '/' + matrix + region + radStr + '.pdf'
        directory = os.path.dirname(fname)
        os.makedirs(directory, exist_ok=True)
        print('savin ',fname)
        fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
        plt.close(fig)

def strainVsTime(): 
  for dirName, case in cases.items():
    print(dirName)
    plotStrainVsTime(dirName,case)
  return

def axesLenVsTime(): 
  cut=20
  for dirName, case in reversed(cases.items()):
    print(dirName)
    fig, ax = plt.subplots(1)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    ax.set_ylabel('$\\frac{L_i}{d}$', rotation=0, size=20, labelpad=10)
    ax.set_xlim([-3,0])
    ax.set_ylim([0,2])
    if False:
      xC=-2.25
      yC=1.75
      h=.15
      w=.5
      x=[]
      y=[]
      for phi in np.linspace(0,2*np.pi,100):
        x.append(xC+w*np.cos(phi))
        y.append(yC+h*np.sin(phi))
      ax.plot(x,y,lw=.5,color='k')
      #ax.plot([xC-.6,xC+.6],[yC,yC], color=colors[2])
      ax.plot([xC+w,xC+.6],[yC,yC], color=colors[2])
      ax.plot([xC-w,xC-.6],[yC,yC], color=colors[2])
      ax.plot([xC,xC],[yC+h,yC+.2], color=colors[0])
      ax.plot([xC,xC],[yC-h,yC-.2], color=colors[0])
      ax.plot(xC,yC,'o',color=colors[1],ms=3)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<case['timeLen']:continue
      #strain = strain[:case['timeLen'], :4]
      strain = strain[cut-case['timeLen']:-1, :4]
      strain=(strain*15/8/np.pi/(256/6)**5)**-.5 / 1.234
      #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count==0: continue
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)
    time = np.array([(i+2+cut-case['timeLen'])*timestep/td for i in range(case['timeLen']-cut-1)])
    #time = np.arange(case['timeLen'])
    for i in range(1, 4):
      ax.plot(time, strainData[31,:,i], color=colors[3-i], alpha=.3)
      ax.plot(time, m[:,i], color=colors[3-i])
      ax.fill_between(time, p10[:,i],p90[:,i],color=colors[3-i],alpha=0.3,edgecolor='none') 
    #fname = 'plots/axesLen_'+dirName+'VsTimeFromStart.pdf'
    fname = 'plots/axesLen_'+dirName+'VsTime.pdf'
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
  fname = 'plots/axesLen_VsStrain.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def aspectRatioVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\frac{2L_3}{L_1+L_2}$', rotation=0, size=20, labelpad=20)
  #ax.set_ylim([0,2])
  ax.set_xlim([-3,0])
  cut=30
  timeLen=100
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/MoInDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        #print(f"File not found: {fname}, skipping.")
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<timeLen+cut:continue
      strain = 2*strain[-timeLen:-1,3]/(strain[-timeLen:-1,2]+strain[-timeLen:-1,1])
      #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5: continue
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)
    time = np.array([(i+2-timeLen)*timestep/td for i in range(timeLen-1)])
    #ax.plot([-3,0], [np.pi/4,np.pi/4], color='k') #above this aspect ratio, the streching is positive
    #ax.plot(time, strainData[31,:,i], color=colors[3-i])#, alpha=.3)
    ax.plot(time, m, color=case['col'])
    ax.fill_between(time,p10,p90,color=case['col'],alpha=0.3,edgecolor='none') 
  fname = 'plots/aspectRatioVsTime.pdf'
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
    time = np.array([(i+2+cut-case['timeLen'])*timestep/td for i in range(case['timeLen']-cut-1)])
    #ax.plot([-3,0], [np.pi/4,np.pi/4], color='k') #above this aspect ratio, the streching is positive
    #ax.plot(time, strainData[31,:,i], color=colors[3-i])#, alpha=.3)
    ax.plot(time, m, color=case['col'], label=f"$\\mathrm{We}={case['WeJfm21']}$")
    #ax.fill_between(time,p10,p90,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.fill_between(time, m-std, m+std, color=case['col'],alpha=0.3,edgecolor='none') 
  ax.legend()
  fname = 'plots/dissVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def areaVsTime(): 
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\frac{A}{\\pi d^2}$', rotation=0, size=22, labelpad=10)
  ax.set_ylim([1,2])
  ax.set_xlim([-3,0])
  cut=20
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if False:
        if strainLen<case['timeLen']:continue
        strain = strain[cut-case['timeLen']:-1,2]/4/np.pi/(256./6*.92)**2
      else:
        if strainLen<2*cut:continue
        strain = strain[cut-case['timeLen']:-1,2]/4/np.pi/(256./6*.92)**2
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5:continue
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)
    time = np.array([(i+2+cut-case['timeLen'])*timestep/td for i in range(case['timeLen']-cut-1)])
    ax.fill_between(time,p10,p90,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.plot(time, m[:], color=case['col'])
  fname = 'plots/areaVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def areaVsTimePad(): 
  normByWeber = False
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70), gridspec_kw={'width_ratios': [1, 1.5**0.5]}, sharey=True)
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
  if normByWeber: 
    ax[0].set_ylabel('$\\left(\\frac{A}{\\pi d^2}-1\\right)\\frac{1}{\\mathrm{We}}$', size=20)
    ax[0].set_ylim([0,.3])
  else: 
    ax[0].set_ylabel('$\\frac{A}{\\pi d^2}$', rotation=0, size=20, labelpad=15)
    ax[0].set_ylim([1,2])
  ax[0].set_xlim([0,2])
  ax[0].set_xlabel('$t/t_d$', rotation=0)
  ax[1].set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax[1].set_xlim([-3,0])
  #ax[1].set_xticks([-2.5,-2,-1.5,-1,-.5,0])
  cut=50
  time = np.array([i*timestep/td for i in range(200)])
  for dirName, case in reversed(cases.items()):
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      if i==0: print(fname)
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen<cut+4:continue
      #strain = strain[:-cut,2]/strain[:-cut,1]**(2/3)/(4*9*np.pi)**(1/3)
      strain = strain[:,2]/strain[:,1]**(2/3)/(4*9*np.pi)**(1/3)
      strain /= strain[1]
      strain = strain[1:1-cut]
      #strain = strain[1:1-cut,2]/4/np.pi/(256./6*.92)**(2)
      #IC 2025 Sep 21: why does the high surface tension case have smaller initial drop?
      #if 'we_02' in dirName: strain/=0.96596215 
      if normByWeber: strain = ( strain - 1 ) / case['WeJfm21']
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((strain,pad_rows))
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5:continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    lbl = rf"$\mathrm{{We}}={case['WeJfm21']:.1f}$"
    ax[0].fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none', label=lbl) 
    ax[0].plot(time, m, color=case['col'])
  time = np.array([(i-199)*timestep/td for i in range(200)])
  for dirName, case in reversed(cases.items()):
    print(dirName)
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<cut+4:continue
      strain = strain[:,2]/strain[:,1]**(2/3)/(4*9*np.pi)**(1/3)
      strain /= strain[1]
      strain = strain[cut:]
      #strain = strain[cut:,2]/4/np.pi/(256./6*.92)**(2)
      #if 'we_02' in dirName: strain/=0.96596215 
      if normByWeber: strain = ( strain - 1 ) / case['WeJfm21']
      pad_rows = np.full(200+cut-strainLen, np.nan)
      strain = np.concatenate((pad_rows,strain))
      count+=1
      allStrain.append(strain)
    print('count',count)
    if count<5:continue
    strainData = np.stack(allStrain, axis=0)
    l = np.nanpercentile(strainData, 25, axis=0)
    u = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    ax[1].fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none', label=f"$\\mathrm{We}={case['WeJfm21']:.1f}$") 
    ax[1].plot(time, m[:], color=case['col'])
  ax[0].legend()
  fname = 'plots/areaVsTimePad.pdf'
  #fig.tight_layout()  
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf')
  return

def survivalVsTime(): 
  fig, ax = plt.subplots(1, 2, figsize=(6.4*100/70,  4.8*50/70))
  fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.2, hspace=0.05)
  for a in ax:
    a.tick_params(which='both', direction='in', top=True, right=True)
    a.set_xlabel('$t/t_d$', rotation=0)
    a.set_xlim([0,6])
  ax[0].set_ylabel('$\\frac{N}{N_0}$', rotation=0, size=20, labelpad=15)
  ax[0].set_ylim([0,1])
  ax[1].set_yscale('log')
  ax[1].set_ylim([1e-3,1])
  time = np.array([i*timestep/td for i in range(200)])
  for dirName, case in (cases.items()):
    allStrain = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/statDrops.txt'
      if i==0: print(fname)
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
    rate = 14.8 * np.exp( -2*3.9/case['WeJfm21'] )
    mask = time < 6
    for a in ax: a.plot(time[mask], np.exp( rate * (1-time[mask]) ), c=case['col'], zorder=-2) 
    lbl=f"$\\mathrm{We}={case['WeJfm21']:.1f}$"
    ax[0].plot(time[mask], m[mask], 'o', c=case['col'], label=lbl, clip_on=False, mfc='None', alpha=.7, zorder=-1) 
    mask = (m > 0) & (time < 6)
    ax[1].plot(time[mask], m[mask], 'o', c=case['col'], label=lbl, clip_on=False, mfc='None', alpha=.7, zorder=-1) 
  ax[1].legend()
  fname = 'plots/survivalVsTime.pdf'
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
    ax.text(xC+w,yC,'$L_3\\mathbf{\\hat l_3}$',color=colors[2], ha='left', va='center')
    ax.text(xC,yC+h+.1,'$L_1\\mathbf{\\hat l_1}$',color=colors[0], ha='center', va='bottom')
    x=[]
    y=[]
    for phi in np.linspace(0,2*np.pi,100):
      x.append(xC+w*np.cos(phi))
      y.append(yC+h*np.sin(phi))
    ax.plot(x,y,lw=.5,color='k')
  if False:
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
    ax.text(xC+w+3*l,yC,'$S_3\\mathbf{\\hat s_3}$',color=colors[2], ha='center', va='center')
    ax.text(xC,yC+w+2*l,'$S_1\\mathbf{\\hat s_1}$',color=colors[0], ha='center', va='center')
  fname = 'plots/strainDiagram.pdf'
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
    #fname = 'plots/back/alignmentVsDissBack'+str(back)+'.pdf'
    fname = f'plots/aspectRatioVsDiss/back{back:03}.pdf'
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
  fig, ax = plt.subplots(1)
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$\\tau/t_d$', rotation=0)
  ax.set_ylabel('$\\langle|\\mathbf{\\hat l_3}(t+\\tau)\\cdot \\mathbf{\\hat s_3}(t)|\\rangle$')#, rotation=0, labelpad=10)
  ax.set_ylim([0.4,.8])
  ax.set_xlim([-1.5,1.5])
  secax = ax.secondary_yaxis(0, transform=ax.transData)
  secax.tick_params(which='both', direction='in', top=True, right=True)
  secax.tick_params(
      labelbottom=False,  # Remove bottom labels
      labeltop=False,     # Remove top labels (secondary axes often default to top)
      labelleft=False,     # Remove top labels (secondary axes often default to top)
      labelright=False,     # Remove top labels (secondary axes often default to top)
      bottom=True,        # Show bottom ticks
      top=True,           # Show top ticks
      direction='inout',  # Extend both in and out
      #length=6            # Length of ticks (adjust as needed)
  )
  cut=20
  for dirName, case in reversed(cases.items()):
    print(dirName)
    timeLen=120
    lagMin=-50
    lagMax=50
    if '_08' in dirName:
      timeLen=80
      lagMin=-30
      lagMax=30
    if '_10' in dirName:
      continue
      timeLen=70
      lagMin=-25
      lagMax=25
    print('avgLen',timeLen-lagMax+lagMin)
    strainLens = []
    listCorr = []
    count=0
    for i in range(1000):
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/strain/SphereR085.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          strain = np.loadtxt(f)
      except FileNotFoundError:
        continue
      strainLen=np.shape(strain)[0]
      strainLens.append(strainLen)
      if strainLen<timeLen+2*cut:continue
      strain = strain[cut:timeLen+cut,:]
      fname='/home/ianto.cannon/drops/boxStrain/'+dirName+case['output']+f'/run_break_{i:03}/MoInDrops.txt'
      with open(fname, encoding = 'utf-8') as f:
        MoI = np.loadtxt(f)
      MoI = MoI[cut:timeLen+cut,:]
      for i in range(len(MoI[:,0])):
        if sum(MoI[i,4:7]**2) > 1.1: print('MoI not unit length',i)
      #inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
      count+=1
      corr = np.zeros(lagMax-lagMin)
      for lag in range(lagMin,lagMax):
        for t in range(-lagMin,timeLen-lagMax):
          #corr[lag-lagMin] += (strain[t,4] * MoI[t+lag,10]) / (timeLen-lagMax+lagMin)
          corr[lag-lagMin] += abs( sum( strain[t,4:7] * MoI[t+lag,10:13] ) ) / (timeLen-lagMax+lagMin)
      listCorr.append(corr)
    print('count',count)
    #strainLens.sort()
    #print('longest strainLen',strainLens[-1])
    #print('5th longest strainLen',strainLens[-5])
    #print('10th longest strainLens',strainLens[-10])
    #if count<10:continue
    ax.plot([case['tsig']/td,case['tsig']/td], [0,1], c=case['col'], alpha=.3)
    allCorr = np.stack(listCorr, axis=0)
    l = np.nanpercentile(allCorr, 45, axis=0)
    u = np.nanpercentile(allCorr, 55, axis=0)
    m = np.nanmean(allCorr, axis=0)
    s = np.nanstd(allCorr, axis=0)/5
    time = np.arange(lagMin,lagMax) * timestep/td 
    #ax.fill_between(time,l,u,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.fill_between(time,m-s,m+s,color=case['col'],alpha=0.3,edgecolor='none') 
    ax.plot(time, m, c=case['col'])
  fname = 'plots/MoIDotStrainVsDelay.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#surPowVsWaveNumber()
#surPowWavVsWaveNumber()
#surPowVsFreg()
#surPowFreqVsFreg()
#energyVsWaveNumber()
#surPowFreqVsWaveNumber()
#areaVsTime()
#areaVsTimePad()
survivalVsTime()
#areaAndForcingVsTime()
#MoIAlignStrainVsTime()
#sanBernado()
#strainVsTime()
#axesLenVsTime()
#axesLenVsStrain()
#aspectRatioVsTime()
#dissVsTime()
#strainDiagram()
#aspecRatioVsDiss()
#MoIDotStrainVsDelay()
