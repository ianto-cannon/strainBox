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
tsig=(rho*diam**3/sigma)**.5
td=diam**(2/3.)

def areaVsTime(): 
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
  fname = 'plots/areaVsTime.pdf'
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
  for i in range(200):
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
  for i in range(200):
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
  for i in range(200):
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
  for i in range(200):
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
  for i in range(200):
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
  for i in range(200):
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
  #colors = ['C0', 'C1', 'C2']
  cmap = plt.get_cmap('plasma')
  colors = [cmap(i) for i in np.linspace(0.8, 0.2, 3)] 
  ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
  ax.set_ylabel('$\\mathbf{s_i\\cdot m_1}$', rotation=0)#, size=18)
  ax.set_ylim([0,1])
  ax.set_xlim([-3,0])
  cut=30
  timeLen=100
  allStrain = []
  allS1 = []
  allS3 = []
  allS2 = []
  count=0
  for i in range(200):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/strainSphereR085.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        strain = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
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
  fname = 'plots/MoIAlignStrainQuartVsTime.pdf'
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

def strainVsTime(): 
  for matrix in [
  'MoInDrops'
  #'strainBoxR021',
  #'strainFarBoxR021',
  #'strainSphereR021',
  #'strainFarSphereR021',
  #'strainFarModesR021',
  #'strainBoxR043',
  #'strainFarBoxR043',
  #'strainSphereR043',
  #'strainFarSphereR043',
  #'strainFarModesR043',
  #'strainBoxR064',
  #'strainFarBoxR064',
  #'strainSphereR064',
  #'strainFarSphereR064',
  #'strainFarModesR064',
  #'strainBoxR085',
  #'strainFarBoxR085',
  #'strainSphereR085',
  #'strainFarSphereR085',
  #'strainFarModesR085',
  #'strainBoxR107',
  #'strainFarBoxR107',
  #'strainSphereR107',
  #'strainFarSphereR107',
  #'strainModes',
  #'strainFarModesR107'
   ]:
    fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$(t-t_b)/t_d$', rotation=0)
    labels = ['$S_1d^{2/3}\\epsilon^{-1/3}$', r'$S_2$', r'$S_3$']
    #colors = ['C0', 'C1', 'C2']
    cmap = plt.get_cmap('plasma')
    colors = [cmap(i) for i in np.linspace(0.8, 0.2, 3)] 
    if matrix=='MoInDrops': 
      ax.set_ylabel('$\\frac{I}{60\\rho\\pi d^5}$', rotation=0, size=20, labelpad=15)
      labels = ['$I_1/I_d$', '$I_2/I_d$', '$I_3/I_d$']
      labels = ['$I_1/I_d$', '$I_2/I_d$', '$I_3/I_d$']
      #ax.set_ylim([0,10])
      xC=-2.25
      yC=8
      #ax.plot([xC-.6,xC+.6],[yC,yC], color=colors[0])
      #ax.plot([xC,xC],[yC-1,yC+1], color=colors[2])
      x=[]
      y=[]
      for phi in np.linspace(0,2*np.pi,100):
        x.append(xC+.5*np.cos(phi))
        y.append(yC+.7*np.sin(phi))
      #ax.plot(x,y, color='k')
    else:
      ax.set_ylabel('$\\frac{Sd^{2/3}}{\\epsilon^{1/3}}$', rotation=0, size=20)
      ax.set_ylim([-.4,.4])
    ax.set_xlim([-3,0])
    cut=30
    timeLen=100
    allStrain = []
    #avgStrain=np.zeros((timeLen,4))
    #avgStrainSq=np.zeros((timeLen,4))
    count=0
    for i in range(200):
      fname = f'/home/ianto.cannon/drops/boxStrain/we_05RBy12/run_break_{i:03}/'+matrix+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          #print('loadin ',fname)
          strain = np.loadtxt(f)
      except FileNotFoundError:
        print(f"File not found: {fname}, skipping.")
        continue
      if matrix=='strainModes':
        strain=strain[4::5,:]
      strainLen=np.shape(strain)[0]
      if strainLen>198:continue
      if strainLen<timeLen+cut:continue
      strain = strain[-timeLen:, :4]
      if matrix=='MoInDrops': strain=(strain*15/8/np.pi/(256/6)**5*3/2)**-.5
      else: strain=strain*(256/3.)**(2/3.)
      #why is factor of 3/2 needed to get I=1 at t=0? It corresponds to R=0.92L/6
      #times=np.shape(strain)[0]
      #if times < timeLen:
      #pad_rows = np.full((200-times, 4), np.nan)
      #strain = np.vstack((pad_rows, strain))
      #else:
      #  strain = strain[-timeLen:]  # trim to last `timeLen` rows
      count+=1
      allStrain.append(strain)
    print('count',count)
    #stDev = np.sqrt(avgStrainSq - avgStrain**2)
    strainData = np.stack(allStrain, axis=0)
    p10 = np.nanpercentile(strainData, 25, axis=0)
    p50 = np.nanpercentile(strainData, 50, axis=0)
    p90 = np.nanpercentile(strainData, 75, axis=0)
    m = np.nanmean(strainData, axis=0)
    std = np.nanstd(strainData, axis=0)
    #time = np.array([(i-199)*timestep/td for i in range(200)])
    time = np.array([(i+1-timeLen)*timestep/td for i in range(timeLen)])
    for i in range(1, 4):
      ax.plot(time[:-1], strainData[1,:-1,i], color=colors[i-1], alpha=.3)
      ax.plot(time, m[:,i], label=labels[i-1], color=colors[i-1])
      ax.fill_between(
        time,p10[:,i],p90[:,i],
        #time,m[:,i]-std[:,i],m[:,i]+std[:,i],
        #avgStrain[:, i] - stDev[:, i],
        #avgStrain[:, i] + stDev[:, i],
        color=colors[i-1],
        alpha=0.3,edgecolor='none') 
    #ax.legend()
    fname = 'plots/Inv' + matrix+'VsTime.pdf'
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
#MoIAlignStrainVsTime()
#sanBernado()
strainVsTime()
