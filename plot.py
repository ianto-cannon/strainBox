import numpy as np
import matplotlib.pyplot as plt

def areaVsTime(): 
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$t$', rotation=0)
  #ax.set_ylabel('$A$', rotation=0)
  fname = 'FSpec.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    df = np.loadtxt(f)
  print('shapedf',np.shape(df))
  #ax.plot(df[1:,0]-df[:-1,0]) to get timestep
  timestep=.05
  #energy injected to drop at each timestep
  surP = np.array([np.sum(df[i,1:]) for i in range(1,len(df[:,0]))])
  surP*=2 #forgot that half of kx is missing in real FFT
  ax.plot(df[1:,0], surP, label='$\\int \\mathbf{\\hat f_\\sigma\\cdot\\hat u}d^3k$')
  fname = 'statDrops.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    df = np.loadtxt(f)
  #ax.plot(df[:,0],df[:,2])#*1e-4*2*np.sqrt(2)/3)
  ax.plot(df[2:,0],(df[2:,2]-df[1:-1,2])*np.sqrt(2)/3,label='$\\sigma dA$')#/timestep)
  ax.legend()
  fname = 'areaVsTime.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def energyVsWaveNumber(): 
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
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
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$kL/2\\pi$', rotation=0)
  ax.set_ylabel('$\\mathbf{\\hat f_\\sigma\\cdot\\hat u^*}$', rotation=0)#, labelpad=10)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([1e0,1e2])
  #ax.set_ylim([1e-2,1e2])
  spec=np.zeros(85)
  count=0
  for i in range(200):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
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
  fname = 'plots/FSpec.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def surPowFreqVsWaveNumber(): 
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$k$', rotation=0)
  ax.set_ylabel('$\\omega$', rotation=0)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([1,129])
  ax.set_ylim([1,98])
  spec=np.zeros((51,85))
  count=0
  for i in range(200):
    fname = f'/home/ianto.cannon/drops/boxStrain/we_05BoxLBy2p5/run_break_{i:03}/FSpec.txt'
    try:
      with open(fname, encoding = 'utf-8') as f:
        print('loadin ',fname)
        df = np.loadtxt(f)
    except FileNotFoundError:
      print(f"File not found: {fname}, skipping.")
      continue
    times=np.shape(df)[0]
    print('times',times)
    if times<=150:continue
    dfre = np.fft.rfft(df[-150:-50,2:], axis=0) / 100 #len(df[50:-50,0])
    count+=1
    spec += (dfre.real-spec)/count
  #dfre = np.fft.rfft(df[:,1:], axis=0) / len(df[:,0])
  #freqs = np.rfft.fftfreq(len(df[:,0]))
  #wavNumb = np.array([i/6 for i in range(len(df[0,:]))])
  #wavNumb = range(1,100)
  #ax.pcolormesh(wavNumb, freqs, dfre, shading='auto')
  #mappable = ax.pcolormesh(np.log(np.abs(dfre)**2))
  mappable = ax.pcolormesh(dfre.real, cmap='bwr', vmin=-4, vmax=4)
  plt.colorbar(mappable, ax=ax, label='Magnitude')
  fname = 'plots/surPowFreqVsWaveNumber.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def MoIAlignStrain(): 
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
  fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.set_xlabel('$t$', rotation=0)
  #ax.set_ylabel('$A$', rotation=0)
  fname = 'MoInDrops.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    MoI = np.loadtxt(f)
  fname = 'strainBox.txt'
  with open(fname, encoding = 'utf-8') as f:
    print('loadin ',fname)
    strain = np.loadtxt(f)
  inProdS3 = [ abs(sum(MoI[i,4:7]*strain[i,10:13])) for i in range(len(MoI[:,0])) ]
  inProdS2 = [ abs(sum(MoI[i,4:7]*strain[i,7:10])) for i in range(len(MoI[:,0])) ]
  inProdS1 = [ abs(sum(MoI[i,4:7]*strain[i,4:7])) for i in range(len(MoI[:,0])) ]
  ax.plot(MoI[:,0], inProdS1, label='$\\mathbf{s_1\\cdot m_1}$')
  ax.plot(MoI[:,0], inProdS2, label='$\\mathbf{s_2\\cdot m_1}$')
  ax.plot(MoI[:,0], inProdS3, label='$\\mathbf{s_3\\cdot m_1}$')
  ax.legend()
  fname = 'MoIAlignStrain.pdf'
  print('savin ',fname)
  fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

def sanBernado(): 
  import matplotlib.image as mpimg
  plt.rcParams.update({"text.usetex": True})
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
  plt.rcdefaults()
  plt.rcParams.update({"text.usetex": True,'font.size' : 12,})
  for matrix in ['MoInDrops','strainModes','strainBox','strainFarModes','strainFarBox','strainLowPassBox','strainSphere','strainFarSphere']:
    fig, ax = plt.subplots(1)#, figsize=[columnWid, .6*columnWid])
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlabel('$t$', rotation=0)
    timelen=100
    avgStrain=np.zeros((timelen,4))
    avgStrainSq=np.zeros((timelen,4))
    count=0
    for i in range(200):
      #fname = f'/Users/iantocannon/Desktop/fsm/drops/boxStrain/we_05/run_break_{i:03}/'+matrix+'.txt'
      fname = f'/home/ianto.cannon/drops/boxStrain/we_05/run_break_{i:03}/'+matrix+'.txt'
      try:
        with open(fname, encoding = 'utf-8') as f:
          print('loadin ',fname)
          strain = np.loadtxt(f)
      except FileNotFoundError:
        print(f"File not found: {fname}, skipping.")
        continue
      times=np.shape(strain)[0]
      #times=len(strain[:,0])
      print('times',times)
      if times<timelen:continue
      if times>190:continue
      count+=1
      avgStrain += (strain[-timelen:,:4]-avgStrain)/count
      avgStrainSq += (strain[-timelen:,:4]**2-avgStrainSq)/count
    print('count',count)
    stDev = np.sqrt(avgStrainSq - avgStrain**2)
    labels = [r'$s_1$', r'$s_2$', r'$s_3$']
    colors = ['C0', 'C1', 'C2']
    for i in range(1, 4):
        ax.plot(strain[-timelen:,i], color=colors[i-1], alpha=.3)
        ax.plot(avgStrain[:, i], label=labels[i-1], color=colors[i-1])
        ax.fill_between(
            np.arange(avgStrain.shape[0]),
            avgStrain[:, i] - stDev[:, i],
            avgStrain[:, i] + stDev[:, i],
            color=colors[i-1],
            alpha=0.3
        ) 
    ax.legend()
    fname = 'plots/' + matrix+'VsTimeDash.pdf'
    print('savin ',fname)
    fig.savefig(fname, bbox_inches='tight', transparent=True, format='pdf', dpi=600)
  return

#surPowVsWaveNumber()
#energyVsWaveNumber()
surPowFreqVsWaveNumber()
#areaVsTime()
#MoIAlignStrain()
#sanBernado()
#strainVsTime()
