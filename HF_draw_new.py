from xrsana import xrs_ComptonProfiles
import matplotlib.pyplot as plt
path_file = '/home/hushiqi/work/xrs_ana/xrsana/resources/data/ComptonProfiles.dat'


cp,bing_en,occ_num,shell_name = xrs_ComptonProfiles.PzProfile("Na",path_file)
plt.figure(figsize=(10, 6))
plt.plot(cp[:, 0], cp[:, 1], 'b-', linewidth=1, label='Total Compton Profile')
plt.plot(cp[:, 0], cp[:, 2], 'r-', linewidth=1, label='Total Compton Profile')
plt.plot(cp[:, 0], cp[:, 3], 'g-', linewidth=1, label='Total Compton Profile')
plt.xlabel('Pz (a.u.)')
plt.ylabel('Compton Profile')
plt.title('Compton Profile for Na')
plt.legend()
plt.show()
