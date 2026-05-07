from xrsana import xrs_fileIO
import matplotlib.pyplot as plt
path_file = '/home/hushiqi/work/xrs_ana/xrsana/resources/data/ComptonProfiles.dat'
element = 'Na'

HF_data, occupation, bindingen, colnames = xrs_fileIO.readbiggsdata(path_file, element)

plt.figure(figsize=(10, 6))
pz = HF_data[:, 0]
total = HF_data[:, 1]
plt.plot(pz, total, 'b-', linewidth=1, label='Total Compton Profile')

for i in range(2, HF_data.shape[1]):
    plt.plot(pz, HF_data[:, i], label=f'Shell {i-1} (Occ: {occupation[i-2]}, BE: {bindingen[i-2]} eV)')

plt.yscale('log')
plt.xlabel('pz (atomic units)', fontsize=12)
plt.ylabel('Compton Profile J(pz)', fontsize=12)
plt.title(f'Compton Profile of  {element}', fontsize=14)
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend(fontsize=11)

plt.tight_layout()

plt.show()