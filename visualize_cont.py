import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

CAND_COLORS = np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f',
                '#ff7f00','#cab2d6','#6a3d9a', '#90ee90', '#9B870C', '#2f4554',
                '#61a0a8', '#d48265', '#c23531'])


# variances = [0.1, 0.01, 0.001, 0]
variances = [0.25, 0.2, 0.1, 0.05, 0]
# variances = [0.8, 0.7, 0.5, 0.3, 0.2, 0.1]
probs = [(i+1)/100. for i in range(100)]
xs = np.arange(0, 1, 0.001)

fig = 'A'

file_template_dis = "log/fig2%s/outfile%s"%(fig, fig)+"_K_5000_memory_%f_BWE1_0.950000_BM1E1_%f_BM2E1_%f_mutationAa_0.000010_N1000.txt"
file_template_cont = "log/fig2%s_cont/outfile%s"%(fig, fig)+"_K_5000_memory_%f_BWE1_0.950000_BM1E1_1.450000_BM2E1_0.450000_mutationAa_0.000010_std%f_N1000.txt"
BIRTH_RATE_W_E1 = 0.95
n_run = 10000.

results = []
interpotaed_results = []
for variance in variances:
    dots = []
    for prob in probs:
        if variance == 0:
            filename = file_template_dis%(prob, BIRTH_RATE_W_E1+0.5, BIRTH_RATE_W_E1-0.5)
        elif variance == 0.001:
            filename = file_template_cont%(prob, 0)
        else:
            filename = file_template_cont%(prob, variance)
        no_ext_freq = float(open(filename).readlines()[1].rstrip().split()[1])/n_run
        dots.append(no_ext_freq)
    results.append(dots)
    cs = CubicSpline(probs, dots)
    interpotaed_results.append(cs(xs))


results = np.array(results)

plt.style.use('bmh')
fig, ax1 = plt.subplots(1, 1)
# ax1.plot(num_data, fighting_back, 'orange', label='fighting back')
for i in range(len(variances)):
    ax1.scatter(probs, results[i], c=CAND_COLORS[i*2+1], s=15, alpha=1)
    ax1.plot(xs, interpotaed_results[i], c=CAND_COLORS[i*2], label='std=%.5f'%variances[i])

ax1.set_xlabel('Phenotypic Memory')
ax1.set_ylabel('Probability of evolutionary rescue')
ax1.grid(True)
ax1.legend()
# plt.title('Simulation results with continuous phenotypic distribution')
plt.show()