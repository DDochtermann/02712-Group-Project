import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.stats import linregress

CAND_COLORS = np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f',
                '#ff7f00','#cab2d6','#6a3d9a', '#90ee90', '#9B870C', '#2f4554',
                '#61a0a8', '#d48265', '#c23531'])

# variances = [0.5, 0.4, 0.3, 0.2, 0.1]
# variances = [0.8, 0.7, 0.5, 0.3, 0.2, 0.1]
# variances = [0.9, 0.8, 0.7, 0.6]
variances = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
probs = [(i+1)/100. for i in range(100)]
xs = np.arange(0, 1, 0.001)

file_template = "log/fig2A/outfileA_K_5000_memory_%f_BWE1_0.950000_BM1E1_%f_BM2E1_%f_mutationAa_0.000010_N1000.txt"
# file_template = "log/fig2B/outfileB_K_5000_memory_%f_BWE1_0.950000_BM1E1_%f_BM2E1_%f_mutationAa_0.000010_N1000.txt"
# file_template = "log/fig2AB/outfileA_K_5000_memory_%f_BWE1_0.950000_BM1E1_%f_BM2E1_%f_mutationAa_0.000010_N1000.txt"
BIRTH_RATE_W_E1 = 0.95
n_run = 10000.

results = []
interpotaed_results = []
for variance in variances:
    dots = []
    for prob in probs:
        filename = file_template%(prob, BIRTH_RATE_W_E1+variance, BIRTH_RATE_W_E1-variance)
        no_ext_freq = float(open(filename).readlines()[1].rstrip().split()[1])/n_run
        dots.append(no_ext_freq)
    results.append(dots)
    cs = CubicSpline(probs, dots)
    interpotaed_results.append(cs(xs))
    # slope, intercept, r_value, p_value, std_err = linregress(probs, dots)
    # interpotaed_results.append(slope*xs+intercept)


results = np.array(results)

plt.style.use('bmh')
fig, ax1 = plt.subplots(1, 1)
# ax1.plot(num_data, fighting_back, 'orange', label='fighting back')
for i in range(len(variances)):
    ax1.scatter(probs, results[i], c=CAND_COLORS[i*2+1], s=15, alpha=1)
    ax1.plot(xs, interpotaed_results[i], c=CAND_COLORS[i*2], label='Variance=%.2f'%variances[i]**2)

# ax1.set_ylim(0, 0.09)
ax1.set_xlabel('Phenotypic Memory')
ax1.set_ylabel('Probability of evolutionary rescue')
ax1.grid(True)
ax1.legend()
plt.show()