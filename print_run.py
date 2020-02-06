# variances = [0.8, 0.7, 0.5, 0.3, 0.2, 0.1]
# variances = [0.95, 0.9, 0.6, 0.4]
variances = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
variances = [0.9, 0.8, 0.7, 0.6]
probs = [(i+1)/100. for i in range(100)]
fig = 'AB'
temp = './base%s_dynamic.o %.1f %.2f\n'

with open('run%s_dynamic2.sh'%fig, 'w') as fw:
    for variance in variances:
        for prob in probs:
            fw.write(temp%(fig, variance, prob))


# variance_es = [0.5]
# # variance_gs = [0.1, 0.01, 0.0001]
# variance_gs = [0.01, 0.0001]
# probs = [(i+1)/100. for i in range(100)]
# fig = 'A'
# temp = './base%s_cont.o %.2f %.1f %.2f\n'

# with open('run%s_cont.sh'%fig, 'w') as fw:
#     for variance_e in variance_es:
#         for variance_g in variance_gs:
#             for prob in probs:
#                 fw.write(temp%(fig, variance_g, variance_e, prob))