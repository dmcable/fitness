import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special
import pdb
import time

if __name__ == '__main__':
    # N = 200
    # epsilon = 0.0001
    # selection = 0
    # distribution = np.zeros(N)
    # distribution[0] = 1
    # num_generations = 100
    # mut_rate = 0.01
    # ncr = np.zeros(N)
    # for i in range(N):
    #     ncr[i] = scipy.special.comb(N, i)
    # for trial in range(num_generations):
    #     new_distribution = np.zeros(N)
    #     for i in range(N):
    #         index = i
    #         for j in range(N):
    #             k = j
    #             prob = (1 - selection) * index * 1.0 / N
    #             prob += mut_rate * (1 - 2 * prob)
    #             new_distribution[k] += ncr[k] * math.pow(prob, k) * math.pow(1 - prob, N - k) * distribution[index]
    #     new_distribution /= sum(new_distribution)
    #     if(np.max(abs(distribution - new_distribution)) < epsilon):
    #         break
    #     distribution = new_distribution / sum(new_distribution)
    #     #print(distribution)
    #     print(np.sum(distribution))
    #     distribution[0] = 0
    #     distribution[N - 1] = 0
    # print(trial)
    # plt.plot(distribution)
    # plt.show()

    num_generations = 1000
    N_init = 10000
    num_samples = 100000
    #sim_distribution = np.zeros(N + 1)
    mut_rate = 1e-8  # place holder
    growth_rate = 1.014
    results = []
    for i in range(num_samples):
        N = N_init
        x = 1
        existence_len = int(np.random.rand() * num_generations)
        for g in range(existence_len):
            prob = mut_rate + x * 1.0 / N
            N = int(N * growth_rate)
            x = np.random.poisson(N * prob)
            if(x == 0 or x >= N):
                break
        if(x > 0 and x <= 100):
            results.append(x)
    plt.hist(results, bins=100)
    plt.show()
#
# start = time.time()
# x = np.random.poisson(1)
# end = time.time()
# print(end - start)
