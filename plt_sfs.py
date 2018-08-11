import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    graphs = []
    labels = []
    for i in range(4):
        sfs = np.load('saved_data/sfs_strat_' + str(i) + '.npy')
        sfs /= np.sum(sfs)
        graphs += plt.plot(np.log(sfs[0:300]))
        labels.append('quartile ' + str(i))
        plt.hold(True)
    plt.legend(graphs, labels)
    plt.xlabel('frequency')
    plt.ylabel('log proportion of genes')
    plt.title('log of site frequency spectrum')
    plt.show()
