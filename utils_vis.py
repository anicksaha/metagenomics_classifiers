import matplotlib.pyplot as plt

def generate_variability_plot(variabilities):
    x = []
    for i in range(len(variabilities)):
        x.append(i)
    y = variabilities

    plt.plot(x, y, '-r', linewidth=1.5, markersize=9, color = 'r')
    plt.xlabel("Index")
    plt.ylabel("Identity of Sequence")
    plt.show()