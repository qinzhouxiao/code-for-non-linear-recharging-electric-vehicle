import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy as np

from matplotlib.ticker import MultipleLocator, FixedLocator, FuncFormatter
###### Locators for Y-axis
# set tickmarks at multiples of 1.
majorLocator = MultipleLocator(1.)
# create custom minor ticklabels at logarithmic positions
ra = np.array([ [n+(1.-np.log10(i))]  for n in range(0,10) for i in [2,3,4,5,6,7,8,9][::-1]]).flatten()*-1.
minorLocator = FixedLocator(ra)
###### Formatter for Y-axis (chose any of the following two)
# show labels as powers of 10 (looks ugly)
# majorFormatter= FuncFormatter(lambda x,p: "{:.1e}".format(10**x) )
# or using MathText (looks nice, but not conform to the rest of the layout)
majorFormatter= FuncFormatter(lambda x,p: r"$10^{"+"{x:d}".format(x=int(x))+r"}$" )


plt.figure(figsize=(20, 10))
locs = ['NY', 'BAY', 'COL', 'FLA', 'NW', 'NE', 'CAL', 'LKS']
for i in range(len(locs)):
    loc = locs[i]
    plt.subplot(241 + i)
    ax = plt.gca()
    # Set the locators
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    # Set formatter if you like to have the ticklabels consistently in power notation
    ax.yaxis.set_major_formatter(majorFormatter)
    T1 = []
    D1 = []
    tot_time = 0
    num = 0
    file = open(loc + '-exp-ours.txt')
    for line in file:
        time, dis = line.split()
        T1.append(math.log10(float(time)))
        #T1.append(float(time))
        tot_time += float(time)
        num += 1
        D1.append(int(dis) / 10000)
    # plt.scatter(D1, T, c='blue', s=5, marker='o', label='ours')
    sns.regplot(x=D1, y=T1, order=3, ci=95,
        #        scatter=False,
                scatter_kws={'s': 3, 'color': 'blue', 'alpha': 0.5, 'edgecolors': 'None'},
                line_kws={'linestyle': '-', 'color': 'blue'},
                label='ours')
    # print(loc, "ours", tot_time, num, tot_time / num)

    T2 = []
    D2 = []
    i = 0
    tot_time1 = 0
    tot_time2 = 0
    num1 = 0
    num2 = 0
    file = open(loc + '-exp-pulse.txt')
    for line in file:
        time, dis = line.split()
        T2.append(float(time))
        tot_time1 += float(time)
        num1 += 1
        D2.append(int(dis) / 10000)
        if not (D1[i] != D2[i]):
            tot_time2 += float(time)
            num2 += 1
        i += 1
    # print(loc, "Pulse-all", tot_time1, num1, tot_time1 / num1)
    # print(loc, "Pulse-solved", tot_time2, num2, tot_time2 / num2)

    TO = []
    DO = []
    TX = []
    DX = []
    # print(len(D1), len(D2))
    for i in range(len(D1)):
        if abs(D1[i] - D2[i]) > 1:
            TX.append(math.log10(T2[i]))
            DX.append(D2[i])
            if math.log10(T2[i]) < 0:
                print(D1[i], D2[i])
        else:
            TO.append(math.log10(T2[i]))
            DO.append(D2[i])
    # plt.scatter(DO, TO, c='red', s=5, marker='o', label='pulse solved')
    plt.scatter(DX, TX, c='black', s=5, marker='o', label='pulse unsolved')
    sns.regplot(x=DO, y=TO, order=3, ci=95,
           #     scatter=False,
                scatter_kws={'s': 3, 'color': 'red', 'alpha': 0.5, 'edgecolors': 'None'},
                line_kws={'linestyle': '-', 'color': 'red'},
                label='pulse solved')

    plt.title(loc)
    plt.xlabel('Distance (km)')
    plt.ylabel('Construct Time (sec)')
    #plt.yscale('symlog')
    #plt.yscale('log')
    plt.legend(loc='lower right')

plt.savefig('exp1.pdf')
plt.savefig('exp1.eps')
plt.show()
