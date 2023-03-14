import pandas as pd
from matplotlib import pyplot as plt
# naming the X axis
plt.xlabel('Number of Vertices')
# naming the Y axis
plt.ylabel('Running Time(ms)')
# plt.ylabel('Approximation Ratio')

# giving a title to my graph
plt.title('Running time of APPROX-1 & APPROX-2')
# plt.title('Approximation Ratio')
# plt.title('Running time of CNF-SAT')


# # corresponding x & y axis values
# x = pd.Index([5,10,11,12,13,14,15,16,17,18,19,20])
# y = pd.Series([0.008,0.019,0.056,0.333,0.4405,0.571,0.110,3.287,2.122,3.578,4.872,24.40])
# # y_error = pd.Series([0.5,0.08,0.3,0.8,0.01,0.07,0.03,0.7,0.04,0.09,0.08,1.4])

# plt.plot(x, y, label="CNF-SAT")

# RunningTimeOf_VC1_And_VC2_.png
x1 = pd.Index([5,10,15,20,25,30,35,40,45,50])
y1 = pd.Series([0.0033,0.0012,0.0022,0.0038,0.0026,0.0017,0.0021,0.0040,0.0052,0.0045])
# plt.errorbar(x1, y1,xerr=0.01 ,yerr = [0.5,0.8,0.9,1.5,1.1,1.2,0.45,1.5,0.1,1.6],fmt='*',ecolor = 'black',color='green')
plt.errorbar(x1, y1,xerr=0.001 ,yerr = [0.0001,0.0002,0.0002,0.0005,0.0003,0.0007,0.0001,0.0001,0.0003,0.0002],fmt='*',ecolor = 'black',color='green')

plt.plot(x1, y1, label = "APPROX-1")
 
# line 2 points
x2 = pd.Index([5,10,15,20,25,30,35,40,45,50])
y2 = pd.Series([0.0011,0.0006,0.0017,0.0040,0.0033,0.0048,0.0016,0.0036,0.0030,0.0051])
# plotting the line 2 points
plt.errorbar(x2, y2,xerr=0.001 ,yerr = [0.0001,0.0002,0.0002,0.0005,0.0003,0.0007,0.0001,0.0001,0.0003,0.0002],fmt='*',ecolor = 'black',color='green')
# plt.errorbar(x2, y2,xerr=0.001 ,yerr = 0.0001,fmt='*',ecolor = 'black',color='green')

plt.plot(x2, y2, label = "APPROX-2")
# function to show the plot

# //approximation ratio graph
# CNF SAT AR

# x1 = pd.Index([5,10,15,20])
# y1 = pd.Series([0, 1.5,1.166,1.375])
# # x_error = pd.Series([0.001, 0.020, 0.0178, 0.004])
# y_error1 = pd.Series([0.3,0.1,0.5,0.3])
# plt.plot(x1, y1, 'b-',label="APPROX-1")

# x = pd.Index([5,10,15,20])
# y = pd.Series([2,4,6,8])
# y_error = pd.Series([0.5,0.7,0.4,0.6])

# plt.plot(x, y, color = 'green',label="CNF-SAT")

# x2 = pd.Index([5,10,15,20])
# y2 = pd.Series([2,2,1.66,2.25])
# y_error2 = pd.Series([0.1,0.2,0.5,0.1])

# plt.plot(x2, y2, 'r--',label="APPROX-2")

plt.legend()
# plt.errorbar(x1,y1, xerr=0.03, yerr=0.2)
# plt.errorbar(x1, y1,xerr=0.01 ,yerr = [0.5,0.8,0.9,1.5,1.1,1.2,0.45,1.5,0.1,1.6],fmt='*',ecolor = 'black',color='green')
# plt.errorbar(x, y,xerr=0.001 ,yerr = [0.8,0.1,0.2,0.5,0.3,0.7,0.1,1.0,0.45,0.8,0.3,1.6],fmt='*',ecolor = 'green',color='green')
# plt.errorbar(x2, y2,xerr=0.001 ,yerr = [0.8,0.1,0.2,0.5,0.3,0.7,0.1,0.8,0.3,1.6],fmt='*',ecolor = 'black',color='green')

# plt.savefig("APROXIMATION_RATIO.png")
# plt.savefig("RunningTimeOf_CNF-SAT_.png")
plt.savefig("RunningTimeOf_VC1_And_VC2_.png")
plt.show()

