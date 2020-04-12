import pandas as pd
import matplotlib.pyplot as plt

dataframe = pd.read_csv("test.csv")
x = dataframe.NumThreads
y1= dataframe.Arr1
y2= dataframe.Arr2
y3= dataframe.Arr3
y4= dataframe.Arr4
y5= dataframe.Arr5
y6= dataframe.Arr6



fig1=plt.figure(1)
plt.plot(x, y1)
fig2=plt.figure(2)
plt.plot(x, y2)
fig3=plt.figure(3)
plt.plot(x, y3)
fig4=plt.figure(4)
plt.plot(x, y4)
fig5=plt.figure(5)
plt.plot(x, y5)
fig6=plt.figure(6)
plt.plot(x, y6)


fig1.suptitle('the relation between run-time and number of threads for N equals 10 ', fontsize=14, fontweight='bold')
fig2.suptitle('the relation between run-time and number of threads for N equals 100', fontsize=14, fontweight='bold')
fig3.suptitle('the relation between run-time and number of threads for N equals 1000', fontsize=14, fontweight='bold')
fig4.suptitle('the relation between run-time and number of threads for N equals 10000', fontsize=14, fontweight='bold')
fig5.suptitle('the relation between run-time and number of threads for N equals 100000', fontsize=14, fontweight='bold')
fig6.suptitle('the relation between run-time and number of threads for N equals 1000000', fontsize=14, fontweight='bold')

ax1 = fig1.add_subplot(111)
ax1.set_xlabel('numThreads')
ax1.set_ylabel('Run-time')
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('numThreads')
ax2.set_ylabel('Run-time')
ax3 = fig3.add_subplot(111)
ax3.set_xlabel('numThreads')
ax3.set_ylabel('Run-time')
ax4 = fig4.add_subplot(111)
ax4.set_xlabel('numThreads')
ax4.set_ylabel('Run-time')
ax5 = fig5.add_subplot(111)
ax5.set_xlabel('numThreads')
ax5.set_ylabel('Run-time')
ax6 = fig6.add_subplot(111)
ax6.set_xlabel('numThreads')
ax6.set_ylabel('Run-time')
fig1.savefig("1.png")
fig2.savefig("2.png")
fig3.savefig("3.png")
fig4.savefig("4.png")
fig5.savefig("5.png")
fig6.savefig("5.png")
plt.show()  