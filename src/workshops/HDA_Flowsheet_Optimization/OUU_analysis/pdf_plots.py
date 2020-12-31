import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
import pandas as pd

run_data = pd.read_pickle("run_data_100.pkl")

print(run_data)

# Plot the cap cost and operating cost
plt.figure(1)
plt.subplot(2, 2, 1)
run_data["arrhenius"].plot.hist()
plt.title("arrhenius factor")

plt.subplot(2, 2, 2)
run_data["act_energy"].plot.hist()
plt.title("activation energy")

plt.subplot(2, 2, 3)
run_data["cap_cost"].plot.hist()
plt.title("Capital Cost")

plt.subplot(2, 2, 4)
run_data["op_cost"].plot.hist()
plt.title("Operating Cost")
plt.show()

plt.figure(2)
plt.subplot(2, 3, 1)
run_data["h101_temp"].plot.hist()
plt.title("H101 Temperature")

plt.subplot(2, 3, 2)
run_data["r101_temp"].plot.hist()
plt.title("R101 Temperature")

plt.subplot(2, 3, 3)
run_data["f101_temp"].plot.hist()
plt.title("F101 Temperature")

plt.subplot(2, 3, 4)
run_data["h102_temp"].plot.hist()
plt.title("H102 Temperature")

plt.subplot(2, 3, 5)
run_data["cond_rr"].plot.hist()
plt.title("C101 Reflux Ratio")

plt.subplot(2, 3, 6)
run_data["cond_br"].plot.hist()
plt.title("C101 Boilup Ratio")

plt.show()
