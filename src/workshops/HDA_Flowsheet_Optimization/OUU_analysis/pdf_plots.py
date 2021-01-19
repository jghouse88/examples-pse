import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
import pandas as pd

run_data = pd.read_pickle("run_data_100_det.pkl")

print(run_data)

# Plot the cap cost and operating cost
plt.rcParams.update({'font.size': 12})
plt.figure(1)
plt.subplot(2, 2, 1)
run_data["arrhenius"].plot.hist(
    color='b',
    range=[run_data["arrhenius"].min(), run_data["arrhenius"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(63000000000.0, color="red", label="Deterministic case")
plt.title("arrhenius factor")
plt.legend()

plt.subplot(2, 2, 2)
run_data["act_energy"].plot.hist(
    color='b',
    range=[run_data["act_energy"].min(), run_data["act_energy"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(217600.0, color="red", label="Deterministic case")
plt.title("activation energy")
plt.legend()

plt.subplot(2, 2, 3)
run_data["cap_cost"].plot.hist(
   color='g',
   range=[29000, 31000],
   rwidth=0.95,
   label="Stochastic case")
plt.axvline(29926.759, color="red", label="Deterministic case")
plt.title("Capital Cost")
plt.legend()

plt.subplot(2, 2, 4)
run_data["op_cost"].plot.hist(
    color='g',
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(408342.352, color="red", label="Deterministic case")
plt.title("Operating Cost")
plt.legend()
plt.show()

plt.figure(2)
plt.subplot(2, 3, 1)
run_data["h101_temp"].plot.hist(
    color='teal',
    range=[run_data["h101_temp"].min(), run_data["h101_temp"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(568.784, color="red", label="Deterministic case")
plt.legend()
plt.title("H101 Temperature")

plt.subplot(2, 3, 2)
run_data["r101_temp"].plot.hist(
    color='teal',
    range=[run_data["r101_temp"].min(), run_data["r101_temp"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(790.019, color="red", label="Deterministic case")
plt.legend()
plt.title("R101 Temperature")

plt.subplot(2, 3, 3)
run_data["f101_temp"].plot.hist(
    color='teal',
    range=[290, 310],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(298.149, color="red", label="Deterministic case")
plt.legend()
plt.title("F101 Temperature")

plt.subplot(2, 3, 4)
run_data["h102_temp"].plot.hist(
    color='teal',
    range=[run_data["h102_temp"].min(), run_data["h102_temp"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(368.785, color="red", label="Deterministic case")
plt.legend()
plt.title("H102 Temperature")

plt.subplot(2, 3, 5)
run_data["cond_rr"].plot.hist(
    color='teal',
    range=[run_data["cond_rr"].min(), run_data["cond_rr"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(0.802, color="red", label="Deterministic case")
plt.legend()
plt.title("C101 Reflux Ratio")

plt.subplot(2, 3, 6)
run_data["cond_br"].plot.hist(
    color='teal',
    range=[run_data["cond_br"].min(), run_data["cond_br"].max()],
    rwidth=0.95,
    label="Stochastic case")
plt.axvline(1.134, color="red", label="Deterministic case")
plt.legend()
plt.title("C101 Boilup Ratio")

plt.show()

failed_cases = run_data[run_data["run_status"] == 0].index.values
color_cycles = ["red", "magenta", "yellow", "orange"]
color_index = 0
plt.figure(3)
plt.subplot(1, 3, 1)
run_data["arrhenius"].plot.hist(
    color='b',
    range=[run_data["arrhenius"].min(), run_data["arrhenius"].max()],
    rwidth=0.95,
    label="Feasible")

for i in failed_cases:
    label = "Inf. sc. #" + str(i)
    color = color_cycles[color_index]
    plt.axvline(
        run_data["arrhenius"][i], color=color, label=label)
    color_index = color_index+1
plt.title("arrhenius factor")
plt.legend()

plt.subplot(1, 3, 2)
color_index = 0
run_data["act_energy"].plot.hist(
    color='b',
    range=[run_data["act_energy"].min(), run_data["act_energy"].max()],
    rwidth=0.95,
    label="Feasible")
for i in failed_cases:
    label = "Inf. sc. #" + str(i)
    color = color_cycles[color_index]
    plt.axvline(
        run_data["act_energy"][i], color=color, label=label)
    color_index = color_index+1
plt.title("activation energy")
plt.legend()

# remove infeasible cases
run_data_feasible = run_data.drop(failed_cases)
print(run_data_feasible)
plt.subplot(1, 3, 3)
sum_op_cost = [sum(run_data_feasible["op_cost"]),
               sum(run_data_feasible["op_cost_det"])]
x_names = ["sum_stochastic", "sum_deterministic"]
plt.bar(x_names, sum_op_cost)
plt.ylim([39200000, 39300000])
plt.title("Operating Cost")
plt.legend()
plt.show()
