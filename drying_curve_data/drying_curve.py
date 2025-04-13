import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Apply global style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['axes.edgecolor'] = '#333333'
mpl.rcParams['grid.alpha'] = 0.3
mpl.rcParams['grid.linestyle'] = '--'


import os
print(os.getcwd())  # Prints the current working directory

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

def drying_curve_fit(filepaths, starting_index, extrapolate = False, 
                     color = None, label = None, linewidth=None, 
                     linestyle=None, marker=None, vline = 0):

    
    all_ys = []
    for filepath in filepaths:
        # read JSON file and parse contents
        with open(filepath, 'r') as file:
            python_obj = json.load(file)

        ys = [line[1] for line in python_obj]
        all_ys += ys

    ys = all_ys[starting_index:]

    #extrapolate last value til the end
    if extrapolate and len(ys) < 30000:
        dif = 30000 - len(ys)
        ys += [ys[-1]]*dif

    #remove high outliers by becoming the average of its two neighbors
    for i in range(len(ys)):
        if ys[i] > ys[0]:
            ys[i] = 0.5*(ys[i-1]+ys[i+1])

    #scale so that it starts at y value of 1
    y0 = ys[0]
    ys = [y / y0 for y in ys]

    #make as many xs as ys
    xs = np.array(range(len(ys)))

    #convert xs to minutes
    xs = [x/60 for x in xs]

    #plot the line
    ax.plot(xs, ys, label=label, color=color, linewidth=linewidth, 
            linestyle=linestyle, marker=marker, alpha=0.9)
    
    #create a vertical dashed line where the slope becomes ~0
    y_index = xs.index(vline)
    plt.plot([vline, vline], [0, ys[y_index] ], linestyle='--', color=color, zorder=2)

    #print the average y value past vline
    print(f"Equilibrium for {label} is {np.average(ys[y_index:])}")


legend_labels = ["Untreated WT", "BslA Activated", "WT NaOH Rehydrated"]
filepaths = [r"json_drying_data\WT_Untreated.json"]
drying_curve_fit(filepaths, starting_index = 29, label = legend_labels[0],
                 color='#F5B041', linewidth=2, vline = 1070)

filepaths = [r"json_drying_data\BslA_NaOH.json"]
drying_curve_fit(filepaths, starting_index = 73, label = legend_labels[1], 
                 color='#5DADE2', linewidth=2, vline = 70)

filepaths = [r"json_drying_data\WT_NaOH.json"]
drying_curve_fit(filepaths, starting_index = 96, label = legend_labels[2],
                  color='#dd9789', linewidth = 2, vline = 90)

# Equilibrium for Untreated WT is 0.1781675639229501
# Equilibrium for BslA Activated is 0.90421924502062
# Equilibrium for WT NaOH Rehydrated is 0.7363704366317478

# Ticks and grid
ax.tick_params(axis='both', which='major', labelsize=10, width=1, length=4)
ax.set_xlim([0, 1200])
ax.set_ylim([0, 1.01])
ax.grid(True, linestyle='--', alpha=0.6)

# Remove top and right spines for cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Tight layout and save (optional)
plt.tight_layout()
#plt.savefig('curvature_plot_test.png', bbox_inches='tight', dpi=300)
plt.savefig('drying_plot.png', bbox_inches='tight', dpi=300)

plt.show()



