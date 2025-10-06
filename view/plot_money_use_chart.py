# plot_money_use_chart.py
#
# Author: Jose Deodoro <deodoro.filho@gmail.com> <jdeoliv@gmu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import json
import matplotlib.pyplot as plt
import glob
from collections import defaultdict

# Pick last run output
fname = sorted(glob.glob('out/*.json'))[-1]
with open(fname, 'r') as f:
    data = json.load(f)

runs_per_slope = data['runs_per_slope']

data = defaultdict(list)

for slope_idx, slope_data in enumerate(runs_per_slope):
    if slope_data:
        # Get the slope value from the first entry
        slope = slope_data[0]['Slope']

        for entry in slope_data[:-1]:
            data[slope].append(entry)

# Average per slope
for slope, entries in data.items():
    # Group by t
    t_to_usingmoney = defaultdict(list)
    for entry in entries:
        t = entry['t']
        t_to_usingmoney[t].append(entry['usingmoney'])

    # Average for each t
    times = sorted(t_to_usingmoney.keys())
    averaged_usingmoney = []
    for t in times:
        lists = t_to_usingmoney[t]
        # Average each position
        avg = [sum(x[i] for x in lists) / len(lists) for i in range(10)]
        averaged_usingmoney.append(avg)

    # Plot
    plt.figure(figsize=(10, 6))
    for i in range(10):
        values = [avg[i] for avg in averaged_usingmoney]
        plt.plot(times, values, label=f'Good {i+1}')

    plt.xlabel('Time (t)')
    plt.ylabel('Using Money')
    print([i['t'] for i in entries])
    plt.title(f'Averaged Evolution of Using Money for Slope {slope} ({len([i for i in entries
              if i['t'] == 1])} runs)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'plot_slope_{slope}.png')
    plt.show()
