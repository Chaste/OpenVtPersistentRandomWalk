#  this will call the executable TestPersistentRandomWalkInputing and run the test, inputting the parameter -num_cells

import os
import json
import csv
import math
import re
import time


# read in the json file test.json

script_dir = os.path.dirname(os.path.abspath(__file__))
json_path = os.path.join(script_dir, "input.json")

print(f"Reading parameters from {json_path}")

with open(json_path, "r") as f:
    data = json.load(f)

num_cells = data['sim']['num_sims']
end_time = data['model']['max_time']
dt = 1.0
sampling_rate = data['sim']['output_per']
model_type = data['model']['model']
output_name = data['sim']['output_name']

print(f"Number of cells: {num_cells}")
print(f"End time: {end_time}")
print(f"Time step: {dt}")
print(f"Model type: {model_type}")
print(f"Output name: {output_name}")

base_runner = f"/home/chaste/build/projects/OpenVtPersistentRandomWalk/test/TestPersistentRandomWalk -num_cells {num_cells} -end_time {end_time} -sampling_rate {sampling_rate} -dt {dt} -model_type {model_type} -output_name {output_name}"

if model_type == "MODEL001":
    movement_parameter = 0.01

    model_1_runner = base_runner + f" -movement_parameter {movement_parameter}"

    print(model_1_runner)
    os.system(model_1_runner)

elif model_type == "MODEL004":
    # velocity_change_rate = data['model']['model_args']['ptime']
    velocity_change_rate = data['sim']['output_per']
    motility_parameter =  data['model']['model_args']['speed']
    bias_parameter = data['model']['model_args']['bias']

    bias_direction_x = data['model']['model_args']['bdir']
    bias_direction_y = data['model']['model_args']['bdir']

    model_4_runner = base_runner + f" -velocity_change_rate {velocity_change_rate} -motility_parameter {motility_parameter} -bias_parameter {bias_parameter} -bias_direction_x {bias_direction_x} -bias_direction_y {bias_direction_y}"

    print(model_4_runner)
    os.system(model_4_runner)


time.sleep(2)  # wait for 2 seconds to ensure all data is processed

# read in the results.viznodes data as a plain text file, line by line, and write it to a clean csv

with open(f"/home/chaste/output/OpenVT/Test01PersistentRandomWalk/{output_name}/results_from_time_0/results.viznodes", 'r') as f:
    lines = f.readlines()

    # create a dictionary to hold the data
    data = {}
    num_sims = 100
    for i in range(num_sims):
        data[i] = {'time': [], 'x_pos': [], 'y_pos': []}

    # process each line
    for line in lines:
        # split the line into columns
        split_values = [val for val in re.split(r'[\t \n]+', line) if val]
        # do something with the columns
        for iter, value in enumerate(split_values):
            cell_id = math.floor((iter - 0.5) / 2)
            if cell_id < 0:
                cell_id = 1
            if iter == 0:
                for i in range(num_sims):
                    data[i]['time'].append((value))
            else:
                if (iter % 2 == 1):
                    data[cell_id]['x_pos'].append((value))
                elif (iter % 2 == 0):
                    data[cell_id]['y_pos'].append((value))

time.sleep(2)  # wait for 2 seconds to ensure all data is processed

# # create a csv file to dump the data:
with open(f"/home/chaste/output/OpenVT/Test01PersistentRandomWalk/{output_name}/results_from_time_0/data.csv", 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['time','id', 'com_1', 'com_2'])
    for cell_id, values in data.items():
        for i in range(len(values['time'])):
            csvwriter.writerow([values['time'][i], cell_id, values['x_pos'][i], values['y_pos'][i]])
