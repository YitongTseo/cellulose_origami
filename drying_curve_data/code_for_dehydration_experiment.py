# TO RUN FROM RASPBERRY PI TERMINAL:
# python3 -m venv tinker_env
# source tinker_env/bin/activate
# run with python.

import time
from tinkerforge.ip_connection import IPConnection, Error
from tinkerforge.bricklet_load_cell_v2 import BrickletLoadCellV2
import json
 
# 20 grams was calibrated as 2,000,000 grams (so an effective shift of 100,000)
CALIBRATION_CONVERSION_FACTOR = 100000

HOST = "localhost"
PORT = 4223
LOAD_CELL_UID = "23B3"
 
ipcon = IPConnection()
load_cell = BrickletLoadCellV2(LOAD_CELL_UID, ipcon)
 
try:
    ipcon.connect(HOST, PORT)
    print("Connected to brickd")
except Error as e:
    print(f"Error connecting to brickd: {e}")


def start_experiment(output_name, wait_time=1, total_time = 60* 60 * 3, save_frequency=60):
    weights_by_time = []
    for idx in range(int(total_time / wait_time)):
        weight = load_cell.get_weight() / CALIBRATION_CONVERSION_FACTOR
        weights_by_time.append((idx * wait_time, weight))
        if idx % save_frequency == 0:
            with open(output_name, "w") as f:
                json.dump(weights_by_time, f)
            print('saving! for second ', idx, ' at weight:', weight)
        time.sleep(wait_time)

 

start_experiment(output_name="03.26.25_WT_exp_v1_baselinewithoutcellulose.json")
