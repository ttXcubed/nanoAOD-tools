# -*- coding: utf-8 -*-

import subprocess
import argparse
import re
import itertools
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', type=Path, required=True)
parser.add_argument('--output_dir', type=Path, help="Output dir, takes input dir if not specified.")
parser.add_argument('--year', type=int, default=2018, help="Year of samples [Default: 2018].")
args = parser.parse_args()

# If output directory is not provided, assume we want the output to be
# alongside the input directory.
if args.output_dir is None:
    args.output_dir = args.input_dir

n_trees_list = [500, 100]
learning_rate_list = [0.01, 0.03, 0.05, 0.08]
max_depth_list = [3, 5, 8]
min_child_weight_list = [.01, .05, .1, .5]

output_dir_log = args.output_dir / 'log_scan_BDT_params'
output_dir_log.mkdir(parents=True, exist_ok=True)
output_dir_log.chmod(0o777)
subprocess.call(['chmod', '-R', 'a+w+x+r', str(output_dir_log)])

with open('scan_BDT_params_condor', 'w+') as condor_f:
    condor_f.write(
        "\n".join([
            'Executable = scan_BDT_params_executable.sh',
            'Should_Transfer_Files = NO',
            'Getenv = True',
            'Requirements = ( OpSysAndVer == "CentOS7" )'
        ]) + '\n'  # Add final newline.
    )

    for trees, lr, maxDepth, minChildW in itertools.product(n_trees_list, learning_rate_list, max_depth_list, min_child_weight_list):
        condor_f.write(
            "\n".join([
                f'arguments = "--input_dir {args.input_dir} --n_trees {trees} --learning_rate {lr} --max_depth {maxDepth} --min_child_weight {minChildW} --output_dir {args.output_dir} --year {args.year}"',
                f'Output = {output_dir_log}/log_trees{trees}_lr_{lr}_maxDepth_{maxDepth}_minChildWeight_{minChildW}.$(Process).out',
                f'Error = {output_dir_log}/log_trees{trees}_lr_{lr}_maxDepth_{maxDepth}_minChildWeight_{minChildW}.$(Process).err',
                f'Log = {output_dir_log}/log_trees{trees}_lr_{lr}_maxDepth_{maxDepth}_minChildWeight_{minChildW}.$(Process).log',
                'queue'
            ]) + '\n'
        )
        