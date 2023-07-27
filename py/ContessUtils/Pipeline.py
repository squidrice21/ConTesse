# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

import argparse
import os
from pathlib import Path
import subprocess


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('spec', help='Path to spec file', type=Path)
    parser.add_argument('exe', help='Path to cpp exe file', type=Path)
    parser.add_argument('case', help='Case name and test type', type=str)
    parser.add_argument('--out', help='Output path', type=Path)
    parser.add_argument('--xml', help='Path to a xml file', type=Path)
    args = parser.parse_args()

    # Run contess
    cmd = [str(args.exe), 'spec=' + str(args.spec.resolve()), '[pipeline]']
    if args.xml:
        cmd += ['-r', 'xml', '-o', str(args.xml)]
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        cmd += ['out=' + str(args.out)]
    subprocess.run(cmd)

    # Check result
    test_type = args.case.split(':')[0]
    case_name = args.case.split(':')[1]
    if len(args.case.split(':')) > 2:
        case_name = args.case.split(':')[1] + '_cam' + args.case.split(':')[2]
    case_name = case_name.replace('.json', '')
    output_dir = os.path.join('./', 'test_out')
    if test_type == 'basic_pipeline':
        output_dir = os.path.join(output_dir, 'basic_pipeline')
    else:
        output_dir = os.path.join(output_dir, 'multicamera_pipeline')
    svg_path = Path(os.path.join(output_dir, f'{case_name}_out.svg'))
    if not svg_path.exists():
        exit(-1)
