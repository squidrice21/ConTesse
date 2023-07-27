# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

'''
How to run this script
python /path/to/py/scripts/run_tests.py /path/to/ \
    spec.json(e.g, test-specs/_spec.json) /path/to/stsContourTessellation.exe
run_tests.py ../contour-tessellation/test-specs/wso_failures_0815.json ./TestsContourTessellation --style ../strokes-matching-build-debug/qviewer/StrokesMatching --style-param ../contour-tessellation/test-specs/style_param/param.json

Preferably run this in a separate folder as it creates a bunch of temporary files
'''

import argparse
import asyncio
import re
from datetime import datetime
from pathlib import Path
from typing import List

import pandas as pd

import ContessUtils.Async as Async
from ContessUtils.TestLauncher import TestSpec, TestReport, run_test
from ContessUtils.Log import logger


async def _run_test(*arg, **kwargs):
    semaphore = kwargs.pop('semaphore')
    # prevent more than certain number things to run at the same time
    async with semaphore:
        return await run_test(*arg, **kwargs)


def _match_lines(pattern, line_list: List):
    result = [re.findall(pattern, l)[0]
              for l in line_list if re.findall(pattern, l)]
    return result


def _make_frames(stats):
    frame = pd.DataFrame(stats)

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    # print(frame)

    sum_stats = []
    inputs = frame['name'].unique()
    for case in inputs:
        case_rows = frame.loc[frame['name'] == case]
        succ_case_rows = frame.loc[(frame['name'] == case) & (
            frame['success'] == True)]
        sum_case = {}
        sum_case['Input'] = case
        sum_case['type'] = case_rows['type'].unique()[0]
        sum_case['Success count'] = succ_case_rows.shape[0]
        sum_case['Input count'] = case_rows.shape[0]
        sum_case['Input faces'] = succ_case_rows['input_faces'].mean()
        sum_case['Control faces'] = succ_case_rows['control_faces'].mean()
        sum_case['Output faces'] = succ_case_rows['output_faces'].mean()
        sum_case['Runtime'] = succ_case_rows['runtime'].mean()
        wso_df = succ_case_rows['wso'].value_counts()
        wso_str = ''
        for key, value in wso_df.items():
            wso_str += f'{key}/{value};'
        sum_case['Trial count'] = wso_str
        sum_stats.append(sum_case)

    # Separate by wso
    sum_wso_stats = []
    for case in inputs:
        wso_l = frame.loc[(frame['name'] == case) & (
            frame['success'] == True)]['wso'].unique()
        for wso in wso_l:
            succ_case_rows = frame.loc[(frame['name'] == case) & (
                frame['success'] == True) & (
                frame['wso'] == wso)]
            case_rows = frame.loc[(frame['name'] == case) & (
                frame['wso'] == wso)]
            sum_case = {}
            sum_case['Input'] = case
            sum_case['type'] = case_rows['type'].unique()[0]
            sum_case['Success count'] = succ_case_rows.shape[0]
            sum_case['Input count'] = case_rows.shape[0]
            sum_case['Input faces'] = succ_case_rows['input_faces'].mean()
            sum_case['Control faces'] = succ_case_rows['control_faces'].mean()
            sum_case['Output faces'] = succ_case_rows['output_faces'].mean()
            sum_case['Runtime'] = succ_case_rows['runtime'].mean()
            wso_df = succ_case_rows['wso'].value_counts()
            wso_str = ''
            for key, value in wso_df.items():
                wso_str += f'{key}/{value};'
            sum_case['Trial count'] = wso_str
            sum_wso_stats.append(sum_case)

    sum_frame = pd.DataFrame(sum_stats)
    sum_frame = sum_frame.sort_values(
        by=['type', 'Input faces'], ascending=True)
    sum_frame = sum_frame.drop(columns=['type'])

    sum_wso_frame = pd.DataFrame(sum_wso_stats)
    sum_wso_frame = sum_wso_frame.sort_values(
        by=['type', 'Input faces'], ascending=True)
    sum_wso_frame = sum_wso_frame.drop(columns=['type'])

    return frame, sum_frame, sum_wso_frame


async def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('spec', help='path to spec file', type=Path)
    parser.add_argument('exe', help='path to cpp exe file', type=Path)
    parser.add_argument('--out', help='Output path', type=Path)
    parser.add_argument(
        '--max-procs', help='Max number of tests to run at the same time', type=int, default=4)
    parser.add_argument('--log-file', help='Path to a log file', type=Path)
    parser.add_argument(
        '--verbose', help='Verbose terminal log', action='store_true')
    parser.add_argument(
        '--latex', help='Path to the table latex file', type=Path)
    args = parser.parse_args()

    if args.log_file:
        args.log_file.resolve().parent.mkdir(parents=True, exists_ok=True)
        logger.add_log_file(args.log_file)

    if args.verbose:
        logger.setLevel('DEBUG')

    semaphore = asyncio.Semaphore(args.max_procs)
    spec = TestSpec(args.spec, args.exe, Path.cwd(), args.out)
    coroutines = []
    for pipeline_case in spec.pipeline_cases:
        coroutines.append(_run_test(spec, TestSpec.PIPELINE,
                          pipeline_case, semaphore=semaphore))
    for multi_camera_case in spec.multi_camera_cases:
        coroutines.append(_run_test(spec, TestSpec.MULTI_CAMERA,
                          multi_camera_case, semaphore=semaphore))

    # once the gather operation is done here all tests are executed
    results: List[TestReport] = await Async.gather(*coroutines)

    # Use log to truly determine success/failure
    for result in results:
        found_input = len(_match_lines(
            r'Input faces: .+', [l for l in result.out if 'Input faces:' in l])) > 0
        found_control = len(_match_lines(
            r'Control mesh faces: .+',
            [l for l in result.out if 'Control mesh faces:' in l])) > 0
        found_output = len(_match_lines(
            r'Output faces: .+', [l for l in result.out if 'Output faces:' in l])) > 0
        if not (found_input and found_control and found_output):
            result.success = False
        else:
            result.success = True

    # you can do further processing here
    num_success = 0
    num_failure = 0
    for result in results:
        if result.success:
            num_success += 1
        else:
            num_failure += 1
    logger.info("Num success = %s", num_success)
    logger.info("Num failure = %s", num_failure)
    logger.info('=========================================')

    stats = []
    for result in results:
        logger.info("name=%s, success=%s", result.name, result.success)

        case = {'name': result.name.replace('.json', ''), 'cam': '000'}
        if ':' in case['name']:
            case['cam'] = case['name'].split(':')[1]
            case['name'] = case['name'].split(':')[0]
            case['type'] = 'multi'
        else:
            case['cam'] = case['name'].split('_')[-1]
            case['name'] = '_'.join(case['name'].split('_')[0:-1])
            case['type'] = 'anim'
        case['success'] = result.success

        # Get mesh statistics
        mesh_open = [l for l in result.out if 'Input mesh' in l][0]
        case['open'] = 'open' in mesh_open
        case['input_faces'] = _match_lines(
            r'Input faces: .+', [l for l in result.out if 'Input faces:' in l])[-1].replace('Input faces: ', '')
        case['input_faces'] = int(case['input_faces'])
        case['control_faces'] = _match_lines(
            r'Control mesh faces: .+',
            [l for l in result.out if 'Control mesh faces:' in l])[-1].replace('Control mesh faces: ', '')
        case['control_faces'] = int(case['control_faces'])
        if result.success:
            case['output_faces'] = _match_lines(
                r'Output faces: .+', [l for l in result.out if 'Output faces:' in l])[-1].replace('Output faces: ', '')
            case['output_faces'] = int(case['output_faces'])
        else:
            case['output_faces'] = 0

        # Get runtime based on the timestamps in the log
        start_str = _match_lines(
            r'^\[[^\]]+\]', [l for l in result.out if '## Perturbation round 1 ##' in l])[0]
        if result.success:
            end_str = _match_lines(
                r'^\[[^\]]+\]', [l for l in result.out if 'Contess remeshing done' in l])[0]
        else:
            end_str = _match_lines(
                r'^\[[^\]]+\]', [l for l in result.out if '[contess]' in l])[-1]
        start_time = datetime.strptime(start_str, '[%Y-%m-%d %H:%M:%S.%f]')
        end_time = datetime.strptime(end_str, '[%Y-%m-%d %H:%M:%S.%f]')
        logger.info(f'Log-based timestamps: {start_time}, {end_time}')
        logger.info(f'Runtime: {end_time - start_time}')

        case['runtime'] = (end_time - start_time).total_seconds()

        # WSO round
        result_lines = _match_lines(r'## WSO round.*', result.out)
        logger.info(f'Run WSO for {len(result_lines)} round(s)')

        case['wso'] = int(len(result_lines))
        stats.append(case)

        if result.success:
            # Output info
            result_lines = _match_lines(r'Wrote to svg:.*', result.out)
            logger.info('\n'.join(result_lines))
        else:
            logger.info(f'Error details: {result.xml_report_path}')

        logger.info('----------------------------------------')

    frame, sum_frame, sum_wso_frame = _make_frames(stats)
    if not args.latex:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(sum_frame)
            # print(sum_wso_frame)
    else:
        with open(str(args.latex).replace('.tex', '_all.csv'), 'w') as f:
            f.write(frame.to_csv(index=False))
        with open(str(args.latex).replace('.tex', '_stats.csv'), 'w') as f:
            f.write(sum_frame.to_csv(index=False))
        with open(str(args.latex).replace('.tex', '_wso_stats.csv'), 'w') as f:
            f.write(sum_wso_frame.to_csv(index=False))
        with open(args.latex, 'w') as f:
            def _header_format(c):
                if ' ' not in c:
                    return c
                else:
                    return '\\makecell{' + c.replace(' ', ' \\\\ ') + '}'

            def _format(x):
                if '_' in x:
                    x = x.replace('_', '\\_')
                return x

            header = [_header_format(c) for c in sum_frame.columns]
            f.write(sum_frame.to_latex(
                header=header, index=False, escape=False, formatters={'Input': _format}))

            header = [_header_format(c) for c in sum_wso_frame.columns]
            f.write(sum_wso_frame.to_latex(
                header=header, index=False, escape=False, formatters={'Input': _format}))

if __name__ == '__main__':
    try:
        Async.run_async_entry_point(main())
    except KeyboardInterrupt:
        logger.warning("Keyboard interrupt. Exiting...")
    except Exception:
        logger.exception("Exception happened. Code failed.")
