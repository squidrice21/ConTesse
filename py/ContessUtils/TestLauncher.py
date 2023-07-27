# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

import json
import os
from pathlib import Path

import ContessUtils.Async as Async
from ContessUtils.Log import logger
from ContessUtils.Exception import ContessException


def unique_name():
    import uuid
    return uuid.uuid4().hex


class TestSpec:
    MULTI_CAMERA = 0
    PIPELINE = 1
    INDIVIDUAL_CAM = ['000', '001', '002', '010', '011', '012', '020', '021', '022', '100', '101', '102',
                      '110', '112', '120', '121', '122', '200', '201', '202', '210', '211', '212', '220', '221', '222']

    def __init__(self, spec_file: Path, cpp_exe: Path, cwd: Path, output: Path):
        if '.json' in os.path.splitext(spec_file):
            try:
                spec = json.loads(spec_file.read_text())
            except Exception as e:
                raise ContessException(
                    f'Could not read file {spec_file}') from e
        else:
            spec = {'pipeline_test_cases': [str(spec_file.name)]}
        if not cpp_exe.exists():
            raise ContessException(f'Could not find cpp executable {cpp_exe}')

        self.spec_root = spec_file.parent.absolute()
        self.spec_json = spec
        self.cpp_exe = cpp_exe

        multi_camera_lines = []
        for file in spec.get('multicamera_test_cases', []):
            if ':' in file:
                multi_camera_lines.append(file)
            else:
                multi_camera_lines += [file + ':' +
                                       cam for cam in self.INDIVIDUAL_CAM]

        self.multi_camera_cases = [Path(file) for file in multi_camera_lines]
        logger.debug('read following multicamera_test_cases:\n %s',
                     self.multi_camera_cases)
        self.pipeline_cases = [Path(file)
                               for file in spec.get('pipeline_test_cases', [])]
        logger.debug('read following pipeline_test_cases:\n %s',
                     self.pipeline_cases)
        self.cwd = cwd
        self.output = output


class TestReport:
    def __init__(self, name, success, out, err, xml_report, xml_report_path):
        self.name = name
        self.success = success
        self.out = out
        self.err = err
        self.xml_report = xml_report
        self.xml_report_path = xml_report_path


async def run_test(spec: TestSpec, test_type, test_json_path: Path):
    if test_type not in [TestSpec.MULTI_CAMERA, TestSpec.PIPELINE]:
        raise ContessException(f'Invalid type_type: {test_type}')
    if not test_json_path.is_absolute():
        test_json_path = spec.spec_root / test_json_path

    tmp_dir: Path = spec.cwd / ".tmp" / unique_name()
    tmp_dir.mkdir(exist_ok=True, parents=True)
    tmp_json = tmp_dir / "spec.json"
    tmp_xml = tmp_dir / "output.xml"

    # create the temp file
    def _write_spec():
        with tmp_json.open("w") as stream:
            tmp_spec = {
                'multicamera_test_cases': [],
                'pipeline_test_cases': [],
            }
            if test_type == TestSpec.MULTI_CAMERA:
                tmp_spec['multicamera_test_cases'].append(
                    test_json_path.as_posix())
            elif test_type == TestSpec.PIPELINE:
                tmp_spec['pipeline_test_cases'].append(
                    test_json_path.as_posix())
            stream.write(json.dumps(tmp_spec, indent=2))
    await Async.run_func_in_bg(_write_spec)

    # start running the test
    case_name = 'multicamera_pipeline:' + \
        os.path.basename(os.path.basename(test_json_path))
    if test_type == TestSpec.PIPELINE:
        case_name = 'basic_pipeline:' + \
            os.path.basename(os.path.basename(test_json_path))
    args = [
        'python3',
        os.path.join(Path(__file__).parent.absolute(), 'Pipeline.py'),
        f'{tmp_json}',
        spec.cpp_exe.absolute(),
        case_name,
        '--xml',
        tmp_xml.as_posix()
    ]
    if spec.output:
        args += ['--out',
                 spec.output.absolute()
                 ]
    code, out, err = await Async.run_subprocess_in_bg(case_name, args, cwd=spec.cwd)

    # collect output
    xml_report = None
    if tmp_xml.exists():
        xml_report = tmp_xml.read_text()

    return TestReport(test_json_path.name, code == 0, out, err, xml_report, tmp_xml.as_posix())
