# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

def on_load():
	import sys
	import os
	import paraview.servermanager

	# To have the __file__ variable when running a paraview script from the GUI
	# we need version 5.8 or higher.
	major_version = paraview.servermanager.vtkSMProxyManager.GetVersionMajor()
	minor_version = paraview.servermanager.vtkSMProxyManager.GetVersionMinor()
	assert major_version >= 5, "Please use paraview version 5.8.0 or higher"
	assert minor_version >= 8, "Please use paraview version 5.8.0 or higher"

	# Get current directory that has scripts in it
	cwd = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
	if cwd not in sys.path: sys.path.append(cwd)

	# Remember this directory
	import ContessParaview
	ContessParaview.Globals.set_scripts_dir(cwd)
	return ContessParaview


# make the ContessParaview package available
on_load()
from ContessParaview import Pipeline
from ContessParaview.Pipeline import Entity, Color, ShowRay

# For test only
#c=Pipeline.SerializedViewer("/Users/hoshyari/code_adobe/npr/contour-tessellation/_build-debug/split_sub_tri.json")
#c=Pipeline.SerializedViewer("/Users/hoshyari/Downloads/torus_fine2_v1.json")
#c=Pipeline.SerializedViewer("/Users/chenxili/project/contour-tessellation-build-debug/test_out/ray_cast_QI/torus_fine2_v1.json")
c=Pipeline.SerializedViewer("/Users/hoshyari/Downloads/bunny_000001_stitched_chain.json")