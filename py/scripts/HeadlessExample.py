# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

# This file will open a file and make some edits and take a screen shot.
# To run it use  pvpython HeadlessExample.py (pvpython can be found in the paraview distribution folder)

# Not a very nice thing do to. But okay for now.
import os
import sys
cwd = os.path.abspath(os.path.dirname(os.path.dirname(__file__))) # add the folder containing this folder: ContessParaview 
if cwd not in sys.path:
    sys.path.append(cwd)


import ContessParaview.Pipeline as Pipeline


def main():

    if len(sys.argv) >= 2:
        output = sys.argv[1]
    else:
        output = '/Users/hoshyari/Desktop/script_image.png'

    c = Pipeline.SerializedViewer("/Users/hoshyari/Downloads/bunny_000001_stitched_chain.json")
    # You can use any other camera
    c.use_contour_camera()
    # You can make other edits here.
    c.set_source_config(Pipeline.Entity.CONTOUR_CAMERA, is_visible=False)
    Pipeline.set_background_color(Pipeline.Color.WHITE)
    # Note that this also works in paraview app itself.
    Pipeline.take_screenshot(output, dims = (2000, 1000))


if __name__ == '__main__':
    main()