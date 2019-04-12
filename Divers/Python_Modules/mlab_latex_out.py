#!/usr/bin/env python
# -*- coding: utf-8 -*-







import os
import os.path
import shutil
import subprocess
import tempfile

import numpy
import PIL.Image as Image

import pylab as pl
from mayavi import mlab


__all__ = ['mlab_imshow_latex']


DPI = 4800


template = '''
\\documentclass[margin=10pt]{standalone}
\\begin{document}
    $ %s $
\\end{document}
'''.strip()


class TemporaryDirectory(object):
    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()

    def __enter__(self):
        return self.tmp_dir

    def __exit__(self, type, value, traceback):
        shutil.rmtree(self.tmp_dir)
        self.tmp_dir = None


def pillow_to_array(img):
    rgba_img_array = numpy.array(img.getdata(), numpy.uint8)
    rgba_img_array.shape = (img.size[1], img.size[0], 4)
    return rgba_img_array


def color_array(rgba_img_array, color):
    if len(color) == 3:
        color = (color[0], color[1], color[2], 255)
    index_array = numpy.where(rgba_img_array == (0, 0, 0, 255))
    rgba_img_array[index_array[0], index_array[1]] = color
    return rgba_img_array


def mlab_imshow_color(img_array, **kwargs):
    """
    Plot a color image with mayavi.mlab.imshow.
    img_array is a ndarray with dim (n, m, 4) and scale (0->255]
    **kwargs is passed onto mayavi.mlab.imshow(..., **kwargs)
    """
    my_lut = pl.c_[img_array.reshape(-1, 4)]
    my_lut_lookup_array = pl.arange(img_array.shape[0] * img_array.shape[1]).reshape(img_array.shape[0],
                                                                                     img_array.shape[1])

    the_imshow = mlab.imshow(my_lut_lookup_array, colormap='binary', **kwargs)      # temporary colormap
    the_imshow.module_manager.scalar_lut_manager.lut.table = my_lut
    mlab.draw()

    return the_imshow


def tex_out(string):
    document = template % (string)
    with TemporaryDirectory() as tmp_dir:
        with open(os.path.join(tmp_dir, 'math.tex'), 'w') as math_file:
            math_file.write(document)
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(['latex', 'math.tex'], cwd=tmp_dir, stdout=devnull, stderr=devnull)
            subprocess.check_call(['dvipng', '-bg', 'Transparent', '-D', str(DPI), '-T', 'tight', '-o', 'math.png',
                                   'math.dvi'], cwd=tmp_dir, stdout=devnull, stderr=devnull)
        math_img = Image.open(os.path.join(tmp_dir, 'math.png'))
    math_img_array = pillow_to_array(math_img)
    return math_img_array


def mlab_imshow_latex(math_string, textcolor=None, **kwargs):
    img_array = tex_out(math_string)
    if textcolor is not None:
        img_array = color_array(img_array, textcolor)
    return mlab_imshow_color(img_array, **kwargs)
