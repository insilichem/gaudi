#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sys
from pychimera import patch_environ, enable_chimera

if __name__ == '__main__':
    patch_environ()
    enable_chimera()
    pytest.main(sys.argv[1:])