.. GaudiMM: Genetic Algorithms with Unrestricted
   Descriptors for Intuitive Molecular Modeling
   
   http://bitbucket.org/insilichem/gaudi
  
   Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
   
        http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Quick usage guide
=================

Running GAUDI jobs is quite easy with :mod:`gaudi.cli.gaudi_run`:

.. code-block:: console

    gaudi run /path/to/some_file.gaudi-input

After the job is completed, you can check the results with :mod:`gaudi.cli.gaudi_view`:

.. code-block:: console

    gaudi view /path/to/some_file.gaudi-output

You can choose between two molecular visualization tools: Chimera itself (using `GaudiView <https://bitbucket.org/insilichem/gaudiview>`_ extension), or our in-house GUI, `GAUDInspect <https://bitbucket.org/insilichem/gaudinspect>`_.

Both tools support lazy-load, filtering and sorting, so choose whichever you prefer.