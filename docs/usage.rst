Quick usage guide
=================

Running GAUDI jobs is quite easy with :mod:`gaudi.cli.gaudi_run`:

.. code-block:: console

    gaudi run /path/to/some_file.gaudi-input

After the job is completed, you can check the results with :mod:`gaudi.cli.gaudi_view`:

.. code-block:: console

    gaudi view /path/to/some_file.gaudi-output

You can choose between two molecular visualization tools: Chimera itself (using `GaudiView <https://bitbucket.org/jrgp/gaudinspect>`_ extension), or our in-house GUI, `GAUDInspect <https://bitbucket.org/jrgp/gaudiview>`_.

Both tools support lazy-load, filtering and sorting, so choose whichever you prefer.