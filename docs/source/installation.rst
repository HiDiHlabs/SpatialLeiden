Installation
============


PyPi and ``pip``
----------------

To install ``spatialleiden`` from `PyPI <https://pypi.org/project/spatialleiden/>`_
using ``pip`` just run

.. code-block:: bash

    pip install spatialleiden


bioconda and ``conda``
----------------------

``spatialleiden`` is available for `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
installations from the `bioconda <https://bioconda.github.io/recipes/spatialleiden/README.html>`_
channel.


.. code-block:: bash

    conda install bioconda::spatialleiden

.. note::

    Of course, it is also possible to use ``mamba`` instead of ``conda``
    to speed up the installation.


From GitHub
-----------

You can install the latest versions directly from
`GitHub <https://github.com/HiDiHlabs/SpatialLeiden>`_. To do so clone the repository
using the ``git clone`` command. Navigate into the downloaded directory and install
using

.. code-block:: bash

    pip install -e .

If you want to install the development version you can install the additional optional
dependencies with

.. code-block:: bash

    pip install -e '.[dev]'
