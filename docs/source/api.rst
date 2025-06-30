API
===

.. currentmodule:: spatialleiden

.. .. toctree::
..    :maxdepth: 2
..    :caption: Modules

..    api/gpu


Multiplex Leiden
----------------

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   leiden
   spatialleiden
   spatialleiden_multimodal
   multiplex_leiden


Multiplex Leiden (GPU)
______________________

The :py:mod:`spatialleiden.gpu` module provides functions to generate clusters using
GPU acceleration.

While the functions in the ``gpu`` module are similar to the functions in the main package,
some differences exist due to different feature-completeness of the underlying
implementations.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   gpu.leiden
   gpu.spatialleiden
   gpu.spatialleiden_multimodal
   gpu.multiplex_leiden


Resolution search
-----------------

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   search_resolution
   search_resolution_latent
   search_resolution_spatial


Utility functions
-----------------

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   distance2connectivity
