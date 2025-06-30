GPU: ``gpu``
==========

.. module:: spatialleiden.gpu

.. currentmodule:: spatialleiden

The :py:mod:`spatialleiden.gpu` module provides functions to generate clusters using
GPU acceleration.

While the functions in the `gpu` module are similar to the functions in the main package,
some differences exist due to different feature-completeness of the underlying
implementations.

.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   gpu.leiden
   gpu.spatialleiden
   gpu.spatialleiden_multimodal
   gpu.multiplex_leiden
