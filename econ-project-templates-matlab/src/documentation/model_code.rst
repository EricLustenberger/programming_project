.. _model_code:


Model code
======================

The directory **src.model_code** contains source files that might differ by model and which are potentially used at various steps of the analysis.

For example, you may have a class that is used both in the :ref:`analysis` and the :ref:`final` steps. Additionally, maybe you have different utility functions in the baseline version and for your robustness check. You can just inherit from the baseline class and override the utility function then.

Schelling example, :file:`move_until_happy.m`:

.. include:: ../model_code/move_until_happy.m
   :start-after: %{
   :end-before: %}