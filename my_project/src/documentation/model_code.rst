.. _model_code:


Model code
======================

The Model code is composed of a total of four different m-files. The :file:`aiyagari_solver.m` contains the solution algorithms used to simulate and solve the :cite:'aiyagari' model. :file:`cash_equivalent.m` and :file:`consumption_equivalent.m` calculate the welfare effects of a policy change and :file:`setup.m`: holds some model specifications, which are endogenous. 

The Aiyagari Model
------------------
The Aiyagari Model, :file:`aiyagari_solver.m`:

.. include:: ../model_code/aiyagari_solver.m
   :start-after: %{
   :end-before: %}

*Solution Algorithms and methods*
-To solve the household problem
Endogenous Gridpoint by :cite:'endgridmethod'
or Fixed-Point method depending on the policy change considered,
to assure convergence 
- Aggregation method
bisection method or gradual updating depending on the policy change considered, to assure 
convergence 



The Cash Equivalent
-------------------
The Cash Equivalent, :file:`cash_equivalent.m`:

.. include:: ../model_code/aiyagari_solver.m
   :start-after: %{
   :end-before: %}


The Consumption Equivalent
--------------------------
The Consumption Equivalent, :file:`consumption_equivalent.m`:

.. include:: ../model_code/consumption_equivalent.m
   :start-after: %{
   :end-before: %}


The Setup
--------------------------
The Setup, :file:`setup.m`:

.. include:: ../model_code/setup.m
   :start-after: %{
   :end-before: %}