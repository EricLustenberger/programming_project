.. raw:: latex 
   
    \pagebreak

.. _introduction:

Introduction
===============

The goal of this Programming project was to use Waf to improve reproducibility and understanding of the code used for my project module. The implementation in Waf is particularly suitable, since the Paper consists of a baseline case and an analysis 
case, which are both compared to show that model specification matters in the context of welfare analysis. The same code is thus used twice with a minor modification in the second case. The following document explains how the code works, such that my results can be replicated with relative ease. The main code is written in Matlab. I therefore used the Matlab template of the course. Documentation on the design rationale, Waf, and more background is at http://hmgaudecker.github.io/econ-project-templates/ . Moreover, the Matlab version of the template uses a modified and translated version of Stachurski's and Sargent's code accompanying their Online Course :cite:`StachurskiSargent13` for Schelling's (1969, :cite:`Schelling69`) segregation model as the running example.

.. _idea_of_Waf:

Waf and the division of tasks
-----------------------------

The core idea is to split up parts of the code, which will be utilized more than once in order to minimize the occurence of errors and to improve efficiency and readability of the code. Input and output are kept apart and when rerunning the code, the output is deleated and replaced by an updated version. Furthermore, the different source files are kept in different folders connected by wscript files, which govern the order and the dependencies of the different files. 

.. literalinclude:: ../wscript

.. _subdirectories:

How the directories are structured:
-----------------------------------

The code of the project will be divided in subdirectories:

    * :ref:`model_specifications`, consisting of assets_baseline.json and endogenous_job_finding_duration.json each containing the parameter values and fixed model specifications used for the analysis. 
    * :ref:`model_code`, consisting of different m-files containing the model code.
    * *Analysis*, is the main file calling the model specifications and parameters to run the analysis. It can be found in this documentation in :ref:`analysis`.
    * *Final*, consists of a graph_equivalents.m file, which creates the plots. It can be found in this documentation in :ref:`final`.
    * *Paper*, consisting of a tex-document, which can be found in this documentation in :ref:`paper`.

.. _projectpaths:

Project paths:
--------------

In the wscript file in the main directory, all project-wide paths are defined, making the definition of paths easier and less problematic if our project is run on different machines. 
The following code is taken from the above-mentioned file.

.. literalinclude:: ../../wscript
    :start-after: out = 'bld'
    :end-before:     # Convert the directories into Waf nodes. 