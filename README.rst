=============================
paper-transfocators-resources
=============================

Select slit opening
===================

Select the aperture to match the desired coherent fraction:

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1cf.png

(f1,f2) trajectories
====================

For a given f1, calculate f2 in order to have the focus at the sample. 

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1traj_comparison1.png

Focal sizes
===========

This are the data from the two methods, each method use its own (f1,f2) trajectories, that are not the same. 

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1size_comparison1.png

To compare the methods I calculates the sizes using the same (f1,f2) trajectory (the one from WOFRY)

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1size_comparison2.png



TODO
====
- Double-check the papameters used in GSMM and WOFRY
- Recalculate using rectangular slit (Gaussian used) and real lenses (ideal used) and undulator coherent mode decomposition
- Calculate all wanted energies
- compute transmission (intensity at center and integrated intensity)
- Start thinking in 2D simulations and compare

DONE
====
- Obtain GSM trajectories (f1,f2) and corresponding size for better comparison
