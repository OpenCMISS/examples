.. _examples-finiteElasticity-HomogeneousPipeAxialExtension:

Axial extension in a homogeneous pipe
-------------------------------------

| **Example created by:** `Jagir Hussan <https://unidirectory.auckland.ac.nz/profile/rjag008>`_
| **Documentation by:** `David Nickerson <http://about.me/david.nickerson/>`_

This example demonstrates the application of an mechanical extension to a homogeneous tissue cylinder model along the axis of the cylinder. The Mooney-Rivlin constitutive law is used via CellML.

This is a Python example. In order to exectute this example, you will need to have the Iron Python bindings available. This :download:`Cmgui command file <viewAxialStretchResults.com>` can be used to visualise the simulation output. The geometric finite element mesh used in this example was created using PyZinc with the script :download:`quadcylinderthreematerials.py`.

The example is encoded in :download:`HomogeneousPipeAxialExtension.py`. Below we describe some of the key aspects of this example.

.. contents::

Example description
+++++++++++++++++++

In this example the The Python script :download:`exfile.py` is used to load the geometry from a file in the Zinc `EX format <http://www.cmiss.org/cmgui/wiki/TheCmguiEXFormatGuideExnodeAndExelemFiles>`_.

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START imports
   :end-before: #DOC-END imports

The :file:`exfile.py` script in needs to be in the current folder or available in your Python environment for the ``import`` on line 1 to succeed. Whereas on line 3 we are specifically adding the location we expect to find the Iron Python bindings using the :envvar:`OPENCMISS_ROOT` environment variable. Using the ``exfile`` module we are able to load the finite element model geometry from the ``EXREGION`` file: :download:`hetrogenouscylinder.exregion`, as shown below:

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START load exfile
   :end-before: #DOC-END load exfile

Throughout the remainder of the Python script you can see the data from the now defined ``exregion`` object used to create the geometric mesh via the standard Iron methods. For example, in this code:

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: # DOC-START define node coordinates
   :end-before: # DOC-END define node coordinates
   
we set the nodal coordinates and derivatives in the models geometric field. There is also code to keep track of nodes located at either end of the vessel, which will be used later when defining boundary conditions.

CellML
******

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START create cellml environment 
   :end-before: #DOC-END create cellml environment

In the above code, we use the standard Iron methods to create the :term:`CellML environment` (line 1) and import the CellML model defining the Mooney-Rivlin constitutive law. Lines 8-13 then set the variables from the model that are :term:`known` and lines 15-20 are the :term:`wanted` variables from the CellML model.

Results
+++++++

.. figure:: doc/results.png
   :align: center
   :width: 40%
   :figwidth: 80%
   
   **Figure:** Results from running this pipe extension simulation. The gold lines show the original, undeformed, cylinder geometry. The coloured lines show the deformed geometry, with the colour varying to show the difference in strain through the wall of the cylinder. The cones represent the three normalised principal strains at material points throughout the tissue volume (red for compression and blue for extension).