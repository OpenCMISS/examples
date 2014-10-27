.. _examples-finiteElasticity-HomogeneousPipeAxialExtension:

Axial extension in a homogeneous pipe
-------------------------------------

| Example created by: `Jagir Hussan <https://unidirectory.auckland.ac.nz/profile/rjag008>`_.
| Documentation by: `David Nickerson <http://about.me/david.nickerson/>`_.

This example demonstrates the application of an mechanical extension to a homogeneous pipe model along the axis of the pipe. The Mooney-Rivlin constitutive law is used via CellML and PyZinc is used to load the geometric model.

This is a Python example. In order to exectute this example, you will need to have both the Zinc and Iron Python bindings available. This :download:`Cmgui command file <viewAxialStretchResults.com>` can be used to visualise the simulation output. The geometric finite element mesh used in this example was created using PyZinc with the script :download:`quadcylinderthreematerials.py`.

The example is encoded in :download:`HomogeneousPipeAxialExtension.py`. Below we describe some of the key aspects of this example.

Example description
+++++++++++++++++++

In a slight twist on the usual Iron Python scripts, in this example we also make use of PyZinc to load the model so we require the ``exfile`` class from PyZinc as well as the standard Iron Python bindings:

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START imports
   :end-before: #DOC-END imports

As shown on line 1, we assume PyZinc is available to your Python environment. Whereas on line 3 we are specifically searching for the Iron Python bindings using the :envvar:`OPENCMISS_ROOT` environment variable. Using the PyZinc ``exfile`` module we are able to load the finite element model geometry from the ``EXREGION`` file: :download:`hetrogenouscylinder.exregion`, as shown below:

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START load exfile
   :end-before: #DOC-END load exfile

Throughout the remainder of the Python script you can see the data from the now defined ``exregion`` object used to create the geometric mesh via the standard Iron methods.

.. literalinclude:: HomogeneousPipeAxialExtension.py
   :language: python
   :linenos:
   :start-after: #DOC-START create cellml environment 
   :end-before: #DOC-END create cellml environment

In the above code, we use the standard Iron methods to create the :term:`CellML environment` (line 1) and import the CellML model defining the Mooney-Rivlin constitutive law. Lines 8-13 then set the variables from the model that are :term:`known` and lines 15-20 are the :term:`wanted` variables from the CellML model.