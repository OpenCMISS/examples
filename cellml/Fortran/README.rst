.. _examples-cellml-fortran:

Fortran introduction to CellML in OpenCMISS-Iron
------------------------------------------------

This example provides an entry level demonstration of creating fields in OpenCMISS-Iron to be defined using CellML models. This demonstration uses Fortran, see :ref:`examples-cellml-python` for the corresponding Python example.

Following the usual Iron practices, you first need to declare the objects to be used in you application. For CellML, we would normally declare:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml type declarations
   :end-before: !DOC-END cellml type declarations
   
which declares a single CellML environment that we will use and the fields...

And then we do some stuff:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml setup environment
   :end-before: !DOC-END cellml setup environment

and on line 3 above you can see that this happens, and then on line 6 that is done.

We then do the actual field mapping:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml define field maps
   :end-before: !DOC-END cellml define field maps
   
and here you can see ...

