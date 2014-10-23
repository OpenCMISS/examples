.. _examples-cellml-fortran:

Fortran introduction to CellML in OpenCMISS-Iron
------------------------------------------------

This example provides an entry level demonstration of creating fields in OpenCMISS-Iron to be defined using CellML models. This demonstration uses Fortran, see :ref:`examples-cellml-python` for the correpsponding Python example.

Following the usual Iron practices, you first need to declare the objects to be used in you application. For CellML, we would normally declare:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :start-after: !DOC-START cellml type declarations
   :end-before: !DOC-END cellml type declarations
   
which declares a single CellML environment that we will use and the fields...