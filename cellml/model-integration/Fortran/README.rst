.. _examples-cellml-fortran:

Fortran introduction to CellML in OpenCMISS
-------------------------------------------

This example provides an entry level demonstration of creating fields in OpenCMISS to be defined using CellML models. This demonstration uses Fortran, see :ref:`examples-cellml-python` for the corresponding Python example. The source code for this example is available here: :download:`src/FortranExample.f90`.

Following the usual OpenCMISS practices, you first need to declare the objects to be used in you application. For CellML, we would normally declare:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml type declarations
   :end-before: !DOC-END cellml type declarations
   
which declares a single CellML environment that we will use and the fields that we will use with CellML. We also fetch the URL of the CellML model we will load from the command line arguments, if it is provided. Otherwise we fall back on the default :download:`Noble 1998 <n98.xml>` model.

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START get cellml model URL
   :end-before: !DOC-END get cellml model URL

We then create our :term:`CellML environment` and import the model defined above. Having imported the model we are able to flag the parameters of interest, namely the stimulus current and the sodium channel conductance, as :term:`known` variables. We also flag all the membrane currents as :term:`wanted` variables in order to have access to their values from OpenCMISS.

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml setup environment
   :end-before: !DOC-END cellml setup environment

Having flagged the desired variables, we finish the creation of the CellML environment. This triggers the instantion of the CellML model into a computable black box for use by solvers in OpenCMISS control loops.

We are then able to specify the actual field mapping, in this case assuming that we are going to perform an electrophysiology simulation with the CellML model so needing to map the membrane potential variable from the model to the OpenCMISS dependent field.

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START cellml define field maps
   :end-before: !DOC-END cellml define field maps
   
Following the maps we create the CellML fields:

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START create CellML fields
   :end-before: !DOC-END create CellML fields
   
We then define the spatial variation for the stimulus current variable flagged as known above. When setting up this spatial distribution we have to take into account the fact that this example application may be run in parallel and therefore only set values relevant to each specific computational node.

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START define CellML stimulus current
   :end-before: !DOC-END define CellML stimulus current

And finally, in a similar manner as the electrical stimulus above, we define the spatial distribution of the sodium channel conductance variable also flagged as known above.

.. literalinclude:: src/FortranExample.f90
   :language: fortran
   :linenos:
   :start-after: !DOC-START define CellML g_Na
   :end-before: !DOC-END define CellML g_Na
   