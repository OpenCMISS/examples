.. _glossary:

========
Glossary
========

.. glossary::
   :sorted:

   CellML environment
   CellML Environment
      The ``CellML Environment`` is the primary object in the OpenCMISS API which manages the use of CellML models in Iron. A single instance of a CellML environment is able to manage multiple CellML models (and multiple instances of the same model, if required). 
      
   known
   Known variables
      CellML models used in OpenCMISS are required to be mathematically complete and suitably constrained. Often for such models to be useful in simulations, the definition of specific variables from the model(s) must be overridden in order to allow their value to be controlled during the simulation by OpenCMISS fields. Such variables are said to be *known*. When addressing variables from CellML models in the OpenCMISS API we use the standard :term:`CellML variable ID`.
      
   wanted
   Wanted variables
      In order for OpenCMISS to have access to retrive the value of variables from CellML models, those variables must be flagged to OpenCMISS. Such variables are said to be *wanted*. When addressing variables from CellML models in the OpenCMISS API we use the standard :term:`CellML variable ID`.
      
   CellML variable ID
   CellML variable name
      In order to uniquely identify variables from CellML models using a text string, OpenCMISS methods use identifiers of the form: ``component_name/variable_name``. Where ``component_name`` is the value of the name attribute of the component in which the desired variable is located and the ``variable_name`` is the value of that variable's name attribute. Variables in the CellML model which are connected can be addressed by any of the relevant CellML variable ID's.
      
   CellML models field
      The field which associates CellML models within a given :term:`CellML environment` to finite element models in OpenCMISS.
      
   CellML parameters field
      The OpenCMISS field which is used to define spatially varying :term:`known` parameters in the CellML models for a given :term:`CellML environment`.
      
   CellML intermediate field
      The OpenCMISS field which is used to contain the intermediate variables for the CellML models for a given :term:`CellML environment`, these are :term:`wanted` variables which are not state variables in the CellML model(s).