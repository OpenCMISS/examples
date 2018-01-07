/*
 * \file
 * \author Chris Bradley
 * \brief This is an example program which sets up a field which uses a more complex mesh using OpenCMISS calls from C.
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss/iron.h"

#define STRING_SIZE 20

#define REGION_USER_NUMBER 1

#define CHECK_ERROR(S) \
  if(Err != CMFE_NO_ERROR) { \
    if(Err == CMFE_ERROR_CONVERTING_POINTER) { \
      fprintf(stderr,"Error: %s: Error converting pointer.\n",(S)); \
    } \
    else if(Err == CMFE_POINTER_IS_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is null.\n",(S)); \
    } \
    else if(Err == CMFE_POINTER_NOT_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is not null.\n",(S)); \
    } \
    else if(Err == CMFE_COULD_NOT_ALLOCATE_POINTER) { \
      fprintf(stderr,"Error: %s: Could not allocate pointer.\n",(S)); \
    } \
    exit(Err); \
  }

int main() 
{
  /* int WorldCoordinateSystemUserNumber;
     int WorldRegionUserNumber; */
  cmfe_CoordinateSystemType WorldCoordinateSystem=(cmfe_CoordinateSystemType)NULL;
  cmfe_RegionType WorldRegion=(cmfe_RegionType)NULL,Region=(cmfe_RegionType)NULL;
  char Label[STRING_SIZE];
  int Err;

  Err = cmfe_CoordinateSystem_Initialise(&WorldCoordinateSystem);
  CHECK_ERROR("Initialising world coordinate system");
  Err = cmfe_Region_Initialise(&WorldRegion);
  CHECK_ERROR("Initialising world region");
  if(cmfe_Initialise(WorldCoordinateSystem,WorldRegion) == CMFE_NO_ERROR)
    {

      Err = cmfe_Region_LabelGet(WorldRegion,STRING_SIZE,Label);
      printf("The world region label is '%s'.\n",Label);

      Err = cmfe_Region_Initialise(&Region);
      Err = cmfe_Region_CreateStart(REGION_USER_NUMBER,WorldRegion,Region);
      Err = cmfe_Region_LabelSet(Region,8,"Testing");
      Err = cmfe_Region_CreateFinish(Region);

      Err = cmfe_Region_LabelGet(Region,STRING_SIZE,Label);	       
      printf("The region label is '%s'.\n",Label);

      Err = cmfe_Region_Finalise(&Region);

      Err = cmfe_Finalise();
    }

  return Err;
}
