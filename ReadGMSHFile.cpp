/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/*
 * ReadGMSHFile.cc -- subroutine of the RWGSurface class constructor
 *
 * homer reid    -- 3/2007 
 * -------------------------------------------------------------------------------
 * This function has been adpated for MC-RayTracing code
 * Francisco Ramirez -- 10/2018
 * -------------------------------------------------------------------------------
 */
/*************************************************************/
/* constants needed in this file only  ***********************/
/*************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#define TYPE_TRIANGLE 2
#define TYPE_POINT    15

#define NODE_START_KEYWORD1 "$NOD"
#define NODE_START_KEYWORD2 "$Nodes"

#define NODE_END_KEYWORD1   "$ENDNOD"
#define NODE_END_KEYWORD2   "$EndNodes"

#define ELM_START_KEYWORD1  "$ELM"
#define ELM_START_KEYWORD2  "$Elements"

#define ELM_END_KEYWORD1    "$ENDELM"
#define ELM_END_KEYWORD2    "$EndElements"

#define FORMAT_LEGACY 0
#define FORMAT_NEW    1

#define MAXREFPTS 100

#include "mcphoton_lib.h"

/*************************************************************/
/* Read vtx and panels from a GMSH .msh file to specify */
/* a surface.                                                */
/* If MeshTag==-1, then all panels are read.          */
/* Otherwise, only panels on the specified physical region   */
/* are read.                                                 */
/*************************************************************/
void Surface::ReadGMSHFile(FILE *MeshFile, const string &FileName){

  //Panel *P;
  char buffer[100];
  int VI[3]; //, Temp;
  int *GMSH2HR;
  int jGMSH;
  int i, j, nv, ne; // , np;
  int LineNum, NumElements, NodeIndex;
  int ElNum, ElType, RegPhys, RegElem, NodeCnt, nConv;

  int NumRefPts;

  int WhichMeshFormat, KeywordFound, nt, nTags, nRead, iDummy, bufPos;

  /* stuff to figure out correct orientation of panel normals */
  //double dRP, dRPMin;
  //double *RP, *RPMin;
  //double CentroidDisplaced[3];
  //int nrp;
  //int RefPntIndices[MAXREFPTS];
 
  /*------------------------------------------------------------*/
  /* Read lines until we hit the keyword indicating the start   */
  /* of the 'nodes' section, then read the number of nodes.     */
  /*------------------------------------------------------------*/
  LineNum=0;
  KeywordFound=0;
  while(KeywordFound==0)
   { 
     if (!fgets(buffer,100,MeshFile))
		 ErrorMsg(FileName + ": failed to find node start keyword");
     LineNum++;
     if( !strncmp(buffer,NODE_START_KEYWORD1,strlen(NODE_START_KEYWORD1)))
      { WhichMeshFormat=FORMAT_LEGACY; KeywordFound=1; }
     else if( !strncmp(buffer,NODE_START_KEYWORD2,strlen(NODE_START_KEYWORD2)))
      { WhichMeshFormat=FORMAT_NEW; KeywordFound=1; }
   };
  if (!fgets(buffer, 100, MeshFile) || !(nVertex = atoi(buffer)))
	  ErrorMsg(FileName + ": invalid number of nodes");
  LineNum++;

  /*------------------------------------------------------------*/
  /*- Read in the vtx (which GMSH calls 'nodes.')          */
  /*- Note that the numbering of the vtx in GMSH does not  */
  /*- necessarily correspond to their ordering in the mesh      */
  /*- file. To remedy this situation, we construct a mapping    */
  /*- between GMSH's vtx indices and our internal vertex   */
  /*- indices, which works like this: The vertex that GMSH      */
  /*- calls 'node 3' is stored in slot GMSH2HR[3] within our    */
  /*- internal vtx array.                                  */ 
  /*------------------------------------------------------------*/
  vtx = new Point3D [nVertex];
  GMSH2HR= new int [2*nVertex];
  for(i=0; i<2*nVertex; i++)
   GMSH2HR[i]=-1;
  for (nv=0; nv<nVertex; nv++)
   { if (!fgets(buffer,100,MeshFile))
	  ErrorMsg(FileName + ":  too few nodes");
     LineNum++;
	 nConv = sscanf(buffer, "%i %le %le %le", &NodeIndex,
		 &vtx[nv].x, &vtx[nv].y, &vtx[nv].z);
          //                vtx+3*nv,vtx+3*nv+1,vtx+3*nv+2); 
	 if (nConv != 4)
		 ErrorMsg(FileName + ": " + to_string(LineNum) +": invalid node specification");
     if (NodeIndex>2*nVertex)
		 ErrorMsg(FileName + ": internal error in ReadGMSHFile");
     GMSH2HR[NodeIndex]=nv;
   }; /* for (nv=0; nv<nVertex; nv++) */
 
  /*------------------------------------------------------------*/
  /*- Eliminate any redundant vtx from the vertex list.   -*/
  /*- Explain me in greater detail please.                     -*/
  /*------------------------------------------------------------*/
  int NumRedundantvtx=0;
  for (i = 0; i < nVertex; i++)
	  for (j = i + 1; j < nVertex; j++) {
		  Vector3D VTest(vtx[i] - vtx[j]);
		  if (VTest.magnitude() < 1.0e-6)
		  {
			  /* remap all references to my node j so that they now refer to my node i*/
			  for (jGMSH = 0; jGMSH < 2 * nVertex; jGMSH++)
				  if (GMSH2HR[jGMSH] == j)
					  GMSH2HR[jGMSH] = i;

			  NumRedundantvtx++;
			  //fprintf(stderr,"\n*\n* redundant nodes found!!(%i,%i)\n*\n",i,j);
		  };
	  }
  /*------------------------------------------------------------*/
  /* Confirm that the next two lines in the file are the        */
  /* end-of-nodes-section keyword and the start-of-elements-section */
  /* keyword, then read the number of elements.                 */
  /*------------------------------------------------------------*/
  if ( !fgets(buffer,100,MeshFile) )
	  ErrorMsg(FileName + ":  bad file format (nodes section not terminated)");
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(buffer,NODE_END_KEYWORD1,strlen(NODE_END_KEYWORD1)))
	  ErrorMsg(FileName + ": " + to_string(LineNum) + ": unexpected keyword");
   }
  else
   { if ( strncmp(buffer,NODE_END_KEYWORD2,strlen(NODE_END_KEYWORD2)))
	  ErrorMsg(FileName + ": " + to_string(LineNum) + ": unexpected keyword");
   };

  if ( !fgets(buffer,100,MeshFile) )
	  ErrorMsg(FileName + ": bad file format (elements section not initiated)");
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(buffer,ELM_START_KEYWORD1,strlen(ELM_START_KEYWORD1)))
	  ErrorMsg(FileName + ": " + to_string(LineNum) + ": unexpected keyword");
   }
  else
   { if ( strncmp(buffer,ELM_START_KEYWORD2,strlen(ELM_START_KEYWORD2)))
	  ErrorMsg(FileName + ": " + to_string(LineNum) + ": unexpected keyword");
   }

  if ( !fgets(buffer,100,MeshFile))
	  ErrorMsg(FileName + ": bad file format (invalid number of elements)");
  LineNum++;
  nConv=sscanf(buffer,"%i",&NumElements);
  if (nConv != 1 || NumElements < 0)
	  ErrorMsg(FileName + ": " + to_string(LineNum) + ": invalid number of elements");
 
  /*------------------------------------------------------------*/
  /*- Now read each line of the elements section.               */ 
  /*------------------------------------------------------------*/
  NumPanels=NumRefPts=0;
  Panels = new Panel [NumElements];
  for (ne=0; ne<NumElements; ne++)
   { 
     if (!fgets(buffer,100,MeshFile))
		 ErrorMsg(FileName + ": too few elements in input file");
     LineNum++;

     if (WhichMeshFormat==FORMAT_LEGACY)
      { 
        nConv=sscanf(buffer,"%i %i %i %i %i %i %i %i",
                     &ElNum,&ElType,&RegPhys,&RegElem,&NodeCnt,VI,VI+1,VI+2);
        if (nConv<5)
			ErrorMsg(FileName + ": " + to_string(LineNum) + ": invalid element specification");
      }
     else
      { 
        nConv=sscanf(buffer,"%i %i %i%n",&ElNum,&ElType,&nTags,&nRead);
        bufPos=nRead;
        if (nConv<3)
			ErrorMsg(FileName + ": " + to_string(LineNum) + ": invalid element specification");

        // read the first 'tag,' which should be the physical region
        if (nTags==0)
         RegPhys=0;
        else 
         sscanf(buffer+bufPos,"%i%n",&RegPhys,&nRead);
        bufPos+=nRead;

        // read any remaining tags 
        for(nt=0; nt<nTags-1; nt++)
         { sscanf(buffer+bufPos,"%i%n",&iDummy,&nRead);
           bufPos+=nRead;
         };

        // finally, read the vertex indices. 
        nConv=sscanf(buffer+bufPos,"%i %i %i",VI,VI+1,VI+2);

        if (ElType==TYPE_TRIANGLE && nConv!=3)
			ErrorMsg(FileName + ": " + to_string(LineNum) + ": invalid element specification");
        else if (ElType==TYPE_POINT && nConv!=1)
				ErrorMsg(FileName + ": " + to_string(LineNum) + ": invalid element specification");
      };
   
     /* we only process elements that are points or triangles */
     switch(ElType)
      { 
        /***************************************************************/
        /* add new reference point to list of reference points *********/
        /***************************************************************/
#if 0
20140921 I am eliminating this, as we have no more use for 
reference points in the msh file, and this code was causing 
problems in mesh files that had large numbers of physical points
defined
        case TYPE_POINT:
          if (NumRefPts==MAXREFPTS)
           ErrExit("%s:%i: too many reference points",FileName,LineNum); 
          RefPntIndices[NumRefPts++]=GMSH2HR[ VI[0] ];
          break;
#endif

        /***************************************************************/
        /* add new triangle to list of panels                          */
        /***************************************************************/
        case TYPE_TRIANGLE:
			if ((MeshTag == -1) || (MeshTag == RegPhys))
			{
				Panels[NumPanels].SetCoord(	vtx[GMSH2HR[VI[0]]],
											vtx[GMSH2HR[VI[1]]],
											vtx[GMSH2HR[VI[2]]]);
			  /*Panels[NumPanels]=NewPanel(vtx, GMSH2HR[ VI[0] ], GMSH2HR[ VI[1] ], 
                                           GMSH2HR[ VI[2] ])*/;
             Panels[NumPanels].Index=NumPanels;
             NumPanels++;
           };
          break;

        default:
          //ErrExit("%s:%i: unknown element type %i",FileName,LineNum,ElType);
          //fprintf(stderr,"%s:%i: warning: ignoring unknown element type %i\n",FileName,LineNum,ElType);
          break;

      }; // switch(ElType)
   }; //for (ne=0; ne<NumElements; ne++)
  free(GMSH2HR);
  fclose(MeshFile);
}