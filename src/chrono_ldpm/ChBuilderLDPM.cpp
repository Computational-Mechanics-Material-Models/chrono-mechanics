// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
//
// =============================================================================
// Authors: Erol Lale
// =============================================================================
// Builder class for LDPM 
//
//  Builder reads node coordinates, element connectivity, facet info from the file created by preprocessor
//
// 
// =============================================================================
#include "chrono/core/ChVector3.h"
#include "chrono/physics/ChSystem.h"
#include "chrono_ldpm/ChBuilderLDPM.h"
#include <initializer_list>

using namespace chrono;
using namespace fea;

namespace chrono {
namespace ldpm {

// ------------------------------------------------------------------
// ChBuilderLDPM
// ------------------------------------------------------------------


void ChBuilderLDPM::read_LDPM_info(std::shared_ptr<ChMesh> my_mesh,  std::shared_ptr<ChMaterialVECT> vect_mat, std::string& LDPM_data_path, 
				std::string& LDPM_GeoName){
	//
	std::vector<std::vector<double>> Mfacet;
	std::vector<std::vector<int>> MtetIDs;
	//std::vector<vertdata> MvertIDs;
	//
	//
	//
	std::string nodesFilename=LDPM_data_path+LDPM_GeoName+"-data-nodes.dat";
    	std::string elemFilename=LDPM_data_path+LDPM_GeoName+"-data-tets.dat";
    	std::string facetFilename=LDPM_data_path+LDPM_GeoName+"-data-facets.dat";    	
    	std::string verticesFilename=LDPM_data_path+LDPM_GeoName+"-data-facetsVertices.dat";
    	//
	std::ifstream nodesFile(nodesFilename);
        std::ifstream elemFile(elemFilename);
	std::ifstream facetFile(facetFilename);	
	std::ifstream verticesFile(verticesFilename);
    ////////////////////////////////////////////////////////////////////////
    //
    // Read node info from Freecad produced "<projectname>-data-nodes.dat" file
	//          ( x, y, z coordinates )
	// And create node instances using ChNodeFEAxyzrot class. 
	// 			( node with 6 deegrees of freedom)
	//
	////////////////////////////////////////////////////////////////////////
	std::vector< std::shared_ptr<ChNodeFEAxyzrot> > LDPMnodeList;	
    if (nodesFile.is_open()) {
       	chrono::ChVector3d pos;
       	double x, y, z;
       	unsigned int idnode=0;
		std::string line;
		
		while (std::getline(nodesFile, line)) {			
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> x >> y >> z;				
				//
				auto mnode= chrono_types::make_shared<ChNodeFEAxyzrot>(chrono::ChFrame<>(chrono::ChVector3d(x, y, z)));
				mnode->SetIndex(idnode);
				LDPMnodeList.push_back(mnode);
                		my_mesh->AddNode(mnode);								
                		++idnode;
			}
		}
    }
    else{
    	throw std::runtime_error("ERROR opening nodes info file: " + std::string(nodesFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
    
    /*
    for(int i=0; i<my_mesh->GetNnodes(); i++){
    
    	auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i));
    	std::cout<< i <<". node: "<<node->Frame().GetPos().x()<<"\t"<<node->Frame().GetPos().y()<<"\t"<<node->Frame().GetPos().z()<<std::endl;
    
    }
	*/

	/*
	auto msection = chrono_types::make_shared<ChSectionLDPM>();
	msection->Set_material(vect_mat);
	msection->SetDrawThickness(2, 2);
	*/	
    ////////////////////////////////////////////////////////////////////////
    //
    // Read element connectivity info from Freecad produced "<projectname>-data-edges.dat" file
	//          ( node1, node2 )
	// And create element instances using ChElementLDPM class. 
	// 			( LDPM beam element with )
	//
	////////////////////////////////////////////////////////////////////////
    if (elemFile.is_open()) {       	
       	int elid;
       	int node1, node2, node3, node4;
       	unsigned int idelem=0;
		std::string line;
    	while (std::getline(elemFile, line)) {
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> node1 >> node2 >> node3 >> node4;				
				auto nodeI = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node1));
				auto nodeJ = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node2));
				auto nodeK = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node3));
				auto nodeL = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node4));
				auto mel = chrono_types::make_shared<ChElementLDPM>();
				mel->SetNodes(nodeI, nodeJ, nodeK, nodeL);
				//mel->SetSection(msection);
				my_mesh->AddElement(mel);
                ++idelem;

			}
		}
    }
    else{
    	throw std::runtime_error("ERROR opening element info file: " + std::string(elemFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
	
	int elementnum=my_mesh->GetNumElements ();
	std::cout<<"elementnum "<<elementnum<<std::endl;
	
	////////////////////////////////////////////////////////////////////////
    //
    // Read coordinate infor of facet vertices from Freecad produced "<projectname>-data-facetsVertices.dat" file
	//       ( vertex_X, vertex_Y, vertex_Z )
	//
	////////////////////////////////////////////////////////////////////////   	
	//double Mvertices[my_mesh->GetNumElements ()*36][3];
	//double Mvertices[2000*36][3];
	std::vector<ChVector3d> Mvertices;
	if (verticesFile.is_open()) {
       	//
       	double x, y, z;
		std::string line;
		while (std::getline(verticesFile, line)) {
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> x >> y >> z ;	
				Mvertices.push_back(ChVector3d(x, y, z));
			}
		}		
    }
    else{
    	throw std::runtime_error("ERROR opening vertices file: " + std::string(verticesFilename) + "\n");
    	exit(EXIT_FAILURE);
    }	
    
    /*
     for(int i=0; i<my_mesh->GetNumElements ()*36; i++){
    
    	
    	std::cout<< i <<". vertice: "<<Mvertices[i].x()<<"\t"<<Mvertices[i].y()<<"\t"<<Mvertices[i].z()<<std::endl;
    
    }
    */
    	
	////////////////////////////////////////////////////////////////////////
    //
    // Read facet info from Freecad produced "<projectname>-data-facets.dat" file
	//       ( Tet IDx IDy IDz Vol pArea cx cy cz px py pz qx qy qz sx sy sz mF )
	//
	////////////////////////////////////////////////////////////////////////  
	if (facetFile.is_open()) {
       	//
		double val;
		//
		double area=0;
		int matFlag=0;
		chrono::ChVector3d  facetC;
		chrono::ChMatrix33<double> nmL;				
		
		//
		std::string line;		
		while (std::getline(facetFile, line)) {
			//std::cout<< line[0] <<std::endl;

			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				std::vector<double> facetvec;
				while(sline>>val){
					facetvec.push_back(val);
				}
												
				// Element IDs 
				int iel=int(facetvec[0]);	
							
				//
				// Get element object and create a section object for it
				//
				auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(iel));
				//
				auto msection = chrono_types::make_shared<ChSectionLDPM>();				
				//
				area=facetvec[5];
				//
				facetC[0]=facetvec[6];
				facetC[1]=facetvec[7];
				facetC[2]=facetvec[8];
				//
				nmL.setZero();
				nmL(0,0)=facetvec[9];
				nmL(0,1)=facetvec[10];
				nmL(0,2)=facetvec[11];
				//
				nmL(1,0)=facetvec[12];
				nmL(1,1)=facetvec[13];
				nmL(1,2)=facetvec[14];
				//
				nmL(2,0)=facetvec[15];
				nmL(2,1)=facetvec[16];
				nmL(2,2)=facetvec[17];

				matFlag=int(facetvec[18]);
				
				msection->Set_material(vect_mat);
				msection->Set_area(area);
				msection->Set_center(facetC);
				msection->Set_facetFrame(nmL);				
				elem->AddFacetI(msection);							
				//
				// Store vertices coordinate info for mass calculation
				//
				elem->AddVertNodeVec(std::initializer_list<ChVector3d>{Mvertices[int(facetvec[1])], Mvertices[int(facetvec[2])], Mvertices[int(facetvec[3])]});
			}	
		}
    } else {
    	throw std::runtime_error("ERROR opening facet info file: " + std::string(facetFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
	Mvertices.clear();
}

void ChBuilderLDPM::read_LDPM_info(std::shared_ptr<ChMesh> my_mesh,  std::vector<std::shared_ptr<ChMaterialVECT>> vect_mat, std::string& LDPM_data_path, 
				std::string& LDPM_GeoName){
	//
	std::vector<std::vector<double>> Mfacet;
	std::vector<std::vector<int>> MtetIDs;
	//std::vector<vertdata> MvertIDs;
	//
	//
	//
	std::string nodesFilename=LDPM_data_path+LDPM_GeoName+"-data-nodes.dat";
    	std::string elemFilename=LDPM_data_path+LDPM_GeoName+"-data-tets.dat";
    	std::string facetFilename=LDPM_data_path+LDPM_GeoName+"-data-facets.dat";    	
    	std::string verticesFilename=LDPM_data_path+LDPM_GeoName+"-data-facetsVertices.dat";
    	//
	std::ifstream nodesFile(nodesFilename);
        std::ifstream elemFile(elemFilename);
	std::ifstream facetFile(facetFilename);	
	std::ifstream verticesFile(verticesFilename);
    ////////////////////////////////////////////////////////////////////////
    //
    // Read node info from Freecad produced "<projectname>-data-nodes.dat" file
	//          ( x, y, z coordinates )
	// And create node instances using ChNodeFEAxyzrot class. 
	// 			( node with 6 deegrees of freedom)
	//
	////////////////////////////////////////////////////////////////////////
	std::vector< std::shared_ptr<ChNodeFEAxyzrot> > LDPMnodeList;	
    if (nodesFile.is_open()) {
       	chrono::ChVector3d pos;
       	double x, y, z;
       	unsigned int idnode=0;
		std::string line;
		
		while (std::getline(nodesFile, line)) {			
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> x >> y >> z;				
				//
				auto mnode= chrono_types::make_shared<ChNodeFEAxyzrot>(chrono::ChFrame<>(chrono::ChVector3d(x, y, z)));
				mnode->SetIndex(idnode);
				LDPMnodeList.push_back(mnode);
                		my_mesh->AddNode(mnode);								
                		++idnode;
			}
		}
    }
    else{
    	throw std::runtime_error("ERROR opening nodes info file: " + std::string(nodesFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
    
    
    ////////////////////////////////////////////////////////////////////////
    //
    // Read element connectivity info from Freecad produced "<projectname>-data-edges.dat" file
	//          ( node1, node2 )
	// And create element instances using ChElementLDPM class. 
	// 			( LDPM beam element with )
	//
	////////////////////////////////////////////////////////////////////////
    if (elemFile.is_open()) {       	
       	int elid;
       	int node1, node2, node3, node4;
       	unsigned int idelem=0;
		std::string line;
    	while (std::getline(elemFile, line)) {
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> node1 >> node2 >> node3 >> node4;				
				auto nodeI = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node1));
				auto nodeJ = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node2));
				auto nodeK = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node3));
				auto nodeL = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(node4));
				auto mel = chrono_types::make_shared<ChElementLDPM>();
				mel->SetNodes(nodeI, nodeJ, nodeK, nodeL);
				//mel->SetSection(msection);
				my_mesh->AddElement(mel);
                ++idelem;

			}
		}
    }
    else{
    	throw std::runtime_error("ERROR opening element info file: " + std::string(elemFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
	
	int elementnum=my_mesh->GetNumElements ();
	std::cout<<"elementnum "<<elementnum<<std::endl;
	
	
	
	
	////////////////////////////////////////////////////////////////////////
    //
    // Read coordinate infor of facet vertices from Freecad produced "<projectname>-data-facetsVertices.dat" file
	//       ( vertex_X, vertex_Y, vertex_Z )
	//
	////////////////////////////////////////////////////////////////////////   	
	//double Mvertices[my_mesh->GetNumElements ()*36][3];
	//double Mvertices[2000*36][3];
	std::vector<ChVector3d> Mvertices;
	if (verticesFile.is_open()) {
       	//
       	double x, y, z;
		std::string line;
		while (std::getline(verticesFile, line)) {
			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				sline>> x >> y >> z ;	
				Mvertices.push_back(ChVector3d(x, y, z));
			}
		}		
    } else {
    	throw std::runtime_error("ERROR opening vertices file: " + std::string(verticesFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
    	
	////////////////////////////////////////////////////////////////////////
    //
    // Read facet info from Freecad produced "<projectname>-data-facets.dat" file
	//       ( Tet IDx IDy IDz Vol pArea cx cy cz px py pz qx qy qz sx sy sz mF )
	//
	////////////////////////////////////////////////////////////////////////  
	if (facetFile.is_open()) {
       	//
		double val;
		//
		double area=0;
		int matFlag=0;
		chrono::ChVector3d  facetC;
		chrono::ChMatrix33<double> nmL;				
		
		//
		std::string line;		
		while (std::getline(facetFile, line)) {
			//std::cout<< line[0] <<std::endl;

			//
			if( line[0]!='/') {				
				std::istringstream sline(line);
				std::vector<double> facetvec;
				while(sline>>val){
					facetvec.push_back(val);
				}
												
				// Element IDs 
				int iel=int(facetvec[0]);	
							
				//
				// Get element object and create a section object for it
				//
				auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(iel));
				
				//
				auto msection = chrono_types::make_shared<ChSectionLDPM>();				
				//
				area=facetvec[5];
				//
				facetC[0]=facetvec[6];
				facetC[1]=facetvec[7];
				facetC[2]=facetvec[8];
				//
				nmL.setZero();
				nmL(0,0)=facetvec[9];
				nmL(0,1)=facetvec[10];
				nmL(0,2)=facetvec[11];
				//
				nmL(1,0)=facetvec[12];
				nmL(1,1)=facetvec[13];
				nmL(1,2)=facetvec[14];
				//
				nmL(2,0)=facetvec[15];
				nmL(2,1)=facetvec[16];
				nmL(2,2)=facetvec[17];

				matFlag=int(facetvec[18]);
				
				msection->Set_material(vect_mat[matFlag]);
				msection->Set_area(area);
				msection->Set_center(facetC);
				msection->Set_facetFrame(nmL);				
				elem->AddFacetI(msection);							
				//
				// Store vertices coordinate info for mass calculation
				//	
				elem->AddVertNodeVec(std::initializer_list<ChVector3d>{Mvertices[int(facetvec[1])], Mvertices[int(facetvec[2])], Mvertices[int(facetvec[3])]});
			}
		}
    } else {
    	throw std::runtime_error("ERROR opening facet info file: " + std::string(facetFilename) + "\n");
    	exit(EXIT_FAILURE);
    }
	Mvertices.clear();
}



}  // end namespace ldpm
}  // end namespace chrono
