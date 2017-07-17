//
//  Initial.cpp
//  FEM
//
//  Created by Sagar Dolas on 24/05/16.
//  Copyright �� 2016 Sagar Dolas. All rights reserved.
//

#include "Initial.hpp"

SIWIR2::Initial::Initial(int argc, const char * argv[]) {
    
    if (argc < 3) {
    	std::cout<<std::endl;
    	std::cout<<"/**************Essential Input Parameters missing********************/ "<<std::endl ;
    	std::cout<<std::endl;
    	std::cout<<"/**************Aborting Program Due to Insufficient Data********************/ "<<std::endl ;
        exit(1) ;

    }
    else if (argc == 3 ){

    	std::cout<<std::endl;
    	std::cout<<"/*********** Refinement Parameter and Input file missing ***********/ "<<std::endl ;
    	std::cout<<"/**** Taking Default Parameter for Refinement and Input file *******/ "<<std::endl ;
    	filename_ ="unit_circle.txt" ;
    	delta = std::stod(argv[1]) ;
    	epsilon = std::stod(argv[2]) ;
    	refinement=  0 ;

    }
    else if (argc == 4){

    	std::cout<<std::endl;
    	std::cout<<"/******** Filename missing ******/"<<std::endl;
    	std::cout<<"/*********Taking Default Parameter for the filename"<<std::endl;
    	filename_ = "unit_circle.txt";
    	delta = std::stod(argv[1]) ;
    	epsilon = std::stod(argv[2]) ;
    	refinement = std::stod(argv[3]) ;
    }
    else {

    	std::cout<<std::endl;
    	std::cout<<"/****** All five parameters provided  ********/"<<std::endl;
    	delta = std::stod(argv[1]) ;
    	epsilon = std::stod(argv[2]) ;
    	refinement = std::stod(argv[3]) ;
    	filename_ = argv[4] ;
    }
}

SIWIR2::Initial::~Initial(){
    
}

void SIWIR2::Initial::filldata(SIWIR2::FEM::DomainData<real_d> &vertex, SIWIR2::FEM::DomainData<real_l> &faces){
    
    std::ifstream file_ ;
    file_.open(filename_) ;

    if (file_.is_open()) {
        
        real_l Nvertex  ;
        file_>>Nvertex ;
        
        numVertex = Nvertex   ;
        
        //Resizing the vector
        vertex.resize(numVertex * SIWIR2::FEM::geometricMap::vertexindex,0.0) ;
        
        std::string tempstring ;
        std::getline(file_, tempstring) ;
        std::getline(file_, tempstring) ;
        
        real_d x = 0.0 , y=0.0 ;
        real_l index = 0 ;
        // Reading all the input data for vertex
        
        for (real_l i = 0 ; i < numVertex; ++i) {
            file_>>index>>x>>y ;
            
            real_l vindex = i * SIWIR2::FEM::geometricMap::vertexindex ;
            vertex[vindex]   = x ;
            vertex[vindex+1] = y ;
        }
        
        real_l Nface ;
        file_>>Nface ;
        
        numFaces = Nface ;
        
        //Resizing the vector
        faces.resize(numFaces * SIWIR2::FEM::geometricMap::faceindex ,0) ;
        
        std::getline(file_, tempstring) ;
        std::getline(file_, tempstring) ;
        real_l v1 =0 , v2 =0 , v3 =0 ;
        
        for (real_l i = 0 ; i < numFaces; ++i) {
            
            file_>>v1>>v2>>v3 ;
            
            real_l findex  = i * SIWIR2::FEM::geometricMap::faceindex ;
            faces[findex] = v1 ;
            faces[findex+1] = v2 ;
            faces[findex+2] = v3 ;
        }
    }
    else {
    	std::cout<<"File not opened"<<std::endl;
    }
    file_.close() ;
}

void SIWIR2::Initial::outOnScreen() {

    std::cout<<"Parameter Specification as follows"<<std::endl;
    std::cout<<"Filename :="<<" "<<filename_<<std::endl;
    std::cout<<"Delta :="<<" "<<delta<<std::endl;
    std::cout<<"Epsilon :="<<" "<<epsilon<<std::endl;
    std::cout<<"Refinement Level :="<<" "<<refinement<<std::endl;
    std::cout<<"Number of Faces:= "<<numFaces<<std::endl;
    std::cout<<"Number of Vertices :="<<numVertex<<std::endl;
    std::cout<<std::endl;

}
