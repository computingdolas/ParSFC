//
//  main.cpp
//  FEM
//
//  Created by Sagar Dolas on 24/05/17.
//  Copyright �� 2016 Sagar Dolas. All rights reserved.
//

// Standard Template Library 
#include <iostream>
#include <fstream>
#include <chrono>

// Local Library 
#include "Initial.hpp"
#include "functions.hpp"
#include "SolverFunction.hpp"
#include "SolverFunction.cpp"
#include "Assemble_Matrix.hpp"
#include "Assemble_Matrix.cpp"
#include "writeMatrix.hpp"
#include "writeMatrix.cpp"
#include "Compresed_Row_Storage.hpp"
#include "Compresed_Row_Storage.cpp"
#include "Grid_Refinement.hpp"
#include "Grid_Refinement.cpp"

using namespace SIWIR2 ;

int main(int argc, const char * argv[]) {
    
    // Declaring the buffer
    FEM::DomainData<real_d> vertex(FEM::dataType::vertex) ;
    FEM::DomainData<real_l> face(FEM::dataType::face) ;
    
    Initial input(argc,argv) ;
    input.filldata(vertex, face) ;
    input.outOnScreen() ;


     // Making the Ksqaure
    // std::cout<<std::endl;
    //  std::cout<<"////////////////////////////////////////"<<std::endl;;
    //  logonScreen("Initialising the ksq");
    //  FEM::DomainData<real_d> ksq(vertex,input.delta);
    // 	std::ofstream out ;
    //  out.open("ksq.txt") ;
    //  printdata(ksq, vertex, out) ;
    //  out.close() ;
    //  logonScreen("Written into ksq.txt") ;
    // std::cout<<"////////////////////////////////////////"<<std::endl;;
    
    //logonScreen("");
    //Delcaring the Refinement buffer
    // SIWIR2::FEM::Refinement<DomainData<real_d>, DomainData<real_l>> R(vertex,input.refinement) ;
    // FEM::DomainData<real_l> facenew(R.refinedGrid(face)) ;

    // std::cout<<std::endl;
    // std::cout<<"////////////////////////////////////////"<<std::endl;;
    // std::cout<<"Domain Data after refinement "<<std::endl;
    // std::cout<<"Number of vertex := "<<vertex.size()/2<<std::endl;
    // std::cout<<"Number of Face   := "<<facenew.size()/3<<std::endl;
    // std::cout<<"////////////////////////////////////////"<<std::endl<<std::endl;
    
    // std::ofstream outfile ; 
    // outfile.open("unit_circle_refine6.txt") ; 
    
    // outfile << vertex.size()/2 <<"\t vertices in the domain" << std::endl ; 
    // outfile << "Index of vertex | x cordinate | y cordiate "<< std::endl ; 

    // // Writing the vertices 
    // size_t index = 0 ; 
   	// for (size_t i = 0 ; i < vertex.size() ; i+=2 ){

   	// 	index  = i / 2 ; 
   	// 	outfile <<"\t"<< index <<"\t"<<vertex[i]<<"\t"<<vertex[i+1]<< std::endl ;
   	// }

    // outfile<< facenew.size() / 3 <<" \t elements in the domain " << std::endl ; 
    // outfile<<"\tindex of vertex 0 | index of vertex 1 | index of vertex 2 "<< std::endl ;
   	// // Writing the faces 
   	// for ( size_t i = 0 ; i < facenew.size() ; i+=3 ){
   	// 	outfile <<"\t"<<facenew[i] <<"\t "<< facenew[i+1] <<"\t"<< facenew[i+2] << std::endl ; 
   	// }

   	// outfile.close() ;     

    //Storing the matrix map
    SIWIR2::FEM::MatrixMap<real_d> A((vertex.size()/2)) ;
    SIWIR2::FEM::MatrixMap<real_d> M((vertex.size()/2)) ;  
	
    std::cout << "Starting Assembly " << std::endl ;     
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    //for (int i =0 ; i < 100 ; ++i){
    	SIWIR2::FEM::Assemble(A,M,vertex,face,(vertex.size()/2),(face.size()/3),input.delta) ;
    //}
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    
    std::cout << "Assembly took =; " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<" microseconds"<< std::endl ; 
    
//     //Writing the data out
//     SIWIR2::FEM::writeMatrix<real_d> w(A) ;
//     SIWIR2::FEM::writeMatrix<real_d> m(M) ;

//     out.open("A.txt") ;
//     w.write_to_file(out) ;
//     out.close() ;
//     logonScreen("Matrix A is written in A.txt") ;

//     out.open("M.txt") ;
//     m.write_to_file(out) ;
//     out.close() ;
//     logonScreen("Matrix M is written in M.txt") ;
// // 
    //Compressed Row Storage Format
    SIWIR2::FEM::Compressed_Row_Storage<real_d, MatrixMap<real_d>> crsA(A) ;
    crsA.StoreCRSFormat();
    SIWIR2::FEM::Compressed_Row_Storage<real_d, MatrixMap<real_d>> crsM(M) ;
    crsM.StoreCRSFormat() ;
    
    // Initialise u for the inverse power iteration
    FEM::DomainData<real_d> u((vertex.size()/2),FEM::dataType::normal) ;
    real_d constant = 1.0 / std::sqrt(u.size()/2) ;
 
    
    // Starting the solver
    SIWIR2::FEM::Solver<Compressed_Row_Storage<real_d, MatrixMap<real_d>>> s(crsA, crsM,u, input.epsilon) ;
    
   // for(int i = 0 ; i < 200 ; ++i){
        
        for (real_l i =0 ; i < u.size(); ++i) {
            u[i] =  1.0  * constant ;
    }
    start = std::chrono::steady_clock::now();
   
        s.InversePowerIteration() ;
   // }
    end = std::chrono::steady_clock::now();

    std::cout << "Solver took =; " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<" microseconds"<< std::endl ; 
// // 
// //     
// //     out.open("Eigenmode.txt") ;
// //     printdata(u, vertex, out) ;
// //     out.close() ;
// //     logonScreen("Written into eigenmode.txt") ;
// //     
// //     out.open("Lamda.txt") ;
// //     printdata(s.lambda, out) ;
// //     out.close() ;
// //     logonScreen("Written into Lambda.txt") ;

    return 0;
}
