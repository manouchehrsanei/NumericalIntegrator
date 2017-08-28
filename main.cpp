
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"

#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"


#include "TPZVTKGeoMesh.h"


#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticprism.h"

#include "tpzarc3d.h"
#include "tpzellipse3d.h"
#include "tpzgeoblend.h"

#include "pzvec.h"
#include "pzquad.h"
#include "tpzintrulep3d.h"





// ************************ (Create a linear Omega description : 1D) ********************************************

TPZGeoMesh *CreateOmega1D(long num_el, REAL a1, REAL a2);




// ************************************** (Create linear meshes: 2D) ********************************************


TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly);

TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly);


// ************************************* (Create linear meshes: 3D) ********************************************


TPZGeoMesh *CreateThreeDHexPrytetGMesh(long nnodesthrhpt, REAL Lx, REAL Ly, REAL Lz);







// Example Functions

REAL f_1D(TPZManVector<REAL,3> &x);


// ******************************************** (main of program) ***********************************************


int main() {

    int order = 10;
    long num_el = 10;
    REAL a1 = 0.0;
    REAL a2 = 1.0;
    TPZGeoMesh * Omega = CreateOmega1D(num_el, a1, a2);
//    TPZGeoMesh * Omega = CreateOmega2D(num_el, a1, a2, b1, b2);
    long n_elements = Omega->NElements();

    REAL NumericalIntegral = 0.0;
    
    for (long iel = 0; iel < n_elements; iel++) { // for a all element in the mesh
    
        TPZGeoEl * element = Omega->Element(iel);
        int itself = element->NSides() - 1;
        TPZIntPoints * IntegrationRule = element->CreateSideIntegrationRule(itself, order);
        int n_points = IntegrationRule ->NPoints();
        
        TPZManVector<REAL,3> x;
        TPZManVector<REAL,3> par_space;
        REAL w, detjac;
        TPZFMatrix<REAL> jac, jacinv, axes;
        
        for (int ip = 0; ip < n_points; ip++) {
            IntegrationRule->Point(ip, par_space, w);
            
            element->X(par_space, x);
            element->Jacobian(par_space, jac, axes, detjac, jacinv);
            
            NumericalIntegral += w * detjac * (1.0);
        }
    }
        
    std::cout << "Numerical integration = " << NumericalIntegral << std::endl;
        
    
    
    return 0;
}




// ++++++++++++++++++++++++++++++++++++++++++++ functions +++++++++++++++++++++++++++++++++++++++++++++++++++

REAL f_1D(TPZManVector<REAL,3> &x){
    REAL x1 = x[0];
    REAL f = (1.0 - x1)*x1*sin(2.0*x1);
    return f;
}


// ********************************* (linear elements) ******************************************************


// ************************************** Create 1D linear meshes ***************************************

TPZGeoMesh *CreateOmega1D(long num_el, REAL a1, REAL a2)
{
    REAL size_el = (a2-a1)/REAL(num_el);
    
    TPZGeoMesh * gmesh_OneDL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("geomesh OneDL"); // geometry name
    
    gmesh_OneDL->SetName(name);
    gmesh_OneDL->SetDimension(geometry_dim);
    
    long num_nodes = num_el + 1; // Number of the nodes
    gmesh_OneDL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
    
    
    const int physical_id = 1; // Define id for material
    
    for (long i = 0 ; i < num_nodes; i++)
    {
        const REAL valElem = i * size_el;
        TPZVec <REAL> coord(3,0.);
        coord[0] = valElem;
        gmesh_OneDL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
        gmesh_OneDL->NodeVec()[i].SetNodeId(i); // The id identification
    }
    
    // Creating linear element and  zero-dimensional boundary element
    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
    long element_id = 0;

    
    // Elements
    
    for (long iel = 0; iel < num_el; iel++)
    {
        const long inod_l = iel;
        const long inod_r = iel + 1;
        Linear_topology[0] = inod_l;
        Linear_topology[1] = inod_r;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDL);

    }
    
    element_id++;
    
    
    gmesh_OneDL->BuildConnectivity(); // Construct mesh neighbor connectivity
    
    
    std::ofstream outgmeshOneDL("geomesh_OneDL.txt");
    gmesh_OneDL->Print(outgmeshOneDL);
    
    std::ofstream vtkgmeshOneDL("geomesh_OneDL.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDL, vtkgmeshOneDL);
    
    
    return gmesh_OneDL;
}



// ************************************** Create 2D triangle meshes ***************************************

TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDTri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    
    std::string name("geomesh TwoDTri"); // geometry name
    gmesh_TwoDTri->SetName(name);
    gmesh_TwoDTri->SetDimension(geometry_dim);
    
    
    gmesh_TwoDTri->NodeVec().Resize(nnodes); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodes);
    
    
    TPZVec<long> Triangle_topology(3);
    TPZVec <long> Linear_topology(2);

    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    long elementid = 0;
    int physical_id = 1;

    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    {
        
        // 0th element
        {
            
            coord[0] = Lx; // x coordinate
            coord[1] = Ly/2; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[0].SetNodeId(3);
            gmesh_TwoDTri->NodeVec()[0].SetCoord(coord);
            Triangle_topology[0] = 0; // index
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[1].SetNodeId(7);
            gmesh_TwoDTri->NodeVec()[1].SetCoord(coord);
            Triangle_topology[1] = 1;
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = Ly/2; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[2].SetNodeId(4);
            gmesh_TwoDTri->NodeVec()[2].SetCoord(coord);
            Triangle_topology[2] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        
        // 1st element
        
        {
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_TwoDTri->NodeVec()[3].SetNodeId(0);
        gmesh_TwoDTri->NodeVec()[3].SetCoord(coord);
        Triangle_topology[0] = 3;
            
        coord[0] = Lx/2; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_TwoDTri->NodeVec()[4].SetNodeId(1);
        gmesh_TwoDTri->NodeVec()[4].SetCoord(coord);
        Triangle_topology[1] = 4;
            
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly/2; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_TwoDTri->NodeVec()[5].SetNodeId(5);
        gmesh_TwoDTri->NodeVec()[5].SetCoord(coord);
        Triangle_topology[2] = 5;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
        
        }
        
        // 2nd element
        
        {

            Triangle_topology[0] = 4;
            
            coord[0] = Lx; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[6].SetNodeId(2);
            gmesh_TwoDTri->NodeVec()[6].SetCoord(coord);
            Triangle_topology[1] = 6;
            
            Triangle_topology[2] = 0;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        
        // 3rd element

        {

            Triangle_topology[0] = 4;
            
            Triangle_topology[1] = 0;
            
            Triangle_topology[2] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        
        // 4th element
        
        {

            Triangle_topology[0] = 5;
            
            Triangle_topology[1] = 2;
            
            Triangle_topology[2] = 1;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        
        
        // 5th element
        
        {

            Triangle_topology[0] = 5;
            
            Triangle_topology[1] = 1;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[7].SetNodeId(8);
            gmesh_TwoDTri->NodeVec()[7].SetCoord(coord);
            Triangle_topology[2] = 7;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        
        
        // 6th element
        
        {
 
            Triangle_topology[0] = 4;
  
            Triangle_topology[1] = 2;

            Triangle_topology[2] = 5;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
        // 7th element
        
        {
 
            Triangle_topology[0] = 0;
            
            coord[0] = Lx; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[8].SetNodeId(6);
            gmesh_TwoDTri->NodeVec()[8].SetCoord(coord);
            Triangle_topology[1] = 8;

            Triangle_topology[2] = 1;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }
    }
   
    // bottom
    {
        {
            
            Linear_topology[0] = 3;

            Linear_topology[1] = 4;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom
            elementid++;
        }
    
        {

            Linear_topology[0] = 4;
            
            Linear_topology[1] = 6;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom
            elementid++;
        }
    }
       // right
    
    {
        {

            Linear_topology[0] = 6;
            
            Linear_topology[1] = 0;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; right
            elementid++;
        }
        
        {

            Linear_topology[0] = 0;

            Linear_topology[1] = 8;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; right
            elementid++;
        }
    }
    // top
    {
        {

            Linear_topology[0] = 8;
            
            Linear_topology[1] = 1;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDTri); // create boundary element; top
            elementid++;
        }
        
        {
            
            Linear_topology[0] = 1;
            
            Linear_topology[1] = 7;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDTri); // create boundary element; top
            elementid++;
        }
    }
      // left
    {
        {

            Linear_topology[0] = 7;

            Linear_topology[1] = 5;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; left
            elementid++;
        }
        
        {

            Linear_topology[0] = 5;
            
            Linear_topology[1] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; left
            elementid++;
        }
    }
    
    // Build the mesh
    gmesh_TwoDTri->BuildConnectivity();
    
    std::ofstream outgmeshTwoDTri("geomesh_TwoDTri.txt");
    gmesh_TwoDTri->Print(outgmeshTwoDTri);
    
    std::ofstream vtkgmeshTwoDTri("geomesh_TwoDTri.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDTri, vtkgmeshTwoDTri);
    return gmesh_TwoDTri;
    
}



// ************************************** Create 2D quadrilateral meshes ***************************************

TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDQuad = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    
    std::string name("geomesh TwoDQuad"); // geometry name
    gmesh_TwoDQuad->SetName(name);
    gmesh_TwoDQuad->SetDimension(geometry_dim);
    
    
    gmesh_TwoDQuad->NodeVec().Resize(nnodesqu); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodesqu);
    
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec<long> Triangle_topology(3);
    TPZVec <long> Linear_topology(2);
    TPZVec<long> point_topology(1);

    TPZVec<REAL> coord(3,0.0);
    
    // Index of element

    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    
    {
        // 0th element
        {
            
            coord[0] = Lx; // x coordinate
            coord[1] = Ly/2; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[0].SetNodeId(3);
            gmesh_TwoDQuad->NodeVec()[0].SetCoord(coord);
            Quadrilateral_topology[0] = 0; // index
            
            coord[0] = Lx; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[1].SetNodeId(9);
            gmesh_TwoDQuad->NodeVec()[1].SetCoord(coord);
            Quadrilateral_topology[1] = 1;
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[2].SetNodeId(8);
            gmesh_TwoDQuad->NodeVec()[2].SetCoord(coord);
            Quadrilateral_topology[2] = 2;
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = Ly/2; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[3].SetNodeId(4);
            gmesh_TwoDQuad->NodeVec()[3].SetCoord(coord);
            Quadrilateral_topology[3] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
            elementid++;
            
        }
    
        // 1st element
        {
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[4].SetNodeId(1);
            gmesh_TwoDQuad->NodeVec()[4].SetCoord(coord);
            Quadrilateral_topology[0] = 4; // index
            

            Quadrilateral_topology[1] = 3;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = Ly/2; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[5].SetNodeId(5);
            gmesh_TwoDQuad->NodeVec()[5].SetCoord(coord);
            Quadrilateral_topology[2] = 5;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[6].SetNodeId(0);
            gmesh_TwoDQuad->NodeVec()[6].SetCoord(coord);
            Quadrilateral_topology[3] = 6;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
            elementid++;
            
        }
  
        // 2nd element
        {
            
            coord[0] = Lx; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[7].SetNodeId(2);
            gmesh_TwoDQuad->NodeVec()[7].SetCoord(coord);
            Quadrilateral_topology[0] = 7; // index
            
            
            Quadrilateral_topology[1] = 0;
            
            Quadrilateral_topology[2] = 3;

            Quadrilateral_topology[3] = 4;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
            elementid++;
            
        }
        // 3rd element
        
        {

            Triangle_topology[0] = 3;
            
            Triangle_topology[1] = 2;
            
            coord[0] = Lx/4; // x coordinate
            coord[1] =0.75*Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[8].SetNodeId(7);
            gmesh_TwoDQuad->NodeVec()[8].SetCoord(coord);
            Triangle_topology[2] = 8;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
            elementid++;
        }
        
        // 4th element
        
        {
            
            Triangle_topology[0] = 5;
            
            Triangle_topology[1] = 8;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[9].SetNodeId(6);
            gmesh_TwoDQuad->NodeVec()[9].SetCoord(coord);
            Triangle_topology[2] = 9;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
            elementid++;
        }
        
        // 5th element
        
        {
            
            Triangle_topology[0] = 3;
            
            Triangle_topology[1] = 8;

            Triangle_topology[2] = 5;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
            elementid++;
        }
        
        // 6th element
        
        {
            
            Triangle_topology[0] = 8;
            
            Triangle_topology[1] = 9;
            
            Triangle_topology[2] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
            elementid++;
        }
        
        // 7th element point
        
        {

            point_topology[0] = 3;
            new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (elementid, point_topology, physical_id, *gmesh_TwoDQuad); //create point element
            elementid++;

        }

    }
        
        // bottom
        {
            {
                
                Linear_topology[0] = 6;
                
                Linear_topology[1] = 4;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom
                elementid++;
            }
            
            {
                
                Linear_topology[0] = 4;
                
                Linear_topology[1] = 7;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom
                elementid++;
            }
        }
        // right
        
        {
            {
                
                Linear_topology[0] = 7;
                
                Linear_topology[1] = 0;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; right
                elementid++;
            }
            
            {
                
                Linear_topology[0] = 0;
                
                Linear_topology[1] = 1;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; right
                elementid++;
            }
        }
        // top
        {
            {
                
                Linear_topology[0] = 1;
                
                Linear_topology[1] = 2;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; top
                elementid++;
            }
            
            {
                
                Linear_topology[0] = 2;
                
                Linear_topology[1] = 9;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; top
                elementid++;
            }
        }
        // left
        {
            {
                
                Linear_topology[0] = 9;
                
                Linear_topology[1] = 5;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; left
                elementid++;
            }
            
            {
                
                Linear_topology[0] = 5;
                
                Linear_topology[1] = 6;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; left
                elementid++;
            }
        }

        
        // Build the mesh
        gmesh_TwoDQuad->BuildConnectivity();
        
    std::ofstream outgmeshTwoDQuad("geomesh_TwoDQuad.txt");
    gmesh_TwoDQuad->Print(outgmeshTwoDQuad);
    
    std::ofstream vtkgmeshTwoDQuad("geomesh_TwoDQuad.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDQuad, vtkgmeshTwoDQuad);

    
        return gmesh_TwoDQuad;
        
    }
   

// ************************************** Create 3D hexahedral and prism meshes ***************************************


TPZGeoMesh *CreateThreeDHexPrytetGMesh(long nnodesthrhpt, REAL Lx, REAL Ly, REAL Lz)
{
    TPZGeoMesh * gmesh_ThreeDHexPrytet = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    
    std::string name("geomesh ThreeDHexPrytet"); // geometry name
    gmesh_ThreeDHexPrytet->SetName(name);
    gmesh_ThreeDHexPrytet->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDHexPrytet->NodeVec().Resize(nnodesthrhpt); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodesthrhpt);
    
    TPZVec<long> Hexahedron_topology(8);
    TPZVec <long> Quadrilateral_topology(4);
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_front = -1; // define id for a material (border in front)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_back = -3; // define id for a material (border back)
    const int bc_left = -4; // define id for a material (border left)
    const int bc_bottom = -5; // define id for a material (border bottom)
    const int bc_top = -6; // define id for a material (border top)
    
    // 0th element
    
    {
        
        // 0th node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDHexPrytet->NodeVec()[0].SetCoord(coord);
        Hexahedron_topology[0] = 0;
        
        // 1st node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDHexPrytet->NodeVec()[1].SetCoord(coord);
        Hexahedron_topology[1] = 1;
        
        // 2nd node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDHexPrytet->NodeVec()[2].SetCoord(coord);
        Hexahedron_topology[2] = 2;
        
        // 3rd node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDHexPrytet->NodeVec()[3].SetCoord(coord);
        Hexahedron_topology[3] = 3;
        
        // 4th node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[4].SetNodeId(4);
        gmesh_ThreeDHexPrytet->NodeVec()[4].SetCoord(coord);
        Hexahedron_topology[4] = 4;
        
        // 5th node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[5].SetNodeId(5);
        gmesh_ThreeDHexPrytet->NodeVec()[5].SetCoord(coord);
        Hexahedron_topology[5] = 5;
        
        // 6th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[6].SetNodeId(6);
        gmesh_ThreeDHexPrytet->NodeVec()[6].SetCoord(coord);
        Hexahedron_topology[6] = 6;
        
        // 7th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[7].SetNodeId(7);
        gmesh_ThreeDHexPrytet->NodeVec()[7].SetCoord(coord);
        Hexahedron_topology[7] = 7;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
        elementid++;
        
    }
    
    // 1st element
    
    {
        
        // 1st node
        Hexahedron_topology[0] = 1;
        
        // 2nd node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[8].SetNodeId(8);
        gmesh_ThreeDHexPrytet->NodeVec()[8].SetCoord(coord);
        Hexahedron_topology[1] = 8;
        
        // 3rd node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[9].SetNodeId(9);
        gmesh_ThreeDHexPrytet->NodeVec()[9].SetCoord(coord);
        Hexahedron_topology[2] = 9;
        
        // 4th node
        Hexahedron_topology[3] = 2;
        
        // 5th node
        Hexahedron_topology[4] = 5;
        
        // 6th node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[10].SetNodeId(10);
        gmesh_ThreeDHexPrytet->NodeVec()[10].SetCoord(coord);
        Hexahedron_topology[5] = 10;
        
        // 7th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[11].SetNodeId(11);
        gmesh_ThreeDHexPrytet->NodeVec()[11].SetCoord(coord);
        Hexahedron_topology[6] = 11;
        
        // 8th node
        Hexahedron_topology[7] = 6;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
        elementid++;
        
    }
    
    // ********************* boundary ****************
    
    // in front
    {
        {
            
            Quadrilateral_topology[0] = 0;
            
            Quadrilateral_topology[1] = 1;
            
            Quadrilateral_topology[2] = 2;
            
            Quadrilateral_topology[3] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPrytet); // create boundary element; in front
            elementid++;
            
        }
        {
            
            Quadrilateral_topology[0] = 1;
            
            Quadrilateral_topology[1] = 8;
            
            Quadrilateral_topology[2] = 9;
            
            Quadrilateral_topology[3] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPrytet); // create boundary element; in front
            elementid++;
            
        }
        
    }
    

    // right
    {
        {
            
            Quadrilateral_topology[0] = 4;
            
            Quadrilateral_topology[1] = 5;
            
            Quadrilateral_topology[2] = 1;
            
            Quadrilateral_topology[3] = 0;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPrytet); // create boundary element; right
            elementid++;
            
        }
        {
            
            Quadrilateral_topology[0] = 5;
            
            Quadrilateral_topology[1] = 10;
            
            Quadrilateral_topology[2] = 8;
            
            Quadrilateral_topology[3] = 1;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPrytet); // create boundary element; right
            elementid++;
            
        }
        
    }
    
    // back
    {
        {
            
            Quadrilateral_topology[0] = 4;
            
            Quadrilateral_topology[1] = 5;
            
            Quadrilateral_topology[2] = 6;
            
            Quadrilateral_topology[3] = 7;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPrytet); // create boundary element; back
            elementid++;
            
        }
        {
            
            Quadrilateral_topology[0] = 5;
            
            Quadrilateral_topology[1] = 10;
            
            Quadrilateral_topology[2] = 11;
            
            Quadrilateral_topology[3] = 6;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPrytet); // create boundary element; back
            elementid++;
            
        }
        
    }
    
    // left
    {
        {
            
            Quadrilateral_topology[0] = 7;
            
            Quadrilateral_topology[1] = 6;
            
            Quadrilateral_topology[2] = 2;
            
            Quadrilateral_topology[3] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPrytet); // create boundary element; left
            elementid++;
            
        }
        {
            
            Quadrilateral_topology[0] = 6;
            
            Quadrilateral_topology[1] = 11;
            
            Quadrilateral_topology[2] = 9;
            
            Quadrilateral_topology[3] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPrytet); // create boundary element; left
            elementid++;
            
        }
        
    }
    
    
    
    // bottom

    {
            
            Quadrilateral_topology[0] = 0;
            
            Quadrilateral_topology[1] = 4;
            
            Quadrilateral_topology[2] = 7;
            
            Quadrilateral_topology[3] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_bottom,*gmesh_ThreeDHexPrytet); // create boundary element; bottom
            elementid++;
            
   }
    
    
    // top
   {
            
            Quadrilateral_topology[0] = 8;
            
            Quadrilateral_topology[1] = 10;
            
            Quadrilateral_topology[2] = 11;
            
            Quadrilateral_topology[3] = 9;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_top,*gmesh_ThreeDHexPrytet); // create boundary element; top
            elementid++;
            
   }

    
    // Build the mesh
    gmesh_ThreeDHexPrytet->BuildConnectivity();
    std::ofstream outgmeshThreeDHexPrytet("geomesh_ThreeDHexPrytet.txt");
    gmesh_ThreeDHexPrytet->Print(outgmeshThreeDHexPrytet);
    
    std::ofstream vtkgmeshThreeDHexPrytet("geomesh_ThreeDHexPrytet.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDHexPrytet, vtkgmeshThreeDHexPrytet);
    
    return gmesh_ThreeDHexPrytet;
    
}

// ----------------------------------------------------------------------------------------
