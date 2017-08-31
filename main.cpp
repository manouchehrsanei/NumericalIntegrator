
#include <iostream>
#include <cmath>
#include <fstream>


#include "pzreal.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"

#include "TPZVTKGeoMesh.h"

#include "TPZGmshReader.h"



// ************************ (Create a linear Omega description : 1D) ********************************************

TPZGeoMesh *CreateOmega1D(long num_el, REAL a1, REAL a2);


// ************************ (Create a linear Omega description : 2D) ********************************************

TPZGeoMesh *CreateOmega2DTri(long nnodes, REAL a1, REAL a2,REAL b1, REAL b2);

TPZGeoMesh *CreateOmega2DQuad(long nnodesqu, REAL a1, REAL a2,REAL b1, REAL b2);

// ************************ (Create a linear Omega description : 3D) ********************************************

TPZGeoMesh *CreateOmega3DHexa(long nnodeshexa, REAL a1, REAL a2,REAL b1, REAL b2, REAL c1, REAL c2);


// ************************************************ (Example Functions) ****************************************

REAL f_1D(TPZManVector<REAL,3> &x);
REAL f_2D(TPZManVector<REAL,3> &x);
REAL f_3D(TPZManVector<REAL,3> &x);


// ************************************************** (From Gmsh) ***********************************************

TPZGeoMesh *CreateOneDCircleGMesh();
TPZGeoMesh *CreateTwoDCircleGMesh();
TPZGeoMesh *CreateTwoDSphereGMesh();
TPZGeoMesh *CreateThreeDSphereGMesh();


void IntegrateOneCirGeometry(TPZGeoMesh * OmegaOneCir, int order);
void IntegrateTwoCirGeometry(TPZGeoMesh * OmegaTwoCir, int order);
void IntegrateTwoSpheGeometry(TPZGeoMesh * OmegaTwoSphe, int order);
void IntegrateThreeSpheGeometry(TPZGeoMesh * OmegaThreeSphe, int order);



// ******************************************** (main of program) **********************************************


int main() {


    
    int order = 10;
    
//    TPZGeoMesh * OmegaOneCir = CreateOneDCircleGMesh();
//    IntegrateOneCirGeometry(OmegaOneCir, order);
    
//    TPZGeoMesh * OmegaTwoCir = CreateTwoDCircleGMesh();
//    IntegrateTwoCirGeometry(OmegaTwoCir, order);
//
//    TPZGeoMesh * OmegaTwoSphe = CreateTwoDSphereGMesh();
//    IntegrateTwoSpheGeometry(OmegaTwoSphe, order);

    TPZGeoMesh * OmegaThreeSphe = CreateThreeDSphereGMesh();
    IntegrateThreeSpheGeometry(OmegaThreeSphe, order);


    
    
    return 0;
}


// ++++++++++++++++++++++++++++++++++++++++++++ functions +++++++++++++++++++++++++++++++++++++++++++++++++++


void IntegrateOneCirGeometry(TPZGeoMesh * OmegaOneCir, int order){
    
    
    //    long num_el = 1;
    //    REAL a1 = 0.0;
    //    REAL a2 = 10.0;
    //    TPZGeoMesh * Omega = CreateOmega1D(num_el, a1, a2);

    long n_elements = OmegaOneCir->NElements();
    REAL NumericalIntegral = 0.0;
    
    for (long iel = 0; iel < n_elements; iel++)
    { // for a all element in the mesh
        
        TPZGeoEl * element = OmegaOneCir->Element(iel);
        int itself = element->NSides() - 1;
        TPZIntPoints * IntegrationRule = element->CreateSideIntegrationRule(itself, order);
        int n_points = IntegrationRule ->NPoints();
        
        TPZManVector<REAL,3> x;
        TPZManVector<REAL,3> par_space;
        REAL w, detjac;
        TPZFMatrix<REAL> jac, jacinv, axes;
        
        for (int ip = 0; ip < n_points; ip++)
        {
            IntegrationRule->Point(ip, par_space, w);
            
            element->X(par_space, x);
            element->Jacobian(par_space, jac, axes, detjac, jacinv);
            
            NumericalIntegral += w * detjac * (1.0);
        }
    }
    
    std::cout << "Numerical integration = " << NumericalIntegral << std::endl;
}

// ----------------------------------------------------------------------------------------


void IntegrateTwoCirGeometry(TPZGeoMesh * OmegaTwoCir, int order){
    
    
//        int order = 20;
//        long nnodes = 9;
//        long num_el = nnodes-1;
//        REAL a1 = 0.0;
//        REAL a2 = 10.0;
//        REAL b1 = 0.0;
//        REAL b2 = 10.0;
//    
//        TPZGeoMesh * Omega = CreateOmega2DTri(nnodes, a1, a2, b1, b2);
    
        long n_elements = OmegaTwoCir->NElements();

        REAL NumericalIntegral = 0.0;
    
        for (long iel = 0; iel < n_elements; iel++)
        { // for a all element in the mesh
    
            TPZGeoEl * element = OmegaTwoCir->Element(iel);
            int itself = element->NSides() - 1;
    
            TPZIntPoints * IntegrationRule = element->CreateSideIntegrationRule(itself, order);
            int n_points = IntegrationRule ->NPoints();
    
            TPZManVector<REAL,3> x;
    
            TPZManVector<REAL,3> par_space;
            REAL w, detjac;
            TPZFMatrix<REAL> jac, jacinv, axes;
    
            for (int ip = 0; ip < n_points; ip++)
            {
                IntegrationRule->Point(ip, par_space, w);
    
                element->X(par_space, x);
                element->Jacobian(par_space, jac, axes, detjac, jacinv);
    
//                NumericalIntegral += w * detjac * (f_2D(x));
                NumericalIntegral += w * detjac * (1);

            }
        }
        
        std::cout << "Numerical integration = " << NumericalIntegral << std::endl;
    
}

// ----------------------------------------------------------------------------------------


void IntegrateTwoSpheGeometry(TPZGeoMesh * OmegaTwoSphe, int order){
    
    //    int order = 15;
    //    long nnodesqu = 9;
    //    long num_el = nnodesqu/2;
    //    REAL a1 = 0.0;
    //    REAL a2 = 10.0;
    //    REAL b1 = 0.0;
    //    REAL b2 = 10.0;
    //
    //    TPZGeoMesh * Omega = CreateOmega2DQuad(nnodesqu, a1, a2, b1, b2);
    
        long n_elements = OmegaTwoSphe->NElements();
    
        REAL NumericalIntegral = 0.0;
    
        for (long iel = 0; iel < n_elements; iel++)
        { // for a all element in the mesh
    
            TPZGeoEl * element = OmegaTwoSphe->Element(iel);
            int itself = element->NSides() - 1;
    
            TPZIntPoints * IntegrationRule = element->CreateSideIntegrationRule(itself, order);
            int n_points = IntegrationRule ->NPoints();
    
            TPZManVector<REAL,3> x;
    
            TPZManVector<REAL,3> par_space;
            REAL w, detjac;
            TPZFMatrix<REAL> jac, jacinv, axes;
    
            for (int ip = 0; ip < n_points; ip++)
            {
                IntegrationRule->Point(ip, par_space, w);
    
                element->X(par_space, x);
                element->Jacobian(par_space, jac, axes, detjac, jacinv);
    
//                NumericalIntegral += w * detjac * (f_3D(x));
                NumericalIntegral += w * detjac * (1);

            }
        }
        
        std::cout << "Numerical integration = " << NumericalIntegral << std::endl;
    
}

// ----------------------------------------------------------------------------------------

void IntegrateThreeSpheGeometry(TPZGeoMesh * OmegaThreeSphe, int order){
    
    //
    //    int order = 20;
    //    long nnodeshexa = 12;
    //    long num_el = nnodeshexa/6;
    //    REAL a1 = 0.0;
    //    REAL a2 = 10.0;
    //    REAL b1 = 0.0;
    //    REAL b2 = 10.0;
    //    REAL c1 = 0.0;
    //    REAL c2 = 10.0;
    //
    //    TPZGeoMesh * Omega = CreateOmega3DHexa(nnodeshexa, a1, a2, b1, b2, c1, c2);
    
        long n_elements = OmegaThreeSphe->NElements();
    
        REAL NumericalIntegral = 0.0;
    
        for (long iel = 0; iel < n_elements; iel++)
        { // for a all element in the mesh
    
            TPZGeoEl * element = OmegaThreeSphe->Element(iel);
            int itself = element->NSides() - 1;
    
            TPZIntPoints * IntegrationRule = element->CreateSideIntegrationRule(itself, order);
            int n_points = IntegrationRule ->NPoints();
    
            TPZManVector<REAL,3> x;
    
            TPZManVector<REAL,3> par_space;
            REAL w, detjac;
            TPZFMatrix<REAL> jac, jacinv, axes;
    
            for (int ip = 0; ip < n_points; ip++)
            {
                IntegrationRule->Point(ip, par_space, w);
    
                element->X(par_space, x);
                element->Jacobian(par_space, jac, axes, detjac, jacinv);
    
//                NumericalIntegral += w * detjac * (f_3D(x));
                NumericalIntegral += w * detjac * (1);

            }
        }
        
        std::cout << "Numerical integration = " << NumericalIntegral << std::endl;

}


// ------------------------------------------ 1D -------------

REAL f_1D(TPZManVector<REAL,3> &x)
{
    REAL x1 = x[0];
    
    REAL f1D = (1.0 - x1)*x1*sin(2.0*x1);
    return f1D;
}

// ------------------------------------------ 2D -------------

REAL f_2D(TPZManVector<REAL,3> &x)
{
    REAL x1 = x[0];
    REAL x2 = x[1];

    REAL f2D = ((1.0 - x1)*x1*sin(2.0*x1))*((1.0 - x2)*x2*sin(2.0*x2));
    return f2D;
}

// ------------------------------------------ 3D -------------

REAL f_3D(TPZManVector<REAL,3> &x)
{
    REAL x1 = x[0];
    REAL x2 = x[1];
    REAL x3 = x[2];

    REAL f3D = ((1.0 - x1)*x1*sin(2.0*x1))*((1.0 - x2)*x2*sin(2.0*x2))*((1.0 - x3)*x3*sin(2.0*x3));
    return f3D;
}




// ********************************* (linear elements) ******************************************************


// ******************************* Create a linear Omega description : 1D ************************************

TPZGeoMesh *CreateOmega1D(long num_el, REAL a1, REAL a2)
{
    REAL size_el = (a2-a1)/REAL(num_el);
    
    TPZGeoMesh * gmesh_OneDL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("geomesh OneD"); // geometry name
    
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
    
    
    std::ofstream outgmeshOneDL("geomesh_OneD.txt");
    gmesh_OneDL->Print(outgmeshOneDL);
    
    std::ofstream vtkgmeshOneDL("geomesh_OneD.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDL, vtkgmeshOneDL);
    
    
    return gmesh_OneDL;
}



// ******************************* Create a linear Omega description : 2D ************************************

TPZGeoMesh *CreateOmega2DTri(long nnodes, REAL a1, REAL a2,REAL b1, REAL b2)
{
    TPZGeoMesh * gmesh_TwoDTri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    REAL Lx = a2-a1;
    REAL Ly = b2-b1;
    
    std::string name("geomesh TwoDTri"); // geometry name
    gmesh_TwoDTri->SetName(name);
    gmesh_TwoDTri->SetDimension(geometry_dim);
    
    
    gmesh_TwoDTri->NodeVec().Resize(nnodes); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodes);
    
    
    TPZVec<long> Triangle_topology(3);

    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    long elementid = 0;
    int physical_id = 1;

    
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
    
    // Build the mesh
    gmesh_TwoDTri->BuildConnectivity();
    
    std::ofstream outgmeshTwoDTri("geomesh_TwoDTri.txt");
    gmesh_TwoDTri->Print(outgmeshTwoDTri);
    
    std::ofstream vtkgmeshTwoDTri("geomesh_TwoDTri.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDTri, vtkgmeshTwoDTri);
    return gmesh_TwoDTri;
    
}



// ******************************* Create a linear Omega description : 2D ************************************

TPZGeoMesh *CreateOmega2DQuad(long nnodesqu, REAL a1, REAL a2,REAL b1, REAL b2)

{
    TPZGeoMesh * gmesh_TwoDQuad = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    REAL Lx = a2-a1;
    REAL Ly = b2-b1;
    

    std::string name("geomesh TwoDQuad"); // geometry name
    gmesh_TwoDQuad->SetName(name);
    gmesh_TwoDQuad->SetDimension(geometry_dim);
    
    
    gmesh_TwoDQuad->NodeVec().Resize(nnodesqu); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodesqu);
    
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec<long> Triangle_topology(3);


    TPZVec<REAL> coord(3,0.0);
    
    // Index of element

    long elementid = 0;
    int physical_id = 1;

    
    
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
            gmesh_TwoDQuad->NodeVec()[1].SetNodeId(8);
            gmesh_TwoDQuad->NodeVec()[1].SetCoord(coord);
            Quadrilateral_topology[1] = 1;
            
            coord[0] = Lx/2; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[2].SetNodeId(7);
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
        // 3d element
        {
            

            Quadrilateral_topology[0] = 3; // index
            
            
            Quadrilateral_topology[1] = 2;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[8].SetNodeId(6);
            gmesh_TwoDQuad->NodeVec()[8].SetCoord(coord);
            Quadrilateral_topology[2] = 8;
            
            Quadrilateral_topology[3] = 5;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
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
   

// ******************************* Create a linear Omega description : 3D ************************************

TPZGeoMesh *CreateOmega3DHexa(long nnodeshexa, REAL a1, REAL a2,REAL b1, REAL b2, REAL c1, REAL c2)
{
    TPZGeoMesh * gmesh_ThreeDHexPrytet = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    REAL Lx = a2-a1;
    REAL Ly = b2-b1;
    REAL Lz = c2-c1;

    std::string name("geomesh ThreeDHexa"); // geometry name
    gmesh_ThreeDHexPrytet->SetName(name);
    gmesh_ThreeDHexPrytet->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDHexPrytet->NodeVec().Resize(nnodeshexa); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodeshexa);
    
    TPZVec<long> Hexahedron_topology(8);
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;

    
    // 0th element
    
    {
        
        // 0th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDHexPrytet->NodeVec()[0].SetCoord(coord);
        Hexahedron_topology[0] = 0;
        
        
        // 1st node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDHexPrytet->NodeVec()[1].SetCoord(coord);
        Hexahedron_topology[1] = 1;
        
        // 2nd node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDHexPrytet->NodeVec()[2].SetCoord(coord);
        Hexahedron_topology[2] = 2;
        
        // 3rd node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDHexPrytet->NodeVec()[3].SetCoord(coord);
        Hexahedron_topology[3] = 3;
        
        // 4th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[4].SetNodeId(4);
        gmesh_ThreeDHexPrytet->NodeVec()[4].SetCoord(coord);
        Hexahedron_topology[4] = 4;
        
        // 5th node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[5].SetNodeId(5);
        gmesh_ThreeDHexPrytet->NodeVec()[5].SetCoord(coord);
        Hexahedron_topology[5] = 5;
        
        // 6th node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[6].SetNodeId(6);
        gmesh_ThreeDHexPrytet->NodeVec()[6].SetCoord(coord);
        Hexahedron_topology[6] = 6;
        
        // 7th node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz/2; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[7].SetNodeId(7);
        gmesh_ThreeDHexPrytet->NodeVec()[7].SetCoord(coord);
        Hexahedron_topology[7] = 7;

        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
        elementid++;
        
    }
    
    // 1st element
    
    {
        
        // oth node
        Hexahedron_topology[0] = 4;
        
        // 1st node
        Hexahedron_topology[1] = 5;
        
        // 2nd node
        Hexahedron_topology[2] = 6;
        
        // 3rd node
        Hexahedron_topology[3] = 7;
        
        // 4th node
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[8].SetNodeId(8);
        gmesh_ThreeDHexPrytet->NodeVec()[8].SetCoord(coord);
        Hexahedron_topology[4] = 8;
        
        // 5th node
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[9].SetNodeId(9);
        gmesh_ThreeDHexPrytet->NodeVec()[9].SetCoord(coord);
        Hexahedron_topology[5] = 9;
        
        // 6th node
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[10].SetNodeId(10);
        gmesh_ThreeDHexPrytet->NodeVec()[10].SetCoord(coord);
        Hexahedron_topology[6] = 10;
        
        
        // 7th node
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexPrytet->NodeVec()[11].SetNodeId(11);
        gmesh_ThreeDHexPrytet->NodeVec()[11].SetCoord(coord);
        Hexahedron_topology[7] = 11;
        
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
        elementid++;
        
    }
    

    
    // Build the mesh
    gmesh_ThreeDHexPrytet->BuildConnectivity();
    std::ofstream outgmeshThreeDHexPrytet("geomesh_ThreeDHexa.txt");
    gmesh_ThreeDHexPrytet->Print(outgmeshThreeDHexPrytet);
    
    std::ofstream vtkgmeshThreeDHexPrytet("geomesh_ThreeDHexa.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDHexPrytet, vtkgmeshThreeDHexPrytet);
    
    return gmesh_ThreeDHexPrytet;
    
}

// ----------------------------------------------------------------------------------------


TPZGeoMesh *CreateOneDCircleGMesh()
{
    
    // Creating geometric mesh
    
    TPZGeoMesh *gmesh_OneDCircle = new TPZGeoMesh();
    
    // Implementation meshes with GMSH
    
    std::string grid = "1DCircle.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    gmesh_OneDCircle = Geometry.GeometricGmshMesh(grid);
    const std::string name("OneDCircle from gmsh script");
    gmesh_OneDCircle->SetName(name);
    
    std::ofstream outgmeshOneDCircle("geomesh_OneDCircle.txt");
    gmesh_OneDCircle->Print(outgmeshOneDCircle);
    
    std::ofstream vtkgmeshOneDCircle("geomesh_OneDCircle.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDCircle, vtkgmeshOneDCircle);
    return gmesh_OneDCircle;
    
}

// ----------------------------------------------------------------------------------------


TPZGeoMesh *CreateTwoDCircleGMesh()
{
    
    // Creating geometric mesh
    
    TPZGeoMesh *gmesh_TwoDCircle = new TPZGeoMesh();
    
    // Implementation meshes with GMSH
    
    std::string grid = "2DCircle.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    gmesh_TwoDCircle = Geometry.GeometricGmshMesh(grid);
    const std::string name("TwoDCircle from gmsh script");
    gmesh_TwoDCircle->SetName(name);
    
    std::ofstream outgmeshTwoDCircle("geomesh_TwoDCircle.txt");
    gmesh_TwoDCircle->Print(outgmeshTwoDCircle);
    
    std::ofstream vtkgmeshTwoDCircle("geomesh_TwoDCircle.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDCircle, vtkgmeshTwoDCircle);
    return gmesh_TwoDCircle;
    
}

// ----------------------------------------------------------------------------------------


TPZGeoMesh *CreateTwoDSphereGMesh()
{
    
    // Creating geometric mesh
    
    TPZGeoMesh *gmesh_TwoDSphere = new TPZGeoMesh();
    
    // Implementation meshes with GMSH
    
    std::string grid = "2DSphere.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    gmesh_TwoDSphere = Geometry.GeometricGmshMesh(grid);
    const std::string name("TwoDSphere from gmsh script");
    gmesh_TwoDSphere->SetName(name);
    
    std::ofstream outgmeshTwoDSphere("geomesh_TwoDSphere.txt");
    gmesh_TwoDSphere->Print(outgmeshTwoDSphere);
    
    std::ofstream vtkgmeshTwoDSphere("geomesh_TwoDSphere.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDSphere, vtkgmeshTwoDSphere);
    return gmesh_TwoDSphere;
    
}


TPZGeoMesh *CreateThreeDSphereGMesh()
{
    
    // Creating geometric mesh
    
    TPZGeoMesh *gmesh_ThreeDSphere = new TPZGeoMesh();
    
    // Implementation meshes with GMSH
    
    std::string grid = "3DSphere.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    gmesh_ThreeDSphere = Geometry.GeometricGmshMesh(grid);
    const std::string name("ThreeDSphere from gmsh script");
    gmesh_ThreeDSphere->SetName(name);
    
    std::ofstream outgmeshThreeDSphere("geomesh_ThreeDSphere.txt");
    gmesh_ThreeDSphere->Print(outgmeshThreeDSphere);
    
    std::ofstream vtkgmeshThreeDSphere("geomesh_ThreeDSphere.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDSphere, vtkgmeshThreeDSphere);
    return gmesh_ThreeDSphere;
    
}
