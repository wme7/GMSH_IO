#ifndef TRIANGLE
#define TRIANGLE

#include <iostream>
#include "Node.hpp"

class Triangle{
    /*
    ================================================================
    % -------------- %
    % Triangle class %
    % -------------- %
        Describes an element by its three vertex ids, the global identification
    number (id) and the element type.

    Impelmentation: MD
    ---------------  
    Input:
        ------
        Id (size_t) -> global indentification number,
        N1 (Node) -> first vertex id,
        N2 (Node) -> second vertex id
        N3 (Node) -> third vertex id,
        ElementType (integer) -> element type.
        Partition (size_t) -> partition id
    ================================================================       
    */

private:
    size_t Id;
    Node<double> N1;
    Node<double> N2;
    Node<double> N3;
    int ElementType;
    size_t Partition;
    double Area;

public:
    Triangle(): Id(0), ElementType(-1), Partition(0), Area(0.){}
    Triangle(size_t Id_in, Node<double> N1_in, Node<double> N2_in, Node<double> N3_in, int ElementType_in, size_t Partition_in)
    {
        Id = Id_in;
        N1 = N1_in;
        N2 = N2_in;
        N3 = N3_in;
        ElementType = ElementType_in;
        Partition = Partition_in;
        set_Area();
    }
    // Member function that do modify the object
    void set_Area() 
    {   // Area = ((ax-cx).*(by-cy)-(bx-cx).*(ay-cy)) / 2  where:
        // ax = N1.get_x;  bx = N2.get_x;  cx = N3.get_x;
        // ay = N1.get_y;  by = N2.get_y;  cy = N3.get_y;
        Area = 0.5*((N1.get_x()-N3.get_x())*(N2.get_y()-N3.get_y())-(N2.get_x()-N3.get_x())*(N1.get_y()-N3.get_y()));
    }
    double testNodeOrder(const bool print)
    {   // If Area is negative, re-order the nodes inside the element
        if (Area<0)
        {   // SWAP nodes: EToV(i,:) = EToV(i,[1 3 2]);
            Node<double> Nt=N2; N2=N3; N3=Nt;
        }
        set_Area(); if(print) std::cout << "element("<< Id <<").Area= " << Area << std::endl;
        return Area;
    }
    // Member functions that do NOT modify the object
    size_t get_Id() const {return Id;}
    size_t get_N1Id() const {return N1.get_Id();}
    size_t get_N2Id() const {return N2.get_Id();}
    size_t get_N3Id() const {return N3.get_Id();}
    double get_N1x() const {return N1.get_x();}
    double get_N1y() const {return N1.get_y();}
    double get_N1z() const {return N1.get_z();}
    double get_N2x() const {return N2.get_x();}
    double get_N2y() const {return N2.get_y();}
    double get_N2z() const {return N2.get_z();}
    double get_N3x() const {return N3.get_x();}
    double get_N3y() const {return N3.get_y();}
    double get_N3z() const {return N3.get_z();}
    int get_FaceType() const {return ElementType;}
    int get_ElementType() const {return ElementType;}
    size_t get_Partition() const {return Partition;}
    void print() const
    {
        std::cout << "Triangle: " << Id << ", ElementType: " << ElementType << ", Partition: " << Partition << ", Area: " << Area << std::endl;
        std::cout << "\t"; N1.print();
        std::cout << "\t"; N2.print();
        std::cout << "\t"; N3.print();
    }
};

#endif
