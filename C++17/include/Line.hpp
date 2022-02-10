#ifndef LINE
#define LINE

#include <iostream>
#include "Node.hpp"

class Line{
    /*
    ================================================================
    % ---------- %
    % Line class %
    % ---------- %
        Describes an edge by the ids of the two vertex, the global identification
    number (id) and the edge type.

    Impelmentation: MD
    ---------------  
    Input:
        ------
        Id (size_t) -> global indentification number,
        N1 (size_t) -> first vertex id,
        N2 (size_t) -> second vertex id
        EdgeType (integer) -> edge type.         
        Proc (size_t) -> proc id
    ================================================================       
    */

private:
    size_t Id;
    Node<double> N1;
    Node<double> N2;
    int EdgeType;
    size_t Partition;

public:
    Line(): Id(0), EdgeType(-1), Partition(0){}
    Line(size_t Id_in, Node<double> N1_in, Node<double> N2_in, int EdgeType_in, size_t Partition_in)
    {
        Id = Id_in;
        N1 = N1_in;
        N2 = N2_in;
        EdgeType = EdgeType_in;
        Partition = Partition_in;
    }
    // Member functions that do not modify this object
    size_t get_Id() const {return Id;}
    size_t get_N1Id() const {return N1.get_Id();}
    size_t get_N2Id() const {return N2.get_Id();}
    double get_N1x() const {return N1.get_x();}
    double get_N1y() const {return N1.get_y();}
    double get_N1z() const {return N1.get_z();}
    double get_N2x() const {return N2.get_x();}
    double get_N2y() const {return N2.get_y();}
    double get_N2z() const {return N2.get_z();}
    int get_EdgeType() const {return EdgeType;}
    size_t get_Partition() const {return Partition;}
    void print() const
    {   
        std::cout << "Line: " << Id << ", EdgeType: " << EdgeType << ", Partition: " << Partition << std::endl;
        std::cout << "\t"; N1.print();
        std::cout << "\t"; N2.print();
    }
};

#endif