#ifndef POINT
#define POINT

#include <iostream>
#include "Node.hpp"

class Point{
    /*
     ================================================================
    % ------------ %
    %  Node class  %
    % ------------ %
        Describes a Point by the ID of a single node, the global identification
    number (id) and the point type.

    Implementation: MD
    ---------------  
    Input:
        ------
        Id (size_t) -> global indentification number,
        N1 (size_t) -> vertex id,
        PointType (integer) -> point type.         
        Proc (size_t) -> proc id
    ================================================================       
    */
private:
    size_t Id;
    Node<double> N1;
    int PointType;
    size_t Partition;

public:
    Point(): Id(0), PointType(-1), Partition(0){}
    Point(size_t Id_in, Node<double> N1_in, int PointType_in, size_t Partition_in)
    {
        Id = Id_in;
        N1 = N1_in;
        PointType = PointType_in;
        Partition = Partition_in;
    }
    // Member functions that do not modify this object
    size_t get_Id() const {return Id;}
    size_t get_N1Id() const {return N1.get_Id();}
    double get_N1x() const {return N1.get_x();}
    double get_N1y() const {return N1.get_y();}
    double get_N1z() const {return N1.get_z();}
    int get_PointType() const {return PointType;}
    size_t get_Partition() const {return Partition;}
    void print() const
    {   
        std::cout << "Point: " << Id << ", PointType: " << PointType << ", Partition: " << Partition << std::endl;
        std::cout << "\t"; N1.print();
    }
};
#endif