#ifndef NODE
#define NODE

#include <iostream>

template<typename T>
class Node{
    /*
     ================================================================
    % ------------ %
    %  Node class  %
    % ------------ %
        Describes a node by its (x,y,z) coordinates and the global number (id)

    Implementation: MD
    ---------------  
    Input:
        ------
           Id (integer) -> global indentification number,
            x (float/double) -> x coordinate,
            y (float/double) -> y coordinate,
            z (float/double) -> z coordinate.
    ================================================================       
    */
    private:
        size_t Id;
        T x, y, z;
    public:
        // Constructor
        Node(): Id(0), x(0.0), y(0.0), z(0.0){}
        Node(size_t Id_in, T x_in): y(0.0), z(0.0) // 1D point constructor
        {
            Id=Id_in; x=x_in; 
        }
        Node(size_t Id_in, T x_in, T y_in): z(0.0) // 2D point constructor
        {
            Id=Id_in; x=x_in; y=y_in;
        }
        Node(size_t Id_in, T x_in, T y_in, T z_in) // 3D point constructor
        {
            Id=Id_in; x=x_in; y=y_in; z=z_in;
        }
        // Copy Constructor
        Node(Node<T> const& other)
        {
            Id=other.Id; x=other.x; y=other.y; z=other.z;
        }
        // Overload "=" operator
        Node& operator=(Node<T> other)
        {
            Id=other.Id; x=other.x; y=other.y; z=other.z;
            return *this;
        }
        // Member functions that do not modify the object
        size_t get_Id() const {return Id;}
        T get_x() const {return x;}
        T get_y() const {return y;}
        T get_z() const {return z;}
        void print() const
        {
            std::cout << "Node\t" << Id << ":\t(" << x << "\t" << y << "\t" << z <<")" << std::endl;
        }
    };
#endif