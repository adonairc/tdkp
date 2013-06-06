#ifndef MOC_GRID_READER_H_
#define MOC_GRID_READER_H_

#include "tdkp/io/BaseGridReader.h"

namespace tdkp {
	
	class MocGridReaderBase;
	class MocGridReader1D;
	class MocGridReader2DTriangle;
	class MocGridReader2DRectangle;
	class MocGridReader3DTetrahedron;
	
	class MocGridReaderBase : public BaseGridReader {
	public:		
		unsigned int get_num_vertices() const { return vertices.size(); }
		unsigned int get_num_elements() const { return elements.size(); }
		unsigned int get_num_elements_in_region(unsigned int region_index_global) const { return elements.size(); }
		unsigned int get_num_regions()  const { return 1; }
		string get_region_name(unsigned int region_index_global) const { return string("caudillo"); }
		string get_region_material_name(unsigned int region_index_global) const { return string("GaAs"); }
		unsigned int get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const { return element_index_region; }
		Element::ElementShape get_element_shape(unsigned int element_index) const { return element_shapes[element_index]; }
		unsigned int get_num_vertices_in_element(unsigned int element_index) const { return elements[element_index].size(); }
		unsigned int get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const { return elements[element_index][local_vertex_index]; }  
		double get_vertex_x_coord(unsigned int vertex_index) const { return vertices[vertex_index].coords[0]; }
		double get_vertex_y_coord(unsigned int vertex_index) const { return vertices[vertex_index].coords[1]; }
		double get_vertex_z_coord(unsigned int vertex_index) const { return vertices[vertex_index].coords[2]; }
		string get_unique_identifier() const { return string("unique"); }
	protected:
		class vertex_coords {
		public:
			vertex_coords() { coords[0] = coords[1] = coords[2] = 0.0; }
			vertex_coords(double x, double y, double z) { coords[0] = x; coords[1] = y; coords[2] = z; }
			double coords[3];
		};
		vector<vector<unsigned int> > elements;
		vector<Element::ElementShape> element_shapes;
		vector<vertex_coords>         vertices;
	};
	
	class MocGridReader1D : public MocGridReaderBase {
	public:
		MocGridReader1D();
		unsigned int get_dimension() const { return 1; }
			
	};
	class MocGridReader2DTriangle : public MocGridReaderBase {
	public:
		MocGridReader2DTriangle();
		unsigned int get_dimension() const { return 2; }			
	};
	class MocGridReader2DRectangle : public MocGridReaderBase {
	public:
		MocGridReader2DRectangle();
		unsigned int get_dimension() const { return 2; }			
	};
	class MocGridReader3DTetrahedron : public MocGridReaderBase {
	public:
		MocGridReader3DTetrahedron();
		unsigned int get_dimension() const { return 3; }			
	};
	
	/** two line segments at -1.5 0 1.5 */
	MocGridReader1D::MocGridReader1D() {
		// elements		
		elements.resize(2);
		elements[0].resize(2);
		elements[0][0] = 0;
		elements[0][1] = 1;
		elements[1].resize(2);
		elements[1][0] = 1;
		elements[1][1] = 2;	
		// shapes
		element_shapes.push_back(Element::line);
		element_shapes.push_back(Element::line);
		// vertices
		vertices.push_back(vertex_coords(-1.5, 0, 0));
		vertices.push_back(vertex_coords(0, 0, 0));
		vertices.push_back(vertex_coords(1.5, 0, 0));
	}
	
	/** two triangles */
	MocGridReader2DTriangle::MocGridReader2DTriangle() {
		// elements
		elements.resize(2);
		elements[0].push_back(0);
		elements[0].push_back(1);
		elements[0].push_back(2);
		elements[1].push_back(3);
		elements[1].push_back(1);
		elements[1].push_back(2);
		// shapes
		element_shapes.push_back(Element::triangle);
		element_shapes.push_back(Element::triangle);
		// vertices
		vertices.push_back(vertex_coords(0, 0, 0));
		vertices.push_back(vertex_coords(0.5, 0, 0));
		vertices.push_back(vertex_coords(0.5, 0.5, 0));
		vertices.push_back(vertex_coords(1.0, 0.5, 0));
		
	}
	
	/** two rectangles */
	MocGridReader2DRectangle::MocGridReader2DRectangle() {
		// elements
		elements.resize(2);
		elements[0].push_back(0);
		elements[0].push_back(1);
		elements[0].push_back(4);
		elements[0].push_back(5);
		elements[1].push_back(1);
		elements[1].push_back(2);
		elements[1].push_back(3);
		elements[1].push_back(4);	
		// element shapes
		element_shapes.push_back(Element::rectangle);
		element_shapes.push_back(Element::rectangle);
		// vertices
		vertices.push_back(vertex_coords(-1.0, 0.0, 0));
		vertices.push_back(vertex_coords(0.5,  0.0, 0));
		vertices.push_back(vertex_coords(1.5,  0.0, 0));
		vertices.push_back(vertex_coords(1.5,  0.6, 0));
		vertices.push_back(vertex_coords(0.5,  0.6, 0));
		vertices.push_back(vertex_coords(-1.0, 0.6, 0));
	}
	
	/** two elements with common face in yz plane */
	MocGridReader3DTetrahedron::MocGridReader3DTetrahedron() {
		// elements
		elements.resize(2);
		elements[0].push_back(0);
		elements[0].push_back(1);
		elements[0].push_back(2);
		elements[0].push_back(4);
		elements[1].push_back(0);
		elements[1].push_back(2);
		elements[1].push_back(3);
		elements[1].push_back(4);
		// element shapes
		element_shapes.push_back(Element::tetrahedron);
		element_shapes.push_back(Element::tetrahedron);
		// vertices
		vertices.push_back(vertex_coords(0,0,0));
		vertices.push_back(vertex_coords(1.0,0.6,0));
		vertices.push_back(vertex_coords(0,1.2,0));
		vertices.push_back(vertex_coords(-1,0.6,0));
		vertices.push_back(vertex_coords(0,0.6,1));
	}
}; // end of namespace

#endif /*MOC_GRID_READER_H_*/
