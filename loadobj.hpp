#include "rapidobj.hpp"
#include "vec3.hpp"
#include "mat44.hpp"

struct Bbox
{
    Vec3f cornerH,cornerL;

};

struct Node {
    Bbox bbox;
   
    std::vector<Vec3f> positions;
    std::vector<Vec3f> normals;
    std::vector<Node> children;
    Node* parent;
   
    bool enabled; // Only leaf nodes store triangles


};


struct TriangleMesh
{
	std::vector<Vec3f> positions;
	std::vector<Vec3f> normals;
    Bbox box;
};




TriangleMesh load_wavefront_obj( char const* aPath,Mat44f aPreTransform)
{

    auto result = rapidobj::ParseFile( aPath );

    //if( result.error )
    //throw Error( "Unable to load OBJ file ’%s’: %s", aPath, result.error.code.message().c_str() );

    rapidobj::Triangulate( result );

    TriangleMesh ret;

    for( auto const& shape : result.shapes )
    {
        for( std::size_t i = 0; i < shape.mesh.indices.size(); ++i )
        {
            auto const& idx = shape.mesh.indices[i];

            ret.positions.emplace_back( Vec3f{
                result.attributes.positions[idx.position_index*3+0],
                result.attributes.positions[idx.position_index*3+1],
                result.attributes.positions[idx.position_index*3+2]
            });

            
            ret.normals.emplace_back( Vec3f{
                result.attributes.normals[idx.normal_index*3+0],
                result.attributes.normals[idx.normal_index*3+1],
                result.attributes.normals[idx.normal_index*3+2]
            });



        }
        printf("hell0 ,%lu",result.attributes.positions.size());

    }
	for (auto& p : ret.positions)
	{
		Vec4f p4{ p.x, p.y, p.z, 1.f };
		Vec4f t = aPreTransform * p4;
		t /= t.w;

		p = Vec3f{ t.x, t.y, t.z };
	}

    return ret;
    

}
