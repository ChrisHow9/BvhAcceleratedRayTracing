
#include "EasyBMP.h"
#include "vec3.hpp"
#include "loadobj.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream> 
#include <iomanip> 


struct Tri
{
    Vec3f v1,v2,v3;
};

struct Light
{
	Vec3f pos,col;
	Vec3f intensity;
};

struct Color
{
    int r,g,b;
};

struct Ray
{
    Vec3f origin, direction;
};


const int numbOfObjects = 2;
TriangleMesh objects[numbOfObjects];
TriangleMesh mesh;

Tri hellotri;
float fov = 60;
int xRes = 20;
int yRes = 20;

float scaler = tan(fov / 2 * M_PI / 180);
float piDiv = M_PI / 180;
float aspectRatio = xRes/ yRes;

static uint64_t numRayTrianglesTests ;
static uint64_t numRayTrianglesIsect ;
static uint64_t numPrimaryRays ;
static uint64_t bBoxum = 0;






void createInitalBoundingBox(TriangleMesh& mesh)
{
    mesh.box.cornerH = mesh.positions[0];
    mesh.box.cornerL = mesh.positions[0];

    for (int i = 1; i < mesh.positions.size(); i++)
    {
        if (mesh.box.cornerH.x < mesh.positions[i].x) { mesh.box.cornerH.x = mesh.positions[i].x; }
        if (mesh.box.cornerH.y < mesh.positions[i].y) { mesh.box.cornerH.y = mesh.positions[i].y; }
        if (mesh.box.cornerH.z < mesh.positions[i].z) { mesh.box.cornerH.z = mesh.positions[i].z; }

        if (mesh.box.cornerL.x > mesh.positions[i].x) { mesh.box.cornerL.x = mesh.positions[i].x; }
        if (mesh.box.cornerL.y > mesh.positions[i].y) { mesh.box.cornerL.y = mesh.positions[i].y; }
        if (mesh.box.cornerL.z > mesh.positions[i].z) { mesh.box.cornerL.z = mesh.positions[i].z; }
    }
}

bool isVertexInsideBox(const Vec3f& vertex, const Bbox& bbox) {
    return 
    (
    bbox.cornerH.x > vertex.x 
    && bbox.cornerH.y > vertex.y 
    && bbox.cornerH.z > vertex.z 
    && bbox.cornerL.x < vertex.x 
    && bbox.cornerL.y < vertex.y 
    && bbox.cornerL.z < vertex.z
    );
}

//to optimise check parent node's node list
bool boxHasVertex(Node& node, Node& parent)
{
    // Check if the parent node has positions

    for (int i = 0; i < parent.positions.size(); i++)
    {
        if (isVertexInsideBox(parent.positions[i], node.bbox )) {return true;}
    }
    return false;
}

//used for initial bounding box

bool boxAddTris(Node& node) {
    bBoxum++;
    bool addedTriangles = false;
    
    
    for (int i = 0; i < mesh.positions.size(); i += 3) 
    {
        if (isVertexInsideBox(mesh.positions[i], node.bbox) ||
            isVertexInsideBox(mesh.positions[i + 1], node.bbox) ||
            isVertexInsideBox(mesh.positions[i + 2], node.bbox)) 
            {
            node.positions.push_back(mesh.positions[i]);
            node.positions.push_back(mesh.positions[i + 1]);
            node.positions.push_back(mesh.positions[i + 2]);

            node.normals.push_back(mesh.normals[i]);
            node.normals.push_back(mesh.normals[i + 1]);
            node.normals.push_back(mesh.normals[i + 2]);
            addedTriangles = true;
            }
    } 
    return addedTriangles;
}


bool boxAddTris(Node& node, Node& parent) {
    bBoxum++;
    bool addedTriangles = false;
    
    for (int i = 0; i < parent.positions.size(); i += 3) 
    {
        if (isVertexInsideBox(parent.positions[i], node.bbox) ||
            isVertexInsideBox(parent.positions[i + 1], node.bbox) ||
            isVertexInsideBox(parent.positions[i + 2], node.bbox)) 
            {
            node.positions.push_back(parent.positions[i]);
            node.positions.push_back(parent.positions[i + 1]);
            node.positions.push_back(parent.positions[i + 2]);

            node.normals.push_back(parent.normals[i]);
            node.normals.push_back(parent.normals[i + 1]);
            node.normals.push_back(parent.normals[i + 2]);
            addedTriangles = true;
            }
    } 
    return addedTriangles;
}


/*  //old implentation of function, more readable than replacement
void divideBbox(Node& node, int i, int subdiv)
{
    //will -z axis cause issue?
    Node node1,node2,node3,node4;
    Vec3f lowest = node.bbox.cornerL;
    Vec3f highest = node.bbox.cornerH;

    Vec3f midBox =(highest + lowest)/2;

    //node.children.clear();

    // dividing the node into 4 sub nodes, a more accurate solution would be 8.
    node1.bbox.cornerH = {midBox.x,highest.y, midBox.z} ;
    node1.bbox.cornerL = node.bbox.cornerL;
    if (boxHasVertex(node1,node))
    {
        boxAddTris(node1,node);
        node.children.emplace_back(node1);
    }
    

    node2.bbox.cornerH = {highest.x,highest.y,midBox.z} ;
    node2.bbox.cornerL = {midBox.x, lowest.y, lowest.z};
    if (boxHasVertex(node2,node))
    {
        boxAddTris(node2,node);
        node.children.emplace_back(node2);
    }

    node3.bbox.cornerH = {midBox.x, highest.y,highest.x};
    node3.bbox.cornerL = {lowest.x, lowest.y, midBox.z};
    if (boxHasVertex(node3,node))
    {
        boxAddTris(node3,node);
        node.children.emplace_back(node3); 
    }

    node4.bbox.cornerH = highest;
    node4.bbox.cornerL = midBox;
    if (boxHasVertex(node4,node))
    {
        boxAddTris(node4,node);
        node.children.emplace_back(node4);
    }


    ret;urn ;

} */

void divideBbox(Node& node, int i, int subdiv) {

    Vec3f lowest = node.bbox.cornerL;
    Vec3f highest = node.bbox.cornerH;
    Vec3f midBox = (highest + lowest) / 2;
    float padding = 1.1;

    for (int i = 0; i < 2; ++i) 
    {
        for (int j = 0; j < 2; ++j) 
        {
            for (int k = 0; k < 2; ++k) 
            {
                Node subNode;
                subNode.parent = &node;
                subNode.bbox.cornerL = {i ? midBox.x : lowest.x,
                                        j ? midBox.y : lowest.y,
                                        k ? midBox.z : lowest.z};

                subNode.bbox.cornerH = {i ? highest.x : midBox.x,
                                        j ? highest.y : midBox.y,
                                        k ? highest.z : midBox.z};

                if (boxHasVertex(subNode,node)) 
                {
                    boxAddTris(subNode,node);;
                    node.children.emplace_back(subNode);
                }
            }
        }
    }

}



void buildBVHTree(Node& node, TriangleMesh& mesh, int subdiv)
{
    // simulate no bounding box- will still add a bounding box check to every ray
    node.enabled = true;
    if (subdiv == -1) {node.enabled = false;}

    createInitalBoundingBox(mesh);
    
    node.bbox = mesh.box;
    boxAddTris(node);
    
    if (subdiv == 0 || subdiv == -1) 
    {
        node.parent = &node; 
        return;
    }
    
    node.parent = &node;
    std::vector<Node*> currentLayerNodes = {&node};

    for (int i = 0; i < subdiv; i++) 
    {
        
        std::vector<Node*> nextLayerNodes;

        for (auto* currentLayerNode : currentLayerNodes) 
        {
            divideBbox(*currentLayerNode,i,subdiv);

            for ( auto& child : currentLayerNode->children) 
            {
                nextLayerNodes.emplace_back(&child); // Add each child individually
            }
        }
        currentLayerNodes = nextLayerNodes;
    }

}


//from https://tavianator.com/2011/ray_box.html
bool boxIntersection(Bbox b, Ray r) {
    double tmin = -INFINITY, tmax = INFINITY;

    for (int i = 0; i < 3; ++i)
    {
        if (r.direction[i] != 0.0)
        {
            double t1 = (b.cornerL[i] - r.origin[i]) / r.direction[i];
            double t2 = (b.cornerH[i] - r.origin[i]) / r.direction[i];

            // Handle NaN cases
            if (std::isnan(t1) || std::isnan(t2)) {return false;}

            tmin = std::max(tmin, std::min(t1, t2));
            tmax = std::min(tmax, std::max(t1, t2));
        } else if (r.origin[i] < b.cornerL[i] || r.origin[i] > b.cornerH[i]) 
        {
            return false;
        }
    }

    return tmax >= tmin && tmax >= 0.0;
}

// Moeller Trumbore intersection algorithm
bool triangleIntersect(Ray ray, Tri tr, Tri vNormals, Vec3f& normal, Vec3f& weights, Vec3f& intersectpoint)
{
    
    numRayTrianglesTests++;


    float  eps = 0.00001;
    Vec3f e1 = tr.v2 - tr.v1;
    Vec3f e2 = tr.v3 - tr.v1;
    //normal = cross_product(e1, e2);
    Vec3f pvec = cross_product(ray.direction, e2);

    float det = dot(e1, pvec);

    if (det < eps) { return false; }

    float inv_det = 1.0 / det;


    Vec3f tvec = ray.origin - tr.v1;
    float u = dot(tvec, pvec) * inv_det;

    weights.x = u;

    if (u < 0.f || u >1.f) { return false; }

    Vec3f qvec = cross_product(tvec, e1);
    float  v = dot(ray.direction, qvec) * inv_det;

    weights.y = v;

    if (v < 0.f || u + v >1.f) { return false; }

    float t = dot(e2, qvec) * inv_det;


    if (t < eps) { return false; }
    intersectpoint = ray.origin + (ray.direction *t); 

    weights.z = 1 - u - v;


    //normal interploation for smooth shading
    
    Vec3f normalv1 = weights.x * vNormals.v2;
    Vec3f normalv2 = weights.y * vNormals.v3;
    Vec3f normalv3 = weights.z * vNormals.v1;

    normal = normalv1 + normalv2 + normalv3;


    numRayTrianglesIsect++;

    return true;

}

bool BVHIntersection(Node& currentNode, Node& intersectingNode,  Ray& r)
{
    if (!boxIntersection(currentNode.bbox, r)) {return false;} // No intersection with this node
    

    if (currentNode.children.empty()) 
    {
        // Update intersecting node
        intersectingNode = currentNode;
         // Leaf node reached
        return true; 
    }

    // Recursively traverse child nodes
    bool found = false;
    
    for (Node& child : currentNode.children) 
    {
        if (BVHIntersection(child, intersectingNode, r)) {found = true;}
    }

    return found;
}



bool objectIntersect(Ray& ray, Vec3f& normal, Vec3f& weights, Vec3f& intersectpoint, Node& rootNode)
{
    Vec3f normalN, weightsN, intersectN;
    bool hit = false;
    
    Node smallestNode;

    //enable bounding box 
    if (rootNode.enabled)
    {
        if (!BVHIntersection(rootNode, smallestNode, ray)){return hit;}

    }
    else
    {
        smallestNode.positions = rootNode.positions; 
        smallestNode.normals = rootNode.normals;
    }

    //for debugging, return true if intersects bounding box

    //return true;

    float recentIntersection = -1000000;
    for (int i = 0; i< smallestNode.positions.size(); i+=3)

    {
        normalN = {0,0,0};
        weightsN = {0,0,0} ;
        intersectN = {0,0,0};
        Tri tri = {smallestNode.positions[i],smallestNode.positions[i+1],smallestNode.positions[i+2]};
        //not best pratice for readility  but will reuse triangle data structure for vertx normals
        Tri vNormals = {smallestNode.normals[i],smallestNode.normals[i+1],smallestNode.normals[i+2]};
        

        /* calulating normals for backface triangle culling, unsure of performance gain
        Vec3f e1 = tri.v2 - tri.v1;
        Vec3f e2 = tri.v3 - tri.v1;
        Vec3f faceNormal = cross_product(e1, e2);
        faceNormal =  normalise(faceNormal);

        //check if ray direction and normal are in the same direction to cull back facing triagnle
        if (dot(ray.direction, faceNormal) > 0) {continue;}
        */

        if (triangleIntersect(ray, tri, vNormals, normalN, weightsN, intersectN)== true &&  intersectN.z >  recentIntersection )
        {
            recentIntersection = intersectN.z;
            normal = normalN;
            weights = weightsN;
            intersectpoint = intersectN;
            hit = true;
                
        }
        
    }

    return hit;

    
}

Color renderPixel(int x, int y, Node& rootNode)
{
    

    //converting from screen to viewport

    float px = (2 * ((x + 0.5) / xRes) - 1) * scaler * aspectRatio;
    float py = (1 - 2 * (y + 0.5) / yRes) * scaler ;
    
    //Create casting ray
    Ray r;
    r.origin = {0,0,0};
    r.direction = Vec3f{px,py,-1} - r.origin;
    r.direction = normalise(r.direction);

    numPrimaryRays++;

    //init values for triangles
    Vec3f normal, weights , intersectPoint;
   
    Light light;
    light.pos = { 0,0.05,0 };
    light.intensity = { 0.9,0.9,0.9};
    light.col = { 100,100,255 };
    Vec3f amb = {0.1,0.1,0.1};
    Color pixelCol;

    for (int i = 0; i < numbOfObjects; i++  )
    {
        normal = {0,0,0};
        weights = {0,0,0};
        intersectPoint = {0,0,0};
        float recentIntersection = 10000;
        


        if (objectIntersect(r, normal, weights, intersectPoint, rootNode) == false)
        {
            //return background color
            pixelCol = {int(0),int(0),int(0)};

            return pixelCol;

        }

        //for debugging intersections of bounding box
        //pixelCol = {int(255),int(0),int(0)}; return pixelCol;

        

        Vec3f lightdir = light.pos - intersectPoint;
        lightdir = normalise(lightdir);
         
        normal = normalise(normal);
        float diffuse= std::max(0.f, dot(lightdir, normal));
        Vec3f viewdir = normalise(r.origin - intersectPoint);
        Vec3f halfway = normalise( lightdir +viewdir);

        float specularangle = std::max(0.f, dot(halfway, normal));
        float specular = pow(specularangle, 64); 
        float specular_intensity = 0.1;
        Vec3f specular_color = light.intensity * specular * specular_intensity;
        Vec3f lightam = (light.intensity * diffuse) + specular_color + amb; 
        
        pixelCol = 
        {
            //pixel will be black if value exceesd 255
            int(255 * clip(lightam.x,0,1)),
            int(255 * clip(lightam.y,0,1)),
            int(255 * clip(lightam.z,0,1))
        };   
        return pixelCol;        
  
        /* for debug
        std::cout << pixelCol.r << '\n';
        std::cout << pixelCol.g << '\n';
        std::cout << pixelCol.b << '\n';
        */

    }
    return{0,0,0};
}

void writeDataToCSV(const char* csvFileName, int bBoxLevel, float buildBvhTime, float renderTime, float totalTime, size_t numTriangles, unsigned long long numPrimaryRays, unsigned long long numRayTrianglesTests, unsigned long long numBoundingBoxes) {
    std::ofstream csvFile(csvFileName, std::ios::app); // Open the CSV file for appending

    if (!csvFile.is_open()) 
    {
        std::cerr << "Error: Unable to open CSV file." << std::endl;
        return;
    }
    csvFile << std::fixed << std::setprecision(4);

    // Write CSV data
    csvFile << bBoxLevel << "," << buildBvhTime << "," << renderTime << "," << totalTime << "," << numTriangles << "," << numPrimaryRays << "," << numRayTrianglesTests << "," << numBoundingBoxes << std::endl;

    csvFile.close(); // Close the CSV file
}


void render(int xRes, int yRes, const char* outputFileName, Mat44f matrix, const char* objFileName, int bBoxLevel)
{

    mesh = load_wavefront_obj(objFileName, matrix);
    Node rootNode;

    //------Build BVH structure-------
    clock_t timeStartBuildBvh= clock();

    buildBVHTree(rootNode, mesh, bBoxLevel);
    objects[0] = mesh;

    clock_t timeEndBuildBvh= clock();
    //---------------------------------


    //------Render Image---------------

    clock_t timeStartRender = clock();
    

    EasyBMP::Image img(xRes, yRes, outputFileName, EasyBMP::RGBColor(255, 255, 0));

    for (int y = 0; y < yRes; y++)
    {
        
        for (int x = 0; x < xRes; x++)
        {
            Color pixel = renderPixel(x, y, rootNode);
            img.SetPixel(x, y, EasyBMP::RGBColor(pixel.r, pixel.g, pixel.b));
            
        }
    }

    img.Write();

    clock_t timeEndRender = clock();

    //----------------------------------


    //-----------Results----------------
    printf("Results of bvh level %d:\n",bBoxLevel);
    printf("Time to build Bvh                           : %04.4f (sec)\n", (float)(timeEndBuildBvh - timeStartBuildBvh) / CLOCKS_PER_SEC);
    printf("Time to render image                        : %04.4f (sec)\n", (float)(timeEndRender - timeStartRender) / CLOCKS_PER_SEC);
    printf("Toatal time                                 : %04.2f (sec)\n", (float)(timeEndRender - timeStartBuildBvh) / CLOCKS_PER_SEC);
    printf("Total number of triangles                   : %lu\n", mesh.positions.size() / 3);
    printf("Total number of primary rays                : %llu\n", numPrimaryRays);
    printf("Total number of ray-triangles tests         : %llu\n", numRayTrianglesTests);
    printf("Total number of bounding boxes              : %llu\n", bBoxum);
    printf("\n");

    writeDataToCSV("results.csv", 
    bBoxLevel, 
    (float)(timeEndBuildBvh - timeStartBuildBvh) / CLOCKS_PER_SEC, 
    (float)(timeEndRender - timeStartRender) / CLOCKS_PER_SEC, 
    (float)(timeEndRender - timeStartBuildBvh) / CLOCKS_PER_SEC, 
    mesh.positions.size() / 3, 
    numPrimaryRays, 
    numRayTrianglesTests, 
    bBoxum);
}


int main()
{
    xRes = 100;  // Adjust resolution as needed
    yRes = 100;  // Adjust resolution as needed
    Mat44f matrix = make_translation({0, -1, -4}) * make_scaling(0.1, 0.1, 0.1) ;
    const char* objFileName = "teapot.obj";
    int bBoxLevelTest = 7;

    
    //test cases -----
    
    for (int i = -1; i< bBoxLevelTest; i++)
    {
        numRayTrianglesTests = 0;
        numPrimaryRays = 0;
        bBoxum = 0; 

        render(xRes, yRes, "test.bmp", matrix, objFileName ,i);
    } 

    //final render ------
    render(xRes, yRes, "test.bmp", matrix, objFileName ,0);

    return 0;
}