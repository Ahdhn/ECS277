//#define HELLO_WORLD
#define CREATE_WRITE_GRID
//#define HANDLING_METADATA
//#define ITERATION
//#define INDEXSPACE_SAMPLERS
//#define GRIDSAMPLER
//#define DUALGRID_SAMPLER
//#define GRADIENT

#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/ChangeBackground.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/math/FiniteDifference.h>
#include <openvdb/math/Operators.h>

#include <iostream>


template <class GridType>
void makeSphere(GridType& grid, float radius, const openvdb::Vec3f&c) {
    //populate the given grid with a narrow-band level set representation of 
    //a sphere. The width of the narrow band is determined by the grid's 
    //background value 
    using ValueT = typename GridType::ValueType;

    //distacne value for the constant region exterior to the narrow band 
    const ValueT outside = grid.background();

    //distance value for the constant region interior to the narrow band 
    const ValueT inside = -outside;

    //use the background value as the width in voxels of the narrow band 
    int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
    //the bounding box of the narrow band is 2*dim voxels on a side 
    int dim = int(radius + padding);

    //get a voxel accessor 
    typename GridType::Accessor accessor = grid.getAccessor();

    //compute the signed distance from the surface of the sphere of each voxel
    //within the bounding box and insert the value into the grid 
    openvdb::Coord ijk;
    int &i = ijk[0], &j = ijk[1], &k = ijk[2];
    for (i = c[0] - dim; i < c[0] + dim; ++i) {
        const float x2 = openvdb::math::Pow2(i - c[0]);
        for (j = c[1] - dim; j < c[1] + dim; ++j) {
            const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
            for (k = c[2] - dim; k < c[2] + dim; ++k) {

                //the distance from the sphere surface in voxels
                const float dist = openvdb::math::Sqrt(x2y2 +
                    openvdb::math::Pow2(k - c[2])) - radius;

                //convert the floating-point distance to the grid's valued type
                ValueT val = ValueT(dist);

                //only insert distance that smaller in magnitude then 
                //the backgroound value 
                if (val < inside || outside < val) {
                    continue;
                }

                //set the distance for voxel (i,j,k)
                accessor.setValue(ijk, val);
            }

        }
    }

}

int main() {

    std::cout << " Running openvdb_2\n" << std::endl;

#ifdef HELLO_WORLD
    //init openvdb
    openvdb::initialize();

    //create an empty floating point grid with background value 0
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();

    std::cout << " Testing rand access" << std::endl;

    //Get an accessor for coordinate-based access to voxels 
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    //define a coordinate with larger signed indices
    openvdb::Coord xyz(1000, -200000000, 30000000);

    //set the voxel value at (1000, -200000000, 30000000) to 1
    accessor.setValue(xyz, 1.0);

    //verify that the voxel at (1000, -200000000, 30000000) is 1
    std::cout << " Grid " << xyz << " = " << accessor.getValue(xyz) << std::endl;

    //rest the coordinate to those of a different voxel 
    xyz.reset(1000, 200000000, -30000000);

    //verify that the voxel at (1000, 200000000, -30000000) is the background 
    //value 0
    std::cout << " Grid " << xyz << " = " << accessor.getValue(xyz) << std::endl;

    //set the voxel value at (1000, 200000000, -30000000) to 2
    accessor.setValue(xyz, 2.0);

    //set the voxels at the two extremes of the available coordinates space
    //for 32-bit signed coordinates these are (-2147483648, -2147483648, -2147483648)
    //and (2147483648, 2147483648, 2147483648)
    //std::cout << "openvdb::Coord::min() = " << openvdb::Coord::min() << std::endl;
    accessor.setValue(openvdb::Coord::min(), 3.0f);
    accessor.setValue(openvdb::Coord::max(), 3.0f);

    std::cout << " Testing sequential access" << std::endl;
    //print all active voxels by means of an iterator 
    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter;
        ++iter) {
        std::cout << " Grid " << iter.getCoord() << " = " << *iter << std::endl;

    }
#endif

#ifdef CREATE_WRITE_GRID
    openvdb::initialize();

    if (false) {
        //create FloatGrid-typed grid with 2 as background value 
        openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(2.0);

        //populate the grid with sparse, narrow-band level set representation of a 
        //spherw with radius 50 voxels, located at (1.5,2,3) in the index space 
        makeSphere(*grid, 50.0, openvdb::Vec3f(1.5, 2, 3));

        //associate some metadata with the grid 
        grid->insertMeta("radius", openvdb::FloatMetadata(50.0));

        //associate a scaling transform with the grid that set the voxel size to 
        //0.5 units in world space 
        grid->setTransform(openvdb::math::Transform::createLinearTransform(/*voxel size*/0.5));

        //identify the grid as level set 
        grid->setGridClass(openvdb::GRID_LEVEL_SET);

        //name the grid "levelsetsphere:
        grid->setName("LevelSetSphere");

        //create a vdb file object 
        openvdb::io::File file("mygrids.vdb");

        //add the grid pointer to container 
        openvdb::GridPtrVec grids;
        grids.push_back(grid);

        //write out the content of the container 
        file.write(grids);
        file.close();
    }
    else {
        //using tools 
        openvdb::FloatGrid::Ptr grid =
            openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
                /*radius*/50.0, /*center*/openvdb::Vec3f(0.0, 0.0, 0.0),
                /*voxel size */0.5, /*width*/8.0);
        
        openvdb::CoordBBox ijk_bb = grid->evalActiveVoxelBoundingBox();
        openvdb::Coord ijk = grid->evalActiveVoxelDim();

        //openvdb::tools::changeBackground(grid->tree(), 
        //    std::numeric_limits<float>::max());

        //associate meta data with the grid 
        grid->insertMeta("radius", openvdb::FloatMetadata(50.0));

        //name the grid 
        grid->setName("LevelSetSphere");

        //create vdb file object 
        openvdb::io::File file("my_sphere.vdb");

        //add the grid pointer to a container 
        openvdb::GridPtrVec grids;
        grids.push_back(grid);

        //write out the content of the container
        file.write(grids);
        file.close();

    }



#endif 

#ifdef HANDLING_METADATA
    openvdb::initialize();

    openvdb::Vec3SGrid::Ptr grid = openvdb::Vec3SGrid::create();

    grid->insertMeta("vector type",
        openvdb::StringMetadata("covariant (gradient)"));
    grid->insertMeta("radius", openvdb::FloatMetadata(50.0));
    grid->insertMeta("center",
        openvdb::Vec3SMetadata(openvdb::Vec3s(10, 15, 10)));

    //overwrites existing value 
    grid->insertMeta("center",
        openvdb::Vec3SMetadata(openvdb::Vec3s(10.5, 15, 30)));

    //can not overwrite a value of type vec3s with one of type float
    //grid->insertMeta("center", openvdb::FloatMetadata(0.0));

    std::string s = grid->metaValue<std::string>("vector type");
    float r = grid->metaValue<float>("radius");

    //error can't read a value of type vec3s as float
    //float center = grid->metaValue<float>("center");

    /***********************/
    //use Grid::beginMeta() to get std::map iterator over the metadata
    for (openvdb::MetaMap::MetaIterator iter = grid->beginMeta();
        iter != grid->endMeta(); ++iter) {
        const std::string& name = iter->first;
        openvdb::Metadata::Ptr value = iter->second;
        std::string valueAsString = value->str();
        std::cout << name << " = " << valueAsString << std::endl;
    }
    /***********************/

    //remove metadata
    grid->removeMeta("vector type");
    grid->removeMeta("center");
    grid->removeMeta("vector type");//OK
#endif

#ifdef ITERATION

    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid =
        openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
            /*radius*/50.0, /*center*/openvdb::Vec3f(1.5, 2, 3),
            /*voxel size */0.5, /*width*/4.0);

    // node iterator 
    for (openvdb::FloatGrid::TreeType::NodeIter iter = grid->tree().beginNode();
        iter; ++iter) {
        switch (iter.getDepth()) {
        case 0: {
            //Root
            openvdb::FloatGrid::TreeType::RootNodeType*node = nullptr;
            iter.getNode(node);
            if (node) {
                //play with the node 
            }
            break;
        }
        case 1: {
            //Level 2 - internal node
            openvdb::FloatGrid::TreeType::RootNodeType::ChildNodeType*node = nullptr;
            iter.getNode(node);
            if (node) {
                //play with the node 
            }
            break;
        }
        case 2: {
            //Level 1 - internal node
            openvdb::FloatGrid::TreeType::RootNodeType::ChildNodeType::ChildNodeType *node = nullptr;
            iter.getNode(node);
            if (node) {
                //play with the node 
            }
            break;

        }
        case 3: {
            //Level 0 - leaf node
            openvdb::FloatGrid::TreeType::LeafNodeType *node = nullptr;
            iter.getNode(node);
            if (node) {
                //play with the node 
            }
            break;
        }
        default:
            break;
        }
    }

    //leaf node iterator 
    //iterate over reference to const leafnodes 
    for (openvdb::FloatGrid::TreeType::LeafCIter iter = grid->tree().cbeginLeaf();
        iter; ++iter) {
        const openvdb::FloatGrid::TreeType::LeafNodeType& leaf = *iter;
        //...
    }

    //iterate over reference to non-const leafnode 
    for (openvdb::FloatGrid::TreeType::LeafIter iter = grid->tree().beginLeaf();
        iter; ++iter) {
        const openvdb::FloatGrid::TreeType::LeafNodeType& leaf = *iter;
        //...
    }

    //iterate over pointer to const leadnodes
    for (openvdb::FloatGrid::TreeType::LeafCIter iter = grid->tree().cbeginLeaf();
        iter; ++iter) {
        const openvdb::FloatGrid::TreeType::LeafNodeType* leaf = iter.getLeaf();
    }

    //iterate over pointer to non-const leadnodes
    for (openvdb::FloatGrid::TreeType::LeafIter iter = grid->tree().beginLeaf();
        iter; ++iter) {
        const openvdb::FloatGrid::TreeType::LeafNodeType* leaf = iter.getLeaf();
    }

    /***********************/
    //Value iterator
    //iterate over all active value     
    openvdb::Vec3SGrid::Ptr grid3s = openvdb::Vec3SGrid::create();
    for (openvdb::Vec3SGrid::ValueOnCIter iter = grid3s->cbeginValueOn();
        iter.test(); ++iter) {
        const openvdb::Vec3f&  value = *iter;
        //print coordinates of all voxels whose vector value has a length >10
         //print the bounding box coordinates of all tiles whose vector len
         //>10
        if (value.length() > 10.0) {
            if (iter.isVoxelValue()) {
                std::cout << iter.getCoord() << std::endl;
            }
            else {
                //if it does not have a voxel value, then it mush be a tile 
                openvdb::CoordBBox bbox;
                iter.getBoundingBox(bbox);
                std::cout << bbox << std::endl;
            }
        }
    }

    //Iterate over and normalize all inactive values 
    for (openvdb::Vec3SGrid::ValueOffIter iter = grid3s->beginValueOff();
        iter.test(); ++iter) {
        openvdb::Vec3f value = *iter;
        value.normalize();
        iter.setValue(value);
    }

    //normalize the (inactive) backgrun values as well
    openvdb::tools::changeBackground(grid3s->tree(), grid3s->background().unit());




#endif

#ifdef INDEXSPACE_SAMPLERS
    const openvdb::FloatGrid grid;

    //choose fractional coordinates in the index space
    const openvdb::Vec3R ijk(10.5, -100.2, 50.3);

    //compute the value of the grid at ijk via nearest-neighbour (zero-order)
    //interpolation
    openvdb::FloatGrid::ValueType v0 = openvdb::tools::PointSampler::sample(grid.tree(), ijk);

    //compute the value via trilinear (first order) interpolation 
    openvdb::FloatGrid::ValueType v1 = openvdb::tools::BoxSampler::sample(grid.tree(), ijk);

    //compute with triquadratic interpolation 
    openvdb::FloatGrid::ValueType v2 = openvdb::tools::QuadraticSampler::sample(grid.tree(), ijk);

    //We can get better performance by passing the accessor 
    openvdb::FloatGrid::ConstAccessor accessor = grid.getConstAccessor();
    openvdb::FloatGrid::ValueType v0p = openvdb::tools::PointSampler::sample(accessor, ijk);
    openvdb::FloatGrid::ValueType v1p = openvdb::tools::BoxSampler::sample(accessor, ijk);
    openvdb::FloatGrid::ValueType v2p = openvdb::tools::QuadraticSampler::sample(accessor, ijk);



#endif

#ifdef GRIDSAMPLER
    const openvdb::FloatGrid grid;

    //instaniate the GridSampler template on the grid type and on a box
    //sampler for thread-safe but uncahed trilinear interpolation 
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>
        sampler(grid);

    //compute the value of the gird at frational coordinates in index space
    openvdb::FloatGrid::ValueType indexValue = 
        sampler.isSample(openvdb::Vec3R(10.5, -100.2, 50.3));

    //compute the value of the grid at a location in world space
    openvdb::FloatGrid::ValueType worldValue =
        sampler.wsSample(openvdb::Vec3R(0.25, 1.4, -1.1));

    //request a value accessor for accelerated access 
    //because value accessors employ a cache, it is important to declare one
    //accessor per thread
    openvdb::FloatGrid::ConstAccessor accessor = grid.getConstAccessor();

    //instantiate the GridSampler template on the accessor type and on 
    //a box sampler for accelerated trilinear interpolation 
    openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor,
        openvdb::tools::BoxSampler>fastSampler(accessor, grid.transform());
    
    //compute the value of the grid at fractional coordinates in the index 
    //space 
    indexValue = fastSampler.isSample(openvdb::Vec3R(10.5, -100.2, 50.3));

    //compute the value of the grid at a location in the world space
    worldValue = fastSampler.wsSample(openvdb::Vec3R(0.25, 1.4, -1.1));
#endif

#ifdef DUALGRID_SAMPLER
    const openvdb::FloatGrid sourceGrid;
    const openvdb::FloatGrid targetGrid;

    //instantiate the DualGridSampler template on the grid type and on a box
    //sampler for thread-sade but uncached trilinear iterpolation 
    openvdb::tools::DualGridSampler<openvdb::FloatGrid,
        openvdb::tools::BoxSampler> sampler(sourceGrid,
            targetGrid.constTransform());
    //compute the value of the source grid at a location in the target grid's 
    //index space 
    openvdb::FloatGrid::ValueType value = sampler(openvdb::Coord(-200,-50,202));

    //request a value accessor for accelerated access to the source grid 
    //because value accessor employ a cache, it is important to declare
    //one accessor per thread
    openvdb::FloatGrid::ConstAccessor accessor = sourceGrid.getConstAccessor();

    //instantiate the DualGridSampler templare on the accessor type and on 
    //a box sampler for accelated trilinear interpolation 
    openvdb::tools::DualGridSampler<openvdb::FloatGrid::ConstAccessor,
        openvdb::tools::BoxSampler> fastSampler(accessor,
            sourceGrid.constTransform(), targetGrid.constTransform());

    //compute the value of the source grid at a location in the target 
    //grid's index space
    value = fastSampler(openvdb::Coord(-200, -50, 202));

#endif 

#ifdef GRADIENT
    //https://people.cs.clemson.edu/~jtessen/cpsc8190/OpenVDB-dpawiki.pdf
    //cretae a sphere level set and accessor 
    openvdb::FloatGrid::Ptr grid = openvdb::tools::createLevelSetSphere<
        openvdb::FloatGrid>(50.0, openvdb::Vec3f(0, 0, 0), 0.5);
    openvdb::FloatGrid::ConstAccessor acc = grid->getConstAccessor();
    openvdb::math::Transform::Ptr transform = grid->transformPtr();
   
    openvdb::math::UniformScaleMap::Ptr map = 
        transform->map<openvdb::math::UniformScaleMap>();

    //create a gradient operator using a uniform scale map and 2nd order 
    //center difference finite difference scheme
    //computation is done in the range space of the map (world space)
    openvdb::math::Gradient<openvdb::math::UniformScaleMap,
        openvdb::math::CD_2ND> gradop;
    //transform a worldspace coordinates (outside the level set) into index
    //space 
    openvdb::Vec3f pos(50.5, 0.5, 0.5);
    openvdb::Coord ijk(transform->worldToIndexNodeCentered(pos));

    //find the closes point on the surface of the level set 
    std::cout << pos << " grad: " << gradop.result((*map), acc, ijk) << std::endl;


#endif 

    return 0;

}
