///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2016 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////

//#include <openvdb/viewer/Viewer.h>
#include "viewer/Viewer.h"
#include <boost/algorithm/string/classification.hpp> // for boost::is_any_of()
#include <boost/algorithm/string/predicate.hpp> // for boost::starts_with()
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <openvdb/tools/GridOperators.h>
#include "openvdb/Grid.h"
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetRebuild.h>

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility> 
#include <exception>
#ifdef DWA_OPENVDB
#include <logging_base/logging.h>
#include <usagetrack.h>
#endif

float cone_angle = 90;
float num_rays = 500;
bool do_AO = false;
bool use_distance = true;
bool use_curv = false;

//********************** STRINGIFY
//http://www.decompile.com/cpp/faq/file_and_line_error_string.htm
#ifndef STRINGIFY
#define STRINGIFY(x) TOSTRING(x)
#define TOSTRING(x) #x
#endif
//******************************************************************************

void
usage(const char* progName, int status)
{
    (status == EXIT_SUCCESS ? std::cout : std::cerr) <<
        "Usage: " << progName << " file.vdb [file.vdb ...] [options]\n" <<
        "Which: displays OpenVDB grids\n" <<
        "Options:\n" <<
        "    -i                 print grid information\n" <<
        "    -h, -help          print this usage message and exit\n" <<
        "    -version           print version information\n" <<
        "\n" <<
        "Controls:\n" <<
        "    Esc                exit\n" <<
        "    -> (Right)         show next grid\n" <<
        "    <- (Left)          show previous grid\n" <<
        "    1                  toggle tree topology view on/off\n" <<
        "    2                  toggle surface view on/off\n" <<
        "    3                  toggle data view on/off\n" <<
        "    G                  (\"geometry\") look at center of geometry\n" <<
        "    H                  (\"home\") look at origin\n" <<
        "    I                  toggle on-screen grid info on/off\n" <<
        "    left mouse         tumble\n" <<
        "    right mouse        pan\n" <<
        "    mouse wheel        zoom\n" <<
        "\n" <<
        "    X + wheel          move right cut plane\n" <<
        "    Shift + X + wheel  move left cut plane\n" <<
        "    Y + wheel          move top cut plane\n" <<
        "    Shift + Y + wheel  move bottom cut plane\n" <<
        "    Z + wheel          move front cut plane\n" <<
        "    Shift + Z + wheel  move back cut plane\n" <<
        "    Ctrl + X + wheel   move both X cut planes\n" <<
        "    Ctrl + Y + wheel   move both Y cut planes\n" <<
        "    Ctrl + Z + wheel   move both Z cut planes\n";
    exit(status);
}


////////////////////////////////////////


int
main(int argc, char *argv[])
{
#ifdef DWA_OPENVDB
    USAGETRACK_report_basic_tool_usage(argc, argv, /*duration=*/0);
    logging_base::configure(argc, argv);
#endif

    const char* progName = argv[0];
    if (const char* ptr = ::strrchr(progName, '/')) progName = ptr + 1;

    

    int status = EXIT_SUCCESS;

    try {
        openvdb::initialize();

        bool printInfo = false, printGLInfo = false, printVersionInfo = false;

        // Parse the command line.
        std::vector<std::string> filenames;
        for (int n = 1; n < argc; ++n) {
            std::string str(argv[n]);           
            if (str[0] != '-') {
               //str = STRINGIFY(INPUT_DIR) + std::string(argv[n]);
                filenames.push_back(STRINGIFY(INPUT_DIR) + std::string("/") + str);
            }
            else if (str == "-i") {
                printInfo = true;
            }
            else if (str == "-d") { // deprecated
                printGLInfo = true;
            }
            else if (str == "-h" || str == "-help" || str == "--help") {
                usage(progName, EXIT_SUCCESS);
            }
            else if (str == "-version" || str == "--version") {
                printVersionInfo = true;
                printGLInfo = true;
            }
            else if (str == "-ao") {
                do_AO = true;
                std::cout << " Ambient Occlusion enabled " << std::endl;
            }
            else if (str == "-ray") {
                num_rays = atoi(argv[n+1]);
                n++;
                std::cout << " num_rays = " << num_rays << std::endl;
            }
            else if (str == "-ang") {
                cone_angle = atoi(argv[n + 1]);
                n++;
                std::cout <<" cone_angle = "<<cone_angle << std::endl;
            }
            else {
                std::cout << str << std::endl;
                usage(progName, EXIT_FAILURE);
            }
        }

        const size_t numFiles = filenames.size();

        if (printVersionInfo) {
            std::cout << "OpenVDB library version: "
                << openvdb::getLibraryVersionString() << "\n";
            std::cout << "OpenVDB file format version: "
                << openvdb::OPENVDB_FILE_VERSION << std::endl;
            // If there are no files to view, don't print the OpenGL version,
            // since that would require opening a viewer window.
            if (numFiles == 0) return EXIT_SUCCESS;
        }
        if (numFiles == 0 && !printGLInfo) usage(progName, EXIT_FAILURE);

        openvdb_viewer::Viewer viewer = openvdb_viewer::init(progName, /*bg=*/false);

        if (printGLInfo) {
            // Now that the viewer window is open, we can get the OpenGL version, if requested.
            if (!printVersionInfo) {
                // Preserve the behavior of the deprecated -d option.
                std::cout << viewer.getVersionString() << std::endl;
            }
            else {
                // Print OpenGL and GLFW versions.
                std::ostringstream ostr;
                ostr << viewer.getVersionString(); // returns comma-separated list of versions
                const std::string s = ostr.str();
                std::vector<std::string> elems;
                boost::split(elems, s, boost::algorithm::is_any_of(","));
                for (size_t i = 0; i < elems.size(); ++i) {
                    boost::trim(elems[i]);
                    // Don't print the OpenVDB library version again.
                    if (!boost::starts_with(elems[i], "OpenVDB:")) {
                        std::cout << elems[i] << std::endl;
                    }
                }
            }
            if (numFiles == 0) return EXIT_SUCCESS;
        }

        openvdb::GridCPtrVec allGrids;          
        
        // Load VDB files.
        std::string indent(numFiles == 1 ? "" : "    ");
        for (size_t n = 0; n < numFiles; ++n) {
            openvdb::io::File file(filenames[n]);
            std::cout << "openning: " << filenames[n] << std::endl;
            file.open();

            openvdb::GridPtrVecPtr grids = file.getGrids();
                        

            if (grids->empty()) {
                OPENVDB_LOG_WARN(filenames[n] << " is empty");
                continue;
            }
            allGrids.insert(allGrids.end(), grids->begin(), grids->end());

          
            if (printInfo) {
                if (numFiles > 1) std::cout << filenames[n] << ":\n";
                for (size_t i = 0; i < grids->size(); ++i) {
                    const std::string name = (*grids)[i]->getName();
                    openvdb::Coord dim = (*grids)[i]->evalActiveVoxelDim();
                    std::cout << indent << (name.empty() ? "<unnamed>" : name)
                        << " (" << dim[0] << " x " << dim[1] << " x " << dim[2]
                        << " voxels)" << std::endl;
                }
            }
        }

      
        
        for (size_t n = 0; n < numFiles; ++n) {
            auto base_grid_ptr = allGrids[n];
                        
                        
            auto grid_ptr = openvdb::gridConstPtrCast<openvdb::FloatGrid>
                (base_grid_ptr);



            /*openvdb::io::File re_file("explosion_rebuilt.vdb");
            auto grid_rebuilt = openvdb::tools::levelSetRebuild(
                *grid_ptr, 
                grid_ptr->getGridClass() == openvdb::GRID_LEVEL_SET ? 0.0 : 0.01,
                5);
            openvdb::GridPtrVec grids_rebuilt;
            grids_rebuilt.push_back(grid_rebuilt);
            re_file.write(grids_rebuilt);
            re_file.close();

            exit(0);*/
           

            for (auto iter = grid_ptr->beginMeta();
                iter != grid_ptr->endMeta(); ++iter) {
                const std::string& name = iter->first;
                openvdb::Metadata::Ptr value = iter->second;
                std::string valueAsString = value->str();
                std::cout << name << " = " << valueAsString << std::endl;
            }


            float bg_value = grid_ptr->background();
            std::cout << "\n grid_ptr->background() = " 
                << bg_value << std::endl;

            std::cout << " grid_ptr->metaCount() = "<< grid_ptr->metaCount() 
                <<std::endl;

            std::cout << " grid_ptr->hasUniformVoxels() = " 
                << grid_ptr->hasUniformVoxels() << std::endl;
            std::cout << " grid_ptr->voxelSize() = "
                << grid_ptr->voxelSize() << std::endl;

            auto tree_ptr = grid_ptr->treePtr();            
            std::cout << std::endl << std::endl;
            std::cout << "tree_ptr->activeLeafVoxelCount() = "
                << tree_ptr->activeLeafVoxelCount() << std::endl;
            std::cout << "tree_ptr->activeTileCount() = "
                << tree_ptr->activeTileCount() << std::endl;
            std::cout << "tree_ptr->activeVoxelCount() = "
                << tree_ptr->activeVoxelCount() << std::endl;
            std::cout << "tree_ptr->getBackgroundValue() = "
                << tree_ptr->getBackgroundValue()->str() << std::endl;
            std::cout << "tree_ptr->nonLeafCount() = "
                << tree_ptr->nonLeafCount() << std::endl;
            std::cout << "tree_ptr->valueType() = " <<
                tree_ptr->valueType() << std::endl;
            std::cout << "tree_ptr->treeDepth() = " <<
                tree_ptr->treeDepth() << std::endl;                      
            std::cout << "grid_ptr->gridType() = " 
                <<grid_ptr->gridType() << std::endl;

            std::cout << std::endl;

            //loop over all leaf nodes 
            /* int num_leaf_node = 0;
            int num_zero_value_leaf = 0;

            openvdb::FloatGrid::ConstAccessor grid_accessor =
                (*grid_ptr).getAccessor();
            openvdb::math::SevenPointStencil<openvdb::FloatGrid>
                stencil(*grid_ptr);
                        
            for(auto leaf_iter =
                grid_ptr->cbeginValueOn(); leaf_iter.test(); ++leaf_iter){
                //aaccess the leaf node and compute its AO
                
                if (abs(leaf_iter.getValue()) <  0.0002) {
                    num_zero_value_leaf++;

                    //Indices of this voxel
                    openvdb::Coord ijk = leaf_iter.getCoord();                  
                    //std::cout << "ijk = " << ijk << std::endl;
                        
                    const std::string ijk_string = ijk.str();

                    //std::cout << "ijk_string = " << ijk_string << std::endl;

                    float voxel_value = grid_accessor.getValue(ijk);
                    //std::cout << " voxel_value = " << voxel_value << std::endl;
                    //value of this voxel
                    //std::cout << " grid_accessor.getValue(ijk)= "
                    //    << grid_accessor.getValue(ijk) << std::endl;                 
                    //std::cout << " leaf_iter.getValue() = " 
                    //    << leaf_iter.getValue() << std::endl;

                    //voxel size of this voxel
                    //std::cout <<" grid_ptr->voxelSize(ijk.asVec3d()) = " 
                    //    << grid_ptr->voxelSize(ijk.asVec3d()) << std::endl;
                    
                    //The world position of the voxel
                    openvdb::Vec3d xyz = grid_ptr->indexToWorld(ijk);
                    //std::cout << " xyz = " << xyz << " -> " << " ijk = " 
                    //    << ijk << std::endl;
                 
                    //compute gradient/nomral of this voxel and normalize it
                    stencil.moveTo(ijk);
                    auto grad_result_cd_2nd =
                        openvdb::math::ISGradient<openvdb::math::CD_2ND>::result(stencil);
                    grad_result_cd_2nd.normalize();
                    //std::cout << " grad_result_cd_2nd= " << grad_result_cd_2nd
                    //    <<" with len = "<< grad_result_cd_2nd.length() 
                    //    <<std::endl;

                    //keeping walking along the normal direction until you fall
                    //into the background grid 
                    auto normal = grad_result_cd_2nd;       
                    float distance_walked = 0;
                    float ao = 0;
                    while(abs(abs(bg_value) - abs(voxel_value)) > 0.0001) {
                        //new world location
                        xyz += normal;
                        distance_walked += 1.0f;

                        //std::cout << " new xyz = " << xyz << std::endl;
                        //new (fractional) index 
                        openvdb::Vec3d new_ijk = grid_ptr->worldToIndex(xyz);
                       // std::cout << " new_ijk = " << new_ijk << std::endl;

                        //use the sampler to the new values from the world
                        //coordinates directly
                        voxel_value = openvdb::tools::BoxSampler::sample
                        (grid_accessor, new_ijk);
                        //std::cout << " voxel_value= " <<
                        //    voxel_value << std::endl;

                        float dist_diff = abs(abs(voxel_value) - distance_walked);

                        if (dist_diff > 0.0001) {
                            //TODO make sure you are not falling in the bg grid
                            ao += dist_diff;
                        }

                        //new grid index (actual voxel)
                        //openvdb::Coord new_ijk_coord(new_ijk.x(), new_ijk.y(),
                        //    new_ijk.z());
                        openvdb::Coord new_ijk_coord;
                         new_ijk_coord.round(new_ijk);
                        //std::cout << " new new_ijk_coord = " << new_ijk_coord 
                        //    << std::endl;
                        ////get the voxel value 
                        //voxel_value = grid_accessor.getValue(new_ijk_coord);
                        //std::cout << " new voxel_val = " << voxel_value
                        //    << std::endl;
                    }

                    //std::cout << std::endl << std::endl;

                   // std::cout << ijk_string << std::endl;

                   

                    AO_map.insert(std::make_pair(ijk_string,ao));
                    
                }
                
                num_leaf_node++;
            }

            std::cout << "\n num_leaf_node = " << num_leaf_node << std::endl;
            std::cout << "\n num_zero_value_leaf = " << num_zero_value_leaf 
                << std::endl;*/
        }

      
        viewer.open();

        viewer.view(allGrids);

        openvdb_viewer::exit();

    }
    catch (const char* s) {
        OPENVDB_LOG_ERROR(progName << ": " << s);
        status = EXIT_FAILURE;
    }
    catch (std::exception& e) {
        OPENVDB_LOG_ERROR(progName << ": " << e.what());
        status = EXIT_FAILURE;
    }
    return status;
}

// Copyright (c) 2012-2016 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
