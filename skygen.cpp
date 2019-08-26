#include "skygen.h"


using std::string;
using std::cout;
using std::getline;
using std::vector;
using std::endl;
namespace filesys = std::experimental::filesystem;

skygen::skygen(const std::string& terrainFile_, const std::string& latLonListFullPath_,
               const std::string &imageRootPath_, const std::string &imageExt_,
               const std::string &imageFolderPrefix_,
               double &width_, double &height_, double &numResamples_, osg::Vec4f &backgroundColor_,
               double &fovy_, double &aspectRatio_, double &zNear_, double &zFar_,
               double &heightAboveTerrainDesired_, double &maxElevationInDatabase_,
               double &LODScale_, double &SmallFeatureCullingPixelSize_)
            :terrainFile(terrainFile_),latLonListFullPath(latLonListFullPath_),
              imageRootPath(imageRootPath_), imageExt(imageExt_),
              imageFolderPrefix(imageFolderPrefix_),
              width(width_), height(height_), numResamples(numResamples_), backgroundColor(backgroundColor_),
              fovy(fovy_), aspectRatio(aspectRatio_), zNear(zNear_), zFar(zFar_),
              heightAboveTerrainDesired(heightAboveTerrainDesired_),
              maxElevationInDatabase(maxElevationInDatabase_),
              LODScale(LODScale_), SmallFeatureCullingPixelSize(SmallFeatureCullingPixelSize_){
}

int skygen::setCameraFacing(const string& cameraFacing_){
    if(cameraFacing_ == "north" || cameraFacing_ == "south" || cameraFacing_ == "east" || cameraFacing_ =="west" || cameraFacing_ == "up" || cameraFacing_ == "down"){
        cameraFacing = cameraFacing_;
        return 0;
    }else{
    cout << "ERROR: Error in skygen::setCameraFacing(): Direction invalid. Please enter a valid direction (\"north\", \"south\", \"east\" or \"west\")";
    return -1;
    }
}

int skygen::convertStartStopPercentToIndex(const int startPercent_,const int stopPercent_){
    size_t len = numAllStrands;

    startIndex = size_t(len*startPercent_/100);
    stopIndex = size_t(len*stopPercent_/100);
    //cout << "(startPercent,stopPercent) = (" << startPercent_ <<", " << stopPercent_<<")"<<endl;
    //cout << "(len,startIndex,stopIndex) = ("<< len<<", "<<startIndex<<", "<<stopIndex<<")"<<endl;

    return 0;
}


int skygen::setActiveStrands(){
    //startIndex = 2381;
    //stopIndex = 363;

    if(startIndex>numAllStrands || startIndex < 0){
        cout << "ERROR: Error in setActiveStrands(): startIndex aceptable range is (0, " << numAllStrands << "). However it was calculated to be = " << startIndex<<endl;
        return -12;
    }
    if(stopIndex>numAllStrands || stopIndex < 0){
        cout << "ERROR: Error in setActiveStrands(): stopIndex aceptable range is (0, " << numAllStrands << "). However it was calculated to be = " << startIndex<<endl;
        return -13;
    }
    if(stopIndex < startIndex){
        cout << "ERROR: Error in setActiveStrands(): stopIndex is less than startIndex. (startIndex, stopIndex) = ("<<startIndex<<", " <<stopIndex<<")"<<endl;
        return -14;
    }
    //startIndex and stopIndex are valid. Set strands to be rendered.
    strands.clear();
    if(cameraFacing == "north" || cameraFacing == "south" || cameraFacing == "up" || cameraFacing == "down"){
        for(size_t ii = startIndex; ii<stopIndex; ++ii){
            strands.push_back(northSouthStrands[ii]);
        }
        if(cameraFacing == "south"){
            //reverse the order of the strands so they're naturally rendered going from north to south
            size_t numStrands = strands.size();
            for(size_t ii = 0; ii < numStrands; ++ii){
                vector<geoPoint>& currStrand = strands[ii];
                std::reverse(currStrand.begin(),currStrand.end());
            }
        }
    }
    if(cameraFacing == "east" || cameraFacing == "west"){
        for(size_t ii = startIndex; ii<stopIndex; ++ii){
            strands.push_back(eastWestStrands[ii]);
        }
        if(cameraFacing == "west"){
            //reverse the order of the strands to they're naturally rendered going from east to west
            size_t numStrands = strands.size();
            for(size_t ii = 0; ii < numStrands; ++ii){
                vector<geoPoint>& currStrand = strands[ii];
                std::reverse(currStrand.begin(),currStrand.end());
            }
        }
    }
    return 0;

}


size_t skygen::getStrandPercent(const size_t& currStrandIndex){
    size_t csPercent =  (currStrandIndex*100)/numAllStrands;
    cout << "numAllStrands = " << numAllStrands <<endl;
    cout << "currStrandPercent = " << csPercent << endl;
    return csPercent;
}

int skygen::preRenderSetup(const std::string &directionMode_){
    int errCode73 = setCameraFacing(directionMode_);
    if(errCode73<0){
        cout << "ERROR: Error in Skygen::preRenderSetup(). Invalid direction."<<endl;
        return -3;
    }
    if(cameraFacing == "north" || cameraFacing == "south" || cameraFacing == "up" || cameraFacing == "down"){
        numAllStrands = northSouthStrands.size();
    }
    if(cameraFacing == "east" || cameraFacing == "west"){
        numAllStrands = eastWestStrands.size();
    }
    return 0;
}

int skygen::renderWrapper(){
    //Likewise, the direction mode along with the start/stop indexes from above
    //determines which vector<geoPoints> we're interested in rendering.
    //So we preprocesses the strands to be in appropriate order for rendering.
    if(subLatLonList.size()>0){
        int errCode83 = setActiveStrands_SubLatLonList();
        if(errCode83 < 0){
            cout << "ERROR: Error in Skygen::setActiveStrands_SubLatLonList. Returned errCode = " <<errCode83<<endl;
            return -103;
        }
    }else{
        int errCode84 = setActiveStrands();
        if(errCode84 < 0){
            cout << "ERROR: Error in Skygen::setActiveStrands. Returned errCode = " <<errCode84<<endl;
            return -104;
        }
    }

    //render all the geoPoints for all the now-active strands
    size_t len = 0;
    size_t numStrands = strands.size();
    qPressedNonBlocking qPressedHandler;
    for(size_t ss = 0; ss < numStrands; ++ss){
        strandIndex = ss;
        vector<geoPoint>& strand = strands[strandIndex];
        len = strand.size();
        if(len>0){
            locIndex = 0;
            lat_deg = strand[0].lat;
            lon_deg = strand[0].lon;
            //be sure directory to save images is set
            size_t currentIndex = startIndex+strandIndex;
            int currStrandRenderPercent = getStrandPercent(currentIndex);
            cout << "Rendering Strand: (Start, Current, End) =  ("<< startIndex
                 << ", "<<currentIndex << ", " <<stopIndex-1 <<")"<<endl;
            if(subLatLonList.size()){

                currentIndex = S_subLatLonList[frameCount];
                if(cameraFacing=="south" || "east"){
                    currStrandRenderPercent = (currentIndex*100)/eastWestStrands.size();
                }else{
                    currStrandRenderPercent = (currentIndex*100)/northSouthStrands.size();
                }
            }

            std::string saveDirRoot = imageRootPath + "/"+imageFolderPrefix+ std::to_string(currStrandRenderPercent);
            std::string saveDir = saveDirRoot + "/strand" + std::to_string(currentIndex);
            if(qPressedHandler.qPressed()){
                if(subLatLonList.size()){
                    cout << "Pressed 'q'.  To continue, begin at strand: " <<startIndex+strandIndex<<endl;
                }else{
                    cout << "Pressed 'q'.  To continue, begin at strand: " <<currentIndex<<endl;
                    return 0;
                }
            }
            if(!filesys::is_directory(saveDir)){
                bool mkDirOk = filesys::create_directories(saveDir);
                if(!mkDirOk){
                    cout << "ERROR: Error in skygen::renderImages(). Failed to create saveDir = "
                         << saveDir <<endl;
                    return -846;
                }
            }
            int dirOk = setImageSaveDir(saveDir);
            if(dirOk<0){
                cout << "ERROR: Error in Skygen::renderImages(). Failed to setImageSaveDir()." <<endl;
                return -843;
            }
            //GOGOGO
            renderFirstImage();
        }
        for(size_t ll = 1; ll<len; ++ll){
            locIndex = ll;
            lat_deg = strand[locIndex].lat;
            lon_deg = strand[locIndex].lon;
            renderNextImage();
        }
    }
    return 0;
}


int skygen::renderImagesAutoAll(const std::string& directionMode_){
    cout << "ERROR: skygen::renderImagesAutoAll(). This method is not yet implemented."<<endl;
    return -1;
}
int skygen::renderImagesStartAtStrand(const std::string& directionMode_, const size_t& startStrand_){
    int errPreRender = preRenderSetup(directionMode_);
    if(errPreRender<0){
        cout << "ERROR: Error in skygen::renderImagesStartAtStrand(). Failed PreRenderSetup."<<endl;
        return -33;
    }
    startIndex = startStrand_;
    stopIndex = numAllStrands;
    int errCodeRenderWrapper = renderWrapper();
    if(errCodeRenderWrapper<0){
        cout << "ERROR: Error in skygen::renderImagesStartAtStrand(). RenderWrapper returned errCode = "<<errCodeRenderWrapper<<endl;
        return -837;
    }
    return 0;
}

int skygen::renderImagesByStrand(const std::string& directionMode_,
                         const size_t& startStrand_,
                         const size_t& stopStrand){
    int errPreRender = preRenderSetup(directionMode_);
    if(errPreRender<0){
        cout << "ERROR: Error in skygen::renderImagesByStrand(). Failed PreRenderSetup."<<endl;
        return -33;
    }
    startIndex = startStrand_;
    stopIndex = stopStrand;
    int errCodeRenderWrapper = renderWrapper();
    if(errCodeRenderWrapper<0){
        cout << "ERROR: Error in skygen::renderImagesByStrand(). RenderWrapper returned errCode = "<<errCodeRenderWrapper<<endl;
        return -837;
    }
    return 0;
}

//public member function renderImages
//This function is primary loop(s) that render and save all the images of interest.
//It starts by
//std::string directionMode_ defines the camera heading. (i.e. "north" "south" "east" or "west")
//const int startPercent_ and stopPercent_ are used to have a bit more
//fine control over which images are rendered... input should be an int between 0 and 100.
int skygen::renderImagesByPercent(const std::string &directionMode_, const int startPercent_, const int stopPercent_){
    int errPreRender = preRenderSetup(directionMode_);
    if(errPreRender<0){
        cout << "ERROR: Error in skygen::renderImagesByPercent(). Failed PreRenderSetup."<<endl;
        return -33;
    }
    //ensure that startPercent_ and stopPercent_ are valid entries
    if(typeid(startPercent_).name() != typeid(int).name()){
        cout << "ERROR: Error in Skygen::renderImages. Invalid startPercent_. Please enter an integer between 0 and 100." <<endl;
        return -4;
    }
    if(typeid(stopPercent_).name() != typeid(int).name()){
        cout << "ERROR: Error in Skygen::renderImages. Invalid stopPercent_. Please enter an integer between 0 and 100." <<endl;
        return -5;
    }
    if(startPercent_ < 0 || startPercent_ > 100){
        cout << "ERROR: Error in Skygen::renderImages. Invalid range for startPercent_. Please enter an integer between 0 and 100." <<endl;
        return -6;
    }
    if(stopPercent_ < 0 || stopPercent_ > 100){
        cout << "ERROR: Error in Skygen::renderImages. Invalid range for stopPercent_. Please enter an integer between 0 and 100." <<endl;
        return -7;
    }
    if(startPercent_ > stopPercent_){
        cout << "ERROR: Error in Skygen::renderImages. startPecent must be smaller than or equal to stopPercent_." <<endl;
        return -8;
    }
    //The stop/start percentages should correspond to a specific index into the
    //vector that contains all the points we want to render.  Set the corresponding
    //startIndex and stopIndex.
    startIndex = size_t(numAllStrands*startPercent_/100);
    stopIndex = size_t(numAllStrands*stopPercent_/100);
    int errCodeRenderWrapper = renderWrapper();
    if(errCodeRenderWrapper<0){
        cout << "ERROR: Error in skygen::renderImagesByPercent(). RenderWrapper returned errCode = "<<errCodeRenderWrapper<<endl;
        return -837;
    }
    return 0;
}

void skygen::updateSaveImageInfo(){
    //update the filename string used by the post draw callback which saves the images to disk
    //imageName = "LatDeg_" + std::to_string(lat_deg) + "_LonDeg_" + std::to_string(lon_deg) + "_Facing_" + cameraFacing + "_" + std::to_string(frameCount);
    //imageName = std::to_string(frameCount) + "_LatDeg_" + std::to_string(lat_deg) + "_LonDeg_" + std::to_string(lon_deg) + "_Facing_" + cameraFacing;
    if(subLatLonList.size()){
        imageName = cameraFacing + "_S_"+std::to_string(S_subLatLonList[frameCount])+ "_ss_" + std::to_string(ss_subLatLonList[frameCount])+ "_LatDeg_" + std::to_string(lat_deg) + "_LonDeg_" + std::to_string(lon_deg) + "_"+std::to_string(frameCount) ;
    }else{
        imageName = cameraFacing + "_S_"+std::to_string(startIndex+strandIndex)+"_ss_"+std::to_string(locIndex)+ "_LatDeg_" + std::to_string(lat_deg) + "_LonDeg_" + std::to_string(lon_deg) + "_"+std::to_string(frameCount) ;
    }
    pcb->setImageName(imageName);
    ++frameCount;
}

int skygen::render(){
    if(true){
        //Render some frames before we start saving images to disk.
        //This lets us render more frames for the first image in a strand
        //as it takes more time to load its terrain data than subsequent images
        //which are nearby locations with camera in the same direction
        for(uint32_t ii = 0; ii<numPreRenderFrames; ++ii){
            viewer->frame();
        }
        //now render frames that will be saved to disk

        for(uint32_t ii = 0; ii<numRenderFrames; ++ii){
            //updateSaveImageInfo();
            //pcb->setTakePhotoTrue();
            //viewer->frame();

            //Render frames until the database pager is done loading, then save the image
            updateSaveImageInfo();
            pcb->setTakePhotoTrue();
            bool takePhotoStatus = pcb->getTakePhotoStatus();
            while(takePhotoStatus){
                //cout << "getTakePhotoStatus() = " << takePhotoStatus <<endl;
                //cout << "pager files not yet loaded = " << viewer->getDatabasePager()->getFileRequestListSize() <<endl;
                viewer->frame();
                takePhotoStatus = pcb->getTakePhotoStatus();
            }
        }
        //possibly render some additional frames after those we save to disk.
        for(uint32_t ii = 0; ii<numPostRenderFrames; ++ii){
            viewer->frame();
        }
    }else{
        //Render frames until the database pager is done loading, then save the image
        updateSaveImageInfo();
        pcb->setTakePhotoTrue();
        bool takePhotoStatus = pcb->getTakePhotoStatus();
        while(takePhotoStatus){
            cout << "getTakePhotoStatus() = " << takePhotoStatus <<endl;
            cout << "pager files not yet loaded = " << viewer->getDatabasePager()->getFileRequestListSize() <<endl;
            viewer->frame();
            takePhotoStatus = pcb->getTakePhotoStatus();
        }
    }



    return 0;
}


int skygen::renderFirstImage(){
    numPreRenderFrames = 100;
    numRenderFrames = 1;
    numPostRenderFrames = 0;
    updateCameraPose();
    viewer->frame();
    render();
    return 0;
}

int skygen::renderNextImage(){
    numPreRenderFrames = 1;
    numRenderFrames = 1;
    numPostRenderFrames = 0;
    updateCameraPose();
    render();
    return 0;
}

//Calculate and set the view matrix
void skygen::updateCameraPose(){
    lat_rad = osg::DegreesToRadians(lat_deg);
    lon_rad = osg::DegreesToRadians(lon_deg);
    ellipsoid.convertLatLongHeightToXYZ(lat_rad, lon_rad, maxElevationInDatabase, x,y,z);
    heightAboveTerrainAtMaxElevInDatabase = osgSim::HeightAboveTerrain::computeHeightAboveTerrain(terrain.get(), osg::Vec3d(x,y,z));
    heightOfTerrain = maxElevationInDatabase - heightAboveTerrainAtMaxElevInDatabase;
    heightOfCameraEye = heightOfTerrain + heightAboveTerrainDesired;
    //Now set the view matrix for the camera given the lat, lon and alt.
    //Set the camera flying level looking north (or whatever direction cameraFacing).
    ellipsoid.computeLocalToWorldTransformFromLatLongHeight(lat_rad, lon_rad, heightOfCameraEye, vm);
    vm.invert(vm);
    vm *= rotateCameraUp90Deg;
    if(cameraFacing == "west"){
        vm*= rotateCameraLeft90Deg;
    }
    if(cameraFacing == "east"){
        vm*= rotateCameraRight90Deg;
    }
    if(cameraFacing == "south"){
        vm*= rotateCameraLeft180Deg;
    }
    if(cameraFacing == "up"){
        vm*= rotateCameraUp90Deg;
    }
    if(cameraFacing == "down"){
        vm*= rotateCameraDown90Deg;
    }
    viewer->getCamera()->setViewMatrix(vm);
}



int skygen::setImageSaveDir(std::string& imageSaveDir_){
    pcb->PhotosetImageRootPath(imageSaveDir_);
    imageSaveDir = imageSaveDir_;
    return 0;

}

/*
 * * * * * * * * Functions for rerenderSubset * * * * * * * * *
 * The next half dozen (or so) functions are utilzed by rerenderSubset()
 *
 */
//rerenderSubset() is used for rendering a sub-list of lat/lon from the original latLonList.
//Typically, this is because some of the images failed rendering checks, such as sky pixels shining through ground
int skygen::renderSubset(const std::string& directionMode_, const std::string& subLatLonList_){
    subLatLonList = subLatLonList_;
    int errCodeRender = renderImagesStartAtStrand(directionMode_, 0);
    if(errCodeRender<0){
        cout << "ERROR: Error in rerenderSubset(): renderImagesStartAtStrand() returned errCode = "
             << errCodeRender <<endl;
        return -581;
    }
    subLatLonList.clear();
    return 0;
}



int skygen::setActiveStrands_SubLatLonList(){
    std::unordered_set<string> latLonSet;
    int gotLatLonSublist = strandMaker.GetLatLonSet(subLatLonList, latLonSet);
    if(gotLatLonSublist<0){
        cout << "ERROR: setActiveStrandsSubLatLonList(): getLatLonSet() errCode = " << gotLatLonSublist<<endl;
        return-431;
    }
    string latLonKey;
    vector<geoPoint> curr2Add;
    vector<int> curr2Addss;
    if(cameraFacing == "north" || cameraFacing == "south" || cameraFacing == "up" || cameraFacing == "down"){
        for(size_t ii = startIndex; ii<stopIndex; ++ii){
            vector<geoPoint>& currStrand = northSouthStrands[ii];
            const size_t currStrandSize = currStrand.size();
            curr2Add.clear();
            curr2Addss.clear();
            for(size_t jj = 0; jj<currStrandSize; jj++){
                geoPoint& currPoint = currStrand[jj];
                int keyErr = strandMaker.GetLatLon2IDKey(currPoint.lat, currPoint.lon, latLonKey);
                if(keyErr<0){
                    cout << "ERROR: setActiveStrands_SubLatLonList(): getLatLon2IDKey() failed with errCode = "
                         << keyErr <<endl;
                    return -432;
                }

                if(latLonSet.count(latLonKey)){
                    curr2Add.push_back(currPoint);
                    S_subLatLonList.push_back(ii);
                    if(cameraFacing == "south"){
                        //ss_subLatLonList.push_back(currStrandSize-jj-1);
                        curr2Addss.push_back(currStrandSize-jj-1);
                    }else{
                        //ss_subLatLonList.push_back(jj);
                        curr2Addss.push_back(jj);
                    }
                }
            }
            if(curr2Add.size()>0){
                strands.push_back(curr2Add);
                if(cameraFacing == "south"){
                    std::reverse(curr2Addss.begin(),curr2Addss.end());
                }
                const size_t sizeCurr2Addss = curr2Addss.size();
                for(size_t jj = 0; jj<sizeCurr2Addss; jj++){
                    ss_subLatLonList.push_back(curr2Addss[jj]);
                }
            }
        }
        if(cameraFacing == "south"){
            //reverse the order of the strands so they're naturally rendered going from north to south
            size_t numStrands = strands.size();
            for(size_t ii = 0; ii < numStrands; ++ii){
                vector<geoPoint>& currStrand = strands[ii];
                std::reverse(currStrand.begin(),currStrand.end());
            }
        }
    }
    if(cameraFacing == "east" || cameraFacing == "west"){
        for(size_t ii = startIndex; ii<stopIndex; ++ii){

            vector<geoPoint>& currStrand = eastWestStrands[ii];
            const size_t currStrandSize = currStrand.size();
            curr2Add.clear();
            curr2Addss.clear();
            for(size_t jj = 0; jj<currStrandSize; jj++){
                geoPoint& currPoint = currStrand[jj];
                int keyErr = strandMaker.GetLatLon2IDKey(currPoint.lat, currPoint.lon, latLonKey);
                if(keyErr<0){
                    cout << "ERROR: setActiveStrands_SubLatLonList(): getLatLon2IDKey() failed with errCode = "
                         << keyErr <<endl;
                    return -432;
                }
                if(latLonSet.count(latLonKey)){
                    curr2Add.push_back(currPoint);
                    S_subLatLonList.push_back(ii);
                    //cout << "(ii,jj,currStrandSize,currStrandSize-jj-1) = (" << ii <<"," <<jj<<","<<currStrandSize<<"," <<currStrandSize-jj-1<<")"<<endl;
                    if(cameraFacing == "west"){
                        //ss_subLatLonList.push_back(currStrandSize-jj-1);

                        curr2Addss.push_back(currStrandSize-jj-1);
                    }else{
                        //ss_subLatLonList.push_back(jj);
                        curr2Addss.push_back(jj);
                    }
                }
            }
            if(curr2Add.size()>0){
                strands.push_back(curr2Add);
                if(cameraFacing == "west"){
                    std::reverse(curr2Addss.begin(),curr2Addss.end());
                }
                const size_t sizeCurr2Addss = curr2Addss.size();
                for(size_t jj = 0; jj<sizeCurr2Addss; jj++){
                    ss_subLatLonList.push_back(curr2Addss[jj]);
                }
            }
        }
        if(cameraFacing == "west"){
            //reverse the order of the strands to they're naturally rendered going from east to west
            size_t numStrands = strands.size();
            for(size_t ii = 0; ii < numStrands; ++ii){
                vector<geoPoint>& currStrand = strands[ii];
                std::reverse(currStrand.begin(),currStrand.end());
            }
        }
    }

    numAllStrands = strands.size();
    stopIndex = numAllStrands;
    return 0;
}




int skygen::init(){
    if(eastWestStrands.size()<1 || northSouthStrands.size()<1){
        int errCode253 = strandMaker.GetStrands(latLonListFullPath, eastWestStrands, northSouthStrands);
        if(errCode253 < 0){
            cout << "ERROR: Error in SKYGEN::Init().  strandMaker.GetStrands() returned errCode = " << errCode253<<endl;
            return -7373;
        }
    }
    if(!traits.valid()){
        traits = new osg::GraphicsContext::Traits;
    }
    const string wStr = " traits.valid() = 0";
    if(!traits.valid()){
        // print
        osg::notify(osg::WARN) << wStr << endl;
        // error
        return -1;
    }

    // set traits properties
    traits->screenNum = 0;
    traits->x = 400;
    traits->y = 40;
    traits->width = width;
    traits->height =height;
    traits->doubleBuffer = true;
    traits->windowDecoration = false;
    traits->vsync = true;
    traits->samples = numResamples;
    gc = osg::GraphicsContext::createGraphicsContext(traits.get());
    terrain =  osgDB::readNodeFile(terrainFile);
    //osg::Texture::getTextureObjectManager(0)->setMaxTexturePoolSize( 64000 );
    if(!viewer.valid()){
        viewer = new osgViewer::Viewer();
    }
    viewer->getDatabasePager()->setDoPreCompile(true);
    viewer->getDatabasePager()->setTargetMaximumNumberOfPageLOD( 10 );
    camera = viewer->getCamera();
    camera->setGraphicsContext(gc);
    camera->setProjectionMatrixAsPerspective(fovy,aspectRatio,zNear,zFar);//90.0,1.0,0.5,200000.0);
    camera->setLODScale(LODScale);//0.1
    camera->setSmallFeatureCullingPixelSize(SmallFeatureCullingPixelSize);//0.9
    camera->setCullingMode(osg::CullSettings::VIEW_FRUSTUM_SIDES_CULLING | osg::CullSettings::SMALL_FEATURE_CULLING);
    camera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
    /*
    bool disableAutoFrustum = true;
    if(disableAutoFrustum){
        traits->screenNum = 0;
        traits->x = 400;
        traits->y = 400;
        traits->width = width;
        traits->height =height;
        traits->doubleBuffer = true;
        traits->windowDecoration = true;
        traits->vsync = true;
        traits->samples = numResamples;
        camera->setLODScale(LODScale);
        camera->setSmallFeatureCullingPixelSize(SmallFeatureCullingPixelSize);
        camera->setCullingMode(osg::CullSettings::VIEW_FRUSTUM_SIDES_CULLING | osg::CullSettings::SMALL_FEATURE_CULLING);
        camera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
    }*/
    //set the viewport
    camera->setViewport(0, 0, width, height);
    camera->setClearColor(backgroundColor);//white
    camera->setCullingActive(true);


    // Adding the elements
    if(!root.valid()){
        root = new osg::Group;
        root->addChild( terrain.get() );
    }

    // Setting up the view
    viewer->setSceneData( root.get() );
    if(!stats.valid()){
        stats = new osgViewer::StatsHandler();
        viewer->addEventHandler(stats.get());
    }

    if(!image.valid()){
        image = new osg::Image;
    }

    if(!pcb.valid()){
        pcb = new skygen::PhotoCallback(image.get(), imageRootPath, imageExt, false, viewer);
    }
    camera->setPostDrawCallback( pcb.get() );
    //init commonly used rotation matrices
    rotateCameraUp90Deg.makeRotate(-osg::PI_2, osg::Vec3f(1.0, 0.0, 0.0));
    rotateCameraDown90Deg.makeRotate(osg::PI_2, osg::Vec3f(1.0, 0.0, 0.0));
    rotateCameraLeft90Deg.makeRotate(-osg::PI_2, osg::Vec3f(0.0, 1.0, 0.0));
    rotateCameraRight90Deg.makeRotate(osg::PI_2, osg::Vec3f(0.0, 1.0, 0.0));
    rotateCameraLeft180Deg.makeRotate(-osg::PI, osg::Vec3f(0.0, 1.0, 0.0));

    frameCount = 0;
    strandIndex = 0;
    locIndex = 0;
    return 0;
}

