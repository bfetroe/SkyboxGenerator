#ifndef SKYGEN_H
#define SKYGEN_H

#include <osg/Config>
#include <osg/Referenced>
#include <osg/ref_ptr>
#include <osg/GraphicsContext>
#include <osg/GLObjects>
#include <osg/Image>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/ReadFile>
#include <osg/PositionAttitudeTransform>
#include <osgSim/HeightAboveTerrain>
#include <osgDB/WriteFile>
#include <iostream>
#include <string>
#include <ostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "StrandMaker.h"
#include <typeinfo>
#include <unordered_set>
#include <experimental/filesystem>  //fs
#include "qPressedNonBlocking.h"

#include <osg/Texture>
#include <osg/Node>
#include <osgDB/DatabasePager>
#include <osgDB/ReadFile>
#include <osgUtil/PrintVisitor>


class skygen
{
public:
    skygen(const std::string& terrainFile_, const std::string& latLonListFullPath_,
           const std::string& imageRootPath_, const std::string& imageExt_,
           const std::string &imageFolderPrefix_,
           double& width_, double& height_, double& numResamples_, osg::Vec4f& backgroundColor_,
           double& fovy_, double& aspectRatio_, double& zNear_, double& zFar_,
           double& heightAboveTerrainDesired_, double& maxElevationInDatabase_,
           double& LODScale_, double& SmallFeatureCullingPixelSize_);
    int init();
    int setCameraFacing(const std::string& cameraFacing_);
    //There are several ways to start the rendering.
    // 1) renderImagesByStrand (where you defind the start and stop strands)
    // 2) renderImagesStartAtStrand (where you define the start strand and it goes until the end)
    // 3) rednerImagesByPercent (where you give the approximate start and stop as percentages)
    //Every one of the aforementioned methods needs to know the directionMode.
    //std::string directionMode_ defines the camera heading. (i.e. "north" "south" "east" "west" or "up")
    //const int startPercent_ and stopPercent_ are used to have a bit more
    //fine control over which images are rendered... input should be an int between 0 and 100.
    //note: given the interger arithmetic involved, the startPercent_ and
    //stopPercent are approximate and can map to a strand that is actually
    //one less than expected.
    //(Ex: setting both start and stop to 100 would render only the last strand,
    // which would be filed in folder .../color99
    int renderImagesByPercent(const std::string& directionMode_, const int startPercent_, const int stopPercent_);
    int renderImagesByStrand(const std::string& directionMode_,
                             const size_t &startStrand_,
                             const size_t &stopStrand);
    int renderImagesStartAtStrand(const std::string& directionMode_, const size_t& startStrand_);
    int renderSubset(const std::string& directionMode_, const std::string& subLatLonList_);

    //NOTE: the AutoAll option is not yet implemented, and currently returns -1;
    /* renderImagesAutoAll() looks in the imageRootPath and finds the
     * last strand that was fully rendered (say strand 1123) and then starts
     * rendering on the next strand (say, 1124).  This can result in rendering
     * frames that have been rendered previously (if, say, strand 1124 contains
     * 1001 locations, but only 999 were rendered previoulsy, it would re-render
     * those 999 frames and then the last 2, then move onto strand 1125.
     * NOTE: As this method only looks at strands from last to first, if there
     * is a fully rendered strand (say 445 and nothing higher) but strands 10,
     * 11,12 are missing, it would start on strand 446.
     */
    int renderImagesAutoAll(const std::string& directionMode_);

    class PhotoCallback : public osg::Camera::DrawCallback
    {
    public:

        skygen::PhotoCallback( osg::Image* img, std::string& imageRootPath, std::string& imageExt, bool takePhoto,
                               osg::ref_ptr<osgViewer::Viewer>& viewer)
        : _image(img), _fileIndex(0), _imageRootPath(imageRootPath), _imageExt(imageExt) , _takePhoto(takePhoto),
        _initialized(false), _viewer(viewer){}
        skygen::~PhotoCallback(){

        }

        void PhotosetImageRootPath(std::string& imageRootPath){
            _imageRootPath = imageRootPath;
        }
        void setImageName(std::string& imageName){
            _imageName = imageName;
        }
        void setImageExt(std::string& imageExt){
            _imageExt = imageExt;
        }
        void setWidth(double& width){
            _width = width;
        }
        void setHeight(double& height){
            _height = height;
        }
        void setTakePhotoTrue(){
            _takePhoto = true;
        }
        void setTakePhotoFalse(){
            _takePhoto = false;
        }
        void setSinglePhotoModeTrue(){
            _singlePhotoMode = true;
        }
        void setSinglePhotoModeFalse(){
            _singlePhotoMode = false;
        }
        void takePhoto(){
            setSinglePhotoModeTrue();
            setTakePhotoTrue();
        }

        bool getTakePhotoStatus(){
            return _takePhoto;
        }


        virtual void operator()( osg::RenderInfo& renderInfo ) const{
            //_takePhoto = false;
            if(_takePhoto == true){
                if(!_initialized){
                    osg::ref_ptr<osg::GraphicsContext> _gc = renderInfo.getState()->getGraphicsContext();
                    if (_gc->getTraits() ){
                        _width = _gc->getTraits()->width;
                        _height = _gc->getTraits()->height;
                        _pixelFormat = (_gc->getTraits()->alpha ? GL_RGBA : GL_RGB);
                        _initialized = true;
                    }else{
                        std::cerr << "ERROR: ERROR IN PHOTOCALLBACK:INIT(): Failed to _gc->getTraits()" <<std::endl;
                    }
                }else{
                    //check if database pager is finished loading
                    unsigned int pagerSize =_viewer->getDatabasePager()->getFileRequestListSize();
                    if(pagerSize == 0){
                        //time to save the rendered image
                        _image->readPixels( 0, 0, _width, _height, _pixelFormat, GL_UNSIGNED_BYTE );
                        std::stringstream _ss;
                        _ss << _imageRootPath << "/" << _imageName <<"_"<<_fileIndex<< "." << _imageExt;
                        _fileIndex++;
                        if ( !osgDB::writeImageFile(*_image, _ss.str()) ){
                            std::cerr << "ERROR: Error in skygen::PhotoCallBack: Failed to writeImageToFile()" <<std::endl;
                        }else{
                            //image saved successfully.
                            _takePhoto=false;
                        }
                    }else{
                        //Not finished loading, so render another frame with the same pose.
                        //We don't currently need to do anything in this else statement, as simply
                        //maintaining the state of _takePhoto=true will result in rendering
                        //another frame at the same location (given this is still in the skygen class)
                    }
                }
            }

        }


    protected:
        osg::ref_ptr<osg::Image> _image;
        //osg::ref_ptr<osg::Image> _image_gray8;
        mutable uint32_t _fileIndex;

    private:
        //photocallback
        mutable bool _initialized;
        mutable GLenum _pixelFormat;
        mutable std::string _imageRootPath;
        std::string _imageName;
        std::string _imageExt;
        mutable double _width;
        mutable double _height;
        mutable bool _takePhoto;
        bool _singlePhotoMode;
        osg::ref_ptr<osgViewer::Viewer> _viewer;

    };
private:
    //StrandMaker stuff
    StrandMaker strandMaker;
    const std::string latLonListFullPath;
    std::vector<std::vector<geoPoint> > eastWestStrands;
    std::vector<std::vector<geoPoint> > northSouthStrands;
    //skygen functions
    std::vector<std::vector<geoPoint> > strands;
    size_t currStrand;
    size_t currPoint;
    void updateCameraPose();
    void updateSaveImageInfo();
    int preRenderSetup(const std::string &directionMode_);
    int renderWrapper();
    int renderFirstImage();
    int renderNextImage();
    int render();
    int convertStartStopPercentToIndex(int startPercent_,int stopPercent_);
    int setActiveStrands();
    int setImageSaveDir(std::string& imageSaveDir_);
    size_t getStrandPercent(const size_t& currStrandIndex);
    int createDir(std::string path);
    int createDirRecur(std::string path);
    int setActiveStrands_SubLatLonList();
    //skygen variables, strings, pointers, etc
    const std::string terrainFile;
    std::string imageRootPath; //= "C:/Nav/swissDEM_30m/skylines/grid_1km";
    std::string imageSaveDir; //internally generated from imageRootPath + currStrandRenderPercent
    std::string imageExt; //= "png";
    std::string imageFolderPrefix; //= "color"
    osg::ref_ptr<osg::GraphicsContext::Traits> traits;
    osg::ref_ptr<osg::GraphicsContext> gc;
    osg::ref_ptr<osg::Group> root;
    osg::ref_ptr<osg::Node> terrain;
    osg::ref_ptr<osgViewer::Viewer> viewer;
    osg::ref_ptr<osg::Camera> camera;
    osg::EllipsoidModel ellipsoid;
    double fovy;
    double aspectRatio;
    double zNear;
    double zFar;
    double width;
    double height;
    double numResamples;
    osg::Vec4f backgroundColor;
    double SmallFeatureCullingPixelSize;
    double LODScale;
    osg::ref_ptr<skygen::PhotoCallback> pcb;
    osg::ref_ptr<osg::Image> image;
    osg::ref_ptr<osgViewer::StatsHandler> stats;
    uint32_t frameCount;
    osg::Matrixd vm;
    osg::Matrixd rotateCameraLeft90Deg;
    osg::Matrixd rotateCameraRight90Deg;
    osg::Matrixd rotateCameraLeft180Deg;
    osg::Matrixd rotateCameraUp90Deg;
    osg::Matrixd rotateCameraDown90Deg;
    double lat_deg;// = 46.8183;
    double lon_deg;// = 8.2275;
    double lat_rad;// = osg::DegreesToRadians(lat_deg);
    double lon_rad;// = osg::DegreesToRadians(lon_deg);
    double heightAboveTerrainDesired;// = 1.8; //meters (if using terrain built in vpb from heightmaps in meters)
    double maxElevationInDatabase;// = 10000.0; // meters
    double x,y,z;
    double heightAboveTerrainAtMaxElevInDatabase;
    double heightOfTerrain;
    double heightOfCameraEye;
    std::string cameraFacing;
    std::string imageName;
    uint32_t numPreRenderFrames;
    uint32_t numRenderFrames;
    uint32_t numPostRenderFrames;
    size_t startIndex;
    size_t stopIndex;
    int currStrandRenderPercent;
    size_t numAllStrands;
    size_t strandIndex;
    size_t locIndex;
    std::string subLatLonList;
    std::vector<size_t> S_subLatLonList;
    std::vector<size_t> ss_subLatLonList;



};





#endif // SKYGEN_H
