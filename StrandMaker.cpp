#include "StrandMaker.h"

#include <iostream>     //cout
#include <vector>
#include <fstream>      //ifstream and ofstream
#include <sstream>
#include <algorithm>    //sort
#include <string>

using std::string;
using std::endl;
using std::cout;
using std::vector;
using std::getline;
//using strandify::strandedGeoPoint;

const double PI = 3.141592653589793;
const double R_SAFETY_FACTOR = 1.5;


template <class point>
void StrandMaker::PrintVector(std::vector<point>& vec){
    size_t vecSize = vec.size();
    if(vecSize>0){
        std::cout << "("<<vec[0].lat<<", "<<vec[0].lon<<")"<<endl;
        point thing;
        for(auto vecItr = vec.begin()+1; vecItr != vec.end(); vecItr++){
            thing = *vecItr;
            std::cout << ", ("<<thing.lat<<", "<<thing.lon<<")";
        }
        std::cout << endl;
    }else{
        std::cout << "Vector is empty." <<endl;
    }
}

template <class point>
int StrandMaker::SaveStrandsToFile(std::vector<std::vector<point> >& strands, const string& fullPathAndFilename){
    std::ofstream myfile(fullPathAndFilename);
    if (myfile.is_open()){
        size_t currStrandIndex = 0;
        for(auto strandItr = strands.begin(); strandItr != strands.end(); strandItr++){
            vector<point>& strand = *strandItr;
            size_t strandSize= strand.size();
            if(strandSize==0){
                myfile << "NaN\n";
                continue;
            }
            myfile << currStrandIndex <<", "<< strand[0].lat << ", " << strand[0].lon<< "\n";
            if(strandSize>1){
                for(size_t index = 1; index < strandSize; index++){
                    myfile << currStrandIndex << ", " << strand[index].lat << ", " << strand[index].lon<< "\n";
                }
            }
            ++currStrandIndex;
        }
      myfile.close();
    }else{
        cout << "Unable to open file";
        return -1;
    }
    return 0;
}


int StrandMaker::GetLatLonList(const std::string& filenameFullPath, std::vector<geoPoint>& latLonList, const bool appendMode){
    if(appendMode==false){
        latLonList.clear();
    }
    std::ifstream inFile;
    string headerLine;
    inFile.open(filenameFullPath);
    if(!inFile.is_open()){
        cout << "ERROR: Failed to open latLonList. Make sure the filename is valid." <<endl;
        return -74;
    }
    std::getline(inFile, headerLine);

    if(headerLine != expectedHeader){
        //This is not a valid input file
        //...or it is formatted differently than expected
        latLonList.clear();
        cout << "ERROR IN GetLatLonList(): first line was not of the "
             << "expected form. Should be = \"FID,Id,POINT_X,POINT_Y\""
             << " but instead was = \"" << headerLine << "\"" << endl;
        return -1;
    }
    string token;
    geoPoint latLonGeoPoint;

    int column = 0;
    cout.precision(std::numeric_limits<double>::max_digits10);
    while(getline(inFile,token,',')) {
       // cout << "Count = " << count <<endl;
        switch(column){
        case 0:
            //cout << "Token 0 = " << token <<endl;
            break;
        case 1:
            //cout << "Token 1 = " << token <<endl;
            break;
        case 2:{
            //cout << "Token 2 = " << token <<endl;
            std::istringstream converter(token);
            converter >> latLonGeoPoint.lon;
            break;}
        case 3:{
            //cout << "Token 3 = " << token <<endl;
            std::istringstream converter(token);
            converter >> latLonGeoPoint.lat;
            latLonList.push_back(latLonGeoPoint);
            column = 0;
            break;}
        }
        column += 1;
    }
    inFile.close();
    return 0;
}

int StrandMaker::GetGridSpacing(std::vector<geoPoint>& fullLatLonList, double& gridSpacing){
    if(fullLatLonList.size() <= 2){
        cout<<"ERROR: Error in GetGridSpacing(). There are less than 2 geopoints in fullLatLonList"<<endl;
        return -1;
    }
    //look at the first 100 geopoints of fullLatLonList to see what the common distance is between geopoints
    size_t testLen;
    if (fullLatLonList.size() < 100){
        testLen = fullLatLonList.size();
    }else{
        testLen = 100;
    }
    double rMin = 1000.0;
    double r = rMin;
    double x1, y1, x2, y2;
    for(size_t i=0; i<testLen; i++){
        //calc euclidan distance (r) between two geopoints and keep track of smallst r thus far
        x1 = fullLatLonList[i].lon;
        x2 = fullLatLonList[i+1].lon;
        y1 = fullLatLonList[i].lat;
        y2 = fullLatLonList[i+1].lat;
        r = sqrt( pow((x2-x1),2.0) + pow((y2-y1),2.0) );
        if(r<rMin){
            rMin = r;
        }
    }
    if(rMin == 1000.0){
        cout << "ERROR: GetEastWestStrands: Unable to calculate realistic distance metric between geopoints."<<endl;
        gridSpacing = 0.0;
        return -2;
    }
    gridSpacing = rMin;
    return 0;
}

/*
int StrandMaker::GetEastWestStrands(std::vector<geoPoint>& fullLatLonList,
                       std::vector<std::vector<geoPoint>>& eastWestStrands){
    if(fullLatLonList.size() <= 2){
        cout<<"ERROR: Error in GetEastWestStrands: There are less than 2 geopoints in fullLatLonList"<<endl;
        return -1;
    }
    double rMin;
    int errCode52 = GetGridSpacing(fullLatLonList, rMin);
    if(errCode52<0){
        cout << "ERROR: Error in GetEastWestStrands::GetGridSpacing. Returned error code = " << errCode52<<endl;
    }
    //rMin is now valid.
    //Separate points into "strands" based on if they're more than rSafetyFactor away from the previous point.
    const double rSafetyFactor = R_SAFETY_FACTOR; //typically R_SAFETY_FACTOR = 1.5
    const double rFactor = rMin * rSafetyFactor;
    const auto lenFullLatLonList = fullLatLonList.size();
    double x1 = fullLatLonList[0].lon;
    double y1 = fullLatLonList[0].lat;
    double x2, y2, rSquared;
    geoPoint point;
    point.lon = x1;
    point.lat = y1;
    size_t currStrand = 0;
    eastWestStrands.clear();
    vector<geoPoint> emptyGeoVector;
    vector<double> angles;
    eastWestStrands.push_back(emptyGeoVector);
    eastWestStrands[currStrand].push_back(point);
    for(size_t i =1; i<=lenFullLatLonList; i++){
        x2 = fullLatLonList[i].lon;
        y2 = fullLatLonList[i].lat;
        rSquared = pow((x2-x1),2.0) + pow((y2-y1),2.0);
        point.lon = x2;
        point.lat = y2;
        if(rSquared<rFactor){
            //this is a continuation of the same strand
            eastWestStrands[currStrand].push_back(point);
            ///////
            double dot12 = x2-x1;
            double mag1 = 1;
            double mag2 = sqrt(pow(x2-x1,2.0)+pow(y2-y1,2.0));
            double angle = dot12 / (mag1*mag2);
            angles.push_back(angle);
            ///////
        }else{
            //need to form a new strand
            eastWestStrands.push_back(emptyGeoVector);
            ++currStrand;
            eastWestStrands[currStrand].push_back(point);
        }

        //save x2 as x1 and repeat
        x1 = x2;
        y1 = y2;


    }
    //print the first few entries of eastWestStrands to ensure this function is working correctly
    const bool printDiag1 = false;
    if(printDiag1 == true){
        size_t angleCount = 0;
        double minAngle = 1.0;
        double angle = 0.0;
        for(size_t ii = 0; ii<eastWestStrands.size(); ii++){
            for(size_t jj = 0; jj<eastWestStrands[ii].size(); jj++){
                point = eastWestStrands[ii][jj];
                if(jj>1){
                    angle = angles[angleCount];
                    if(angle < minAngle){
                        minAngle = angle;
                    }
                    ++angleCount;
                }
            }
        }
        cout << "minAngle = " << minAngle << endl;
    }
    return 0;
}

*/

int StrandMaker::GetBoundingBox(std::vector<geoPoint>& latLonList,
                   double& latMin,
                   double& latMax,
                   double& lonMin,
                   double& lonMax){
    size_t listSize = latLonList.size();
    geoPoint point;
    double lat, lon;
    if(listSize==1){
        point = latLonList[0];
        lat = point.lat;
        lon = point.lon;
        latMin = lat;
        latMax = lat;
        lonMin = lon;
        lonMax = lon;
        return 1;
    }
    if(listSize==0){
        //latLonList is empty
        latMin = 0.0;
        latMax = 0.0;
        lonMin = 0.0;
        lonMax = 0.0;
        cout << "ERROR: Error in GetBoundingBox(). LatLonList is empty." <<endl;
        return -1;
    }
    //listSize is >= 2
    point = latLonList[0];
    lat = point.lat;
    lon = point.lon;
    latMin = lat;
    latMax = lat;
    lonMin = lon;
    lonMax = lon;
    for(size_t ii = 1; ii<latLonList.size(); ii++){
        point = latLonList[ii];
        lat = point.lat;
        lon = point.lon;
        if(lat > latMax){
            latMax = lat;
        }
        if(lat < latMin){
            latMin = lat;
        }
        if(lon > lonMax){
            lonMax = lon;
        }
        if(lon < lonMin){
            lonMin = lon;
        }
    }
    //we've gone through the whole latLonList and the max and min lat and lon.
    return 0;
}

int StrandMaker::CreateStrandedGrid(std::vector<geoPoint>& latLonList,
                       std::vector<std::vector<std::vector<strandedGeoPoint> > >& strandedGrid,
                       std::vector< std::vector<strandedGeoPoint> >& strandedList,
                       std::vector<size_t>& strandByID){
    //Note to developers: interaction with the strandedGrid is of the form: strandedGrid[lonCell][latCell]
    strandedGrid.clear();
    strandedList.clear();
    //for each geoPoint in latLonList add it to strandedGrid
    //but first a few checks...
    //first we must determine the appropriate spacing between points
    double rMin = -1.0;
    int errCode43 = GetGridSpacing(latLonList, rMin);
    if(errCode43<0){
        cout << "ERROR: Error in CreateStrandedGrid::GetGridSpacing(). Returned error code = " << errCode43<<endl;
        return -19886;
    }
    if(rMin<0){
        cout << "ERROR: Error in CreateStrandedGrid::GetGridSpacing(). rMin is negative: rMin = "<< rMin <<endl;
        return -12654;
    }
    //rMin is now set to whatever grid spacing we found
    //now let's find the bounding box that contains all the geopoints in latLonList
    double latMin, latMax, lonMin, lonMax;
    int errCode93 = GetBoundingBox(latLonList, latMin, latMax, lonMin, lonMax);
    if(errCode93<0){
        cout << "ERROR: Error in CreateStrandedGrid::GetBoundingBox(). Returned error code = " << errCode93<<endl;
        return -23432;
    }
    //the desired cell size is 1.5rMin. Choose bigger lenCell if memory is a problem.
    //Smaller than roughly 1.5rMin is not advized as you have no guarantee that
    //there will be geoPoints in the self or neighboring cells.
    //So, given the min/max lat and lon, create the grid of cells with lenCell spacing.
    const double lenCell = R_SAFETY_FACTOR * rMin; //typically R_SAFETY_FACTOR=1.5
    double boundingBoxSizeDegreesLat = latMax-latMin;
    double boundingBoxSizeDegreesLon = lonMax-lonMin;
    size_t gridSizeLat = 2 + size_t(ceil(boundingBoxSizeDegreesLat/lenCell));
    size_t gridSizeLon = 2 + size_t(ceil(boundingBoxSizeDegreesLon/lenCell));
    vector<vector<strandedGeoPoint> > gridCols(gridSizeLat);
    //fill the empty strandedGrid with empty cells
    for(size_t ii = 0; ii < gridSizeLon; ii++){
        strandedGrid.push_back(gridCols);
    }
    //fill the grid cells from latLonList
    size_t latCell, lonCell;
    double lat, lon;
    size_t numStrands = 0;
    strandedGeoPoint sgPoint;
    std::vector<strandedGeoPoint> sgStrand;
    for(auto point : latLonList){
        //add point to appropriate cell
        lat = point.lat;
        lon = point.lon;
        if(lat < latMin || lon < lonMin || lat > latMax || lon > lonMax){
            cout << "ERROR: Error in CreateStrandedGrid: Found geopoint outside of bounding box." <<endl;
            cout << "(latMin, lat, latMax) = (" << latMin << ", " << lat << ", " << latMax << ")" << endl;
            cout << "(lonMin, lon, lonMax) = (" << lonMin << ", " << lon << ", " << lonMax << ")" << endl;
            return -3;
        }
        latCell = 1 + size_t(floor( (lat-latMin) / lenCell ));
        lonCell = 1 + size_t(floor( (lon-lonMin) / lenCell ));
        sgPoint.lat = lat;
        sgPoint.lon = lon;
        sgPoint.id = numStrands;
        strandedGrid[lonCell][latCell].push_back(sgPoint);
        sgStrand.clear();
        sgStrand.push_back(sgPoint);
        strandedList.push_back(sgStrand);
        strandByID.push_back(numStrands);
        ++numStrands;
    }
    size_t strandedListSize = strandedList.size();
    if(strandedListSize != numStrands){
        cout << "ERROR: Error in CreatedStrandedGrid: After instantiation, strandedList.size() != numStrands." <<endl;
        cout << "strandedList.size() = " << strandedListSize <<endl;
        cout << "numStrands = " << numStrands <<endl;
        return -9823;
    }
    //finished instantiation of grid of strandedGeoPoints
    return 0;
}


template <class point>
const double StrandMaker::getDistanceBetweenPoints(point& p1, point& p2){
    const double x1 = p1.lon;
    const double x2 = p2.lon;
    const double y1 = p1.lat;
    const double y2 = p2.lat;
    return sqrt( pow((x2-x1),2.0) + pow((y2-y1),2.0) );
}

//returns angle in radians that the line between p1 and p2 makes with east
template <class point>
inline const double StrandMaker::getAngleRadians(point& p1, point& p2){
    return atan2(p2.lat - p1.lat, p2.lon - p1.lon);
}

//returns angle in radians that the line between p1 and p2 makes with east
template <class point>
const double StrandMaker::getAngleDegrees(point& p1, point& p2){
    const double angleRad = getAngleRadians(p1,p2);
    const double angleDeg = angleRad * 180.0 /PI;
    return angleDeg;
}

//IsAngleInside(targetAngle, angle1, angle2)
//Returns true is targetAngle is between the acute angle between angle1 and 2.
//Note: doens't work if angle1 and angle2 are exactly 180degrees apart
bool StrandMaker::IsAngleInside_0to360(double target, double angle1, double angle2)
{
    // make the angle from angle1 to angle2 to be <= 180 degrees
    double diffAngle21 = angle2 - angle1;
    double rAngle = fmod((fmod(diffAngle21,360.0) + 360.0), 360.0);
    if (rAngle >= 180){
        std::swap(angle1, angle2);
    }
    // check if it passes through zero
    if (angle1 <= angle2){
        return target >= angle1 && target <= angle2;
    }else{
        return target >= angle1 || target <= angle2;
    }
}

bool StrandMaker::IsAngleInside_Neg180to180(double target, double angle1, double angle2){
    double t = target;
    double a1 = angle1;
    double a2 = angle2;
    //convert all angles to within [0.0,360.0)
    while(t<0.0){t += 360.0;}
    while(a1<0.0){a1 += 360.0;}
    while(a2<0.0){a2 += 360.0;}
    return IsAngleInside_0to360(t,a1,a2);
}

//first, we add the neighbor (and other points in its strand) to the selfStrand
//then we update the neighbor and other points in its strand to new strand numbers
//then we update the strandByID for the moved points
//finally we remove the neighbor point (and others) from the old strand
int StrandMaker::JoinStrands(std::vector<strandedGeoPoint>& selfStrand,
                std::vector<strandedGeoPoint>& neighborStrand,
                std::vector<size_t>& strandByID){
    //go through all the points in the neighborStrand and
    //add them to the selfStrand
    //also update those moved points strandByID
    if(selfStrand.empty()){
        cout << "ERROR: JoinStrands(): selfStrand is empty." << endl;
        return -5253;
    }
    if(neighborStrand.empty()){
        cout << "ERROR: JoinStrands(): neighborStrand is empty." << endl;
        return -5344;
    }
    if(strandByID.empty()){
        cout << "ERROR: JoinStrands(): strandByID is empty." << endl;
        return -5277;
    }
    strandedGeoPoint buddyPoint;
    const size_t selfStrandID = selfStrand[0].id;
    for(auto npItr = neighborStrand.begin(); npItr != neighborStrand.end(); npItr++){
        buddyPoint = *npItr;
        selfStrand.push_back(buddyPoint);
        strandByID[buddyPoint.id] = selfStrandID;
    }
    //all done... so we're done with the neighborStrand... forever, muahahaHAHAA
    neighborStrand.clear();
    return 0;
}


bool CompareStrandedGeoPointByLat(strandedGeoPoint const & p1, strandedGeoPoint const & p2) {
    // return "true" if "p1" is ordered before "p2", for example:
    return p1.lat < p2.lat;
}

bool CompareStrandedGeoPointByLon(strandedGeoPoint const & p1, strandedGeoPoint const & p2) {
    // return "true" if "p1" is ordered before "p2", for example:
    return p1.lon < p2.lon;
}

//all points have been strandified and reside in strandedList[strandNumber][pointInStrand]
//the strands need to be sorted and placed into finalStrands
int StrandMaker::SortAndFinalizeStrands(std::vector<std::vector<strandedGeoPoint> >& jumbledStrands,
                           std::vector<std::vector<geoPoint> >& finalStrands,
                           const std::string& directionMode){
    finalStrands.clear();
    size_t numJumbledStrands = jumbledStrands.size();
    if(numJumbledStrands <=2){
        cout << "ERROR: Error in SortAndFinalizeStrands: jumbledPoints.size() <= 2. Currently size = " << numJumbledStrands <<endl;
        jumbledStrands.clear();
        return -5343;
    }
    std::vector<geoPoint> geoStrand;
    geoPoint point;
    if(directionMode == "north" || directionMode =="south"){
        for(size_t jStrandIndex = 0; jStrandIndex < numJumbledStrands; ++jStrandIndex){
            //for each jumbled strand, sort and add to finalStrands
            std::vector<strandedGeoPoint>& jStrand = jumbledStrands[jStrandIndex];
            if(jStrand.size()>0){
                std::sort(jStrand.begin(), jStrand.end(), CompareStrandedGeoPointByLat);
                geoStrand.clear();
                for(auto sgPoint : jStrand){
                    point.lat = sgPoint.lat;
                    point.lon = sgPoint.lon;
                    geoStrand.push_back(point);
                }
                finalStrands.push_back(geoStrand);
            }
        }
    }else if(directionMode == "east" || directionMode =="west"){
        for(size_t jStrandIndex = 0; jStrandIndex < numJumbledStrands; ++jStrandIndex){
            //for each jumbled strand, sort and add to finalStrands
            std::vector<strandedGeoPoint>& jStrand = jumbledStrands[jStrandIndex];
            if(jStrand.size()>0){
                //NOTE: CompareBy...LON vs LAT
                std::sort(jStrand.begin(), jStrand.end(), CompareStrandedGeoPointByLon);
                geoStrand.clear();
                for(auto sgPoint : jStrand){
                    point.lat = sgPoint.lat;
                    point.lon = sgPoint.lon;
                    geoStrand.push_back(point);
                }
                finalStrands.push_back(geoStrand);
            }
        }
    }else{
        cout << "ERROR: Error in StrandMaker::SortAndFinalizeStrands. Invalid directionMode. Please enter north, south, east or west."<<endl;
        return -77;
    }
    return 0;
}


int StrandMaker::Strandify(std::vector<std::vector<std::vector<strandedGeoPoint> > >& strandedGrid,
                   std::vector<std::vector<geoPoint> >& outputStrands,
                   std::vector< std::vector<strandedGeoPoint> >& strandedList,
                   std::vector<geoPoint>& latLonList,
                   std::vector<size_t>& strandByID,
                   const std::string& directionMode){

    //measured clockwise from east: [0.0,360.0) example: 90.0 is north (default setting)
    double ANGLE_DESIRED_DEGREES_temp = 90.0;
    if(directionMode == "east" || directionMode == "west"){
        ANGLE_DESIRED_DEGREES_temp = 0.0;
    }
    const double ANGLE_DESIRED_DEGREES = ANGLE_DESIRED_DEGREES_temp;
    const double ANGLE_ACCEPTABLE_WINDOW_DEGREES = 15.0;
    const double ANGLE_DESIRED_RADIANS = ANGLE_DESIRED_DEGREES * PI/180.0;
    const double ANGLE_ACCEPTABLE_WINDOW_RADIANS = ANGLE_ACCEPTABLE_WINDOW_DEGREES * PI/180.0;
    const double ANGLE_PLUS_WINDOW_RADIANS = ANGLE_DESIRED_RADIANS + ANGLE_ACCEPTABLE_WINDOW_RADIANS;
    const double ANGLE_MINUS_WINDOW_RADIANS = ANGLE_DESIRED_RADIANS - ANGLE_ACCEPTABLE_WINDOW_RADIANS;
    const double ANGLE_PLUS_WINDOW_DEGREES = ANGLE_DESIRED_DEGREES + ANGLE_ACCEPTABLE_WINDOW_DEGREES;
    const double ANGLE_MINUS_WINDOW_DEGREES = ANGLE_DESIRED_DEGREES - ANGLE_ACCEPTABLE_WINDOW_DEGREES;
    //Note to developers: interaction with the strandedGrid is of the form: strandedGrid[lonCell][latCell]
    outputStrands.clear();
    size_t numCellsLon = strandedGrid.size();
    if(numCellsLon <= 2){
        //strandedGrid is empty
        cout << "ERROR: Error in StrandifyNorth(): strandedGrid.size() <= 2" <<endl;
        cout << "strandedGrid.size() = numCellsLon = " << numCellsLon <<endl;
        cout << "Note that only non-boarder cells contain useful data, so if gridsize<2 then there is no useful data."<<endl;
        return -154;
    }
    size_t numCellsLat = strandedGrid[0].size();
    for(auto strand : strandedGrid){
        if(strand.size() != numCellsLat){
            //strandedGrid should be a big rectangle of cells.  ...but it's not.
            cout << "ERROR: Error in StrandifyNorth(): strandedGrid is not rectangular."<<endl;
            return -248;
        }
    }
    //find the radius R within which neighbor points should be considered for strandification.
    double rMin;
    int errCode523 = GetGridSpacing(latLonList, rMin);
    if(errCode523<0){
        cout << "ERROR: Error in GetEastWestStrands::GetGridSpacing. Returned error code = " << errCode523<<endl;
    }
    const double R_OF_INTEREST = R_SAFETY_FACTOR * rMin;
    //And here we go...
    for(size_t selfCellIndexLon = 1; selfCellIndexLon < numCellsLon-1; ++selfCellIndexLon){
        for(size_t selfCellIndexLat = 1; selfCellIndexLat < numCellsLat-1; ++selfCellIndexLat){
            vector<strandedGeoPoint>& selfCell = strandedGrid[selfCellIndexLon][selfCellIndexLat];
            for(auto selfPointItr=selfCell.begin(); selfPointItr!=selfCell.end(); ++selfPointItr){
                strandedGeoPoint& selfPoint = *selfPointItr;
                size_t selfStrandIndex = strandByID[selfPoint.id];
                //for all the neighbor cells of this strandedGeoPoint selfPoint...
                //look at all the neighbor points and see if they pass the tests for strandification
                for(size_t ioiLon = selfCellIndexLon-1; ioiLon <= 1+selfCellIndexLon; ++ioiLon){
                    for(size_t ioiLat = selfCellIndexLat-1; ioiLat <= 1+selfCellIndexLat; ++ioiLat){
                        vector<strandedGeoPoint>& neighborCell = strandedGrid[ioiLon][ioiLat];
                        for(auto neighborPointItr = neighborCell.begin(); neighborPointItr != neighborCell.end(); ++neighborPointItr){
                            strandedGeoPoint& neighborPoint = *neighborPointItr;
                            size_t neighborStrandIndex = strandByID[neighborPoint.id];
                            //we don't want to repeatedly add the same (self)point to a strand
                            if(selfStrandIndex==neighborStrandIndex){
                                continue;
                            }
                            //if neightborPoint is outside radius R_OF_INTEREST, break to beginning of forloop
                            const double r = getDistanceBetweenPoints(selfPoint, neighborPoint);
                            if(r > R_OF_INTEREST){
                                continue;
                            }
                            //if potentialAngle > desiredAngle, break to beginning of forloop
                            const double angleDeg = getAngleDegrees(selfPoint, neighborPoint);
                            if(!IsAngleInside_Neg180to180(angleDeg, ANGLE_PLUS_WINDOW_DEGREES, ANGLE_MINUS_WINDOW_DEGREES)){
                                continue;
                            }
                            //All checks complete, join the strands...
                            std::vector<strandedGeoPoint>& neighborStrand = strandedList[neighborStrandIndex];
                            std::vector<strandedGeoPoint>& selfStrand = strandedList[selfStrandIndex];
                            int errCode53652 = JoinStrands(selfStrand, neighborStrand, strandByID);
                            if(errCode53652 < 0 ){
                                cout << "ERROR: Error in StrandifyNorth()::JoinStrands(). Returned errCode = " << errCode53652 <<endl;
                                return -8732;
                            }
                        }
                    }
                }
            }
        }
    }
    //all strands are joined, but in jumbled order.. sort them and copy to final output strands
    //we're done with the strandedGrid at this point, so we can clear it out
    int errCode8875 = SortAndFinalizeStrands(strandedList, outputStrands, directionMode);
    if(errCode8875){
        cout<<"ERROR: Failed to SortAndFinalizeStrands: Returned errCode = " << errCode8875<<endl;
    }
    //test that all strands are "reasonble".
    //namely, all angles and distances between
    //strand points are within reasonable parameters.
    for(auto strandItr = outputStrands.begin(); strandItr != outputStrands.end(); strandItr++){
        std::vector<geoPoint>& currStrand = *strandItr;
        if(currStrand.empty()){
            cout << "ERROR: Failed strandification test: outputStrands contains an empty strand." <<endl;
            return -9243;
        }
        size_t currStrandSize = currStrand.size();
        if(currStrandSize < 2){
            //all tests that follow require at lest two points.
            continue;
        }
        geoPoint p1 = currStrand[0];
        //for all neighboring points in the strand, make sure they're within radius r and angle window
        for(size_t index = 1; index < currStrandSize; ++index ){
            geoPoint p2 = currStrand[index];
            //if neightborPoint is outside radius R_OF_INTEREST, break

            const double r = getDistanceBetweenPoints(p1, p2);
            if(r > R_OF_INTEREST){
                cout << "ERROR: Failed strandification test: adjacent strand points are too far apart." <<endl;
                return -3352;
            }
            //if potentialAngle > desiredAngle, break to beginning of forloop
            const double angleDeg = getAngleDegrees(p1, p2);
            if(!IsAngleInside_Neg180to180(angleDeg, ANGLE_PLUS_WINDOW_DEGREES, ANGLE_MINUS_WINDOW_DEGREES)){
                cout << "ERROR: Failed strandification test: adjacent strand points are to great an angle from the desired direction." <<endl;
                return -3353;
            }
            //points pass the tests... move on to the next pair and test them...
            p1 = p2;
        }
    }
    return 0;
}

int StrandMaker::GetNorthSouthStrands(const string& latlonFilenameFullPath,
                                      std::vector<std::vector<geoPoint> >& northSouthStrands){
    //read in the points that need to be strandified.
    vector<geoPoint> latLonList;
    int errCode = GetLatLonList(latlonFilenameFullPath,latLonList,false);
    if(errCode < 0){
        cout << "Failed to GetLatLonList(). ErrorCode = " << errCode << endl;
        return -1;
    }
    //cout<<"Creating Grid for Strand-ization... "<<endl;
    std::vector<std::vector<std::vector<strandedGeoPoint> > > strandedGrid;
    std::vector< std::vector<strandedGeoPoint> > strandedList;
    std::vector<size_t> strandByID;
    int errCode5 = CreateStrandedGrid(latLonList, strandedGrid, strandedList, strandByID);
    if(errCode5 < 0){
        cout << "ERROR: Failed to CreateStrandedGrid(). Returned error code = " << errCode5 <<endl;
        return -3;
    }

    const std::string directionModeNorth = "north";
    int errCode6 = Strandify(strandedGrid, northSouthStrands, strandedList, latLonList, strandByID,directionModeNorth);
    if(errCode6 < 0){
        cout << "ERROR: Strandify() returned error code = " << errCode6 <<endl;
        return -4;
    }
    //cout << "Strandifcation Complete."<<endl;
    return 0;
}

int StrandMaker::GetEastWestStrands(const string& latlonFilenameFullPath,
               std::vector<std::vector<geoPoint> >& eastWestStrands){
    //read in the points that need to be strandified.
    vector<geoPoint> latLonList;
    int errCode = GetLatLonList(latlonFilenameFullPath,latLonList,false);
    if(errCode < 0){
        cout << "Failed to GetLatLonList(). ErrorCode = " << errCode << endl;
        return -1;
    }
    //cout<<"Creating Grid for Strand-ization... "<<endl;
    std::vector<std::vector<std::vector<strandedGeoPoint> > > strandedGrid;
    std::vector< std::vector<strandedGeoPoint> > strandedList;
    std::vector<size_t> strandByID;
    int errCode3 = CreateStrandedGrid(latLonList, strandedGrid, strandedList, strandByID);
    if(errCode3 < 0){
        cout << "ERROR: Failed to CreateStrandedGrid(). Returned error code = " << errCode3 <<endl;
        return -3;
    }

    //cout <<"Starting to Strandify()..."<<endl;
    const std::string directionModeEast = "east";
    int errCode4 = Strandify(strandedGrid, eastWestStrands, strandedList, latLonList, strandByID,directionModeEast);
    if(errCode4 < 0){
        cout << "ERROR: Strandify() returned error code = " << errCode4 <<endl;
        return -4;
    }
    //cout << "Strandifcation Complete."<<endl;
    return 0;
}

int StrandMaker::GetStrands(const string& latlonFilenameFullPath,
               std::vector<std::vector<geoPoint> >& eastWestStrands,
               std::vector<std::vector<geoPoint> >& northSouthStrands){
    //read in the points that need to be strandified.
    vector<geoPoint> latLonList;
    int errCode = GetLatLonList(latlonFilenameFullPath,latLonList,false);
    if(errCode < 0){
        cout << "Failed to GetLatLonList(). ErrorCode = " << errCode << endl;
        return -1;
    }
    //cout<<"Creating Grid for Strand-ization... "<<endl;
    std::vector<std::vector<std::vector<strandedGeoPoint> > > strandedGrid;
    std::vector< std::vector<strandedGeoPoint> > strandedList;
    std::vector<size_t> strandByID;
    int errCode3 = CreateStrandedGrid(latLonList, strandedGrid, strandedList, strandByID);
    if(errCode3 < 0){
        cout << "ERROR: Failed to CreateStrandedGrid(). Returned error code = " << errCode3 <<endl;
        return -3;
    }

    //cout <<"Starting to Strandify()..."<<endl;
    const std::string directionModeEast = "east";
    int errCode4 = Strandify(strandedGrid, eastWestStrands, strandedList, latLonList, strandByID,directionModeEast);
    if(errCode4 < 0){
        cout << "ERROR: Strandify() returned error code = " << errCode4 <<endl;
        return -4;
    }

    strandedGrid.clear();
    strandedList.clear();
    strandByID.clear();
    int errCode5 = CreateStrandedGrid(latLonList, strandedGrid, strandedList, strandByID);
    if(errCode5 < 0){
        cout << "ERROR: Failed to CreateStrandedGrid(). Returned error code = " << errCode5 <<endl;
        return -3;
    }

    const std::string directionModeNorth = "north";
    int errCode6 = Strandify(strandedGrid, northSouthStrands, strandedList, latLonList, strandByID,directionModeNorth);
    if(errCode6 < 0){
        cout << "ERROR: Strandify() returned error code = " << errCode6 <<endl;
        return -4;
    }
    //cout << "Strandifcation Complete."<<endl;
    return 0;
}



//////////////////////////////
/////For rerender strands/////
//////////////////////////////

inline int StrandMaker::getNumEntriesLatLonList(const std::string& latLonListFullPath_, size_t& numEntries){

    std::ifstream inFile;
    inFile.open(latLonListFullPath_);
    if(!inFile.is_open()){
        cout << "ERROR: Failed to open latLonList. Make sure the filename is valid." <<endl;
        return -74;
    }
    string headerLine;
    std::getline(inFile, headerLine);
    if(headerLine != expectedHeader){
        //This is not a valid input file
        //...or it is formatted differently than expected
        cout << "ERROR IN GetLatLonList(): first line was not of the "
             << "expected form. Should be = \"FID,Id,POINT_X,POINT_Y\""
             << " but instead was = \"" << headerLine << "\"" << endl;
        inFile.close();
        return -471;
    }
    //count lines and then reserve that number of entries in the map
    string token;
    numEntries = 0;
    while(getline(inFile,token)){
        ++numEntries;
    }
    inFile.close();
    return 0;
}

int StrandMaker::GetLatLonSet(const std::string& latLonListFullPath, std::unordered_set<string>& latLonSet){
    //count lines and then reserve that number of entries in the map
    string token;
    size_t numLines = 0;
    int errCountNumLines = getNumEntriesLatLonList(latLonListFullPath, numLines);
    if(errCountNumLines < 0){
        return -771;
    }
    cout << "Exected number of images = " << numLines << endl;
    latLonSet.clear();
    latLonSet.reserve(numLines);
    //now go through the file again, this time adding pairs to the map
    std::ifstream inFile;
    string headerLine;
    inFile.open(latLonListFullPath);
    if(!inFile.is_open()){
        cout << "ERROR: Failed to open latLonList. Make sure the filename is valid." <<endl;
        return -74;
    }
    std::getline(inFile, headerLine);
    if(headerLine != expectedHeader){
        //This is not a valid input file
        //...or it is formatted differently than expected
        latLonSet.clear();
        cout << "ERROR IN GetLatLonList(): first line was not of the "
             << "expected form. Should be = \"FID,Id,POINT_X,POINT_Y\""
             << " but instead was = \"" << headerLine << "\"" << endl;
        return -72;
    }
    double lat,lon;
    std::istringstream converter;
    string latLonStr;
    string tokens;
    cout.precision(std::numeric_limits<double>::max_digits10);
    while(getline(inFile,tokens)) {
        converter.clear();
        converter.str(tokens);
        int column = 0;
        while(getline(converter,token,',')){
            switch(column){
            case 0:{
                //cout << "Token 0 = " << token <<endl;
                break;}
            case 1:{
                //cout << "Token 1 = " << token <<endl;
                break;}
            case 2:{
                //cout << "Token 2 = " << token <<endl;
                lon = std::stod(token);
                break;}
            case 3:{
                //cout << "Token 3 = " << token <<endl;
                lat = std::stod(token);
                if(GetLatLon2IDKey(lat, lon, latLonStr)){
                    latLonSet.emplace(latLonStr);
                }else{
                    cout << "ERROR: Error in getMapLatLon2ID(): getLatLon2IDKey() failed." <<endl;
                    return -75;
                }
                break;}
            }
            column += 1;
        }
    }
    inFile.close();
    return 0;
}


inline bool StrandMaker::GetLatLon2IDKey(const double& lat, const double& lon, string& latLon2IDKey){
    string normLon, normLat;
    const int errCodeLon = normLatLonDeg(lon, normLon);
    if(errCodeLon<0){
        cout << ": ERROR IN getLatLon2IDKey(): normLatLonDeg() returned errCodeLon = "
             << errCodeLon <<endl;
        return false;
    }
    const int errCodeLat = normLatLonDeg(lat, normLat);
    if(errCodeLat<0){
        cout << ": ERROR IN getLatLon2IDKey(): normLatLonDeg() returned errCodeLat = "
             << errCodeLat <<endl;
        return false;
    }
    const string commaStr =  ",";
    latLon2IDKey = normLon + commaStr + normLat;
    return true;
}


inline int StrandMaker::normLatLonDeg(const double& latLonDeg, std::string& normalizedStr){
    const size_t latLonDeg_NumDigitsBeforeDecimal = 3;
    const size_t latLonDeg_NumDigitsAfterDecimal = 6;
    if(latLonDeg < 0.0){
        cout << ": ERROR: Error in normLatLonDeg(): attempted to normalize negative latLonDeg. "
             << "Currently this function is only capable of normalizing latLonDeg's that are between"
             <<  " 0.0 and 360.0." <<endl
             << "latLonDeg = " << latLonDeg <<endl;
        return -1;
    }
    if(latLonDeg > 360.0){
        cout <<": ERROR: Error in normLatLonDeg(): attempted to normalize negative latLonDeg. "
             << "Currently this function is only capable of normalizing latLonDeg's that are between"
             <<  " 0.0 and 360.0." <<endl
             << "latLonDeg = " << latLonDeg <<endl;
        return -2;
    }
    normalizedStr.clear();
    const string zeros(latLonDeg_NumDigitsBeforeDecimal+latLonDeg_NumDigitsAfterDecimal, '0');
    const string latLonDegBeforeNorm = std::to_string(latLonDeg);
    const size_t latLonDegBeforeNormSize = latLonDegBeforeNorm.size();
    if(latLonDegBeforeNormSize < 1){
        cout <<": ERROR: Error in normLatLonDeg(): latLonDegBeforeNorm has size < 1."<<endl;
        return -3;
    }
    const string decimalStr = ".";
    const size_t decPos = latLonDegBeforeNorm.find(decimalStr);
    if(decPos==std::string::npos){
        //there is no decimal, this looks just like an int
        //ex: 46   ---> 046.000000
        //but first check to make sure it's not longer than expected
        if(latLonDegBeforeNormSize > latLonDeg_NumDigitsBeforeDecimal){
            cout << ": ERROR: Error in normLatLonDeg(): encounterd latLonDeg string with more than "
                 << latLonDeg_NumDigitsBeforeDecimal << " digits before the decimal." <<endl
                 << "latLonDeg string before normalization = "
                 << latLonDegBeforeNorm <<endl;
            return -44;
        }
        normalizedStr = zeros.substr(0,latLonDeg_NumDigitsBeforeDecimal-latLonDegBeforeNormSize)
                + latLonDegBeforeNorm + decimalStr
                + zeros.substr(0,latLonDeg_NumDigitsAfterDecimal);
    }else{
        //ex:46.435 would yield decPos = 2.
        const size_t lenAfterDecimal = latLonDegBeforeNormSize-decPos-1;
        if(lenAfterDecimal>latLonDeg_NumDigitsAfterDecimal){
            cout << ": ERROR: Error in rgb2gray(): encounterd latDeg string with more than "
                 << latLonDeg_NumDigitsAfterDecimal << " digits after the decimal." <<endl
                 << "latLonDeg string before normalization = "
                 << latLonDegBeforeNorm <<endl;
            return -45;
        }
        normalizedStr = zeros.substr(0,latLonDeg_NumDigitsBeforeDecimal-decPos)
                + latLonDegBeforeNorm
                + zeros.substr(0,latLonDeg_NumDigitsAfterDecimal-lenAfterDecimal);
    }
    return 0;
}

