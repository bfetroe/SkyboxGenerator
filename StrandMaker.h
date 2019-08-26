#ifndef STRANDMAKER_H
#define STRANDMAKER_H

#include <iostream>     //cout
#include <vector>
#include <fstream>      //ifstream and ofstream
#include <sstream>
#include <algorithm>    //sort
#include <string>
#include <unordered_set>


#ifndef GEOPOINT
#define GEOPOINT
struct geoPoint{
    double lat = 0.0;
    double lon = 0.0;
};
#endif

struct strandedGeoPoint {
    double lat = 0.0;
    double lon =0.0;
    size_t id = 0;
};

class StrandMaker
{
public:
    int GetEastWestStrands(const std::string& latlonFilenameFullPath,
                           std::vector<std::vector<geoPoint> >& eastWestStrands);
    int GetNorthSouthStrands(const std::string& latlonFilenameFullPath,
                             std::vector<std::vector<geoPoint> >& northSouthStrands);
    int GetStrands(const std::string& latlonFilenameFullPath,
                   std::vector<std::vector<geoPoint> >& eastWestStrands,
                   std::vector<std::vector<geoPoint> >& northSouthStrands);

    int GetLatLonSet(const std::string& latLonListFullPath, std::unordered_set<std::string>& latLonSet);
    inline bool GetLatLon2IDKey(const double& lat, const double& lon, std::string& latLon2IDKey);

    template <class point>
    int SaveStrandsToFile(std::vector<std::vector<point> >& strands, const std::string& fullPathAndFilename);

    template <class point>
    void PrintVector(std::vector<point>& vec);



private:

    int GetLatLonList(const std::string& filenameFullPath, std::vector<geoPoint>& latLonList, const bool appendMode);
    int GetGridSpacing(std::vector<geoPoint>& fullLatLonList, double& gridSpacing);
    int GetBoundingBox(std::vector<geoPoint>& latLonList,
                       double& latMin,
                       double& latMax,
                       double& lonMin,
                       double& lonMax);

    int CreateStrandedGrid(std::vector<geoPoint>& latLonList,
                           std::vector<std::vector<std::vector<strandedGeoPoint> > >& strandedGrid,
                           std::vector< std::vector<strandedGeoPoint> >& strandedList,
                           std::vector<size_t>& strandByID);

    template <class point>
    inline const double getAngleRadians(point& p1, point& p2);
    template <class point>
    const double getDistanceBetweenPoints(point& p1, point& p2);
    template <class point>
    const double getAngleDegrees(point& p1, point& p2);
    bool IsAngleInside_0to360(double target, double angle1, double angle2);
    bool IsAngleInside_Neg180to180(double target, double angle1, double angle2);
    int JoinStrands(std::vector<strandedGeoPoint>& selfStrand,
                    std::vector<strandedGeoPoint>& neighborStrand,
                    std::vector<size_t>& strandByID);
    //bool CompareStrandedGeoPointByLat(strandedGeoPoint& p1, strandedGeoPoint& p2);
    int SortAndFinalizeStrands(std::vector<std::vector<strandedGeoPoint> >& jumbledStrands,
                               std::vector<std::vector<geoPoint> >& finalStrands,
                               const std::string &directionMode);
    int Strandify(std::vector<std::vector<std::vector<strandedGeoPoint> > >& strandedGrid,
                       std::vector<std::vector<geoPoint> >& strands,
                       std::vector< std::vector<strandedGeoPoint> >& strandedList,
                       std::vector<geoPoint>& latLonList,
                       std::vector<size_t>& strandByID,
                       const std::string& directionMode);

    inline int getNumEntriesLatLonList(const std::string& latLonListFullPath_, size_t& numEntries);
    inline int normLatLonDeg(const double& latLonDeg, std::string& normalizedStr);

    const std::string expectedHeader = "FID,Id,POINT_X,POINT_Y";
};

#endif // STRANDMAKER_H
