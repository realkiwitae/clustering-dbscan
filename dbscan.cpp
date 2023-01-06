#include "dbscan.h"
#include <chrono>
#include <unistd.h>
#include <iostream>
#include "utils/cvutils.h"

using namespace std;
using namespace std::chrono;

#define FPS 10
#define SLEEP_TIME 1000000/FPS
high_resolution_clock::time_point start;
high_resolution_clock::time_point start_move;
#define NOW high_resolution_clock::now()
#define TIME duration_cast<duration<double>>(NOW - start).count()
#define TIME_MOVE duration_cast<duration<double>>(NOW - start_move).count()
#define POINT_RADIUS 10

#define W 1000
#define H 1000
#define NB_COLORS 200

class DBSCANPoint:public ClusterPoint{
    public:
    bool bcore = false;
    std::vector<DBSCANPoint*> neighbours;
    DBSCANPoint():ClusterPoint(){}
    DBSCANPoint(double x, double y):ClusterPoint(x,y){}
};

bool bActiveDraw = false;
int nbClusters = 0;
Point mouse;
std::vector<DBSCANPoint> spots = {DBSCANPoint(500,500)};
std::vector<cv::Scalar> colors;

Point makeRandomSpot(Point c, double radius){
    double angle = (double)std::rand() / RAND_MAX * M_PI * 2;
    double x = c.x + radius*std::cos(angle);
    double y = c.y + radius*std::sin(angle);
    return Point(x,y);
}

double normalRandom(double mean, double sigma){
    double v1 = (double)std::rand() / RAND_MAX;
    double v2 = (double)std::rand() / RAND_MAX;
    return sigma * std::cos(2*M_PI*v2) * std::sqrt(-2*std::log(v1))+mean;
}

void update();
void runClustering();
bool collide(Point a, Point b, double radius);
void onMouse(int event, int x, int y, int, void*);
bool cost(Point a, Point b){
    return norm2(a,b) <= 6*POINT_RADIUS;
}
void DBSCAN(std::vector<DBSCANPoint>& spots ,int& nbClusters, int minPts, FnBPtr costfn); 

void startRenderer(){
    randomColors(colors, NB_COLORS);
    while (true)
    {
        start = NOW;
        newFrame(cv::Size(W,H),cv::Scalar());
        update();
        int us = SLEEP_TIME - TIME*1e6;
        if(us>0)usleep(us);
    }
}

void update(){
    DBSCAN(spots ,nbClusters,3, cost);

    for(DBSCANPoint p : spots){
        drawCircle(p.x,p.y,POINT_RADIUS, colors[p.clusterID%colors.size()]);
    }
    drawCircle(mouse.x,mouse.y,POINT_RADIUS/2,cv::Scalar(255,0,0));
    cv::setMouseCallback("dbscan",onMouse,0);
    show("dbscan",1);
}

void runClustering(){

}

bool collide(Point a, Point b, double radius){
    return norm2(a,b) <= 2*radius;
}

void onMouse(int event, int x, int y, int, void*){
    mouse.x = x;
    mouse.y = y;
	if ( event == cv::EVENT_LBUTTONDOWN ){
        bActiveDraw = ! bActiveDraw;
    }else if(event == cv::EVENT_MOUSEMOVE && bActiveDraw){
        Point p = makeRandomSpot(Point(x,y),normalRandom(0,0));
        bool badd = true;
        for(Point a : spots){
            if(collide(p,a,POINT_RADIUS)){
                badd = false;
                break;
            }
        }
        if(badd)spots.push_back(DBSCANPoint(p.x,p.y));
    }

}



void DBSCAN(std::vector<DBSCANPoint>& spots , int& nbClusters, int minPts, FnBPtr costfn){
    for(auto& p : spots){
        p.clusterID = 0;
        p.bcore = false;
        p.neighbours.clear();
    }
    for(int i = 0; i < (int)spots.size();i++){
        for(int j = i+1; j < (int)spots.size();j++){
            if(costfn(spots[i],spots[j])){
                spots[i].neighbours.push_back(&spots[j]);
                spots[j].neighbours.push_back(&spots[i]);                
            }
        }
    }
    for(auto& p : spots){
        p.bcore = (int)p.neighbours.size() >= minPts;
    }

    std::vector<DBSCANPoint*> visited = {};
    int cidx = 0;
    
    for(auto& s : spots){
        if(!s.bcore)continue;
        if(std::find(visited.begin(),visited.end(),&s) != visited.end())continue;
        cidx++;
        std::vector<DBSCANPoint*> queue = {&s};
        while(queue.size()>0){
            DBSCANPoint* p = queue.back();
            queue.pop_back();
            if(std::find(visited.begin(),visited.end(),p) != visited.end())continue;
            p->clusterID = cidx;
            visited.push_back(p);
            for(DBSCANPoint* n : p->neighbours){
                queue.push_back(n);
            }
        }
    }

    nbClusters = cidx;
}
