#pragma once

#include <ostream>
#include <iterator>
#include <vector>
#include <string>
#include <tuple>
#include <climits>


#include <iostream> 
#include <set>
#include <algorithm>

namespace ISPDParser {

struct TwoPin;
struct RPoint;
struct Point;
class Net;

struct Point {

public:
    int x, y, z;

    Point(void) : x(0), y(0), z(0) {}
    Point(int x, int y) : x(x), y(y), z(0) {}
    Point(int x, int y, int z) : x(x), y(y), z(z) {}
    Point(const Point &p) : x(p.x), y(p.y), z(p.z) {}
    ~Point() {}

    int dist_to(Point to){
        return abs(x - to.x) + abs(y - to.y);
    }
};
// Route Point
struct RPoint {

public:
    int x, y, z;
    bool hori;

    RPoint() { z = 0; }
    RPoint(int x, int y, bool h) : x(x), y(y), hori(h) { z = 0; }
    RPoint(int x, int y, int z, bool h) : x(x), y(y), z(z), hori(h) {}
    RPoint(const RPoint &p) : x(p.x), y(p.y), z(p.z), hori(p.hori) {}
};

struct TwoPin
{
    Point from, to;
    int wlen;
    bool ripup;
    bool overflow;
    Net *parNet;
    bool reroute;
    std::vector<RPoint> path;
    std::vector<Point> corr_path;

    TwoPin()
    {
        overflow = false;
        ripup = false;
        reroute = true;
    }

    TwoPin(const TwoPin &t)
    {
        overflow = t.overflow;
        from = t.from;
        to = t.to;
        ripup = false;
        reroute = true;
        parNet = t.parNet;
    }

    TwoPin(const Point &f, const Point &t){
        from = f;
        to = t;
        overflow = false; 
        ripup = false;
        reroute = true;
    }

    int HPWL()
    {
        return abs(from.x - to.x) + abs(from.y - to.y);
    }
    void create_path(){
        path.clear();
        for(int i = 0; i < corr_path.size() - 1; ++i){
            int x1 = corr_path[i].x, y1 = corr_path[i].y, z1 = corr_path[i].z;
            int x2 = corr_path[i + 1].x, y2 = corr_path[i + 1].y, z2 = corr_path[i + 1].z;
            //go right
            if(x1 + 1 == x2 && y1 == y2){
                path.emplace_back(x1, y1, true);
            }
            //go left
            else if(x1 == x2 + 1 && y1 == y2){
                path.emplace_back(x2, y2, true);
            }
            //go up
            else if(x1 == x2 && y1 + 1 == y2){
                path.emplace_back(x1, y1, false);
            }
            //go down
            else if(x1 == x2 && y1 == y2 + 1){
                path.emplace_back(x2, y2, false);
            }
            else {
                std::cerr<<"error in creating path"<<x1<<", "<<y1<<" -> "<<x2<<", "<<y2<<std::endl;
                exit(-1);
            }
        }
    }
    bool is_same_row(){return from.y == to.y;}
    bool is_same_col(){return from.x == to.x;}
};
class Net {

public:

    Net(const std::string &name, int id, int numPins, int minimumWidth)
              : name(name), id(id), numPins(numPins), minimumWidth(minimumWidth)  {}

    std::string name;
    int id;
    int numPins;
    int minimumWidth;
    std::vector<std::tuple<int, int, int>> pins;

    // Grid coordinates
    std::vector<Point> pin2D;
    std::vector<Point> pin3D;

    // Two pin nets
    std::vector<TwoPin> twopin;


    void decompose(){
        std::vector<int> dist(pin2D.size(), INT32_MAX);
        std::vector<int> parent(pin2D.size(), -1);
        std::vector<bool> isVisit(pin2D.size(), false);
        std::set<int> remainPin;
        for(int i = 0; i < pin2D.size(); i++){remainPin.insert(i);}

        isVisit[0] = true;
        remainPin.erase(0);
        for(int i = 1; i < pin2D.size(); i++){
            dist[i] = pin2D[0].dist_to(pin2D[i]);
            parent[i] = 0;
        }
        
        while(!remainPin.empty()){
            int minDist = INT32_MAX;
            int minPin = -1;
            for(auto& p: remainPin){
                if(dist[p] < minDist){
                    minDist = dist[p];
                    minPin = p;
                }
            }
            
            isVisit[minPin] = true;
            remainPin.erase(minPin);
            for(int i = 0; i < pin2D.size(); i++){
                if(!isVisit[i] && pin2D[minPin].dist_to(pin2D[i]) < dist[i]){
                    dist[i] = pin2D[minPin].dist_to(pin2D[i]);
                    parent[i] = minPin;
                }
            }
        }


        /*setting the Net.twopin*/
        for(int i = 0; i < pin2D.size(); ++i){
            if(parent[i] == -1){continue;}
            TwoPin t_pin(pin3D[i], pin3D[parent[i]]);
            twopin.push_back(t_pin);
        }
    }

    void set_net(Net* net){
        for(int i = 0; i < pin2D.size() - 1; i++){
            twopin[i].parNet = net;
        }
    }
    void sort_twopin() {
        std::sort(twopin.begin(), twopin.end(), 
            [](TwoPin a, TwoPin b) {
                return a.HPWL() < b.HPWL();
            });
    }
};

class CapacityAdj {

public:
    std::tuple<int, int, int> grid1;
    std::tuple<int, int, int> grid2;
    int reducedCapacityLevel;

};

class ispdData {

public:
    ispdData() = default;
    ~ispdData() {
        for (Net *net : nets)
            delete net;
        for (CapacityAdj *capacityAdj : capacityAdjs)
            delete capacityAdj;
    }

    int numXGrid;
    int numYGrid;
    int numLayer;

    std::vector<int> verticalCapacity;
    std::vector<int> horizontalCapacity;
    std::vector<int> minimumWidth;
    std::vector<int> minimumSpacing;
    std::vector<int> viaSpacing;
    
    int lowerLeftX;
    int lowerLeftY;
    int tileWidth;
    int tileHeight;

    int numNet;
    std::vector<Net *> nets;

    int numCapacityAdj;
    std::vector<CapacityAdj *> capacityAdjs;

    friend std::ostream& operator<<(std::ostream &os, const ispdData &data) {
        os << "grid " << data.numXGrid << " " << data.numYGrid << " " << data.numLayer;
        os << "\nvertical capacity ";
        std::copy(data.verticalCapacity.begin(), data.verticalCapacity.end(), 
                std::ostream_iterator<int>(os, " "));
        os << "\nhorizontal capacity ";
        std::copy(data.horizontalCapacity.begin(), data.horizontalCapacity.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nminimum width  ";
        std::copy(data.minimumWidth.begin(), data.minimumWidth.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nminimum spacing  ";
        std::copy(data.minimumSpacing.begin(), data.minimumSpacing.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nvia spacing ";
        std::copy(data.viaSpacing.begin(), data.viaSpacing.end(),
            std::ostream_iterator<int>(os, " "));
        os << "\n" << data.lowerLeftX << " " << data.lowerLeftY << " " 
                    << data.tileWidth << " " << data.tileHeight << "\n";
        os << "num net " << data.numNet << "\n";
        for (const auto *net : data.nets) {
            os << net->name << " " << net->id << " " << net->numPins << " " << net->minimumWidth << "\n";
            for (const auto &coord : net->pins)
            os << std::get<0>(coord) << " " << std::get<1>(coord) << " " << std::get<2>(coord) << "\n";
        }
        os << data.numCapacityAdj << "\n";
        for (const auto *adj : data.capacityAdjs)
            os << std::get<0>(adj->grid1) << " " << std::get<1>(adj->grid1) << " " << std::get<2>(adj->grid1) << " "
            << std::get<0>(adj->grid2) << " " << std::get<1>(adj->grid2) << " " << std::get<2>(adj->grid2) << " "
            << adj->reducedCapacityLevel << "\n";
        return os;
    }

};

static ispdData* parse(std::istream &is) {
    ispdData *data = new ispdData();
    std::string keyword;
    is >> keyword >> data->numXGrid >> data->numYGrid >> data->numLayer;
    is >> keyword >> keyword;
    data->verticalCapacity.clear();
    int val;
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->verticalCapacity.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->horizontalCapacity.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->horizontalCapacity.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->minimumWidth.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->minimumWidth.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->minimumSpacing.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->minimumSpacing.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->viaSpacing.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->viaSpacing.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> data->lowerLeftX >> data->lowerLeftY >> data->tileWidth >> data->tileHeight;
    is >> keyword >> keyword >>data->numNet;
    data->nets.clear();
    for (int i = 0; i < data->numNet; i++) {
        std::string net_name;
        int id, num_pins, min_width;
        is >> net_name >> id >> num_pins >> min_width;
        Net *net = new Net{net_name, id, num_pins, min_width};
        for (int j = 0; j < num_pins; j++) {
            int x, y, z;
            is >> x >> y >> z;
            net->pins.push_back(std::make_tuple(x, y, z));
        }
        data->nets.push_back(net);
    }
    is >> data->numCapacityAdj;
    data->capacityAdjs.clear();
    for (int i = 0; i < data->numCapacityAdj; i++) {
        int x1, y1, z1, x2, y2, z2, reduced_capacity_level;
        is >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> reduced_capacity_level;
        CapacityAdj *adj = new CapacityAdj{{x1, y1, z1}, {x2, y2, z2}, reduced_capacity_level};
        data->capacityAdjs.push_back(adj);
    }
    return data;
}


struct twopinInfo {
    TwoPin* twopin;
    int of;
    twopinInfo();
    twopinInfo(TwoPin* _twopin, int _of){
        twopin = _twopin;
        of = _of;
    }
    bool operator<(const twopinInfo& other) const {
        return of < other.of; // 小于关系
    }
};



}
