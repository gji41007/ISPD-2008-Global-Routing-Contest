#include "ispdData.h"
#include "LayerAssignment.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <queue>

#include <cmath>


int OF(ISPDParser::Point a, ISPDParser::Point b, int minwidth, int minSpace, 
           std::vector<std::vector<int>>& vertCap, std::vector<std::vector<int>>& horiCap,
           std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    int x1 = a.x, y1 = a.y;
    int x2 = b.x, y2 = b.y;
    // if(x2==323){
    //     std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
    //     std::cout<<vertCap.size()<<vertCap[0].size()<<std::endl;
    //     std::cout<<horiCap.size()<<horiCap[0].size()<<std::endl;
    //     std::cout<<vertCurr.size()<<vertCurr[0].size()<<std::endl;
    //     std::cout<<horiCurr.size()<<horiCurr[0].size()<<std::endl;
    // }
    //go right
    if(x1 + 1 == x2 && y1 == y2){
        return std::max((horiCurr[x1][y1] + minwidth + minSpace - horiCap[x1][y1]), 0);
    }
    //go left
    else if(x1 == x2 + 1 && y1 == y2){
        return std::max((horiCurr[x2][y2] + minwidth + minSpace - horiCap[x2][y2]), 0);
    }
    //go up
    else if(x1 == x2 && y1 + 1 == y2){
        return std::max((vertCurr[x1][y1] + minwidth + minSpace - vertCap[x1][y1]), 0);
    }
    //go down
    else if(x1 == x2 && y1 == y2 + 1){
        return std::max((vertCurr[x2][y2] + minwidth + minSpace - vertCap[x2][y2]), 0);
    }
    else{
        std::cerr<<"error in OF: not adjecnt grid\n";
        exit(-1);
    }
}
double cost(ISPDParser::Point a, ISPDParser::Point b, int minwidth, int minSpace, 
           std::vector<std::vector<int>>& vertCap, std::vector<std::vector<int>>& horiCap,
           std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){

    /* param */
    double h1 = 1;
    double h2 = 150;
    double k = -0.3;

    int x1 = a.x, y1 = a.y;
    int x2 = b.x, y2 = b.y;
    int d, cap;
    //go right
    if(x1 + 1 == x2 && y1 == y2){
        d = horiCurr[x1][y1] + minwidth + minSpace;
        cap = horiCap[x1][y1];
    }
    //go left
    else if(x1 == x2 + 1 && y1 == y2){
        d = horiCurr[x2][y2] + minwidth + minSpace;
        cap = horiCap[x2][y2];
    }
    //go up
    else if(x1 == x2 && y1 + 1 == y2){
        d = horiCurr[x1][y1] + minwidth + minSpace;
        cap = horiCap[x1][y1];
    }
    //go down
    else if(x1 == x2 && y1 == y2 + 1){
        d = horiCurr[x2][y2] + minwidth + minSpace;
        cap = horiCap[x2][y2];
    }
    else{
        std::cerr<<"error in cost(): not adjecnt grid\n";
        exit(-1);
    }


    double exponent = -k * (d - cap); 
    double denominator = 1 + std::exp(exponent); 
    return h1 + (h2 / denominator);

}
int currOF(ISPDParser::Point a, ISPDParser::Point b, 
           std::vector<std::vector<int>>& vertCap, std::vector<std::vector<int>>& horiCap,
           std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    int x1 = a.x, y1 = a.y;
    int x2 = b.x, y2 = b.y;
    //go right
    if(x1 + 1 == x2 && y1 == y2){
        return std::max((horiCurr[x1][y1] - horiCap[x1][y1]), 0);
    }
    //go left
    else if(x1 == x2 + 1 && y1 == y2){
        return std::max((horiCurr[x2][y2] - horiCap[x2][y2]), 0);
    }
    //go up
    else if(x1 == x2 && y1 + 1 == y2){
        return std::max((vertCurr[x1][y1] - vertCap[x1][y1]), 0);
    }
    //go down
    else if(x1 == x2 && y1 == y2 + 1){
        return std::max((vertCurr[x2][y2] - vertCap[x2][y2]), 0);
    }
    else{
        std::cerr<<"error in currOF: not adjecnt grid\n";
        exit(-1);
    }
}
void fill_cap(ISPDParser::Net* net, std::vector<ISPDParser::Point>& corr_path, int minSpace, std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    if(corr_path.size() <= 1){return;}
    for(int i = 0; i < corr_path.size() - 1; ++i){
        int x1 = corr_path[i].x, y1 = corr_path[i].y;
        int x2 = corr_path[i+1].x, y2 = corr_path[i+1].y;
        // std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
        //go right
        if(x1 + 1 == x2 && y1 == y2){
            horiCurr[x1][y1] += (net->minimumWidth + minSpace);
        }
        //go left
        else if(x1 == x2 + 1 && y1 == y2){
            horiCurr[x2][y2] += (net->minimumWidth + minSpace);
        }
        //go up
        else if(x1 == x2 && y1 + 1 == y2){
            vertCurr[x1][y1] += (net->minimumWidth + minSpace);
        }
        //go down
        else if(x1 == x2 && y1 == y2 + 1){
            vertCurr[x2][y2] += (net->minimumWidth + minSpace);
        }
        else{
            std::cerr<<"error in fill_cap: not adjecnt grid\n";
            exit(-1);
        }
    }
    
}
void remove_cap(ISPDParser::Net* net, std::vector<ISPDParser::Point>& corr_path, int minSpace, std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    if(corr_path.size() <= 1){return;}
    for(int i = 0; i < corr_path.size() - 1; ++i){
        int x1 = corr_path[i].x, y1 = corr_path[i].y;
        int x2 = corr_path[i+1].x, y2 = corr_path[i+1].y;
        // std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
        //go right
        if(x1 + 1 == x2 && y1 == y2){
            horiCurr[x1][y1] -= (net->minimumWidth + minSpace);
        }
        //go left
        else if(x1 == x2 + 1 && y1 == y2){
            horiCurr[x2][y2] -= (net->minimumWidth + minSpace);
        }
        //go up
        else if(x1 == x2 && y1 + 1 == y2){
            vertCurr[x1][y1] -= (net->minimumWidth + minSpace);
        }
        //go down
        else if(x1 == x2 && y1 == y2 + 1){
            vertCurr[x2][y2] -= (net->minimumWidth + minSpace);
        }
        else{
            std::cerr<<"error in remove_cap: not adjecnt grid\n";
            exit(-1);
        }
    }
    
}
std::vector<ISPDParser::Point> pattern_route(ISPDParser::TwoPin twopin, int minSpace, std::vector<std::vector<int>>& vertCap, std::vector<std::vector<int>>& horiCap,
                                              std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    std::vector<ISPDParser::Point> shortest_path;
    int minSumOF = INT32_MAX;
    if(twopin.HPWL() == 0){return shortest_path;}
    int x1 = twopin.from.x, y1 = twopin.from.y;
    int x2 = twopin.to.x, y2 = twopin.to.y;
    //same column
    if(x1 == x2){
        if(y1 < y2){
            for(int i = y1; i <= y2; ++i){
                shortest_path.push_back(ISPDParser::Point(x1, i));
            }
        }
        else{
            for(int i = y1; i >= y2; --i){
                shortest_path.push_back(ISPDParser::Point(x1, i));
            }
        }
        return shortest_path;
    }
    //same row
    else if(y1 == y2){
        if(x1 < x2){
            for(int i = x1; i <= x2; ++i){
                shortest_path.push_back(ISPDParser::Point(i, y1));
            }
        }
        else{
            for(int i = x1; i >= x2; --i){
                shortest_path.push_back(ISPDParser::Point(i, y1));
            }
        }
        return shortest_path;
    }


    /* upper L shape*/
    int pattSumCost = 0;
    if(y1 < y2){
        for(int i = y1; i < y2; ++i){
            shortest_path.push_back(ISPDParser::Point(x1, i));
            pattSumCost += OF(ISPDParser::Point(x1, i), ISPDParser::Point(x1, i+1), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
        }
        shortest_path.push_back(ISPDParser::Point(x1, y2));
        if(x1 < x2){
            for(int i = x1 + 1; i <= x2; ++i){
                shortest_path.push_back(ISPDParser::Point(i, y2));
                pattSumCost += OF(ISPDParser::Point(i-1, y2), ISPDParser::Point(i, y2), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        else{
            for(int i = x1 - 1; i >= x2; --i){
                shortest_path.push_back(ISPDParser::Point(i, y2));
                pattSumCost += OF(ISPDParser::Point(i+1, y2), ISPDParser::Point(i, y2), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }
    else{
        if(x1 < x2){
            for(int i = x1; i <= x2 - 1; ++i){
                shortest_path.push_back(ISPDParser::Point(i, y1));
                pattSumCost += OF(ISPDParser::Point(i, y1), ISPDParser::Point(i+1, y1), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        else{
            for(int i = x1; i >= x2 + 1; --i){
                shortest_path.push_back(ISPDParser::Point(i, y1));
                pattSumCost += OF(ISPDParser::Point(i, y1), ISPDParser::Point(i-1, y1), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        shortest_path.push_back(ISPDParser::Point(x2, y1));
        for(int i = y1 - 1; i >= y2; --i){
            shortest_path.push_back(ISPDParser::Point(x2, i));
            pattSumCost += OF(ISPDParser::Point(x2, i+1), ISPDParser::Point(x2, i), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
        }
    }
    minSumOF = pattSumCost;
    // std::cout<<pattSumCost<<" "<<minSumOF<<" ";

    /* lower L shape*/
    pattSumCost = 0;
    if(y1 < y2){
        if(x1 < x2){
            for(int i = x1; i <= x2 - 1; ++i){
                // shortest_path.push_back(ISPDParser::Point(i, y1));
                pattSumCost += OF(ISPDParser::Point(i, y1), ISPDParser::Point(i+1, y1), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        else{
            for(int i = x1; i >= x2 + 1; --i){
                // shortest_path.push_back(ISPDParser::Point(i, y1));
                pattSumCost += OF(ISPDParser::Point(i, y1), ISPDParser::Point(i-1, y1), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        // shortest_path.push_back(ISPDParser::Point(x2, y1));
        for(int i = y1 + 1; i <= y2; ++i){
            // shortest_path.push_back(ISPDParser::Point(x2, i));
            pattSumCost += OF(ISPDParser::Point(x2, i-1), ISPDParser::Point(x2, i), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
        }
    }
    else{
        for(int i = y1; i > y2; --i){
            // shortest_path.push_back(ISPDParser::Point(x1, i));
            pattSumCost += OF(ISPDParser::Point(x1, i-1), ISPDParser::Point(x1, i), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
        }
        // shortest_path.push_back(ISPDParser::Point(x1, y2));
        if(x1 < x2){
            for(int i = x1 + 1; i <= x2; ++i){
                // shortest_path.push_back(ISPDParser::Point(i, y2));
                pattSumCost += OF(ISPDParser::Point(i-1, y2), ISPDParser::Point(i, y2), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
        else{
            for(int i = x1 - 1; i >= x2; --i){
                // shortest_path.push_back(ISPDParser::Point(i, y2));
                pattSumCost += OF(ISPDParser::Point(i+1, y2), ISPDParser::Point(i, y2), twopin.parNet->minimumWidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }
    if(pattSumCost < minSumOF){
        minSumOF = pattSumCost;
        shortest_path.clear();
        if(y1 < y2){
            if(x1 < x2){
                for(int i = x1; i <= x2 - 1; ++i){
                    shortest_path.push_back(ISPDParser::Point(i, y1));
                }
            }
            else{
                for(int i = x1; i >= x2 + 1; --i){
                    shortest_path.push_back(ISPDParser::Point(i, y1));
                }
            }
            shortest_path.push_back(ISPDParser::Point(x2, y1));
            for(int i = y1 + 1; i <= y2; ++i){
                shortest_path.push_back(ISPDParser::Point(x2, i));
            }
        }
        else{
            for(int i = y1; i > y2; --i){
                shortest_path.push_back(ISPDParser::Point(x1, i));
            }
            shortest_path.push_back(ISPDParser::Point(x1, y2));
            if(x1 < x2){
                for(int i = x1 + 1; i <= x2; ++i){
                    shortest_path.push_back(ISPDParser::Point(i, y2));
                }
            }
            else{
                for(int i = x1 - 1; i >= x2; --i){
                    shortest_path.push_back(ISPDParser::Point(i, y2));
                }
            }
        }
    }
    return shortest_path;
}
std::vector<std::vector<std::vector<ISPDParser::Point>>> unilateral_mono_route(int x1, int y1, int x2, int y2, ISPDParser::Point source, std::vector<std::vector<int>>& costMap, int minwidth, int minSpace, 
                           std::vector<std::vector<int>>& vertCap, std::vector<std::vector<int>>& horiCap, std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr,
                           std::vector<std::vector<ISPDParser::Point>>& pVert, std::vector<std::vector<ISPDParser::Point>>& pHori){
    int s_x = source.x, s_y = source.y;
    /* vert mono route*/
    std::vector<std::vector<int>> costVert(
        costMap.size(),
        std::vector<int>(costMap[0].size(), 0)
    );
    // std::vector<std::vector<ISPDParser::Point>> pVert(
    //     costMap.size(),
    //     std::vector<ISPDParser::Point>(costMap[0].size(), ISPDParser::Point(-1, -1))
    // );
    // std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
    pVert[s_x][s_y] = ISPDParser::Point(s_x, s_y);
    for(int i = source.x - 1; i >= x1; --i){
        costVert[i][s_y] = (costVert[i+1][s_y] + OF(ISPDParser::Point(i, s_y), ISPDParser::Point(i+1, s_y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));
        pVert[i][s_y] = ISPDParser::Point(i+1, s_y);
    }
    for(int i = source.x + 1; i <= x2; ++i){
        costVert[i][s_y] = (costVert[i-1][s_y] + OF(ISPDParser::Point(i-1, s_y), ISPDParser::Point(i, s_y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));
        pVert[i][s_y] = ISPDParser::Point(i-1, s_y);
    }
    
    // std::cout<<"test"<<std::endl;
    //Top
    for(int y = s_y + 1; y <= y2; ++y){
        for(int x = x1; x <= x2; ++x){
            pVert[x][y] = ISPDParser::Point(x, y-1);
        }
    }
    for(int y = s_y + 1; y <= y2; ++y){
        std::vector<int> rowCost(x2 - x1 + 1, INT32_MAX);
        std::vector<bool> isVisit(x2 - x1 + 1, false); 
        //setting initial cost with verical edge  
        for(int x = x1; x <= x2; ++x){rowCost[x-x1] = (costVert[x][y-1] + OF(ISPDParser::Point(x, y), ISPDParser::Point(x, y-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));}
        
    // std::cout<<x1<<" "<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<s_x<<" "<<s_y<<std::endl;
    // std::cout<<x2 - x1 + 1<<" "<<y<<std::endl;
        for(int _ = 0; _ < x2 - x1 + 1; ++_){
            int minCost = INT32_MAX;
            int minX = -1;
            for(int x = x1; x <= x2; ++x){
                if(!isVisit[x-x1] && rowCost[x-x1] < minCost){
                    minCost = rowCost[x-x1];
                    minX = x;
                }
            }
            if(minX == -1){std::cout<<"--"<<std::endl;}
            isVisit[minX-x1] = true;
            costVert[minX][y] = rowCost[minX-x1];
            //adjust left and right
            if(minX - 1 >= x1 && !isVisit[minX -x1 - 1] && rowCost[minX-x1-1] > rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX-1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pVert[minX-1][y] = ISPDParser::Point(minX, y);
                rowCost[minX-x1-1] = rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX-1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
            if(minX + 1 <= x2 && !isVisit[minX -x1 + 1] && rowCost[minX-x1+1] > rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX+1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pVert[minX+1][y] = ISPDParser::Point(minX, y);
                rowCost[minX-x1+1] = rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX+1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }
    // std::cout<<std::endl;
    // std::cout<<"test"<<std::endl;
   //Bottom
    for(int y = s_y - 1; y >= y1; --y){
        for(int x = x1; x <= x2; ++x){
            pVert[x][y] = ISPDParser::Point(x, y+1);
        }
    }
    for(int y = s_y - 1; y >= y1; --y){
        std::vector<int> rowCost(x2 - x1 + 1, INT32_MAX);
        std::vector<bool> isVisit(x2 - x1 + 1, false); 
        //setting initial cost with verical edge  
        for(int x = x1; x <= x2; ++x){rowCost[x-x1] = (costVert[x][y+1] + OF(ISPDParser::Point(x, y), ISPDParser::Point(x, y+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));}
    
        for(int _ = 0; _ < x2 - x1 + 1; ++_){
            int minCost = INT32_MAX;
            int minX = -1;
            for(int x = x1; x <= x2; ++x){
                if(!isVisit[x-x1] && rowCost[x-x1] < minCost){
                    minCost = rowCost[x-x1];
                    minX = x;
                }
            }
            isVisit[minX-x1] = true;
            costVert[minX][y] = rowCost[minX-x1];
            //adjust left and right
            if(minX - 1 >= x1 && !isVisit[minX -x1 - 1] && rowCost[minX-x1-1] > rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX-1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pVert[minX-1][y] = ISPDParser::Point(minX, y);
                rowCost[minX-x1-1] = rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX-1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
            if(minX + 1 <= x2 && !isVisit[minX -x1 + 1] && rowCost[minX-x1+1] > rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX+1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pVert[minX+1][y] = ISPDParser::Point(minX, y);
                rowCost[minX-x1+1] = rowCost[minX-x1] + OF(ISPDParser::Point(minX, y), ISPDParser::Point(minX+1, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }

    /* hori monotonic route */
    std::vector<std::vector<int>> costHori(
        costMap.size(),
        std::vector<int>(costMap[0].size(), 0)
    );
    // std::vector<std::vector<ISPDParser::Point>> pHori(
    //     costMap.size(),
    //     std::vector<ISPDParser::Point>(costMap[0].size(), ISPDParser::Point(-1, -1))
    // );
    pHori[s_x][s_y] = ISPDParser::Point(s_x, s_y);
    for(int i = s_y - 1; i >= y1; --i){
        costHori[s_x][i] = (costHori[s_x][i+1] + OF(ISPDParser::Point(s_x, i), ISPDParser::Point(s_x, i+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));
        pHori[s_x][i] = ISPDParser::Point(s_x, i+1);
    }
    for(int i = s_y + 1; i <= y2; ++i){
        costHori[s_x][i] = (costHori[s_x][i-1] + OF(ISPDParser::Point(s_x, i), ISPDParser::Point(s_x, i-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));
        pHori[s_x][i] = ISPDParser::Point(s_x, i-1);
    }
    //Right
    for(int x = s_x + 1; x <= x2; ++x){
        for(int y = y1; y <= y2; ++y){
            pHori[x][y] = ISPDParser::Point(x-1, y);
        }
    }
    for(int x = s_x + 1; x <= x2; ++x){
        std::vector<int> colCost(y2 - y1 + 1, INT32_MAX);
        std::vector<bool> isVisit(y2 - y1 + 1, false); 
        //setting initial cost with verical edge  
        for(int y = y1; y <= y2; ++y){colCost[y-y1] = (costHori[x-1][y] + OF(ISPDParser::Point(x-1, y), ISPDParser::Point(x, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));}
    
        for(int _ = 0; _ < y2 - y1 + 1; ++_){
            int minCost = INT32_MAX;
            int minY = -1;
            for(int y = y1; y <= y2; ++y){
                if(!isVisit[y-y1] && colCost[y-y1] < minCost){
                    minCost = colCost[y-y1];
                    minY = y;
                }
            }
            isVisit[minY-y1] = true;
            costHori[x][minY] = colCost[minY-y1];
            //adjust top and bottom 
            if(minY - 1 >= y1 && !isVisit[minY -y1 - 1] && colCost[minY-y1-1] > colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pHori[x][minY-1] = ISPDParser::Point(x, minY);
                colCost[minY-y1-1] = colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
            if(minY + 1 <= y2 && !isVisit[minY -y1 + 1] && colCost[minY-y1+1] > colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pHori[x][minY+1] = ISPDParser::Point(x, minY);
                colCost[minY-y1+1] = colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }
   //Left
    for(int x = s_x - 1; x >= x1; --x){
        for(int y = y1; y <= y2; ++y){
            pHori[x][y] = ISPDParser::Point(x+1, y);
        }
    }
    for(int x = s_x - 1; x >= x1; --x){
        std::vector<int> colCost(y2 - y1 + 1, INT32_MAX);
        std::vector<bool> isVisit(y2 - y1 + 1, false); 
        //setting initial cost with verical edge  
        for(int y = y1; y <= y2; ++y){colCost[y-y1] = (costHori[x+1][y] + OF(ISPDParser::Point(x+1, y), ISPDParser::Point(x, y), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr));}
    
        for(int _ = 0; _ < y2 - y1 + 1; ++_){
            int minCost = INT32_MAX;
            int minY = -1;
            for(int y = y1; y <= y2; ++y){
                if(!isVisit[y-y1] && colCost[y-y1] < minCost){
                    minCost = colCost[y-y1];
                    minY = y;
                }
            }
            isVisit[minY-y1] = true;
            costHori[x][minY] = colCost[minY-y1];
            //adjust top and bottom 
            if(minY - 1 >= y1 && !isVisit[minY -y1 - 1] && colCost[minY-y1-1] > colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pHori[x][minY-1] = ISPDParser::Point(x, minY);
                colCost[minY-y1-1] = colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY-1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
            if(minY + 1 <= y2 && !isVisit[minY -y1 + 1] && colCost[minY-y1+1] > colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr)){
                pHori[x][minY+1] = ISPDParser::Point(x, minY);
                colCost[minY-y1+1] = colCost[minY-y1] + OF(ISPDParser::Point(x, minY), ISPDParser::Point(x, minY+1), minwidth, minSpace, vertCap, horiCap, vertCurr, horiCurr);
            }
        }
    }





    /* select path from hori or vert */
    std::vector<std::vector<std::vector<ISPDParser::Point>>> bestPath(
        costMap.size(),
        std::vector<std::vector<ISPDParser::Point>>( 
            costMap[0].size(), 
            std::vector<ISPDParser::Point>() 
        )
    );
    // for(int x =  52; x <= 56; ++x){
    //     for(int y = 5; y <= 15; ++y){
    //         std::cout<<pVert[x][y].x<<","<<pVert[x][y].y<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"--------------------"<<std::endl;

    for(int x = x1; x <= x2; ++x){
        for(int y = y1; y <= y2; ++y){
            std::vector<ISPDParser::Point> tPath;
            int _x = x, _y = y;
            if(costHori[x][y] < costVert[x][y]){
                costMap[x][y] =  costHori[x][y];
                while(pHori[_x][_y].x != _x || pHori[_x][_y].y != _y){
                    // std::cout<<s_x<<" "<<s_y<<std::endl;
                    // std::cout<<pHori[_x][_y].x<<" "<<pHori[_x][_y].y<<" "<<x<<" "<<std::endl;
                    tPath.push_back(ISPDParser::Point(_x, _y));
                    int t_x = _x, t_y = _y;
                    _x = pHori[t_x][t_y].x;
                    _y = pHori[t_x][t_y].y;
                }
            }
            else{
                costMap[x][y] =  costVert[x][y];
                while(pVert[_x][_y].x != _x || pVert[_x][_y].y != _y){
                    // std::cout<<pVert[_x][_y].x<<" "<<pVert[_x][_y].y<<" "<<x<<" "<<std::endl;
                    tPath.push_back(ISPDParser::Point(_x, _y));
                    int t_x = _x, t_y = _y;
                    _x = pVert[t_x][t_y].x;
                    _y = pVert[t_x][t_y].y;
                }
            }
            tPath.push_back(ISPDParser::Point(s_x, s_y));
            bestPath[x][y] = tPath;
        }
    }

    return bestPath;
}

int calOF(ISPDParser::Net* net, std::vector<ISPDParser::Point>& corr_path, std::vector<std::vector<int>>& vertCap, 
          std::vector<std::vector<int>>& horiCap, std::vector<std::vector<int>>& vertCurr, std::vector<std::vector<int>>& horiCurr){
    int totalOF = 0;
    if(corr_path.size() <= 1){return totalOF;}
    for(int i = 0; i < corr_path.size() - 1; ++i){
        int x1 = corr_path[i].x, y1 = corr_path[i].y;
        int x2 = corr_path[i+1].x, y2 = corr_path[i+1].y;
        totalOF += currOF(ISPDParser::Point(x1, y1), ISPDParser::Point(x2, y2), vertCap, horiCap, vertCurr, horiCurr);
    }
          
    return totalOF;
        
}
bool compareTwoPinHPWL(ISPDParser::TwoPin* a, ISPDParser::TwoPin* b) {
    return a->HPWL() < b->HPWL();
}

int main(int argc, char **argv) {
    assert(argc >= 3 && "Usage: ./router <inputFile> <outputFile>");
    std::ifstream fp(argv[1]);
    assert(fp.is_open() && "Failed to open input file");
    ISPDParser::ispdData *ispdData = ISPDParser::parse(fp);
    fp.close();

    // std::cout << *ispdData << std::endl;

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ispdData->nets.begin(), ispdData->nets.end(), [&](ISPDParser::Net *net) {

        for (auto &pin : net->pins) {

            int x = (std::get<0>(pin) - ispdData->lowerLeftX) / ispdData->tileWidth;
            int y = (std::get<1>(pin) - ispdData->lowerLeftY) / ispdData->tileHeight;
            int z = std::get<2>(pin) - 1;

            if (std::any_of(net->pin3D.begin(), net->pin3D.end(), [x, y, z](const auto &pin) {
                return pin.x == x && pin.y == y && pin.z == z;
            })) continue;
            net->pin3D.emplace_back(x, y, z);

            if (std::any_of(net->pin2D.begin(), net->pin2D.end(), [x, y](const auto &pin) { 
                return pin.x == x && pin.y == y;
            })) continue;
            net->pin2D.emplace_back(x, y);

        }

        return net->pin3D.size() > 1000 || net->pin2D.size() <= 1;

    }), ispdData->nets.end());
    ispdData->numNet = ispdData->nets.size();



    /* Constuct 2D graph*/
    int sumvertCap = 0, sumhoriCap = 0;
    for(int i = 0; i < ispdData->numLayer; ++i){sumvertCap += ispdData->verticalCapacity[i];sumhoriCap += ispdData->horizontalCapacity[i];}
    
    std::vector<std::vector<int>> vertCap(ispdData->numXGrid, std::vector<int>(ispdData->numYGrid - 1, sumvertCap));
    std::vector<std::vector<int>> horiCap(ispdData->numXGrid - 1, std::vector<int>(ispdData->numYGrid, sumhoriCap));

    std::vector<std::vector<int>> vertCurr(ispdData->numXGrid, std::vector<int>(ispdData->numYGrid - 1, 0));
    std::vector<std::vector<int>> horiCurr(ispdData->numXGrid - 1, std::vector<int>(ispdData->numYGrid, 0));

    for(auto& capAdj: ispdData->capacityAdjs){
        int x1 = std::get<0>(capAdj->grid1), y1 = std::get<1>(capAdj->grid1), z1 = std::get<2>(capAdj->grid1);
        int x2 = std::get<0>(capAdj->grid2), y2 = std::get<1>(capAdj->grid2), z2 = std::get<2>(capAdj->grid2);

        if(x1 + 1 == x2 && y1 == y2){
            horiCap[x1][y1] -= ispdData->horizontalCapacity[z1 - 1];
            horiCap[x1][y1] += capAdj->reducedCapacityLevel;
        }
        else if(x1 == x2 && y1 + 1 == y2){
            vertCap[x1][y1] -= ispdData->verticalCapacity[z1 - 1];
            vertCap[x1][y1] += capAdj->reducedCapacityLevel;
            
        }
        else{
            std::cerr<<"error in adjusting cap\n";
            exit(-1);
        }
    }


    /* decompsition, sort two pin net*/
    std::vector<ISPDParser::TwoPin*> twopins;
    for(auto& net: ispdData->nets){
        net->decompose();
        net->set_net(net);
        // net->sort_twopin();
        for(int i = 0; i < net->twopin.size(); ++i){
            twopins.push_back(&net->twopin[i]);
        }
    }


    std::sort(twopins.begin(), twopins.end(), compareTwoPinHPWL);

    /* pattern route */
    for(auto& ptrTwopin: twopins){
        ptrTwopin->corr_path = pattern_route((*ptrTwopin), ispdData->minimumSpacing[0], vertCap, horiCap,vertCurr, horiCurr);
        fill_cap(ptrTwopin->parNet, ptrTwopin->corr_path, ispdData->minimumSpacing[0], vertCurr, horiCurr);
    }
    // for(auto& net: ispdData->nets){
    //     for(int i = 0; i < net->twopin.size(); ++i){
    //         net->twopin[i].corr_path = pattern_route(net->twopin[i], ispdData->minimumSpacing[0], vertCap, horiCap,vertCurr, horiCurr);
    //         fill_cap(net, net->twopin[i].corr_path, ispdData->minimumSpacing[0], vertCurr, horiCurr);
    //     }
    // }


    /* HUM */
    twopins.clear();
    for(auto& net: ispdData->nets){
        for(int i = 0; i < net->twopin.size(); ++i){
            if(calOF(net, net->twopin[i].corr_path, vertCap, horiCap, vertCurr, horiCurr) > 0){
                twopins.push_back(&net->twopin[i]);
            }
        }
    }
    std::sort(twopins.begin(), twopins.end(), compareTwoPinHPWL);
    std::cout<<"iter "<<": "<<twopins.size()<<" pair twopins"<<std::endl;


    std::priority_queue<ISPDParser::twopinInfo> twopinHeap;
    for(auto& net: ispdData->nets){
        for(int i = 0; i < net->twopin.size(); ++i){
            int of = calOF(net, net->twopin[i].corr_path, vertCap, horiCap, vertCurr, horiCurr);
            twopinHeap.push(ISPDParser::twopinInfo(&net->twopin[i], of));
        }
    }


  
    int count = 0;
    // while (!twopinHeap.empty()) {
    for (int i = 0; i < twopins.size(); ++i) {
        ++count;
        if(count%100 == 0){std::cout<<count<<" twopin\n";}
        auto ptrTwopin = twopins[i];

        ISPDParser::TwoPin topTwopin = *(ptrTwopin);
        ISPDParser::Point s = topTwopin.from;
        ISPDParser::Point t = topTwopin.to;

        int x1 = std::min(s.x, t.x), x2 = std::max(s.x, t.x);
        int y1 = std::min(s.y, t.y), y2 = std::max(s.y, t.y);
        int boxSize = 20;
        x1 = std::max(x1 - boxSize, 0);
        y1 = std::max(y1 - boxSize, 0);
        x2 = std::min(x2 + boxSize, ispdData->numXGrid - 1);
        y2 = std::min(y2 + boxSize, ispdData->numYGrid - 1);
        std::vector<std::vector<int>> suCost(
            ispdData->numXGrid,
            std::vector<int>(ispdData->numYGrid, 0)
        );
        std::vector<std::vector<int>> utCost(
            ispdData->numXGrid,
            std::vector<int>(ispdData->numYGrid, 0)
        );
        std::vector<std::vector<ISPDParser::Point>> suVert(
            ispdData->numXGrid,
            std::vector<ISPDParser::Point>(ispdData->numYGrid, ISPDParser::Point(-1, -1))
        );
        std::vector<std::vector<ISPDParser::Point>> utVert(
            ispdData->numXGrid,
            std::vector<ISPDParser::Point>(ispdData->numYGrid, ISPDParser::Point(-1, -1))
        );
        std::vector<std::vector<ISPDParser::Point>> suHori(
            ispdData->numXGrid,
            std::vector<ISPDParser::Point>(ispdData->numYGrid, ISPDParser::Point(-1, -1))
        );
        std::vector<std::vector<ISPDParser::Point>> utHori(
            ispdData->numXGrid,
            std::vector<ISPDParser::Point>(ispdData->numYGrid, ISPDParser::Point(-1, -1))
        );
        // std::cout<<s.x<<","<<s.y<<"->"<<t.x<<","<<t.y<<" "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
        std::vector<std::vector<std::vector<ISPDParser::Point>>> s_u = unilateral_mono_route(x1, y1, x2, y2, s, suCost, topTwopin.parNet->minimumWidth, ispdData->minimumSpacing[0], vertCap, horiCap, vertCurr, horiCurr, suVert, suHori);
        std::vector<std::vector<std::vector<ISPDParser::Point>>> u_t = unilateral_mono_route(x1, y1, x2, y2, t, utCost, topTwopin.parNet->minimumWidth, ispdData->minimumSpacing[0], vertCap, horiCap, vertCurr, horiCurr, utVert, utHori);
        


        int minCost = INT32_MAX;
        ISPDParser::Point minU(-1, -1);
        for(int x = x1; x <= x2; ++x){
            for(int y = y1; y <= y2; ++y){
                if(suCost[x][y] + utCost[x][y] < minCost){
                    minCost = suCost[x][y] + utCost[x][y];
                    minU = ISPDParser::Point(x, y);
                }
            }  
        }  
        // std::cout<<minU.x<<" "<<minU.y<<std::endl;
        // std::cout<<s.x<<" "<<s.y<<" to "<<t.x<<" "<<t.y<<std::endl;
        // std::vector<ISPDParser::Point> suPath = s_u[minU.x][minU.y];
        // for(auto& p: s_u[minU.x][minU.y]){
        //     std::cout<<p.x<<" "<<p.y<<std::endl;
        // }
        std::reverse(s_u[minU.x][minU.y].begin(), s_u[minU.x][minU.y].end());
        s_u[minU.x][minU.y].pop_back();
        s_u[minU.x][minU.y].insert(s_u[minU.x][minU.y].end(), u_t[minU.x][minU.y].begin(), u_t[minU.x][minU.y].end());

        remove_cap(ptrTwopin->parNet, ptrTwopin->corr_path, ispdData->minimumSpacing[0], vertCurr, horiCurr);
        ptrTwopin->corr_path = s_u[minU.x][minU.y];
        fill_cap(ptrTwopin->parNet, ptrTwopin->corr_path, ispdData->minimumSpacing[0], vertCurr, horiCurr);
        // for(auto& p: s_u[minU.x][minU.y]){
        //     std::cout<<p.x<<" "<<p.y<<std::endl;
        // }
        // for(auto& p: u_t[minU.x][minU.y]){
        //     std::cout<<p.x<<" "<<p.y<<std::endl;
        // }
        // std::cout<<t.x<<" "<<t.y<<std::endl;
        // break;
    }

    /* construct RPoint path in each twopin*/
    for(auto& net: ispdData->nets){
        for(auto& pinpair: net->twopin){
            pinpair.create_path();
        }
    }
    LayerAssignment::Graph graph;
    graph.initialLA(*ispdData, 1);
    graph.convertGRtoLA(*ispdData, true);
    graph.COLA(true);

    // Output result
    graph.output3Dresult("3ds1.txt");
  


    // Describe the usage of the given layer assignment algorithm
    // Only works for the given input file "3d.txt"
    if (std::string(argv[1]).find("3d.txt") != std::string::npos) {

        ISPDParser::Net *net = ispdData->nets[0];

        // Decompose multi-pin nets into two-pin nets
        // Since there are only 2 pins in the given net, this step is trivial
        net->twopin.push_back(ISPDParser::TwoPin());
        ISPDParser::TwoPin &twoPin = net->twopin.back();
        twoPin.from = net->pin3D[0];
        twoPin.to   = net->pin3D[1];


        

        // Assume the two pin net is routed
        // The following code is to assign routing paths to the two-pin net
        // The routing path is a sequence of routing segments
        // For a horizontal segment, the start point is the left grid coordinate
        // For a vertical segment, the start point is the bottom grid coordinate
        // Please check the figures in https://www.ispd.cc/contests/08/3d.pdf
        twoPin.parNet = net;
        twoPin.path.emplace_back(0, 0, true);
        twoPin.path.emplace_back(1, 0, false);
        twoPin.path.emplace_back(0, 1, true);
        twoPin.path.emplace_back(0, 1, false);
        twoPin.path.emplace_back(0, 2, true);
        twoPin.path.emplace_back(1, 2, true);
        twoPin.path.emplace_back(2, 1, false);
        twoPin.path.emplace_back(2, 0, false);

        // Assign routing layers to the two-pin net
        LayerAssignment::Graph graph;
        graph.initialLA(*ispdData, 1);
        graph.convertGRtoLA(*ispdData, true);
        graph.COLA(true);

        // Output result
        graph.output3Dresult("3ds1.txt");
    }

    delete ispdData;
    return 0;
}
//./router Benchmarks/adaptec1.gr test
//./router 3d.txt test