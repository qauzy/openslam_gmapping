#ifndef SMMAP_H
#define SMMAP_H
#include <gmapping/grid/map.h>
#include <gmapping/grid/harray2d.h>
#include <gmapping/utils/point.h>
#define SIGHT_INC 1

namespace GMapping {
//    地图单元类型为PointAccumulator
struct PointAccumulator{
	typedef point<float> FloatPoint;
	/* before 
	PointAccumulator(int i=-1): acc(0,0), n(0), visits(0){assert(i==-1);}
	*/
	/*after begin*/
	PointAccumulator(): acc(0,0), n(0), visits(0){}
	PointAccumulator(int i): acc(0,0), n(0), visits(0){assert(i==-1);}
	/*after end*/
        inline void update(bool value, const Point& p=Point(0,0));
	inline Point mean() const {return 1./n*Point(acc.x, acc.y);}
	inline operator double() const { return visits?(double)n*SIGHT_INC/(double)visits:-1; }
	inline void add(const PointAccumulator& p) {acc=acc+p.acc; n+=p.n; visits+=p.visits; }
	static const PointAccumulator& Unknown();
	static PointAccumulator* unknown_ptr;
//    成员变量acc用于累积点坐标
	FloatPoint acc;
//    变量n和visits分别用于记录累积的次数和访问的次数
	int n, visits;
	inline double entropy() const;
};
//    存储结构为HierachicalArray2D<PointAccumulator>
void PointAccumulator::update(bool value, const Point& p){
        //如果点的访问次数大于0,则更新点的访问次数
	if (value) {
		acc.x+= static_cast<float>(p.x); //更新点的x坐标
		acc.y+= static_cast<float>(p.y); //更新点的y坐标
		n++;                             //更新点的访问次数
		visits+=SIGHT_INC;               //更新点的访问次数的总和  SIGHT_INC是一个访问次数,默认值为1
	} else
		visits++;                       //更新点的访问次数的总和  SIGHT_INC是一个访问次数,默认值为1
}

double PointAccumulator::entropy() const{
	if (!visits)
		return -log(.5);
	if (n==visits || n==0)
		return 0;
	double x=(double)n*SIGHT_INC/(double)visits;
	return -( x*log(x)+ (1-x)*log(1-x) );
}


typedef Map<PointAccumulator,HierarchicalArray2D<PointAccumulator> > ScanMatcherMap;

};

#endif 
