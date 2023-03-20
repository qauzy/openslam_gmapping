#ifndef HARRAY2D_H
#define HARRAY2D_H
#include <set>
#include <gmapping/utils/point.h>
#include <gmapping/utils/autoptr.h>
#include "gmapping/grid/array2d.h"

namespace GMapping {
//    在HierarchicalArray2D中把地图分为若干个patch，其拆分的粒度由m_patchMagnitude来描述。
//    在cell函数中首先通过成员函数patchIndexes计算给定的坐标所属的patch的索引。 在patchIndexes
//    函数中将给定的坐标值右移了m_patchMagnitude位，这相当于除以了2n。在计算机中除法的计算是很慢的，
//    所以这里通过移位的方式替代除法操作， 但是这一技巧只能用于整型的数据，对于浮点数不适用。

//    Array2D是一个二维数组，HierarchicalArray2D是一个Array2D的二维数组。 相当于将地图先分割成
//    分辨率比较低的网格，只有当粒子运动到该网格时，才真正分配这个网格的内存。 网格的内存对应着分辨率
//    高的真实的地图，用PointAccumulator来进行计数，记录激光是否通过该点，从而判断改点的状态：占据、
//    空闲、未知。
//    -1 代表栅格状态未知
//    0 代表栅格是空闲的，代表可通过区域
//    100 代表栅格是占用的，代表障碍物，不可通过
template <class Cell>
class HierarchicalArray2D: public Array2D<autoptr< Array2D<Cell> > >{
	public:
		typedef std::set< point<int>, pointcomparator<int> > PointSet;
		HierarchicalArray2D(int xsize, int ysize, int patchMagnitude=5);
		HierarchicalArray2D(const HierarchicalArray2D& hg);
		HierarchicalArray2D& operator=(const HierarchicalArray2D& hg);
		virtual ~HierarchicalArray2D(){}
		void resize(int ixmin, int iymin, int ixmax, int iymax);
		inline int getPatchSize() const {return m_patchMagnitude;}
		inline int getPatchMagnitude() const {return m_patchMagnitude;}
		
		inline const Cell& cell(int x, int y) const;
		inline Cell& cell(int x, int y);
		inline bool isAllocated(int x, int y) const;
		inline AccessibilityState cellState(int x, int y) const ;
		inline IntPoint patchIndexes(int x, int y) const;
		
		inline const Cell& cell(const IntPoint& p) const { return cell(p.x,p.y); }
		inline Cell& cell(const IntPoint& p) { return cell(p.x,p.y); }
		inline bool isAllocated(const IntPoint& p) const { return isAllocated(p.x,p.y);}
		inline AccessibilityState cellState(const IntPoint& p) const { return cellState(p.x,p.y); }
		inline IntPoint patchIndexes(const IntPoint& p) const { return patchIndexes(p.x,p.y);}
		
		inline void setActiveArea(const PointSet&, bool patchCoords=false);
		const PointSet& getActiveArea() const {return m_activeArea; }
		inline void allocActiveArea();
	protected:
		virtual Array2D<Cell> * createPatch(const IntPoint& p) const;
		PointSet m_activeArea;
		int m_patchMagnitude;
		int m_patchSize;
};

template <class Cell>
HierarchicalArray2D<Cell>::HierarchicalArray2D(int xsize, int ysize, int patchMagnitude) 
  :Array2D<autoptr< Array2D<Cell> > >::Array2D((xsize>>patchMagnitude), (ysize>>patchMagnitude)){
	m_patchMagnitude=patchMagnitude;
	m_patchSize=1<<m_patchMagnitude;
}

template <class Cell>
HierarchicalArray2D<Cell>::HierarchicalArray2D(const HierarchicalArray2D& hg)
  :Array2D<autoptr< Array2D<Cell> > >::Array2D((hg.m_xsize>>hg.m_patchMagnitude), (hg.m_ysize>>hg.m_patchMagnitude))  // added by cyrill: if you have a resize error, check this again
{
	this->m_xsize=hg.m_xsize;
	this->m_ysize=hg.m_ysize;
	this->m_cells=new autoptr< Array2D<Cell> >*[this->m_xsize];
	for (int x=0; x<this->m_xsize; x++){
		this->m_cells[x]=new autoptr< Array2D<Cell> >[this->m_ysize];
		for (int y=0; y<this->m_ysize; y++)
			this->m_cells[x][y]=hg.m_cells[x][y];
	}
	this->m_patchMagnitude=hg.m_patchMagnitude;
	this->m_patchSize=hg.m_patchSize;
}

template <class Cell>
void HierarchicalArray2D<Cell>::resize(int xmin, int ymin, int xmax, int ymax){
	int xsize=xmax-xmin;
	int ysize=ymax-ymin;
	autoptr< Array2D<Cell> > ** newcells=new autoptr< Array2D<Cell> > *[xsize];
	for (int x=0; x<xsize; x++){
		newcells[x]=new autoptr< Array2D<Cell> >[ysize];
		for (int y=0; y<ysize; y++){
			newcells[x][y]=autoptr< Array2D<Cell> >(0);
		}
	}
	int dx= xmin < 0 ? 0 : xmin;
	int dy= ymin < 0 ? 0 : ymin;
	int Dx=xmax<this->m_xsize?xmax:this->m_xsize;
	int Dy=ymax<this->m_ysize?ymax:this->m_ysize;
	for (int x=dx; x<Dx; x++){
		for (int y=dy; y<Dy; y++){
			newcells[x-xmin][y-ymin]=this->m_cells[x][y];
		}
		delete [] this->m_cells[x];
	}
	delete [] this->m_cells;
	this->m_cells=newcells;
	this->m_xsize=xsize;
	this->m_ysize=ysize; 
}

template <class Cell>
HierarchicalArray2D<Cell>& HierarchicalArray2D<Cell>::operator=(const HierarchicalArray2D& hg){
//	Array2D<autoptr< Array2D<Cell> > >::operator=(hg);
	if (this->m_xsize!=hg.m_xsize || this->m_ysize!=hg.m_ysize){
		for (int i=0; i<this->m_xsize; i++)
			delete [] this->m_cells[i];
		delete [] this->m_cells;
		this->m_xsize=hg.m_xsize;
		this->m_ysize=hg.m_ysize;
		this->m_cells=new autoptr< Array2D<Cell> >*[this->m_xsize];
		for (int i=0; i<this->m_xsize; i++)
			this->m_cells[i]=new autoptr< Array2D<Cell> > [this->m_ysize];
	}
	for (int x=0; x<this->m_xsize; x++)
		for (int y=0; y<this->m_ysize; y++)
			this->m_cells[x][y]=hg.m_cells[x][y];
	
	m_activeArea.clear();
	m_patchMagnitude=hg.m_patchMagnitude;
	m_patchSize=hg.m_patchSize;
	return *this;
}


template <class Cell>
void HierarchicalArray2D<Cell>::setActiveArea(const typename HierarchicalArray2D<Cell>::PointSet& aa, bool patchCoords){
	m_activeArea.clear();
	for (PointSet::const_iterator it= aa.begin(); it!=aa.end(); it++){
		IntPoint p;
		if (patchCoords)
			p=*it;
		else
			p=patchIndexes(*it);
		m_activeArea.insert(p);
	}
}

template <class Cell>
Array2D<Cell>* HierarchicalArray2D<Cell>::createPatch(const IntPoint& ) const{
	return new Array2D<Cell>(1<<m_patchMagnitude, 1<<m_patchMagnitude);
}


template <class Cell>
AccessibilityState  HierarchicalArray2D<Cell>::cellState(int x, int y) const {
	if (this->isInside(patchIndexes(x,y)))
		if(isAllocated(x,y))
			return (AccessibilityState)((int)Inside|(int)Allocated);
		else
			return Inside;
	return Outside;
}

template <class Cell>
void HierarchicalArray2D<Cell>::allocActiveArea(){
	for (PointSet::const_iterator it= m_activeArea.begin(); it!=m_activeArea.end(); it++){
		const autoptr< Array2D<Cell> >& ptr=this->m_cells[it->x][it->y];
		Array2D<Cell>* patch=0;
		if (!ptr){
			patch=createPatch(*it);
		} else{	
			patch=new Array2D<Cell>(*ptr);
		}
		this->m_cells[it->x][it->y]=autoptr< Array2D<Cell> >(patch);
	}
}

template <class Cell>
bool HierarchicalArray2D<Cell>::isAllocated(int x, int y) const{
	IntPoint c=patchIndexes(x,y);
	autoptr< Array2D<Cell> >& ptr=this->m_cells[c.x][c.y];
	return (ptr != 0);
}

template <class Cell>
IntPoint HierarchicalArray2D<Cell>::patchIndexes(int x, int y) const{
	if (x>=0 && y>=0)
		return IntPoint(x>>m_patchMagnitude, y>>m_patchMagnitude);
	return IntPoint(-1, -1);
}

template <class Cell>
Cell& HierarchicalArray2D<Cell>::cell(int x, int y){
	IntPoint c=patchIndexes(x,y);
	assert(this->isInside(c.x, c.y));
	if (!this->m_cells[c.x][c.y]){
		Array2D<Cell>* patch=createPatch(IntPoint(x,y));
		this->m_cells[c.x][c.y]=autoptr< Array2D<Cell> >(patch);
		//cerr << "!!! FATAL: your dick is going to fall down" << endl;
	}
	autoptr< Array2D<Cell> >& ptr=this->m_cells[c.x][c.y];
	return (*ptr).cell(IntPoint(x-(c.x<<m_patchMagnitude),y-(c.y<<m_patchMagnitude)));
}

template <class Cell>
const Cell& HierarchicalArray2D<Cell>::cell(int x, int y) const{
	assert(isAllocated(x,y));
	IntPoint c=patchIndexes(x,y);
	const autoptr< Array2D<Cell> >& ptr=this->m_cells[c.x][c.y];
	return (*ptr).cell(IntPoint(x-(c.x<<m_patchMagnitude),y-(c.y<<m_patchMagnitude)));
}

};

#endif

