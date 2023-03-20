#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>
#include <gmapping/sensor/sensor_range/sensor_range_export.h>

namespace GMapping{
// 是激光传感器的封装， 它描述了扫描光束的的各种物理特性，继承自类Sensor。
class SENSOR_RANGE_EXPORT RangeSensor: public Sensor{
	friend class Configuration;
	friend class CarmenConfiguration;
	friend class CarmenWrapper;
	public:
	// 在RangeSensor内部嵌套定义了一个结构体Beam，描述了扫描光束的物理特性，包括相对于传感器坐标系的位置、量程、光束角的正余弦值。
		struct Beam{
			OrientedPoint pose;	//pose relative to the center of the sensor
			double span;	//spam=0 indicates a line-like beam
			double maxRange;	//maximum range of the sensor
			double s,c;		//sinus and cosinus of the beam (optimization);
		};	
		RangeSensor(std::string name);
		RangeSensor(std::string name, unsigned int beams, double res, const OrientedPoint& position=OrientedPoint(0,0,0), double span=0, double maxrange=89.0);
		inline const std::vector<Beam>& beams() const {return m_beams;}
		inline std::vector<Beam>& beams() {return m_beams;}
		inline OrientedPoint getPose() const {return m_pose;}
		// 用于更新扫描光束属性
		void updateBeamsLookup();
		bool newFormat;
	protected:
	// RangeSensor定义了两个成员变量，其中m_pose用于记录传感器的位姿，m_beams则记录了扫描光束的信息。同时提供了set-get函数用于访问这些成员变量。
		OrientedPoint m_pose;
		std::vector<Beam> m_beams;
};

};

#endif
