#include "gmapping/gridfastslam/motionmodel.h"
#include <gmapping/utils/stat.h>
#include <iostream>

#define MotionModelConditioningLinearCovariance 0.01
#define MotionModelConditioningAngularCovariance 0.001

namespace GMapping {



OrientedPoint 
MotionModel::drawFromMotion (const OrientedPoint& p, double linearMove, double angularMove) const{
	OrientedPoint n(p);
	double lm=linearMove  + fabs( linearMove ) * sampleGaussian( srr ) + fabs( angularMove ) * sampleGaussian( str );
	double am=angularMove + fabs( linearMove ) * sampleGaussian( srt ) + fabs( angularMove ) * sampleGaussian( stt );
	n.x+=lm*cos(n.theta+.5*am);
	n.y+=lm*sin(n.theta+.5*am);
	n.theta+=am;
	n.theta=atan2(sin(n.theta), cos(n.theta));
	return n;
}
//drawFromMotion()函数将所有粒子在odom坐标的基础上加入高斯白噪声，即所谓的放狗过程
//p是上一时刻粒子的位姿，pnew则是里程计记录的最新位姿，而pold则是建图引擎记录的上一时刻的位姿。
//遍历每个粒子，调用该函数之后，就完成了对各个粒子的位姿预估。
OrientedPoint 
MotionModel::drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const{
	double sxy=0.3*srr;
	OrientedPoint delta=absoluteDifference(pnew, pold);// 计算位姿增量
	OrientedPoint noisypoint(delta);// 里程计增量高斯噪声
        // 该处使用基于里程计的采样算法，与书上公式不相同，是工程化后的公式。

        //sample函數是形參作爲方差，均值爲0的高斯分佈。sample函數是數值分析所近似生成的高斯分佈，
        // 具體的函數實現在stat.cpp中。當前位姿與上一次位姿做差，計算做累計角度偏差和位移偏差。利
        // 用激光雷達測得距離做得分處理
	noisypoint.x+=sampleGaussian(srr*fabs(delta.x)+str*fabs(delta.theta)+sxy*fabs(delta.y));
	noisypoint.y+=sampleGaussian(srr*fabs(delta.y)+str*fabs(delta.theta)+sxy*fabs(delta.x));
	noisypoint.theta+=sampleGaussian(stt*fabs(delta.theta)+srt*sqrt(delta.x*delta.x+delta.y*delta.y));
        // 角度归一计算
	noisypoint.theta=fmod(noisypoint.theta, 2*M_PI);
	if (noisypoint.theta>M_PI)//超过180度 // 角度归一计算
		noisypoint.theta-=2*M_PI;//取-180到0度
        // 上次估算位姿加上添加噪声后的里程计增量，输出估计后的位姿。
	return absoluteSum(p,noisypoint);
}


/*
OrientedPoint 
MotionModel::drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const{
	
	//compute the three stps needed for perfectly matching the two poses if the noise is absent
	
	OrientedPoint delta=pnew-pold;
	double aoffset=atan2(delta.y, delta.x);
	double alpha1=aoffset-pold.theta;
	alpha1=atan2(sin(alpha1), cos(alpha1));
	double rho=sqrt(delta*delta);
	double alpha2=pnew.theta-aoffset;
	alpha2=atan2(sin(alpha2), cos(alpha2));
	
	OrientedPoint pret=drawFromMotion(p, 0, alpha1);
	pret=drawFromMotion(pret, rho, 0);
	pret=drawFromMotion(pret, 0, alpha2);
	return pret;
}
*/


Covariance3 MotionModel::gaussianApproximation(const OrientedPoint& pnew, const OrientedPoint& pold) const{
	OrientedPoint delta=absoluteDifference(pnew,pold);
	double linearMove=sqrt(delta.x*delta.x+delta.y*delta.y);
	double angularMove=fabs(delta.x);
	double s11=srr*srr*linearMove*linearMove;
	double s22=stt*stt*angularMove*angularMove;
	double s12=str*angularMove*srt*linearMove;
	Covariance3 cov;
	double s=sin(pold.theta),c=cos(pold.theta);
	cov.xx=c*c*s11+MotionModelConditioningLinearCovariance;
	cov.yy=s*s*s11+MotionModelConditioningLinearCovariance;
	cov.tt=s22+MotionModelConditioningAngularCovariance;
	cov.xy=s*c*s11;
	cov.xt=c*s12;
	cov.yt=s*s12;
	return cov;
}

};

