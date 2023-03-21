#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/gridslamprocessor.h"

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping {

    const double m_distanceThresholdCheck = 20;

    using namespace std;

    GridSlamProcessor::GridSlamProcessor() : m_infoStream(cout) {

        period_ = 5.0;
        m_obsSigmaGain = 1;
        m_resampleThreshold = 0.5;
        m_minimumScore = 0.;
    }

    GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor &gsp)
            : last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout) {

        period_ = 5.0;

        m_obsSigmaGain = gsp.m_obsSigmaGain;
        m_resampleThreshold = gsp.m_resampleThreshold;
        m_minimumScore = gsp.m_minimumScore;

        m_beams = gsp.m_beams;
        m_indexes = gsp.m_indexes;
        m_motionModel = gsp.m_motionModel;
        m_resampleThreshold = gsp.m_resampleThreshold;
        m_matcher = gsp.m_matcher;

        m_count = gsp.m_count;
        m_readingCount = gsp.m_readingCount;
        m_lastPartPose = gsp.m_lastPartPose;
        m_pose = gsp.m_pose;
        m_odoPose = gsp.m_odoPose;
        m_linearDistance = gsp.m_linearDistance;
        m_angularDistance = gsp.m_angularDistance;
        m_neff = gsp.m_neff;

        cerr << "FILTER COPY CONSTRUCTOR" << endl;
        cerr << "m_odoPose=" << m_odoPose.x << " " << m_odoPose.y << " " << m_odoPose.theta << endl;
        cerr << "m_lastPartPose=" << m_lastPartPose.x << " " << m_lastPartPose.y << " " << m_lastPartPose.theta << endl;
        cerr << "m_linearDistance=" << m_linearDistance << endl;
        cerr << "m_angularDistance=" << m_linearDistance << endl;


        m_xmin = gsp.m_xmin;
        m_ymin = gsp.m_ymin;
        m_xmax = gsp.m_xmax;
        m_ymax = gsp.m_ymax;
        m_delta = gsp.m_delta;

        m_regScore = gsp.m_regScore;
        m_critScore = gsp.m_critScore;
        m_maxMove = gsp.m_maxMove;

        m_linearThresholdDistance = gsp.m_linearThresholdDistance;
        m_angularThresholdDistance = gsp.m_angularThresholdDistance;
        m_obsSigmaGain = gsp.m_obsSigmaGain;

#ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ <<  ": trajectories copy.... ";
#endif
        TNodeVector v = gsp.getTrajectories();
        for (unsigned int i = 0; i < v.size(); i++) {
            m_particles[i].node = v[i];
        }
#ifdef MAP_CONSISTENCY_CHECK
        cerr <<  "end" << endl;
#endif


        cerr << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ...";
        updateTreeWeights(false);
        cerr << ".done!" << endl;
    }

    GridSlamProcessor::GridSlamProcessor(std::ostream &infoS) : m_infoStream(infoS) {
        period_ = 5.0;
        m_obsSigmaGain = 1;
        m_resampleThreshold = 0.5;
        m_minimumScore = 0.;

    }

    GridSlamProcessor *GridSlamProcessor::clone() const {
# ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": performing preclone_fit_test" << endl;
        typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
        PointerMap pmap;
        for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
          const ScanMatcherMap& m1(it->map);
          const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
           for (int x=0; x<h1.getXSize(); x++){
            for (int y=0; y<h1.getYSize(); y++){
              const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
              if (a1.m_reference){
            PointerMap::iterator f=pmap.find(a1.m_reference);
            if (f==pmap.end())
              pmap.insert(make_pair(a1.m_reference, 1));
            else
              f->second++;
              }
            }
          }
        }
        cerr << __func__ <<  ": Number of allocated chunks" << pmap.size() << endl;
        for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
          assert(it->first->shares==(unsigned int)it->second);

        cerr << __func__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
        GridSlamProcessor *cloned = new GridSlamProcessor(*this);

# ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ <<  ": trajectories end" << endl;
        cerr << __func__ << ": performing afterclone_fit_test" << endl;
        ParticleVector::const_iterator jt=cloned->m_particles.begin();
        for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
          const ScanMatcherMap& m1(it->map);
          const ScanMatcherMap& m2(jt->map);
          const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
          const HierarchicalArray2D<PointAccumulator>& h2(m2.storage());
          jt++;
           for (int x=0; x<h1.getXSize(); x++){
            for (int y=0; y<h1.getYSize(); y++){
              const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
              const autoptr< Array2D<PointAccumulator> >& a2(h2.m_cells[x][y]);
              assert(a1.m_reference==a2.m_reference);
              assert((!a1.m_reference) || !(a1.m_reference->shares%2));
            }
          }
        }
        cerr << __func__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
        return cloned;
    }

    GridSlamProcessor::~GridSlamProcessor() {
        cerr << __func__ << ": Start" << endl;
        cerr << __func__ << ": Deleting tree" << endl;
        for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
#ifdef TREE_CONSISTENCY_CHECK
            TNode* node=it->node;
            while(node)
          node=node->parent;
            cerr << "@" << endl;
#endif
            if (it->node)
                delete it->node;
            //cout << "l=" << it->weight<< endl;
        }

# ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": performing predestruction_fit_test" << endl;
        typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
        PointerMap pmap;
        for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
          const ScanMatcherMap& m1(it->map);
          const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
          for (int x=0; x<h1.getXSize(); x++){
        for (int y=0; y<h1.getYSize(); y++){
          const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
          if (a1.m_reference){
            PointerMap::iterator f=pmap.find(a1.m_reference);
            if (f==pmap.end())
              pmap.insert(make_pair(a1.m_reference, 1));
            else
              f->second++;
          }
        }
          }
        }
        cerr << __func__ << ": Number of allocated chunks" << pmap.size() << endl;
        for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
          assert(it->first->shares>=(unsigned int)it->second);
        cerr << __func__ << ": SUCCESS, the error is somewhere else" << endl;
# endif
    }


    void GridSlamProcessor::setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt,
                                                  double aopt,
                                                  int iterations, double likelihoodSigma, double likelihoodGain,
                                                  unsigned int likelihoodSkip) {
        m_obsSigmaGain = likelihoodGain;
        m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma,
                                        likelihoodSkip);
        if (m_infoStream)
            m_infoStream << " -maxUrange " << urange
                         << " -maxUrange " << range
                         << " -sigma     " << sigma
                         << " -kernelSize " << kernsize
                         << " -lstep " << lopt
                         << " -lobsGain " << m_obsSigmaGain
                         << " -astep " << aopt << endl;


    }

    void GridSlamProcessor::setMotionModelParameters
            (double srr, double srt, double str, double stt) {
        m_motionModel.srr = srr;
        m_motionModel.srt = srt;
        m_motionModel.str = str;
        m_motionModel.stt = stt;

        if (m_infoStream)
            m_infoStream << " -srr " << srr << " -srt " << srt
                         << " -str " << str << " -stt " << stt << endl;

    }

    void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold) {
        m_linearThresholdDistance = linear;
        m_angularThresholdDistance = angular;
        m_resampleThreshold = resampleThreshold;
        if (m_infoStream)
            m_infoStream << " -linearUpdate " << linear
                         << " -angularUpdate " << angular
                         << " -resampleThreshold " << m_resampleThreshold << endl;
    }

    //HERE STARTS THE BEEF

    GridSlamProcessor::Particle::Particle(const ScanMatcherMap &m) :
            map(m), pose(0, 0, 0), weight(0), weightSum(0), gweight(0), previousIndex(0) {
        node = 0;
    }


    void GridSlamProcessor::setSensorMap(const SensorMap &smap) {

        /*
          Construct the angle table for the sensor

          FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot
        */
//        目前GMapping只支持一种"FLASER"的雷达，它假设雷达对机器人前面进行扫描，并且安装在机器人的中心(For now detect the readings of only the front laser,
//        and assume its pose is in the center of the robot)。下面的代码片段用于获取传感器"FLASER"，如果map中没有相应的传感器就报错。
        SensorMap::const_iterator laser_it = smap.find(std::string("FLASER"));
        if (laser_it == smap.end()) {
            cerr << "Attempting to load the new carmen log format" << endl;
            laser_it = smap.find(std::string("ROBOTLASER1"));
            assert(laser_it != smap.end());
        }
        const RangeSensor *rangeSensor = dynamic_cast<const RangeSensor *>((laser_it->second));
        assert(rangeSensor && rangeSensor->beams().size());

//        然后我们获取传感器的波束数量用m_beams记录，并用angles记录各个波束对应的角度。最后根据波束数量、波束角度和传感器位置配置扫描匹配器的参数。
        m_beams = static_cast<unsigned int>(rangeSensor->beams().size());
        double *angles = new double[rangeSensor->beams().size()];
        for (unsigned int i = 0; i < m_beams; i++) {
            angles[i] = rangeSensor->beams()[i].pose.theta;
        }
        m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
        delete[] angles;
    }
//    在ROS封装的初始化的最后阶段，调用了建图引擎的init函数完成了整个系统的初始化操作。如下面的代码片段所示，该函数有7个参数。其中size表示粒子
//    数量， xmin、ymin、xmax、ymax一起描述了地图尺寸，它们分别描述了地图中x和y的最小值和最大值。delta则是地图的分辨率或者说比例尺，即一个栅格
//    所对应的实际地图尺寸。 initialPose则是机器人的起始位置。
    void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta,
                                 OrientedPoint initialPose) {
        m_xmin = xmin;
        m_ymin = ymin;
        m_xmax = xmax;
        m_ymax = ymax;
        m_delta = delta;
        if (m_infoStream)
            m_infoStream
                    << " -xmin " << m_xmin
                    << " -xmax " << m_xmax
                    << " -ymin " << m_ymin
                    << " -ymax " << m_ymax
                    << " -delta " << m_delta
                    << " -particles " << size << endl;


        m_particles.clear();
        TNode *node = new TNode(initialPose, 0, 0, 0);
        ScanMatcherMap lmap(Point(xmin + xmax, ymin + ymax) * .5, xmax - xmin, ymax - ymin, delta);
        for (unsigned int i = 0; i < size; i++) {
            m_particles.push_back(Particle(lmap));
            m_particles.back().pose = initialPose;
            m_particles.back().previousPose = initialPose;
            m_particles.back().setWeight(0);
            m_particles.back().previousIndex = 0;

            // this is not needed
            //		m_particles.back().node=new TNode(initialPose, 0, node, 0);

            // we use the root directly
            m_particles.back().node = node;
        }
        m_neff = (double) size;
        m_count = 0;
        m_readingCount = 0;
        m_linearDistance = m_angularDistance = 0;
    }

    void GridSlamProcessor::processTruePos(const OdometryReading &o) {
        const OdometrySensor *os = dynamic_cast<const OdometrySensor *>(o.getSensor());
        if (os && os->isIdeal() && m_outputStream) {
            m_outputStream << setiosflags(ios::fixed) << setprecision(3);
            m_outputStream << "SIMULATOR_POS " << o.getPose().x << " " << o.getPose().y << " ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " << o.getTime()
                           << endl;
        }
    }
//    具体来说，"processScan"函数的实现流程如下：
//
//    获取激光雷达的扫描数据，并将其转换为一系列测量点的坐标。
//    使用"scanmatcher"算法对机器人当前的位置进行估计。
//    使用得到的机器人位置和测量点的坐标来更新地图信息。
//    根据需要，调整地图的分辨率和大小等参数。
//    返回更新后的地图信息。

    // processScan函数功能：通过激光数据来辅助机器人定位，进而获取更准确的机器人位置坐标；当获取到机器人的位置坐标后，我们就可以根据机器人的位置坐标及激光雷达数据更新地图了。
    bool GridSlamProcessor::processScan(const RangeReading &reading, int adaptParticles) {

        /**retireve the position from the reading, and compute the odometry*/
        OrientedPoint relPose = reading.getPose();//获取激光数据的位姿
        //m_count是一个整数，用于记录运行的次数
        //如果m_count为0,则执行下面的代码   意思是一开始没有位姿信息，就把第一次获得的激光数据的位姿作为起始位姿，
        //然后后面以此为起始点，不断处理数据，不断更新位姿   这个if只执行一次；
        if (!m_count) {
            m_lastPartPose = m_odoPose = relPose; //将获取到的激光数据的位姿赋值给m_odoPose
        }

        //write the state of the reading and update all the particles using the motion model
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
            OrientedPoint &pose(it->pose); //获取粒子的位姿

            //drawFromMotion()函数将所有粒子在odom坐标的基础上加入高斯白噪声，即所谓的放狗过程
            //m_particles[i].pose是上一时刻粒子的位姿，relPose则是里程计(应该是激光雷达)记录的最新位姿，而m_odoPose则是建图引擎记录的上一时刻的位姿。
            //遍历每个粒子，调用该函数之后，就完成了对各个粒子的位姿预估。

            //这里的MotionModel是一个十分粗糙的运动模型，只是简单的矢量加减运算。 相比于《Probabilistic Robotics》中提到的速度模型和里程计模型而言，
            // 有很多方面都没有考虑，精度上可能有折扣。 但是计算简单，计算量小。有时为了提高系统的实时性，这种以牺牲精度来换取计算效率的方式是很常用的手段。
            // 模型不必十分精确，够用就行。
            pose = m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose); //------->得到粒子新的位姿
        }

        // update the output file
        if (m_outputStream.is_open()) {
            m_outputStream << setiosflags(ios::fixed) << setprecision(6);
            m_outputStream << "ODOM ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(3) << m_odoPose.x << " " << m_odoPose.y << " ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(6) << m_odoPose.theta << " ";
            m_outputStream << reading.getTime();
            m_outputStream << endl;
        }
        if (m_outputStream.is_open()) {
            m_outputStream << setiosflags(ios::fixed) << setprecision(6);
            m_outputStream << "ODO_UPDATE " << m_particles.size() << " ";
            for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
                OrientedPoint &pose(it->pose);
                m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
                m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it->weight << " ";
            }
            m_outputStream << reading.getTime();
            m_outputStream << endl;
        }

        //invoke the callback
        onOdometryUpdate();


        // accumulate the robot translation and rotation
        //计算机器人的平移和旋转量
        OrientedPoint move = relPose - m_odoPose; //当前获取的激光位姿位姿-上一时刻位姿

        move.theta = atan2(sin(move.theta), cos(move.theta));

        //sqrt()函数是求平方根，即求平移量  此处*是点的重载运算符，即(move.x*move.x+move.y*move.y),括号里的*是正常的乘法，别搞混了
        m_linearDistance += sqrt(move * move);
        //fabs()函数是求绝对值，即求旋转量
        m_angularDistance += fabs(move.theta);

        // if the robot jumps throw a warning
        if (m_linearDistance > m_distanceThresholdCheck) {
            cerr << "***********************************************************************" << endl;
            cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
            cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
            cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y
                 << " " << m_odoPose.theta << endl;
            cerr << "New Odometry Pose (reported from observation)= " << relPose.x << " " << relPose.y
                 << " " << relPose.theta << endl;
            cerr << "***********************************************************************" << endl;
            cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
            cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
            cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
            cerr << "***********************************************************************" << endl;
        }

        m_odoPose = relPose;//更新，将当前获取的激光位姿位姿赋值给m_odoPose

        bool processed = false; //processed是一个bool型变量，用于记录是否已经处理完毕了激光数据  如果处理完毕了，则为true，否则为false

        // process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
        //判断机器人是否移动足够的距离或旋转足够的角度，仅当条件满足时才处理激光数据
        if (!m_count //如果是第一帧，则不处理激光数据
            || m_linearDistance >= m_linearThresholdDistance //如果移动距离超过了阈值，则处理激光数据
            || m_angularDistance >= m_angularThresholdDistance//如果旋转角度超过了阈值，则处理激光数据
            || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_)) { //如果时间间隔超过了阈值，则处理激光数据
            last_update_time_ = reading.getTime();

            if (m_outputStream.is_open()) {
                m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                m_outputStream << "FRAME " << m_readingCount;
                m_outputStream << " " << m_linearDistance;
                m_outputStream << " " << m_angularDistance << endl;
            }

            if (m_infoStream)
                m_infoStream << "update frame " << m_readingCount << endl
                             << "update ld=" << m_linearDistance << " ad=" << m_angularDistance << endl;


            cerr << "Laser Pose= " << reading.getPose().x << " " << reading.getPose().y
                 << " " << reading.getPose().theta << endl;


            //this is for converting the reading in a scan-matcher feedable form
            assert(reading.size() == m_beams);
            double *plainReading = new double[m_beams];
            for (unsigned int i = 0; i < m_beams; i++) {
                plainReading[i] = reading[i];
            }
            m_infoStream << "m_count " << m_count << endl;
            //复制激光雷达数据，对复制的数据进行处理（真正的处理过程）
            RangeReading *reading_copy =
                    new RangeReading(reading.size(),
                                     &(reading[0]),
                                     static_cast<const RangeSensor *>(reading.getSensor()),
                                     reading.getTime());
            //如果不是第一帧，则处理激光数据
            if (m_count > 0) {
                //2.1进行扫描匹配 ，将激光数据传入scanMatch函数中，得到匹配得分
                scanMatch(plainReading);

                if (m_outputStream.is_open()) {
                    m_outputStream << "LASER_READING " << reading.size() << " ";
                    m_outputStream << setiosflags(ios::fixed) << setprecision(2);
                    for (RangeReading::const_iterator b = reading.begin(); b != reading.end(); b++) {
                        m_outputStream << *b << " ";
                    }
                    OrientedPoint p = reading.getPose();
                    m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                    m_outputStream << p.x << " " << p.y << " " << p.theta << " " << reading.getTime() << endl;
                    m_outputStream << "SM_UPDATE " << m_particles.size() << " ";
                    for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++) {
                        const OrientedPoint &pose = it->pose;
                        m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
                        m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it->weight
                                       << " ";
                    }
                    m_outputStream << endl;
                }
                //应该是个空函数
                onScanmatchUpdate();
                //由于scanmatch之中对粒子的权重进行了更新，那么这时候 各个粒子的轨迹上的累计权重都需要重新计算
                updateTreeWeights(false);

                if (m_infoStream) {
                    m_infoStream << "neff= " << m_neff << endl;
                }
                if (m_outputStream.is_open()) {
                    m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                    m_outputStream << "NEFF " << m_neff << endl;
                }
                resample(plainReading, adaptParticles, reading_copy);

                //如果是第一帧，我们需要对每个粒子进行初始化。
                //如果是第一帧激光数据，直接计算ActiveArea，此时机器人的位姿都是（0，0，0）
            } else {
                m_infoStream << "Registering First Scan" << endl;
                for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
                    m_matcher.invalidateActiveArea();
                    m_matcher.computeActiveArea(it->map, it->pose, plainReading);
                    m_matcher.registerScan(it->map, it->pose, plainReading);

                    // cyr: not needed anymore, particles refer to the root in the beginning!
                    // 这是轨迹的根，没有父节点
                    TNode *node = new TNode(it->pose, 0., it->node, 0);
                    //node->reading=0;
                    node->reading = reading_copy;
                    it->node = node;

                }
            }
            //		cerr  << "Tree: normalizing, resetting and propagating weights at the end..." ;
            updateTreeWeights(false);
            //		cerr  << ".done!" <<endl;

            delete[] plainReading;
            m_lastPartPose = m_odoPose; //update the past pose for the next iteration
            m_linearDistance = 0;
            m_angularDistance = 0;
            m_count++;
            processed = true;

            //keep ready for the next step
            for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
                it->previousPose = it->pose;
            }

        }
        if (m_outputStream.is_open())
            m_outputStream << flush;
        m_readingCount++;
        return processed;
    }


    std::ofstream &GridSlamProcessor::outputStream() {
        return m_outputStream;
    }

    std::ostream &GridSlamProcessor::infoStream() {
        return m_infoStream;
    }


    int GridSlamProcessor::getBestParticleIndex() const {
        unsigned int bi = 0;
        double bw = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < m_particles.size(); i++)
            if (bw < m_particles[i].weightSum) {
                bw = m_particles[i].weightSum;
                bi = i;
            }
        return (int) bi;
    }

    void GridSlamProcessor::onScanmatchUpdate() {}

    void GridSlamProcessor::onResampleUpdate() {}

    void GridSlamProcessor::onOdometryUpdate() {}


};// end namespace




