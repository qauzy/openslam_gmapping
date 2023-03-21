
#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif

/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/

//scanMatch()函数采用NDT算法对当前坐标的激光束和每个粒子各自维护的map进行匹配，对粒子坐标进行前后左右左转右转微调，再给出匹配得分
//以传感器的数据作为输入，在函数体中遍历建图引擎的粒子集合
inline void GridSlamProcessor::scanMatch(const double* plainReading){
  // sample a new pose from each scan in the reference
  
  double sumScore=0;
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
    OrientedPoint corrected;
    double score, l, s;
    /* 计算最优的粒子
    optimize 调用了 score 这个函数 （计算粒子得分）
    在score 函数里，首先计算障碍物的坐标phit，然后将phit转换成网格坐标iPhit
    计算光束上与障碍物相邻的非障碍物网格坐标pfree,pfrree由phit沿激光束方向移动一个网格的距离得到，将pfree转换成网格坐标ipfree（增量，并不是实际值）
    在iphit 及其附近8个（m_kernelSize:default=1）栅格（pr,对应自由栅格为pf）搜索最优可能是障碍物的栅格。
    最优准则： pr 大于某一阈值，pf小于该阈值，且pr栅格的phit的平均坐标与phit的距离bestMu最小。
    得分计算： s +=exp(-1.0/m_gaussianSigma*bestMu*besMu)  参考NDT算法,距离越大，分数越小，分数的较大值集中在距离最小值处，符合正态分布模型
    至此 score 函数结束并返回粒子（currentPose）得分，然后回到optimize函数
    optimize 干的事就是 currentPose 的位姿进行微调，前、后、左、右、左转、右转 共6次，然后选取得分最高的位姿，返回最终的得分
    */
    //
    //   optomize的中心思想就是给一个初始位姿结合激光数据和地图求出最优位姿。由于这个初始位姿可能不够准确，但是我们要怎么找倒较准确的位姿（最优位姿）呢？
    //   在初始位姿的基础上，我们给他的x,y,thelta方向上依次给一个增量，计算每一次的位姿得分（得分表示在该位姿下激光束和地图的匹配程度），看看是不是比不
    //   给增量之前得分高，得分高意味着位姿更准确了。这样迭代求出最优位姿。

    //optimize 函数返回了bestScore，同时计算了corrected，即新位姿


    /*
    匹配器的函数optimize是一种爬山算法,用于寻找x^(i)t=argmaxxp(x|m(i)t−1,zt,x′(i)t)，完成对建议分布的采样
    其中corrected是一个引用方式的传参，该函数最终退出之后，corrected所记录的则是x^(i)t。 
    此外，该函数返回的是一个匹配度，记录在局部变量score中。建图引擎的成员变量m_minimumScore是一个匹配度阈值，当超过该阈值时，我们认为匹配成功，使用x^(i)t作为更新样本。 
    否则，我们仍然使用传统的运动模型下的预测状态作为更新样本
    */

    score=m_matcher.optimize(corrected, it->map, it->pose, plainReading); 
    //    it->pose=corrected;
    //判断得分是否符合要求
    //当计算出的分数比较靠谱后，就用corrected值更新当前粒子的位姿
    if (score>m_minimumScore){
      it->pose=corrected;
    } else {
	if (m_infoStream){
	  m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l <<std::endl;
	  m_infoStream << "lp:" << m_lastPartPose.x << " "  << m_lastPartPose.y << " "<< m_lastPartPose.theta <<std::endl;
	  m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " "<< m_odoPose.theta <<std::endl;
	}
    }
    /*   likelihoodAndSocre 作用是计算粒子的权重和（l），如果出现匹配失败，则 l=noHit     */
    //粒子的最优位姿计算了之后，重新计算粒子的权重(相当于粒子滤波器中的观测步骤，计算p(z|x,m))，粒子的权重由粒子的似然来表示。
    //后续会根据粒子的权重进行重采样。
    m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);
    sumScore+=score;
    it->weight+=l;
    it->weightSum+=l;

    /* 计算可活动区域
    computeActiveArea 用于计算每个粒子相应的位姿所扫描到的区域
    计算过程首先考虑了每个粒子的扫描范围会不会超过子地图的大小，如果会，则resize地图的大小
    然后定义了一个activeArea 用于设置可活动区域，调用了gridLine() 函数,这个函数如何实现的，
    请参考百度文库那篇介绍。
    */
    //set up the selective copy of the active area
    //by detaching the areas that will be updated
    m_matcher.invalidateActiveArea();
    m_matcher.computeActiveArea(it->map, it->pose, plainReading);
  }
  if (m_infoStream)
    m_infoStream << "Average Scan Matching Score=" << sumScore/m_particles.size() << std::endl;	
}

inline void GridSlamProcessor::normalize(){
  //normalize the log m_weights
  double gain=1./(m_obsSigmaGain*m_particles.size());
  double lmax= -std::numeric_limits<double>::max();
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
    lmax=it->weight>lmax?it->weight:lmax;
  }
  //cout << "!!!!!!!!!!! maxwaight= "<< lmax << endl;
  
  m_weights.clear();
  double wcum=0;
  m_neff=0;
  for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
    m_weights.push_back(exp(gain*(it->weight-lmax)));
    wcum+=m_weights.back();
    //cout << "l=" << it->weight<< endl;
  }
  
  m_neff=0;
  for (std::vector<double>::iterator it=m_weights.begin(); it!=m_weights.end(); it++){
    *it=*it/wcum;
    double w=*it;
    m_neff+=w*w;
  }
  m_neff=1./m_neff;
  
}

inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading){
  
  bool hasResampled = false;
  
  TNodeVector oldGeneration;
  for (unsigned int i=0; i<m_particles.size(); i++){
    oldGeneration.push_back(m_particles[i].node);
  }
  
  if (m_neff<m_resampleThreshold*m_particles.size()){		
    
    if (m_infoStream)
      m_infoStream  << "*************RESAMPLE***************" << std::endl;
    
    uniform_resampler<double, double> resampler;
    m_indexes=resampler.resampleIndexes(m_weights, adaptSize);
    
    if (m_outputStream.is_open()){
      m_outputStream << "RESAMPLE "<< m_indexes.size() << " ";
      for (std::vector<unsigned int>::const_iterator it=m_indexes.begin(); it!=m_indexes.end(); it++){
	m_outputStream << *it <<  " ";
      }
      m_outputStream << std::endl;
    }
    
    onResampleUpdate();
    //BEGIN: BUILDING TREE
    ParticleVector temp;
    unsigned int j=0;
    std::vector<unsigned int> deletedParticles;  		//this is for deleteing the particles which have been resampled away.
    
    //		cerr << "Existing Nodes:" ;
    for (unsigned int i=0; i<m_indexes.size(); i++){
      //			cerr << " " << m_indexes[i];
      while(j<m_indexes[i]){
	deletedParticles.push_back(j);
	j++;
			}
      if (j==m_indexes[i])
	j++;
      Particle & p=m_particles[m_indexes[i]];
      TNode* node=0;
      TNode* oldNode=oldGeneration[m_indexes[i]];
      //			cerr << i << "->" << m_indexes[i] << "B("<<oldNode->childs <<") ";
      node=new	TNode(p.pose, 0, oldNode, 0);
      //node->reading=0;
      node->reading=reading;
      //			cerr << "A("<<node->parent->childs <<") " <<endl;
      
      temp.push_back(p);
      temp.back().node=node;
      temp.back().previousIndex=m_indexes[i];
    }
    while(j<m_indexes.size()){
      deletedParticles.push_back(j);
      j++;
    }
    //		cerr << endl;
    std::cerr <<  "Deleting Nodes:";
    for (unsigned int i=0; i<deletedParticles.size(); i++){
      std::cerr <<" " << deletedParticles[i];
      delete m_particles[deletedParticles[i]].node;
      m_particles[deletedParticles[i]].node=0;
    }
    std::cerr  << " Done" <<std::endl;
    
    //END: BUILDING TREE
    std::cerr << "Deleting old particles..." ;
    m_particles.clear();
    std::cerr << "Done" << std::endl;
    std::cerr << "Copying Particles and  Registering  scans...";
    for (ParticleVector::iterator it=temp.begin(); it!=temp.end(); it++){
      it->setWeight(0);
      m_matcher.invalidateActiveArea();
      m_matcher.registerScan(it->map, it->pose, plainReading);
      m_particles.push_back(*it);
    }
    std::cerr  << " Done" <<std::endl;
    hasResampled = true;
  } else {
    int index=0;
    std::cerr << "Registering Scans:";
    TNodeVector::iterator node_it=oldGeneration.begin();
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      //create a new node in the particle tree and add it to the old tree
      //BEGIN: BUILDING TREE  
      TNode* node=0;
      node=new TNode(it->pose, 0.0, *node_it, 0);
      
      //node->reading=0;
      node->reading=reading;
      it->node=node;

      //END: BUILDING TREE
      m_matcher.invalidateActiveArea();
      m_matcher.registerScan(it->map, it->pose, plainReading);
      it->previousIndex=index;
      index++;
      node_it++;
      
    }
    std::cerr  << "Done" <<std::endl;
    
  }
  //END: BUILDING TREE
  
  return hasResampled;
}
