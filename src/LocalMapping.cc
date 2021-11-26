/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2020 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/


#include "LocalMapping.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Optimizer.h"
#include "Converter.h"
#include "Config.h"

#include<mutex>
#include<chrono>

namespace ORB_SLAM3
{

LocalMapping::LocalMapping(System* pSys, Atlas *pAtlas, const float bMonocular, bool bInertial, const string &_strSeqName):
    mpSystem(pSys), mbMonocular(bMonocular), mbInertial(bInertial), mbResetRequested(false), mbResetRequestedActiveMap(false), mbFinishRequested(false), mbFinished(true), mpAtlas(pAtlas), bInitializing(false),
    mbAbortBA(false), mbStopped(false), mbStopRequested(false), mbNotStop(false), mbAcceptKeyFrames(true),
    mbNewInit(false), mIdxInit(0), mScale(1.0), mInitSect(0), mbNotBA1(true), mbNotBA2(true), infoInertial(Eigen::MatrixXd::Zero(9,9))
{
    /*
     * mbStopRequested：    外部线程调用，为true，表示外部线程请求停止 local mapping
     * mbStopped：          为true表示可以并终止localmapping 线程
     * mbNotStop：          true，表示不要停止 localmapping 线程，因为要插入关键帧了。需要和 mbStopped 结合使用
     * mbAcceptKeyFrames：  true，允许接受关键帧。tracking 和local mapping 之间的关键帧调度
     * mbAbortBA：          是否流产BA优化的标志位
     * mbFinishRequested：  请求终止当前线程的标志。注意只是请求，不一定终止。终止要看 mbFinished
     * mbResetRequested：   请求当前线程复位的标志。true，表示一直请求复位，但复位还未完成；表示复位完成为false
     * mbFinished：         判断最终LocalMapping::Run() 是否完成的标志。
     */
    mnMatchesInliers = 0;

    mbBadImu = false;

    mTinit = 0.f;

    mNumLM = 0;
    mNumKFCulling=0;

#ifdef REGISTER_TIMES
    nLBA_exec = 0;
    nLBA_abort = 0;
#endif

}

// 设置回环检测线程句柄
void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

// 设置跟踪线程句柄
void LocalMapping::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}

// 线程主函数
void LocalMapping::Run()
{

    // 标记状态，表示当前run函数正在运行，尚未结束
    mbFinished = false;
    // 主循环
    while(1)
    {
        // Tracking will see that Local Mapping is busy
        // Step 1 告诉Tracking，LocalMapping正处于繁忙状态，请不要给我发送关键帧打扰我
        // LocalMapping线程处理的关键帧都是Tracking线程发过来的
        SetAcceptKeyFrames(false);

        // Check if there are keyframes in the queue
        // 等待处理的关键帧列表不为空 并且imu正常
        if(CheckNewKeyFrames() && !mbBadImu)
        {

#ifdef REGISTER_TIMES
            double timeLBA_ms = 0;
            double timeKFCulling_ms = 0;

            std::chrono::steady_clock::time_point time_StartProcessKF = std::chrono::steady_clock::now();
#endif
            // BoW conversion and insertion in Map
            // Step 2 处理列表中的关键帧，包括计算BoW、更新观测、描述子、共视图，插入到地图等
            ProcessNewKeyFrame();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndProcessKF = std::chrono::steady_clock::now();

            double timeProcessKF = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndProcessKF - time_StartProcessKF).count();
            vdKFInsert_ms.push_back(timeProcessKF);
#endif
            // Check recent MapPoints
            // Step 3 根据地图点的观测情况剔除质量不好的地图点
            MapPointCulling();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCulling = std::chrono::steady_clock::now();

            double timeMPCulling = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMPCulling - time_EndProcessKF).count();
            vdMPCulling_ms.push_back(timeMPCulling);
#endif
            // Triangulate new MapPoints
            // Step 4 当前关键帧与相邻关键帧通过三角化产生新的地图点，使得跟踪更稳
            CreateNewMapPoints();

            mbAbortBA = false; // 注意orbslam2中放在了函数SearchInNeighbors（用到了mbAbortBA）后面，应该放这里更合适
			// 已经处理完队列中的最后的一个关键帧
            if(!CheckNewKeyFrames())
            {
                // Find more matches in neighbor keyframes and fuse point duplications
                //  Step 5 检查并融合当前关键帧与相邻关键帧帧（两级相邻）中重复的地图点
                // 先完成相邻关键帧与当前关键帧的地图点的融合（在相邻关键帧中查找当前关键帧的地图点），
                // 再完成当前关键帧与相邻关键帧的地图点的融合（在当前关键帧中查找当前相邻关键帧的地图点）
                SearchInNeighbors();
            }

#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCreation = std::chrono::steady_clock::now();

            double timeMPCreation = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMPCreation - time_EndMPCulling).count();
            vdMPCreation_ms.push_back(timeMPCreation);
#endif

            bool b_doneLBA = false;
            int num_FixedKF_BA = 0;
            int num_OptKF_BA = 0;
            int num_MPs_BA = 0;
            int num_edges_BA = 0;
			// 已经处理完队列中的最后的一个关键帧，并且闭环检测没有请求停止LocalMapping
            if(!CheckNewKeyFrames() && !stopRequested())
            {
                // 当前地图中关键帧数目大于2个
                if(mpAtlas->KeyFramesInMap()>2)
                {
                    // Step 6.1 处于IMU模式并且当前关键帧所在的地图已经完成IMU初始化
                    if(mbInertial && mpCurrentKeyFrame->GetMap()->isImuInitialized())
                    {
                        // 计算上一关键帧到当前关键帧相机光心的距离 + 上上关键帧到上一关键帧相机光心的距离
                        float dist = cv::norm(mpCurrentKeyFrame->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->GetCameraCenter()) +
                                cv::norm(mpCurrentKeyFrame->mPrevKF->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->mPrevKF->GetCameraCenter());
                        // 如果距离大于5厘米，记录当前KF和上一KF时间戳的差，累加到mTinit
                        if(dist>0.05)
                            mTinit += mpCurrentKeyFrame->mTimeStamp - mpCurrentKeyFrame->mPrevKF->mTimeStamp;
                        // 当前关键帧所在的地图尚未完成IMU BA2（IMU第三阶段初始化）   
                        if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2())
                        {
                            // 如果累计时间差小于10s 并且 距离小于2厘米，认为运动幅度太小，不足以初始化IMU，将mbBadImu设置为true
                            if((mTinit<10.f) && (dist<0.02))
                            {
                                cout << "Not enough motion for initializing. Reseting..." << endl;
                                unique_lock<mutex> lock(mMutexReset);
                                mbResetRequestedActiveMap = true;
                                mpMapToReset = mpCurrentKeyFrame->GetMap();
                                mbBadImu = true;    //在跟踪线程里会重置当前活跃地图
                            }
                        }
                        // 判断成功跟踪匹配的点数是否足够多
                        // 条件---------1.1、跟踪成功的内点数目大于75-----1.2、并且是单目--或--2.1、跟踪成功的内点数目大于100-----2.2、并且不是单目             
                        bool bLarge = ((mpTracker->GetMatchesInliers()>75)&&mbMonocular)||((mpTracker->GetMatchesInliers()>100)&&!mbMonocular);
                        // 局部地图+IMU一起优化，优化关键帧位姿、地图点、IMU参数
                        Optimizer::LocalInertialBA(mpCurrentKeyFrame, &mbAbortBA, mpCurrentKeyFrame->GetMap(),num_FixedKF_BA,num_OptKF_BA,num_MPs_BA,num_edges_BA, bLarge, !mpCurrentKeyFrame->GetMap()->GetIniertialBA2());
                        b_doneLBA = true;
                    }
                    else
                    {
						// Step 6.2 不是IMU模式或者当前关键帧所在的地图还未完成IMU初始化
						// 局部地图BA，不包括IMU数据
						// 注意这里的第二个参数是按地址传递的,当这里的 mbAbortBA 状态发生变化时，能够及时执行/停止BA
                        // 局部地图优化，不包括IMU信息。优化关键帧位姿、地图点
                        Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA, mpCurrentKeyFrame->GetMap(),num_FixedKF_BA,num_OptKF_BA,num_MPs_BA,num_edges_BA);
                        b_doneLBA = true;
                    }

                }
#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndLBA = std::chrono::steady_clock::now();

                if(b_doneLBA)
                {
                    timeLBA_ms = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLBA - time_EndMPCreation).count();
                    vdLBASync_ms.push_back(timeLBA_ms);

                    nLBA_exec += 1;
                    if(mbAbortBA)
                    {
                        nLBA_abort += 1;
                    }
                    vnLBA_edges.push_back(num_edges_BA);
                    vnLBA_KFopt.push_back(num_OptKF_BA);
                    vnLBA_KFfixed.push_back(num_FixedKF_BA);
                    vnLBA_MPs.push_back(num_MPs_BA);
                }

#endif

                // Initialize IMU here
                // Step 7 当前关键帧所在地图未完成IMU初始化（第一阶段）
                if(!mpCurrentKeyFrame->GetMap()->isImuInitialized() && mbInertial)
                {
                    // 在函数InitializeIMU里设置IMU成功初始化标志 SetImuInitialized
                    // IMU第一阶段初始化
                    if (mbMonocular)
                        InitializeIMU(1e2, 1e10, true);
                    else
                        InitializeIMU(1e2, 1e5, true);
                }


                // Check redundant local Keyframes
                // 跟踪中关键帧插入条件比较松，交给LocalMapping线程的关键帧会比较密，这里再删除冗余
                // Step 8 检测并剔除当前帧相邻的关键帧中冗余的关键帧
                // 冗余的判定：该关键帧的90%的地图点可以被其它关键帧观测到
                KeyFrameCulling();

#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndKFCulling = std::chrono::steady_clock::now();

                timeKFCulling_ms = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndKFCulling - time_EndLBA).count();
                vdKFCullingSync_ms.push_back(timeKFCulling_ms);
#endif
				// Step 9 如果距离IMU第一阶段初始化成功累计时间差小于100s，进行VIBA
                if ((mTinit<100.0f) && mbInertial)
                {
                    // Step 9.1 根据条件判断是否进行VIBA1（IMU第二阶段初始化）
                    // 条件：1、当前关键帧所在的地图还未完成IMU初始化---并且--------2、正常跟踪状态----------
                    if(mpCurrentKeyFrame->GetMap()->isImuInitialized() && mpTracker->mState==Tracking::OK) // Enter here everytime local-mapping is called
                    {
                        // 当前关键帧所在的地图还未完成VIBA 1
                        if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA1()){
                            // 如果累计时间差大于5s，开始VIBA1（IMU第二阶段初始化）
                            if (mTinit>5.0f)
                            {
                                cout << "start VIBA 1" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA1();
                                if (mbMonocular)
                                    InitializeIMU(1.f, 1e5, true); // 1.f, 1e5
                                else
                                    InitializeIMU(1.f, 1e5, true); // 1.f, 1e5

                                cout << "end VIBA 1" << endl;
                            }
                        }
                        //else if (mbNotBA2){
                        // Step 9.2 根据条件判断是否进行VIBA2（IMU第三阶段初始化）
                        // 当前关键帧所在的地图还未完成VIBA 2
                        else if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2()){
                            // 如果累计时间差大于15s，开始VIBA2（IMU第三阶段初始化）
                            if (mTinit>15.0f){ // 15.0f
                                cout << "start VIBA 2" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA2();
                                if (mbMonocular)
                                    InitializeIMU(0.f, 0.f, true); // 0.f, 0.f
                                else
                                    InitializeIMU(0.f, 0.f, true);

                                cout << "end VIBA 2" << endl;
                            }
                        }

                        // scale refinement
                        // Step 9.3 在关键帧小于100时，会在满足一定时间间隔后多次进行尺度、重力方向优化
                        if (((mpAtlas->KeyFramesInMap())<=100) &&
                                ((mTinit>25.0f && mTinit<25.5f)||       
                                (mTinit>35.0f && mTinit<35.5f)||
                                (mTinit>45.0f && mTinit<45.5f)||
                                (mTinit>55.0f && mTinit<55.5f)||
                                (mTinit>65.0f && mTinit<65.5f)||
                                (mTinit>75.0f && mTinit<75.5f))){
                            cout << "start scale ref" << endl;
                            if (mbMonocular)
                                // 使用了所有关键帧，但只优化尺度和重力方向
                                ScaleRefinement();
                            cout << "end scale ref" << endl;
                        }
                    }
                }
            }

#ifdef REGISTER_TIMES
            vdLBA_ms.push_back(timeLBA_ms);
            vdKFCulling_ms.push_back(timeKFCulling_ms);
#endif

			// Step 10 将当前帧加入到闭环检测队列中
            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);


#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndLocalMap = std::chrono::steady_clock::now();

            double timeLocalMap = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLocalMap - time_StartProcessKF).count();
            vdLMTotal_ms.push_back(timeLocalMap);
#endif
        }
        else if(Stop() && !mbBadImu) // 当要终止当前线程的时候
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                // cout << "LM: usleep if is stopped" << endl;
				// 如果还没有结束利索,那么等等它
                usleep(3000);
            }
            // 然后确定终止了就跳出这个线程的主循环
            if(CheckFinish())
                break;
        }

        // 查看是否有复位线程的请求
        ResetIfRequested();

        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(true);

        // 如果当前线程已经结束了就跳出主循环
        if(CheckFinish())
            break;

        // cout << "LM: normal usleep" << endl;
        usleep(3000);
    }

    // 设置线程已经终止
    SetFinish();
}

// 插入关键帧,由外部（Tracking）线程调用;这里只是插入到列表中,等待线程主函数对其进行处理
void LocalMapping::InsertKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexNewKFs);
    // 将关键帧插入到列表中
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA=true;
}

// 查看列表中是否有等待被插入的关键帧,
bool LocalMapping::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexNewKFs);
    return(!mlNewKeyFrames.empty());
}

/**
 * @brief 处理列表中的关键帧，包括计算BoW、更新观测、描述子、共视图，插入到地图等
 * 
 */
void LocalMapping::ProcessNewKeyFrame()
{
    // Step 1：从缓冲队列中取出一帧关键帧
    // 该关键帧队列是Tracking线程向LocalMapping中插入的关键帧组成
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        // 取出列表中最前面的关键帧，作为当前要处理的关键帧
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        // 取出最前面的关键帧后，在原来的列表里删掉该关键帧
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    // Step 2：计算该关键帧特征点的Bow信息
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    // Step 3：当前处理关键帧中有效的地图点，更新normal，描述子等信息
    // TrackLocalMap中和当前帧新匹配上的地图点和当前关键帧进行关联绑定
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    // 对当前处理的这个关键帧中的所有的地图点展开遍历
    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    // 如果地图点不是来自当前帧的观测，为当前地图点添加观测
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    // 获得该点的平均观测方向和观测距离范围
                    pMP->UpdateNormalAndDepth();
                    // 更新地图点的最佳描述子
                    pMP->ComputeDistinctiveDescriptors();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    // 如果当前帧中已经包含了这个地图点,但是这个地图点中却没有包含这个关键帧的信息
                    // 这些地图点可能来自双目或RGBD跟踪过程中新生成的地图点，或者是CreateNewMapPoints 中通过三角化产生
                    // 将上述地图点放入mlpRecentAddedMapPoints，等待后续MapPointCulling函数的检验
                    mlpRecentAddedMapPoints.push_back(pMP); 
                }
            }
        }
    }

    // Update links in the Covisibility Graph
    // Step 4：更新关键帧间的连接关系（共视图）
    mpCurrentKeyFrame->UpdateConnections();

    // Insert Keyframe in Map
    // Step 5：将该关键帧插入到地图中
    mpAtlas->AddKeyFrame(mpCurrentKeyFrame);
}

/**
 * @brief 处理新的关键帧，使队列为空，注意这里只是处理了关键帧，并没有生成MP
 */
void LocalMapping::EmptyQueue()
{
    while(CheckNewKeyFrames())
        ProcessNewKeyFrame();
}
/**
 * @brief 检查新增地图点，根据地图点的观测情况剔除质量不好的新增的地图点
 * mlpRecentAddedMapPoints：存储新增的地图点，这里是要删除其中不靠谱的
 */
void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    // Step 1：根据相机类型设置不同的观测阈值
    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3; // MODIFICATION_STEREO_IMU here 3
    const int cnThObs = nThObs;

    int borrar = mlpRecentAddedMapPoints.size();
	// Step 2：遍历检查的新添加的MapPoints
    while(lit!=mlpRecentAddedMapPoints.end())
    {
        MapPoint* pMP = *lit;

        if(pMP->isBad())
            // Step 2.1：已经是坏点的MapPoints直接从检查链表中删除
            lit = mlpRecentAddedMapPoints.erase(lit);
        else if(pMP->GetFoundRatio()<0.25f)
        {
            // Step 2.2：跟踪到该MapPoint的Frame数相比预计可观测到该MapPoint的Frame数的比例小于25%，删除
            // (mnFound/mnVisible） < 25%
            // mnFound ：地图点被多少帧（包括普通帧）看到，次数越多越好
            // mnVisible：地图点应该被看到的次数
            // (mnFound/mnVisible）：对于大FOV镜头这个比例会高，对于窄FOV镜头这个比例会低
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=2 && pMP->Observations()<=cnThObs)
        {
            // Step 2.3：从该点建立开始，到现在已经过了不小于2个关键帧
            // 但是观测到该点的关键帧数却不超过cnThObs帧，那么删除该点
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=3)
            // Step 2.4：从建立该点开始，已经过了3个关键帧而没有被剔除，则认为是质量高的点
            // 因此没有SetBadFlag()，仅从队列中删除，放弃继续对该MapPoint的检测
            lit = mlpRecentAddedMapPoints.erase(lit);
        else
        {
            lit++;
            borrar--;
        }
    }
    //cout << "erase MP: " << borrar << endl;
}

/**
 * @brief 用当前关键帧与相邻关键帧通过三角化产生新的地图点，使得跟踪更稳
 * 
 */
void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    // nn表示搜索最佳共视关键帧的数目
    // 不同传感器下要求不一样,单目的时候需要有更多的具有较好共视关系的关键帧来建立地图
    int nn = 10;
    // For stereo inertial case
    if(mbMonocular)
        nn=20;
	// Step 1：在当前关键帧的共视关键帧中找到共视程度最高的nn帧相邻关键帧
    vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    // imu模式下在附近添加更多的帧进来
    if (mbInertial)
    {
        KeyFrame* pKF = mpCurrentKeyFrame;
        int count=0;
        // 在总数不够且上一关键帧存在，且添加的帧没有超过总数时
        while((vpNeighKFs.size()<=nn)&&(pKF->mPrevKF)&&(count++<nn))
        {
            vector<KeyFrame*>::iterator it = std::find(vpNeighKFs.begin(), vpNeighKFs.end(), pKF->mPrevKF);
            if(it==vpNeighKFs.end())
                vpNeighKFs.push_back(pKF->mPrevKF);
            pKF = pKF->mPrevKF;
        }
    }

    float th = 0.6f;
	// 特征点匹配配置 最佳距离 < 0.6*次佳距离，比较苛刻了。不检查旋转
    ORBmatcher matcher(th,false);

    // 取出当前帧从世界坐标系到相机坐标系的变换矩阵
    auto Rcw1 = mpCurrentKeyFrame->GetRotation_();
    auto Rwc1 = Rcw1.t();
    auto tcw1 = mpCurrentKeyFrame->GetTranslation_();
    cv::Matx44f Tcw1{Rcw1(0,0),Rcw1(0,1),Rcw1(0,2),tcw1(0),
                     Rcw1(1,0),Rcw1(1,1),Rcw1(1,2),tcw1(1),
                     Rcw1(2,0),Rcw1(2,1),Rcw1(2,2),tcw1(2),
                     0.f,0.f,0.f,1.f};

    // 得到当前关键帧（左目）光心在世界坐标系中的坐标、内参
    auto Ow1 = mpCurrentKeyFrame->GetCameraCenter_();

    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    // 用于后面的点深度的验证;这里的1.5是经验值
    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactor;

    // Search matches with epipolar restriction and triangulate
    // Step 2：遍历相邻关键帧vpNeighKFs
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {
        // 下面的过程会比较耗费时间,因此如果有新的关键帧需要处理的话,就先去处理新的关键帧吧
        if(i>0 && CheckNewKeyFrames())// && (mnMatchesInliers>50))
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        GeometricCamera* pCamera1 = mpCurrentKeyFrame->mpCamera, *pCamera2 = pKF2->mpCamera;

        // Check first that baseline is not too short
        // 邻接的关键帧光心在世界坐标系中的坐标
        auto Ow2 = pKF2->GetCameraCenter_();
        // 基线向量，两个关键帧间的相机位移
        auto vBaseline = Ow2-Ow1;
        // 基线长度
        const float baseline = cv::norm(vBaseline);

        // Step 3：判断相机运动的基线是不是足够长
        if(!mbMonocular)
        {
            // 如果是双目相机，关键帧间距小于本身的基线时不生成3D点
            // 因为太短的基线下能够恢复的地图点不稳定
            if(baseline<pKF2->mb)
            continue;
        }
        else
        {
            // 单目相机情况
            // 邻接关键帧的场景深度中值
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            // baseline与景深的比例
            const float ratioBaselineDepth = baseline/medianDepthKF2;
            // 如果特别远(比例特别小)，基线太短恢复3D点不准，那么跳过当前邻接的关键帧，不生成3D点
            if(ratioBaselineDepth<0.01)
                continue;
        }

        // Compute Fundamental Matrix
        // Step 4：根据两个关键帧的位姿计算它们之间的基础矩阵
        auto F12 = ComputeF12_(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        // Step 5：通过BoW对两关键帧的未匹配的特征点快速匹配，用极线约束抑制离群点，生成新的匹配点对
        vector<pair<size_t,size_t> > vMatchedIndices;
        // imu相关，非imu时为false
        bool bCoarse = mbInertial &&
                ((!mpCurrentKeyFrame->GetMap()->GetIniertialBA2() && mpCurrentKeyFrame->GetMap()->GetIniertialBA1())||
                 mpTracker->mState==Tracking::RECENTLY_LOST);
        // 通过极线约束的方式找到匹配点（且该点还没有成为MP，注意非单目已经生成的MP这里直接跳过不做匹配，所以最后并不会覆盖掉特征点对应的MP）
        matcher.SearchForTriangulation_(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices,false,bCoarse);

        auto Rcw2 = pKF2->GetRotation_();
        auto Rwc2 = Rcw2.t();
        auto tcw2 = pKF2->GetTranslation_();
        cv::Matx44f Tcw2{Rcw2(0,0),Rcw2(0,1),Rcw2(0,2),tcw2(0),
                         Rcw2(1,0),Rcw2(1,1),Rcw2(1,2),tcw2(1),
                         Rcw2(2,0),Rcw2(2,1),Rcw2(2,2),tcw2(2),
                         0.f,0.f,0.f,1.f};

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        // Step 6：对每对匹配通过三角化生成3D点,和 Triangulate函数差不多
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            // Step 6.1：取出匹配特征点

            // 当前匹配对在当前关键帧中的索引
            const int &idx1 = vMatchedIndices[ikp].first;
            
            // 当前匹配对在邻接关键帧中的索引
            const int &idx2 = vMatchedIndices[ikp].second;

            // 当前匹配在当前关键帧中的特征点
            const cv::KeyPoint &kp1 = (mpCurrentKeyFrame -> NLeft == -1) ? mpCurrentKeyFrame->mvKeysUn[idx1]
                                                                         : (idx1 < mpCurrentKeyFrame -> NLeft) ? mpCurrentKeyFrame -> mvKeys[idx1]
                                                                                                               : mpCurrentKeyFrame -> mvKeysRight[idx1 - mpCurrentKeyFrame -> NLeft];
            // mvuRight中存放着极限校准后双目特征点在右目对应的像素横坐标，如果不是基线校准的双目或者没有找到匹配点，其值将为-1（或者rgbd）
            const float kp1_ur=mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = (!mpCurrentKeyFrame->mpCamera2 && kp1_ur>=0);
            // 查看点idx1是否为右目的点
            const bool bRight1 = (mpCurrentKeyFrame -> NLeft == -1 || idx1 < mpCurrentKeyFrame -> NLeft) ? false
                                                                               : true;
            // 当前匹配在邻接关键帧中的特征点
            const cv::KeyPoint &kp2 = (pKF2 -> NLeft == -1) ? pKF2->mvKeysUn[idx2]
                                                            : (idx2 < pKF2 -> NLeft) ? pKF2 -> mvKeys[idx2]
                                                                                     : pKF2 -> mvKeysRight[idx2 - pKF2 -> NLeft];
            // mvuRight中存放着双目的深度值，如果不是双目，其值将为-1
            // mvuRight中存放着极限校准后双目特征点在右目对应的像素横坐标，如果不是基线校准的双目或者没有找到匹配点，其值将为-1（或者rgbd）
            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = (!pKF2->mpCamera2 && kp2_ur>=0);
            // 查看点idx2是否为右目的点
            const bool bRight2 = (pKF2 -> NLeft == -1 || idx2 < pKF2 -> NLeft) ? false
                                                                               : true;
            // 当目前为左右目时，确定两个点所在相机之间的位姿关系
            if(mpCurrentKeyFrame->mpCamera2 && pKF2->mpCamera2){
                if(bRight1 && bRight2){
                    Rcw1 = mpCurrentKeyFrame->GetRightRotation_();
                    Rwc1 = Rcw1.t();
                    tcw1 = mpCurrentKeyFrame->GetRightTranslation_();
                    Tcw1 = mpCurrentKeyFrame->GetRightPose_();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter_();

                    Rcw2 = pKF2->GetRightRotation_();
                    Rwc2 = Rcw2.t();
                    tcw2 = pKF2->GetRightTranslation_();
                    Tcw2 = pKF2->GetRightPose_();
                    Ow2 = pKF2->GetRightCameraCenter_();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera2;
                }
                else if(bRight1 && !bRight2){
                    Rcw1 = mpCurrentKeyFrame->GetRightRotation_();
                    Rwc1 = Rcw1.t();
                    tcw1 = mpCurrentKeyFrame->GetRightTranslation_();
                    Tcw1 = mpCurrentKeyFrame->GetRightPose_();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter_();

                    Rcw2 = pKF2->GetRotation_();
                    Rwc2 = Rcw2.t();
                    tcw2 = pKF2->GetTranslation_();
                    Tcw2 = pKF2->GetPose_();
                    Ow2 = pKF2->GetCameraCenter_();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera;
                }
                else if(!bRight1 && bRight2){
                    Rcw1 = mpCurrentKeyFrame->GetRotation_();
                    Rwc1 = Rcw1.t();
                    tcw1 = mpCurrentKeyFrame->GetTranslation_();
                    Tcw1 = mpCurrentKeyFrame->GetPose_();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter_();

                    Rcw2 = pKF2->GetRightRotation_();
                    Rwc2 = Rcw2.t();
                    tcw2 = pKF2->GetRightTranslation_();
                    Tcw2 = pKF2->GetRightPose_();
                    Ow2 = pKF2->GetRightCameraCenter_();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera2;
                }
                else{
                    Rcw1 = mpCurrentKeyFrame->GetRotation_();
                    Rwc1 = Rcw1.t();
                    tcw1 = mpCurrentKeyFrame->GetTranslation_();
                    Tcw1 = mpCurrentKeyFrame->GetPose_();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter_();

                    Rcw2 = pKF2->GetRotation_();
                    Rwc2 = Rcw2.t();
                    tcw2 = pKF2->GetTranslation_();
                    Tcw2 = pKF2->GetPose_();
                    Ow2 = pKF2->GetCameraCenter_();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera;
                }
            }

            // Check parallax between rays
			// Step 6.2：利用匹配点反投影得到视差角
			// 特征点反投影,其实得到的是在各自相机坐标系下的一个非归一化的方向向量,和这个点的反投影射线重合
            auto xn1 = pCamera1->unprojectMat_(kp1.pt);
            auto xn2 = pCamera2->unprojectMat_(kp2.pt);
            // 由相机坐标系转到世界坐标系(得到的是那条反投影射线的一个同向向量在世界坐标系下的表示,还是只能够表示方向)，得到视差角余弦值
            auto ray1 = Rwc1*xn1;
            auto ray2 = Rwc2*xn2;
            // 这个就是求向量之间角度公式
            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            // 加1是为了让cosParallaxStereo随便初始化为一个很大的值
            float cosParallaxStereo = cosParallaxRays+1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            // Step 6.3：对于双目，利用双目得到视差角；单目相机没有特殊操作
            if(bStereo1)
                // 传感器是双目相机,并且当前的关键帧的这个点有对应的深度
                // 假设是平行的双目相机，计算出两个相机观察这个点的时候的视差角;
                // ? 感觉直接使用向量夹角的方式计算会准确一些啊（双目的时候），那么为什么不直接使用那个呢？
                // 回答：因为双目深度值、基线是更可靠的，比特征匹配再三角化出来的稳
                cosParallaxStereo1 = cos(2*atan2(mpCurrentKeyFrame->mb/2,mpCurrentKeyFrame->mvDepth[idx1]));
            else if(bStereo2)
                //传感器是双目相机,并且邻接的关键帧的这个点有对应的深度，和上面一样操作
                cosParallaxStereo2 = cos(2*atan2(pKF2->mb/2,pKF2->mvDepth[idx2]));
            
            // 得到双目观测的视差角
            cosParallaxStereo = min(cosParallaxStereo1,cosParallaxStereo2);
			// Step 6.4：三角化恢复3D点
            cv::Matx31f x3D;
            bool bEstimated = false;
            // cosParallaxRays>0 && (bStereo1 || bStereo2 || cosParallaxRays<0.9998)表明视线角正常
            // cosParallaxRays<cosParallaxStereo表明前后帧视线角比双目视线角大，所以用前后帧三角化而来，反之使用双目的，如果没有双目则跳过
            // 视差角度小时用三角法恢复3D点，视差角大时（离相机近）用双目恢复3D点（双目以及深度有效）
            if(cosParallaxRays<cosParallaxStereo && cosParallaxRays>0 && (bStereo1 || bStereo2 ||
               (cosParallaxRays<0.9998 && mbInertial) || (cosParallaxRays<0.9998 && !mbInertial)))
            {
                // Linear Triangulation Method
                // 见Initializer.cpp的Triangulate函数,A矩阵构建的方式类似，不同的是乘的反对称矩阵那个是像素坐标构成的，而这个是归一化坐标构成的
                // Pc = Tcw*Pw， 此处Tcw行数默认为3，因为不会计算第4行所以没有去掉，
                // 左右两面乘Pc的反对称矩阵 [Pc]x * Tcw *Pw = 0 构成了A矩阵，中间涉及一个尺度a，因为都是归一化平面，但右面是0所以直接可以约掉不影响最后的尺度
                //  0 -1 y    Tcw.row(0)     -Tcw.row(1) + y*Tcw.row(2)
                //  1 0 -x *  Tcw.row(1)  =   Tcw.row(0) - x*Tcw.row(2)
                // -y x  0    Tcw.row(2)    x*Tcw.row(1) - y*Tcw.row(0)
                // 发现上述矩阵线性相关，所以取前两维，两个点构成了4行的矩阵，就是如下的操作，求出的是4维的结果[X,Y,Z,A]，所以需要除以最后一维使之为1，就成了[X,Y,Z,1]这种齐次形式
                cv::Matx14f A_r0 = xn1(0) * Tcw1.row(2) - Tcw1.row(0);
                cv::Matx14f A_r1 = xn1(1) * Tcw1.row(2) - Tcw1.row(1);
                cv::Matx14f A_r2 = xn2(0) * Tcw2.row(2) - Tcw2.row(0);
                cv::Matx14f A_r3 = xn2(1) * Tcw2.row(2) - Tcw2.row(1);
                cv::Matx44f A{A_r0(0), A_r0(1), A_r0(2), A_r0(3),
                              A_r1(0), A_r1(1), A_r1(2), A_r1(3),
                              A_r2(0), A_r2(1), A_r2(2), A_r2(3),
                              A_r3(0), A_r3(1), A_r3(2), A_r3(3)};

                cv::Matx44f u,vt;
                cv::Matx41f w;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                cv::Matx41f x3D_h = vt.row(3).t();

                if(x3D_h(3)==0)
                    continue;
                // 归一化成为齐次坐标,然后提取前面三个维度作为欧式坐标
                // Euclidean coordinates
                // x3D = x3D_h.get_minor<3,1>(0,0) / x3D_h(3);
                x3D = cv::Matx31f(x3D_h.get_minor<3,1>(0,0)(0) / x3D_h(3), x3D_h.get_minor<3,1>(0,0)(1) / x3D_h(3), x3D_h.get_minor<3,1>(0,0)(2) / x3D_h(3));

                bEstimated = true;

            }
            else if(bStereo1 && cosParallaxStereo1<cosParallaxStereo2)
            {
                // 如果是双目，用视差角更大的那个双目信息来恢复，直接用已知3D点反投影了
                x3D = mpCurrentKeyFrame->UnprojectStereo_(idx1);
                bEstimated = true;
            }
            else if(bStereo2 && cosParallaxStereo2<cosParallaxStereo1)
            {
                x3D = pKF2->UnprojectStereo_(idx2);
                bEstimated = true;
            }
            else
            {
                continue; //No stereo and very low parallax
            }
            // 为方便后续计算，转换成为了行向量
            cv::Matx13f x3Dt = x3D.t();

            if(!bEstimated) continue;
            //Check triangulation in front of cameras
            // Step 6.5：检测生成的3D点是否在相机前方,不在的话就放弃这个点
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
            // Step 6.6：计算3D点在当前关键帧下的重投影误差
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1(1);
            const float invz1 = 1.0/z1;

            if(!bStereo1)
            {
				// 单目情况下
                cv::Point2f uv1 = pCamera1->project(cv::Point3f(x1,y1,z1));
                float errX1 = uv1.x - kp1.pt.x;
                float errY1 = uv1.y - kp1.pt.y;
                // 假设测量有一个像素的偏差，2自由度卡方检验阈值是5.991
                if((errX1*errX1+errY1*errY1)>5.991*sigmaSquare1)
                    continue;

            }
            else
            {
                // 双目情况
                float u1 = fx1*x1*invz1+cx1;
                // 根据视差公式计算假想的右目坐标
                float u1_r = u1 - mpCurrentKeyFrame->mbf*invz1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                // 自由度为3，卡方检验阈值是7.8
                if((errX1*errX1+errY1*errY1+errX1_r*errX1_r)>7.8*sigmaSquare1)
                    continue;
            }

            //Check reprojection error in second keyframe
            // 计算3D点在另一个关键帧下的重投影误差，操作同上
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2(1);
            const float invz2 = 1.0/z2;
            if(!bStereo2)
            {
                cv::Point2f uv2 = pCamera2->project(cv::Point3f(x2,y2,z2));
                float errX2 = uv2.x - kp2.pt.x;
                float errY2 = uv2.y - kp2.pt.y;
                if((errX2*errX2+errY2*errY2)>5.991*sigmaSquare2)
                    continue;
            }
            else
            {
                float u2 = fx2*x2*invz2+cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf*invz2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if((errX2*errX2+errY2*errY2+errX2_r*errX2_r)>7.8*sigmaSquare2)
                    continue;
            }

            //Check scale consistency
            // Step 6.7：检查尺度连续性

            // 世界坐标系下，3D点与相机间的向量，方向由相机指向3D点
            auto normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            auto normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            if(mbFarPoints && (dist1>=mThFarPoints||dist2>=mThFarPoints)) // MODIFICATION
                continue;
            // ratioDist是不考虑金字塔尺度下的距离比例
            const float ratioDist = dist2/dist1;
            // 金字塔尺度因子的比例
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave]/pKF2->mvScaleFactors[kp2.octave];

            // 距离的比例和图像金字塔的比例不应该差太多，否则就跳过
            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;

            // Triangulation is succesfull
			// Step 6.8：三角化生成3D点成功，构造成MapPoint
            cv::Mat x3D_(x3D);
            MapPoint* pMP = new MapPoint(x3D_,mpCurrentKeyFrame,mpAtlas->GetCurrentMap());

            // Step 6.9：为该MapPoint添加属性：
            // a.观测到该MapPoint的关键帧
            pMP->AddObservation(mpCurrentKeyFrame,idx1);            
            pMP->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP,idx1);
            pKF2->AddMapPoint(pMP,idx2);

            // b.该MapPoint的描述子
            pMP->ComputeDistinctiveDescriptors();

            // c.该MapPoint的平均观测方向和深度范围
            pMP->UpdateNormalAndDepth();

            mpAtlas->AddMapPoint(pMP);
            // Step 6.10：将新产生的点放入检测队列
            // 这些MapPoints都会经过MapPointCulling函数的检验
            mlpRecentAddedMapPoints.push_back(pMP);
        }
    }
}

/**
 * @brief 检查并融合当前关键帧与相邻帧（两级相邻）重复的MapPoints
 * 
 */
void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    // Step 1：获得当前关键帧在共视图中权重排名前nn的邻接关键帧
    // 开始之前先定义几个概念
    // 当前关键帧的邻接关键帧，称为一级相邻关键帧，也就是邻居
    // 与一级相邻关键帧相邻的关键帧，称为二级相邻关键帧，也就是邻居的邻居

    // 单目情况要10个邻接关键帧，双目或者RGBD则要20个
    int nn = 10;
    if(mbMonocular)
        nn=20;

    // 和当前关键帧相邻的关键帧，也就是一级相邻关键帧
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    
    // Step 2：存储一级相邻关键帧及其二级相邻关键帧
    vector<KeyFrame*> vpTargetKFs;
    // 开始对所有候选的一级关键帧展开遍历：
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        // 没有和当前帧进行过融合的操作
        if(pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        // 加入一级相邻关键帧    
        vpTargetKFs.push_back(pKFi);
        // 标记已经加入
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;
    }

    // Add some covisible of covisible
    // Extend to some second neighbors if abort is not requested
	// 以一级相邻关键帧的共视关系最好的5个相邻关键帧 作为二级相邻关键帧
    for(int i=0, imax=vpTargetKFs.size(); i<imax; i++)
    {
        const vector<KeyFrame*> vpSecondNeighKFs = vpTargetKFs[i]->GetBestCovisibilityKeyFrames(20);
        // 遍历得到的二级相邻关键帧
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            // 当然这个二级相邻关键帧要求没有和当前关键帧发生融合,并且这个二级相邻关键帧也不是当前关键帧
            if(pKFi2->isBad() || pKFi2->mnFuseTargetForKF==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            // 存入二级相邻关键帧    
            vpTargetKFs.push_back(pKFi2);
            pKFi2->mnFuseTargetForKF=mpCurrentKeyFrame->mnId;
        }
        if (mbAbortBA)
            break;
    }

    // Extend to temporal neighbors
    if(mbInertial)
    {
        KeyFrame* pKFi = mpCurrentKeyFrame->mPrevKF;
        while(vpTargetKFs.size()<20 && pKFi)
        {
            if(pKFi->isBad() || pKFi->mnFuseTargetForKF==mpCurrentKeyFrame->mnId)
            {
                pKFi = pKFi->mPrevKF;
                continue;
            }
            vpTargetKFs.push_back(pKFi);
            pKFi->mnFuseTargetForKF=mpCurrentKeyFrame->mnId;
            pKFi = pKFi->mPrevKF;
        }
    }

    // Search matches by projection from current KF in target KFs
    // 使用默认参数, 最优和次优比例0.6,匹配时检查特征点的旋转
    ORBmatcher matcher;

    // Step 3：将当前帧的地图点分别与一级二级相邻关键帧地图点进行融合 -- 正向
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        // 将地图点投影到关键帧中进行匹配和融合；融合策略如下
        // 1.如果地图点能匹配关键帧的特征点，并且该点有对应的地图点，那么选择观测数目多的替换两个地图点
        // 2.如果地图点能匹配关键帧的特征点，并且该点没有对应的地图点，那么为该点添加该投影地图点
        // 注意这个时候对地图点融合的操作是立即生效的
        matcher.Fuse(pKFi,vpMapPointMatches);
        if(pKFi->NLeft != -1) matcher.Fuse(pKFi,vpMapPointMatches,true);
    }

    if (mbAbortBA)
        return;

    // Search matches by projection from target KFs in current KF
    // Step 4：将一级二级相邻关键帧地图点分别与当前关键帧地图点进行融合 -- 反向
    // 用于进行存储要融合的一级邻接和二级邻接关键帧所有MapPoints的集合
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapPointMatches.size());
    
    //  Step 4.1：遍历每一个一级邻接和二级邻接关键帧，收集他们的地图点存储到 vpFuseCandidates
    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        // 遍历当前一级邻接和二级邻接关键帧中所有的MapPoints,找出需要进行融合的并且加入到集合中
        for(vector<MapPoint*>::iterator vitMP=vpMapPointsKFi.begin(), vendMP=vpMapPointsKFi.end(); vitMP!=vendMP; vitMP++)
        {
            MapPoint* pMP = *vitMP;
            if(!pMP)
                continue;
            
            // 如果地图点是坏点，或者已经加进集合vpFuseCandidates，跳过
            if(pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;

            // 加入集合，并标记已经加入
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }
    // Step 4.2：进行地图点投影融合,和正向融合操作是完全相同的
    // 不同的是正向操作是"每个关键帧和当前关键帧的地图点进行融合",而这里的是"当前关键帧和所有邻接关键帧的地图点进行融合"
    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);
    if(mpCurrentKeyFrame->NLeft != -1) matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates,true);


    // Update points
    // Step 5：更新当前帧地图点的描述子、深度、观测主方向等属性
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
        MapPoint* pMP=vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                // 在所有找到pMP的关键帧中，获得最佳的描述子
                pMP->ComputeDistinctiveDescriptors();

                // 更新平均观测方向和观测距离
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    // Step 6：更新当前帧的MapPoints后更新与其它帧的连接关系
    // 更新covisibility图
    mpCurrentKeyFrame->UpdateConnections();
}

/**
 * @brief 根据两关键帧的姿态计算两个关键帧之间的基本矩阵
 * @param  pKF1 关键帧1
 * @param  pKF2 关键帧2
 * @return      基本矩阵
 */
cv::Mat LocalMapping::ComputeF12(KeyFrame *&pKF1, KeyFrame *&pKF2)
{
    // 先构造两帧之间的R12,t12
    cv::Mat R1w = pKF1->GetRotation();
    cv::Mat t1w = pKF1->GetTranslation();
    cv::Mat R2w = pKF2->GetRotation();
    cv::Mat t2w = pKF2->GetTranslation();

    cv::Mat R12 = R1w*R2w.t();
    cv::Mat t12 = -R1w*R2w.t()*t2w+t1w;

    // 得到 t12 的反对称矩阵
    cv::Mat t12x = SkewSymmetricMatrix(t12);

    const cv::Mat &K1 = pKF1->mpCamera->toK();
    const cv::Mat &K2 = pKF2->mpCamera->toK();


    return K1.t().inv()*t12x*R12*K2.inv();
}

cv::Matx33f LocalMapping::ComputeF12_(KeyFrame *&pKF1, KeyFrame *&pKF2)
{
    auto R1w = pKF1->GetRotation_();
    auto t1w = pKF1->GetTranslation_();
    auto R2w = pKF2->GetRotation_();
    auto t2w = pKF2->GetTranslation_();

    auto R12 = R1w*R2w.t();
    auto t12 = -R1w*R2w.t()*t2w+t1w;

    auto t12x = SkewSymmetricMatrix_(t12);

    const auto &K1 = pKF1->mpCamera->toK_();
    const auto &K2 = pKF2->mpCamera->toK_();


    return K1.t().inv()*t12x*R12*K2.inv();
}

// 外部线程调用,请求停止当前线程的工作; 其实是回环检测线程调用,来避免在进行全局优化的过程中局部建图线程添加新的关键帧
void LocalMapping::RequestStop()
{
    unique_lock<mutex> lock(mMutexStop);
    mbStopRequested = true;
    unique_lock<mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

// 检查是否要把当前的局部建图线程停止工作,运行的时候要检查是否有终止请求,如果有就执行. 由run函数调用
bool LocalMapping::Stop()
{
    unique_lock<mutex> lock(mMutexStop);
    // 如果当前线程还没有准备停止,但是已经有终止请求了,那么就准备停止当前线程
    if(mbStopRequested && !mbNotStop)
    {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

// 检查mbStopped是否为true，为true表示可以并终止localmapping 线程
bool LocalMapping::isStopped()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
}

// 求外部线程调用，为true，表示外部线程请求停止 local mapping
bool LocalMapping::stopRequested()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopRequested;
}

// 释放当前还在缓冲区中的关键帧指针
void LocalMapping::Release()
{
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);
    if(mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

/**
 * @brief 查看是否接收关键帧，也就是当前线程是否在处理数据，当然tracking线程也不会全看这个值，他会根据队列阻塞情况
 */
bool LocalMapping::AcceptKeyFrames()
{
    unique_lock<mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

// 设置"允许接受关键帧"的状态标志
void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    unique_lock<mutex> lock(mMutexAccept);
    mbAcceptKeyFrames=flag;
}

/**
 * @brief 如果不让它暂停，即使发出了暂停信号也不暂停
 */
bool LocalMapping::SetNotStop(bool flag)
{
    unique_lock<mutex> lock(mMutexStop);

    if(flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

/**
 * @brief 放弃这次BA
 */
void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

/**
 * @brief 检测当前关键帧在共视图中的关键帧，根据地图点在共视图中的冗余程度剔除该共视关键帧
 * 冗余关键帧的判定：90%以上的地图点能被其他关键帧（至少3个）观测到
 */
void LocalMapping::KeyFrameCulling()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
	  // 该函数里变量层层深入，这里列一下：
    // mpCurrentKeyFrame：当前关键帧，本程序就是判断它是否需要删除
    // pKF： mpCurrentKeyFrame的某一个共视关键帧
    // vpMapPoints：pKF对应的所有地图点
    // pMP：vpMapPoints中的某个地图点
    // observations：所有能观测到pMP的关键帧
    // pKFi：observations中的某个关键帧
    // scaleLeveli：pKFi的金字塔尺度
    // scaleLevel：pKF的金字塔尺度
    const int Nd = 21; // MODIFICATION_STEREO_IMU 20 This should be the same than that one from LIBA
    mpCurrentKeyFrame->UpdateBestCovisibles();		// 更新共视关系
    // 1. 根据Covisibility Graph提取当前帧的共视关键帧
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    float redundant_th;
    // 非IMU模式时
    if(!mbInertial)
        redundant_th = 0.9;
    else if (mbMonocular)	// 单目+imu 时
        redundant_th = 0.9;
    else	// 双目+imu时
        redundant_th = 0.5;

    const bool bInitImu = mpAtlas->isImuInitialized();
    int count=0;

    // Compoute last KF from optimizable window:
    unsigned int last_ID;
    if (mbInertial)
    {
        int count = 0;
        KeyFrame* aux_KF = mpCurrentKeyFrame;
        // 找到第前21个关键帧的关键帧id
        while(count<Nd && aux_KF->mPrevKF)
        {
            aux_KF = aux_KF->mPrevKF;
            count++;
        }
        last_ID = aux_KF->mnId;
    }


	// 对所有的共视关键帧进行遍历
    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        count++;
        KeyFrame* pKF = *vit;
        // 如果是地图里第1个关键帧或者是标记为坏帧，则跳过
        if((pKF->mnId==pKF->GetMap()->GetInitKFid()) || pKF->isBad())
            continue;
        // Step 2：提取每个共视关键帧的地图点
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        // 记录某个点被观测次数，后面并未使用
        int nObs = 3;
        // 观测次数阈值，默认为3
        const int thObs=nObs;
        // 记录冗余观测点的数目
        int nRedundantObservations=0;
        int nMPs=0;
        // Step 3：遍历该共视关键帧的所有地图点，判断是否90%以上的地图点能被其它至少3个关键帧（同样或者更低层级）观测到
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        // 对于双目，仅考虑近处（不超过基线的40倍 ）的地图点
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    // pMP->Observations() 是观测到该地图点的相机总数目（单目1，双目2）
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = (pKF -> NLeft == -1) ? pKF->mvKeysUn[i].octave
                                                                     : (i < pKF -> NLeft) ? pKF -> mvKeys[i].octave
                                                                                          : pKF -> mvKeysRight[i].octave;
                        const map<KeyFrame*, tuple<int,int>> observations = pMP->GetObservations();
                        int nObs=0;
                        // 遍历观测到该地图点的关键帧
                        for(map<KeyFrame*, tuple<int,int>>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            tuple<int,int> indexes = mit->second;
                            int leftIndex = get<0>(indexes), rightIndex = get<1>(indexes);
                            int scaleLeveli = -1;
                            if(pKFi -> NLeft == -1)
                                scaleLeveli = pKFi->mvKeysUn[leftIndex].octave;
                            else {
                                if (leftIndex != -1) {
                                    scaleLeveli = pKFi->mvKeys[leftIndex].octave;
                                }
                                if (rightIndex != -1) {
                                    int rightLevel = pKFi->mvKeysRight[rightIndex - pKFi->NLeft].octave;
                                    scaleLeveli = (scaleLeveli == -1 || scaleLeveli > rightLevel) ? rightLevel
                                                                                                  : scaleLeveli;
                                }
                            }
							// 尺度约束：为什么pKF 尺度+1 要大于等于 pKFi 尺度？
							// 回答：因为同样或更低金字塔层级的地图点更准确
                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                // 已经找到3个满足条件的关键帧，就停止不找了
                                if(nObs>thObs)
                                    break;
                            }
                        }
                        // 地图点至少被3个关键帧观测到，就记录为冗余点，更新冗余点计数数目
                        if(nObs>thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }
		// Step 4：该关键帧90%以上的有效地图点被判断为冗余的，则删除该关键帧
        if(nRedundantObservations>redundant_th*nMPs)
        {
            // imu模式下需要更改前后关键帧的连续性，且预积分要叠加起来
            if (mbInertial)
            {
                // 关键帧少于Nd个，跳过不删
                if (mpAtlas->KeyFramesInMap()<=Nd)
                    continue;
                // 关键帧与当前关键帧id差一个，跳过不删
                if(pKF->mnId>(mpCurrentKeyFrame->mnId-2))
                    continue;
                // 关键帧具有前后关键帧
                if(pKF->mPrevKF && pKF->mNextKF)
                {
                    const float t = pKF->mNextKF->mTimeStamp-pKF->mPrevKF->mTimeStamp;
                    // 下面两个括号里的内容一模一样
                    // imu初始化了，且距当前帧的ID超过21，且前后两个关键帧时间间隔小于3s
                    // 或者时间间隔小于0.5s
                    if((bInitImu && (pKF->mnId<last_ID) && t<3.) || (t<0.5))
                    {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    }
                    // 没经过imu初始化的第三阶段，且关键帧与其前一个关键帧的距离小于0.02m，且前后两个关键帧时间间隔小于3s
                    else if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2() && (cv::norm(pKF->GetImuPosition()-pKF->mPrevKF->GetImuPosition())<0.02) && (t<3))
                    {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    }
                }
            }
            // 非imu就没那么多事儿了，直接干掉
            else
            {
                pKF->SetBadFlag();
            }
        }
        // 遍历共视关键帧个数超过一定，就不弄了
        if((count > 20 && mbAbortBA) || count>100) // MODIFICATION originally 20 for mbabortBA check just 10 keyframes
        {
            break;
        }
    }
}

/**
 * @brief 返回反对称矩阵
 * @param v 三维向量
 * @return v的反对称矩阵
 */
cv::Mat LocalMapping::SkewSymmetricMatrix(const cv::Mat &v)
{
    return (cv::Mat_<float>(3,3) <<             0, -v.at<float>(2), v.at<float>(1),
            v.at<float>(2),               0,-v.at<float>(0),
            -v.at<float>(1),  v.at<float>(0),              0);
}

cv::Matx33f LocalMapping::SkewSymmetricMatrix_(const cv::Matx31f &v)
{
    cv::Matx33f skew{0.f, -v(2), v(1),
                     v(2), 0.f, -v(0),
                     -v(1), v(0), 0.f};

    return skew;
}

void LocalMapping::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        cout << "LM: Map reset recieved" << endl;
        mbResetRequested = true;
    }
    cout << "LM: Map reset, waiting..." << endl;
    // 一直等到局部建图线程响应之后才可以退出
    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetRequested)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Map reset, Done!!!" << endl;
}

/**
 * @brief 接收重置当前地图的信号
 */
void LocalMapping::RequestResetActiveMap(Map* pMap)
{
    {
        unique_lock<mutex> lock(mMutexReset);
        cout << "LM: Active map reset recieved" << endl;
        mbResetRequestedActiveMap = true;
        mpMapToReset = pMap;
    }
    cout << "LM: Active map reset, waiting..." << endl;

    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetRequestedActiveMap)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Active map reset, Done!!!" << endl;
}
// 检查是否有复位线程的请求
void LocalMapping::ResetIfRequested()
{
    bool executed_reset = false;
    {
        unique_lock<mutex> lock(mMutexReset);
    	// 执行复位操作:清空关键帧缓冲区,清空待cull的地图点缓冲
        if(mbResetRequested)
        {
            executed_reset = true;

            cout << "LM: Reseting Atlas in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();
        	// 恢复为false表示复位过程完成
            mbResetRequested=false;
            mbResetRequestedActiveMap = false;

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu=false;

            mIdxInit=0;

            cout << "LM: End reseting Local Mapping..." << endl;
        }

        if(mbResetRequestedActiveMap) {
            executed_reset = true;
            cout << "LM: Reseting current map in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu=false;

            mbResetRequestedActiveMap = false;
            cout << "LM: End reseting Local Mapping..." << endl;
        }
    }
    if(executed_reset)
        cout << "LM: Reset free the mutex" << endl;

}
// 请求终止当前线程
void LocalMapping::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

/**
 * @brief 查看完成信号，跳出while循环
 */
bool LocalMapping::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

/**
 * @brief 设置当前线程已经真正地结束了
 */
void LocalMapping::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;    
    unique_lock<mutex> lock2(mMutexStop);
    mbStopped = true;
}

/**
 * @brief 当前线程的run函数是否已经终止
 */
bool LocalMapping::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}

/** 
 * @brief imu初始化
 * @param priorG 陀螺仪偏置的信息矩阵系数，主动设置时一般bInit为true，也就是只优化最后一帧的偏置，这个数会作为计算信息矩阵时使用
 * @param priorA 加速度计偏置的信息矩阵系数
 * @param bFIBA 是否做BA优化
 */
void LocalMapping::InitializeIMU(float priorG, float priorA, bool bFIBA)
{
    // 1. 将所有关键帧放入列表及向量里，且查看是否满足初始化条件
    if (mbResetRequested)
        return;

    float minTime;
    int nMinKF;
    // 从时间及帧数上限制初始化，不满足下面的不进行初始化
    if (mbMonocular)
    {
        minTime = 2.0;
        nMinKF = 10;
    }
    else
    {
        minTime = 1.0;
        nMinKF = 10;
    }

    // 当前地图大于10帧才进行初始化
    if(mpAtlas->KeyFramesInMap()<nMinKF)
        return;

    // Retrieve all keyframe in temporal order
    // 按照顺序存放目前地图里的关键帧，顺序按照前后顺序来，包括当前关键帧
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while(pKF->mPrevKF)
    {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    // 以相同内容再构建一个vector
    vector<KeyFrame*> vpKF(lpKF.begin(),lpKF.end());

    // TODO 跟上面重复？
    if(vpKF.size()<nMinKF)
        return;

    mFirstTs=vpKF.front()->mTimeStamp;
    if(mpCurrentKeyFrame->mTimeStamp-mFirstTs<minTime)
        return;
    
    // 正在做IMU的初始化，在tracking里面使用，如果为true，暂不添加关键帧
    bInitializing = true;

    // 先处理新关键帧，防止堆积且保证数据量充足
    while(CheckNewKeyFrames())
    {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }
    // 2. 正式IMU初始化
    const int N = vpKF.size();
    IMU::Bias b(0,0,0,0,0,0);

    // Compute and KF velocities mRwg estimation
    // 在IMU没有初始化情况下
    if (!mpCurrentKeyFrame->GetMap()->isImuInitialized())
    {
        cv::Mat cvRwg;
        cv::Mat dirG = cv::Mat::zeros(3,1,CV_32F);
        for(vector<KeyFrame*>::iterator itKF = vpKF.begin(); itKF!=vpKF.end(); itKF++)
        {
            if (!(*itKF)->mpImuPreintegrated)
                continue;
            if (!(*itKF)->mPrevKF)
                continue;
            // 初始化时关于速度的预积分定义Ri.t()*(s*Vj - s*Vi - Rwg*g*tij)
            dirG -= (*itKF)->mPrevKF->GetImuRotation()*(*itKF)->mpImuPreintegrated->GetUpdatedDeltaVelocity();
            // 求取实际的速度，位移/时间
            cv::Mat _vel = ((*itKF)->GetImuPosition() - (*itKF)->mPrevKF->GetImuPosition())/(*itKF)->mpImuPreintegrated->dT;
            (*itKF)->SetVelocity(_vel);
            (*itKF)->mPrevKF->SetVelocity(_vel);
        }
        // dirG = sV1 - sVn + n*Rwg*g*t
        // 归一化
        dirG = dirG/cv::norm(dirG);
        // 原本的重力方向
        cv::Mat gI = (cv::Mat_<float>(3,1) << 0.0f, 0.0f, -1.0f);
        // 求速度方向与重力方向的角轴
        cv::Mat v = gI.cross(dirG);
        // 求角轴模长
        const float nv = cv::norm(v);
        // 求转角大小
        const float cosg = gI.dot(dirG);
        const float ang = acos(cosg);
        // 先计算旋转向量，在除去角轴大小
        cv::Mat vzg = v * ang / nv;
        // 获得重力方向到当前速度方向的旋转向量
        cvRwg = IMU::ExpSO3(vzg);
        mRwg = Converter::toMatrix3d(cvRwg);
        mTinit = mpCurrentKeyFrame->mTimeStamp-mFirstTs;
    }
    else
    {
        mRwg = Eigen::Matrix3d::Identity();
        mbg = Converter::toVector3d(mpCurrentKeyFrame->GetGyroBias());
        mba = Converter::toVector3d(mpCurrentKeyFrame->GetAccBias());
    }

    mScale=1.0;

    // 暂时没发现在别的地方出现过
    mInitTime = mpTracker->mLastFrame.mTimeStamp-vpKF.front()->mTimeStamp;

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    // 计算残差及偏置差，优化尺度重力方向及速度偏置，偏置先验为0，双目时不优化尺度
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale, mbg, mba, mbMonocular, infoInertial, false, false, priorG, priorA);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    /*cout << "scale after inertial-only optimization: " << mScale << endl;
    cout << "bg after inertial-only optimization: " << mbg << endl;
    cout << "ba after inertial-only optimization: " << mba << endl;*/

    // 尺度太小的话初始化认为失败
    if (mScale<1e-1)
    {
        cout << "scale too small" << endl;
        bInitializing=false;
        return;
    }

    // 到此时为止，前面做的东西没有改变map

    // Before this line we are not changing the map

    unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    // 尺度变化超过设定值，或者非单目时（无论带不带imu，但这个函数只在带imu时才执行，所以这个可以理解为双目imu）
    if ((fabs(mScale-1.f)>0.00001)||!mbMonocular)
    {
        // 恢复重力方向与尺度信息
        mpAtlas->GetCurrentMap()->ApplyScaledRotation(Converter::toCvMat(mRwg).t(),mScale,true);
        // 更新普通帧的位姿，主要是当前帧与上一帧
        mpTracker->UpdateFrameIMU(mScale,vpKF[0]->GetImuBias(),mpCurrentKeyFrame);
    }
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    // Check if initialization OK
    // 即使初始化成功后面还会执行这个函数重新初始化
    // 在之前没有初始化成功情况下（此时刚刚初始化成功）对每一帧都标记，后面的kf全部都在tracking里面标记为true
    // 也就是初始化之前的那些关键帧即使有imu信息也不算
    if (!mpAtlas->isImuInitialized())
        for(int i=0;i<N;i++)
        {
            KeyFrame* pKF2 = vpKF[i];
            pKF2->bImu = true;
        }

    /*cout << "Before GIBA: " << endl;
    cout << "ba: " << mpCurrentKeyFrame->GetAccBias() << endl;
    cout << "bg: " << mpCurrentKeyFrame->GetGyroBias() << endl;*/

    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    // 代码里都为true
    if (bFIBA)
    {
        // 承接上一步纯imu优化，按照之前的结果更新了尺度信息及适应重力方向，所以要结合地图进行一次视觉加imu的全局优化，这次带了MP等信息
        if (priorA!=0.f)
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, 0, NULL, true, priorG, priorA);
        else
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, 0, NULL, false);
    }

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    // If initialization is OK
    mpTracker->UpdateFrameIMU(1.0,vpKF[0]->GetImuBias(),mpCurrentKeyFrame);
    if (!mpAtlas->isImuInitialized())
    {
        cout << "IMU in Map " << mpAtlas->GetCurrentMap()->GetId() << " is initialized" << endl;
        // ! 重要！标记初始化成功
        mpAtlas->SetImuInitialized();
        mpTracker->t0IMU = mpTracker->mCurrentFrame.mTimeStamp;
        mpCurrentKeyFrame->bImu = true;
    }


    mbNewInit=true;  // 没用到
    mnKFs=vpKF.size();  // 没用到
    mIdxInit++;  // 没用到

    // 里面存放的是来到Local线程，但还未处理的关键帧
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
    {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    mpTracker->mState=Tracking::OK;
    bInitializing = false;


    /*cout << "After GIBA: " << endl;
    cout << "ba: " << mpCurrentKeyFrame->GetAccBias() << endl;
    cout << "bg: " << mpCurrentKeyFrame->GetGyroBias() << endl;
    double t_inertial_only = std::chrono::duration_cast<std::chrono::duration<double> >(t1 - t0).count();
    double t_update = std::chrono::duration_cast<std::chrono::duration<double> >(t3 - t2).count();
    double t_viba = std::chrono::duration_cast<std::chrono::duration<double> >(t5 - t4).count();
    cout << t_inertial_only << ", " << t_update << ", " << t_viba << endl;*/

    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}

/**
 * @brief 通过BA优化进行尺度更新，关键帧小于100，使用了所有关键帧的信息，但只优化尺度和重力方向。每10s在这里的时间段内时多次进行尺度更新
 */
void LocalMapping::ScaleRefinement()
{
    // Minimum number of keyframes to compute a solution
    // Minimum time (seconds) between first and last keyframe to compute a solution. Make the difference between monocular and stereo
    // unique_lock<mutex> lock0(mMutexImuInit);
    if (mbResetRequested)
        return;

    // Retrieve all keyframes in temporal order
    // 1. 检索所有的关键帧（当前地图）
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while(pKF->mPrevKF)
    {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    vector<KeyFrame*> vpKF(lpKF.begin(),lpKF.end());
    // 加入新添加的帧
    while(CheckNewKeyFrames())
    {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    const int N = vpKF.size();
    // 2. 更新旋转与尺度

    // 待优化变量的初值
    mRwg = Eigen::Matrix3d::Identity();
    mScale=1.0;

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    // 优化重力方向与尺度
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (mScale<1e-1) // 1e-1
    {
        cout << "scale too small" << endl;
        bInitializing=false;
        return;
    }

    // Before this line we are not changing the map
    // 3. 开始更新地图
    unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    // 3.1 如果尺度更新较多，或是在双目imu情况下更新地图
    if ((fabs(mScale-1.f)>0.00001)||!mbMonocular)
    {
        mpAtlas->GetCurrentMap()->ApplyScaledRotation(Converter::toCvMat(mRwg).t(),mScale,true);
        mpTracker->UpdateFrameIMU(mScale,mpCurrentKeyFrame->GetImuBias(),mpCurrentKeyFrame);
    }
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    // 3.2 优化的这段时间新进来的kf全部清空不要
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
    {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    double t_inertial_only = std::chrono::duration_cast<std::chrono::duration<double> >(t1 - t0).count();

    // To perform pose-inertial opt w.r.t. last keyframe
    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}

/**
 * @brief 返回是否正在做IMU的初始化，在tracking里面使用，如果为true，暂不添加关键帧
 */
bool LocalMapping::IsInitializing()
{
    return bInitializing;
}

/**
 * @brief 获取当前关键帧的时间戳，System::GetTimeFromIMUInit()中调用
 */
double LocalMapping::GetCurrKFTime()
{

    if (mpCurrentKeyFrame)
    {
        return mpCurrentKeyFrame->mTimeStamp;
    }
    else
        return 0.0;
}

/**
 * @brief 获取当前关键帧
 */
KeyFrame *LocalMapping::GetCurrKF()
{
    return mpCurrentKeyFrame;
}

} // namespace ORB_SLAM3
