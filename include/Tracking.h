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


#ifndef TRACKING_H
#define TRACKING_H

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>
#include <opencv2/video/tracking.hpp>

#include"Viewer.h"
#include"FrameDrawer.h"
#include"Atlas.h"
#include"LocalMapping.h"
#include"LoopClosing.h"
#include"Frame.h"
#include "ORBVocabulary.h"
#include"KeyFrameDatabase.h"
#include"ORBextractor.h"
#include "Initializer.h"
#include "MapDrawer.h"
#include "System.h"
#include "ImuTypes.h"

#include "GeometricCamera.h"

#include <mutex>
#include <unordered_set>

namespace ORB_SLAM3
{

class Viewer;
class FrameDrawer;
class Atlas;
class LocalMapping;
class LoopClosing;
class System;

class Tracking
{  

public:
    Tracking(System* pSys, ORBVocabulary* pVoc, FrameDrawer* pFrameDrawer, MapDrawer* pMapDrawer, Atlas* pAtlas,
             KeyFrameDatabase* pKFDB, const string &strSettingPath, const int sensor, const string &_nameSeq=std::string());

    ~Tracking();

    // Parse the config file
    // 提取配置文件数据
    bool ParseCamParamFile(cv::FileStorage &fSettings);
    bool ParseORBParamFile(cv::FileStorage &fSettings);
    bool ParseIMUParamFile(cv::FileStorage &fSettings);

    // Preprocess the input and call Track(). Extract features and performs stereo matching.
    // 输入图像输出位姿Tcw
    cv::Mat GrabImageStereo(const cv::Mat &imRectLeft,const cv::Mat &imRectRight, const double &timestamp, string filename);
    cv::Mat GrabImageRGBD(const cv::Mat &imRGB,const cv::Mat &imD, const double &timestamp, string filename);
    cv::Mat GrabImageMonocular(const cv::Mat &im, const double &timestamp, string filename);
    // cv::Mat GrabImageImuMonocular(const cv::Mat &im, const double &timestamp);
    // 放置IMU数据
    void GrabImuData(const IMU::Point &imuMeasurement);
    // 设置线程指针
    void SetLocalMapper(LocalMapping* pLocalMapper);
    void SetLoopClosing(LoopClosing* pLoopClosing);
    void SetViewer(Viewer* pViewer);
    void SetStepByStep(bool bSet);

    // Load new settings
    // The focal lenght should be similar or scale prediction will fail when projecting points
    // 更换新的标定参数，未使用
    void ChangeCalibration(const string &strSettingPath);

    // Use this function if you have deactivated local mapping and you only want to localize the camera.
    // 设置是否仅定位模式还是SLAM模式
    void InformOnlyTracking(const bool &flag);
    // localmapping中更新了关键帧的位姿后，更新普通帧的位姿，通过IMU积分更新速度。localmapping中初始化imu中使用
    void UpdateFrameIMU(const float s, const IMU::Bias &b, KeyFrame* pCurrentKeyFrame);
    KeyFrame* GetLastKeyFrame()
    {
        return mpLastKeyFrame;
    }
    // 新建地图
    void CreateMapInAtlas();
    std::mutex mMutexTracks;

    // 更新数据集
    void NewDataset();
    // 获得数据集总数
    int GetNumberDataset();
    // 获取匹配内点总数
    int GetMatchesInliers();
public:

    // Tracking states
    // 跟踪状态
    enum eTrackingState{
        SYSTEM_NOT_READY=-1,            //系统没有准备好的状态,一般就是在启动后加载配置文件和词典文件时候的状态
        NO_IMAGES_YET=0,                //当前无图像
        NOT_INITIALIZED=1,              //有图像但是没有完成初始化
        OK=2,                           //正常跟踪状态
        RECENTLY_LOST=3,                //IMU模式：当前地图中的KF>10,且丢失时间<5秒。纯视觉模式：没有该状态
        LOST=4,                         //IMU模式：当前帧跟丢超过5s。纯视觉模式：重定位失败
        OK_KLT=5                        //未使用
    };

    eTrackingState mState;
    eTrackingState mLastProcessedState;

    // Input sensor
    int mSensor;

    // Current Frame
    Frame mCurrentFrame;
    Frame mLastFrame; //跟踪成功后，保存当前帧数据

    cv::Mat mImGray;

    // Initialization Variables (Monocular)
    std::vector<int> mvIniLastMatches;
    std::vector<int> mvIniMatches;
    std::vector<cv::Point2f> mvbPrevMatched;
    std::vector<cv::Point3f> mvIniP3D;
    Frame mInitialFrame;

    // Lists used to recover the full camera trajectory at the end of the execution.
    // Basically we store the reference keyframe for each frame and its relative transformation
    // 代码结束后保存位姿用的列表
    list<cv::Mat> mlRelativeFramePoses;
    list<KeyFrame*> mlpReferences;
    list<double> mlFrameTimes;
    list<bool> mlbLost;

    // frames with estimated pose
    int mTrackedFr;
    bool mbStep;

    // True if local mapping is deactivated and we are performing only localization
    // true表示仅定位模式，此时局部建图线程和闭环线程关闭
    bool mbOnlyTracking;

    void Reset(bool bLocMap = false);
    void ResetActiveMap(bool bLocMap = false);

    float mMeanTrack;
    bool mbInitWith3KFs;
    double t0; // time-stamp of first read frame
    double t0vis; // time-stamp of first inserted keyframe
    double t0IMU; // time-stamp of IMU initialization

    // 获取局部地图点
    vector<MapPoint*> GetLocalMapMPS();

    bool mbWriteStats;

#ifdef REGISTER_TIMES
    void LocalMapStats2File();
    void TrackStats2File();
    void PrintTimeStats();

    vector<double> vdRectStereo_ms;
    vector<double> vdORBExtract_ms;
    vector<double> vdStereoMatch_ms;
    vector<double> vdIMUInteg_ms;
    vector<double> vdPosePred_ms;
    vector<double> vdLMTrack_ms;
    vector<double> vdNewKF_ms;
    vector<double> vdTrackTotal_ms;

    vector<double> vdUpdatedLM_ms;
    vector<double> vdSearchLP_ms;
    vector<double> vdPoseOpt_ms;
#endif

    vector<int> vnKeyFramesLM;
    vector<int> vnMapPointsLM;

protected:

    // Main tracking function. It is independent of the input sensor.
    void Track();   //主要的跟踪函数

    // Map initialization for stereo and RGB-D
    void StereoInitialization();    //双目初始化

    // Map initialization for monocular
    void MonocularInitialization(); //单目初始化
    void CreateNewMapPoints();      //创建新的地图点
    cv::Mat ComputeF12(KeyFrame *&pKF1, KeyFrame *&pKF2);
    void CreateInitialMapMonocular();   //单目模式下创建初始化地图

    void CheckReplacedInLastFrame();//检查上一帧中的地图点是否需要被替换
    bool TrackReferenceKeyFrame();  // 参考关键帧跟踪
    void UpdateLastFrame();         //更新上一帧位姿，在上一帧中生成临时地图点
    bool TrackWithMotionModel();    //恒速模型跟踪
    bool PredictStateIMU();         //用IMU预测位姿

    bool Relocalization();          //重定位

    void UpdateLocalMap();          //更新局部地图
    void UpdateLocalPoints();       //更新局部地图点
    void UpdateLocalKeyFrames();    //更新局部地图里的关键帧

    bool TrackLocalMap();           //局部地图跟踪
    bool TrackLocalMap_old();       //未定义
    void SearchLocalPoints();       //搜索局部地图点

    bool NeedNewKeyFrame();         //判断是否需要插入关键帧
    void CreateNewKeyFrame();       //创建关键帧

    // Perform preintegration from last frame
    void PreintegrateIMU();         //IMU预积分

    // Reset IMU biases and compute frame velocity
    // 重置并重新计算IMU偏置。未使用
    void ComputeGyroBias(const vector<Frame*> &vpFs, float &bwx,  float &bwy, float &bwz);
    void ComputeVelocitiesAccBias(const vector<Frame*> &vpFs, float &bax,  float &bay, float &baz);


    bool mbMapUpdated;

    // Imu preintegration from last frame
    IMU::Preintegrated *mpImuPreintegratedFromLastKF;

    // Queue of IMU measurements between frames
    std::list<IMU::Point> mlQueueImuData;   //存放两帧之间的IMU数据

    // Vector of IMU measurements from previous to current frame (to be filled by PreintegrateIMU)
    std::vector<IMU::Point> mvImuFromLastFrame;
    std::mutex mMutexImuQueue;

    // Imu calibration parameters
    IMU::Calib *mpImuCalib;     //IMU标定参数

    // Last Bias Estimation (at keyframe creation)
    IMU::Bias mLastBias;

    // In case of performing only localization, this flag is true when there are no matches to
    // points in the map. Still tracking will continue if there are enough matches with temporal points.
    // In that case we are doing visual odometry. The system will try to do relocalization to recover
    // "zero-drift" localization to the map.
    bool mbVO;

    //Other Thread Pointers
    // 其他线程的指针
    LocalMapping* mpLocalMapper;
    LoopClosing* mpLoopClosing;

    //ORB特征提取器
    ORBextractor* mpORBextractorLeft, *mpORBextractorRight;
    ORBextractor* mpIniORBextractor;

    //BoW
    ORBVocabulary* mpORBVocabulary;
    KeyFrameDatabase* mpKeyFrameDB;

    // Initalization (only for monocular)
    Initializer* mpInitializer;
    bool mbSetInit;

    //Local Map
    KeyFrame* mpReferenceKF;
    std::vector<KeyFrame*> mvpLocalKeyFrames;
    std::vector<MapPoint*> mvpLocalMapPoints;
    
    // System
    System* mpSystem;
    
    //Drawers
    Viewer* mpViewer;
    FrameDrawer* mpFrameDrawer;
    MapDrawer* mpMapDrawer;
    bool bStepByStep;

    //Atlas
    Atlas* mpAtlas;

    //Calibration matrix
    cv::Mat mK;
    cv::Mat mDistCoef;
    float mbf;

    //New KeyFrame rules (according to fps)
    int mMinFrames;
    int mMaxFrames;

    int mnFirstImuFrameId;
    // 经过多少帧后可以重置IMU，一般设置为和帧率相同，对应的时间是1s
    int mnFramesToResetIMU;

    // Threshold close/far points
    // Points seen as close by the stereo/RGBD sensor are considered reliable
    // and inserted from just one frame. Far points requiere a match in two keyframes.
    float mThDepth;             //近远点阈值，基线的倍数

    // For RGB-D inputs only. For some datasets (e.g. TUM) the depthmap values are scaled.
    float mDepthMapFactor;      // RGB-D尺度缩放因子

    //Current matches in frame
    int mnMatchesInliers;       //当前帧匹配内点数

    //Last Frame, KeyFrame and Relocalisation Info
    KeyFrame* mpLastKeyFrame;
    unsigned int mnLastKeyFrameId;
    unsigned int mnLastRelocFrameId;
    double mTimeStampLost;
    double time_recently_lost;
    double time_recently_lost_visual;

    unsigned int mnFirstFrameId;
    unsigned int mnInitialFrameId;
    unsigned int mnLastInitFrameId;

    bool mbCreatedMap;


    //Motion Model
    cv::Mat mVelocity;      // 恒速模型的速度。通过位姿增量获得或者IMU积分得到

    //Color order (true RGB, false BGR, ignored if grayscale)
    bool mbRGB;

    list<MapPoint*> mlpTemporalPoints;  // 临时地图点

    //int nMapChangeIndex;

    int mnNumDataset;

    ofstream f_track_stats;

    ofstream f_track_times;
    double mTime_PreIntIMU;
    double mTime_PosePred;
    double mTime_LocalMapTrack;
    double mTime_NewKF_Dec;

    GeometricCamera* mpCamera, *mpCamera2;  //相机类

    int initID, lastID;

    cv::Mat mTlr;

public:
    cv::Mat mImRight;
};

} //namespace ORB_SLAM

#endif // TRACKING_H
