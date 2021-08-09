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


#include "LoopClosing.h"

#include "Sim3Solver.h"
#include "Converter.h"
#include "Optimizer.h"
#include "ORBmatcher.h"
#include "G2oTypes.h"

#include<mutex>
#include<thread>


namespace ORB_SLAM3
{

LoopClosing::LoopClosing(Atlas *pAtlas, KeyFrameDatabase *pDB, ORBVocabulary *pVoc, const bool bFixScale):
    mbResetRequested(false), mbResetActiveMapRequested(false), mbFinishRequested(false), mbFinished(true), mpAtlas(pAtlas),
    mpKeyFrameDB(pDB), mpORBVocabulary(pVoc), mpMatchedKF(NULL), mLastLoopKFid(0), mbRunningGBA(false), mbFinishedGBA(true),
    mbStopGBA(false), mpThreadGBA(NULL), mbFixScale(bFixScale), mnFullBAIdx(0), mnLoopNumCoincidences(0), mnMergeNumCoincidences(0),
    mbLoopDetected(false), mbMergeDetected(false), mnLoopNumNotFound(0), mnMergeNumNotFound(0)
{
    // 连续性阈值
    mnCovisibilityConsistencyTh = 3;
    mpLastCurrentKF = static_cast<KeyFrame*>(NULL);
}

// 设置追踪线程句柄
void LoopClosing::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}
// 设置局部建图线程的句柄
void LoopClosing::SetLocalMapper(LocalMapping *pLocalMapper)
{
    mpLocalMapper=pLocalMapper;
}

// 回环线程主函数
void LoopClosing::Run()
{
    mbFinished =false;

    // 线程主循环
    while(1)
    {
        // Check if there are keyframes in the queue
        // Loopclosing中的关键帧是LocalMapping发送过来的，LocalMapping是Tracking中发过来的
        // 在LocalMapping中通过 InsertKeyFrame 将关键帧插入闭环检测队列mlpLoopKeyFrameQueue
        // Step 1 查看闭环检测队列mlpLoopKeyFrameQueue中有没有关键帧进来
        if(CheckNewKeyFrames())
        {
            if(mpLastCurrentKF)
            {
                mpLastCurrentKF->mvpLoopCandKFs.clear();
                mpLastCurrentKF->mvpMergeCandKFs.clear();
            }
#ifdef REGISTER_TIMES
            timeDetectBoW = 0;
            std::chrono::steady_clock::time_point time_StartDetectBoW = std::chrono::steady_clock::now();
#endif
            bool bDetected = NewDetectCommonRegions();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndDetectBoW = std::chrono::steady_clock::now();
            double timeDetect = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndDetectBoW - time_StartDetectBoW).count();
            double timeDetectSE3 = timeDetect - timeDetectBoW;

            if(timeDetectBoW > 0)
            {
                vTimeBoW_ms.push_back(timeDetectBoW);
            }
            vTimeSE3_ms.push_back(timeDetectSE3);
            vTimePRTotal_ms.push_back(timeDetect);
#endif

            if(bDetected)
            {
                if(mbMergeDetected)
                {
                    if ((mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO) &&
                        (!mpCurrentKF->GetMap()->isImuInitialized()))
                    {
                        cout << "IMU is not initilized, merge is aborted" << endl;
                    }
                    else
                    {
                        Verbose::PrintMess("*Merged detected", Verbose::VERBOSITY_QUIET);
                        Verbose::PrintMess("Number of KFs in the current map: " + to_string(mpCurrentKF->GetMap()->KeyFramesInMap()), Verbose::VERBOSITY_DEBUG);
                        cv::Mat mTmw = mpMergeMatchedKF->GetPose();
                        g2o::Sim3 gSmw2(Converter::toMatrix3d(mTmw.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTmw.rowRange(0, 3).col(3)),1.0);
                        cv::Mat mTcw = mpCurrentKF->GetPose();
                        g2o::Sim3 gScw1(Converter::toMatrix3d(mTcw.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTcw.rowRange(0, 3).col(3)),1.0);
                        g2o::Sim3 gSw2c = mg2oMergeSlw.inverse();
                        g2o::Sim3 gSw1m = mg2oMergeSlw;

                        mSold_new = (gSw2c * gScw1);

                        if(mpCurrentKF->GetMap()->IsInertial() && mpMergeMatchedKF->GetMap()->IsInertial())
                        {
                            if(mSold_new.scale()<0.90||mSold_new.scale()>1.1){
                                mpMergeLastCurrentKF->SetErase();
                                mpMergeMatchedKF->SetErase();
                                mnMergeNumCoincidences = 0;
                                mvpMergeMatchedMPs.clear();
                                mvpMergeMPs.clear();
                                mnMergeNumNotFound = 0;
                                mbMergeDetected = false;
                                Verbose::PrintMess("scale bad estimated. Abort merging", Verbose::VERBOSITY_NORMAL);
                                continue;
                            }
                            // If inertial, force only yaw
                            if ((mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO) &&
                                   mpCurrentKF->GetMap()->GetIniertialBA1()) // TODO, maybe with GetIniertialBA1
                            {
                                Eigen::Vector3d phi = LogSO3(mSold_new.rotation().toRotationMatrix());
                                phi(0)=0;
                                phi(1)=0;
                                mSold_new = g2o::Sim3(ExpSO3(phi),mSold_new.translation(),1.0);
                            }
                        }

                        mg2oMergeSmw = gSmw2 * gSw2c * gScw1;

                        mg2oMergeScw = mg2oMergeSlw;

#ifdef REGISTER_TIMES
                        std::chrono::steady_clock::time_point time_StartMerge = std::chrono::steady_clock::now();
#endif
                        if (mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO)
                            MergeLocal2();
                        else
                            MergeLocal();
#ifdef REGISTER_TIMES
                        std::chrono::steady_clock::time_point time_EndMerge = std::chrono::steady_clock::now();
                        double timeMerge = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMerge - time_StartMerge).count();
                        vTimeMergeTotal_ms.push_back(timeMerge);
#endif
                    }

                    vdPR_CurrentTime.push_back(mpCurrentKF->mTimeStamp);
                    vdPR_MatchedTime.push_back(mpMergeMatchedKF->mTimeStamp);
                    vnPR_TypeRecogn.push_back(1);

                    // Reset all variables
                    mpMergeLastCurrentKF->SetErase();
                    mpMergeMatchedKF->SetErase();
                    mnMergeNumCoincidences = 0;
                    mvpMergeMatchedMPs.clear();
                    mvpMergeMPs.clear();
                    mnMergeNumNotFound = 0;
                    mbMergeDetected = false;

                    if(mbLoopDetected)
                    {
                        // Reset Loop variables
                        mpLoopLastCurrentKF->SetErase();
                        mpLoopMatchedKF->SetErase();
                        mnLoopNumCoincidences = 0;
                        mvpLoopMatchedMPs.clear();
                        mvpLoopMPs.clear();
                        mnLoopNumNotFound = 0;
                        mbLoopDetected = false;
                    }

                }

                if(mbLoopDetected)
                {
                    vdPR_CurrentTime.push_back(mpCurrentKF->mTimeStamp);
                    vdPR_MatchedTime.push_back(mpLoopMatchedKF->mTimeStamp);
                    vnPR_TypeRecogn.push_back(0);


                    Verbose::PrintMess("*Loop detected", Verbose::VERBOSITY_QUIET);

                    mg2oLoopScw = mg2oLoopSlw; //*mvg2oSim3LoopTcw[nCurrentIndex];
                    if(mpCurrentKF->GetMap()->IsInertial())
                    {
                        cv::Mat Twc = mpCurrentKF->GetPoseInverse();
                        g2o::Sim3 g2oTwc(Converter::toMatrix3d(Twc.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(Twc.rowRange(0, 3).col(3)),1.0);
                        g2o::Sim3 g2oSww_new = g2oTwc*mg2oLoopScw;

                        Eigen::Vector3d phi = LogSO3(g2oSww_new.rotation().toRotationMatrix());
                        //cout << "tw2w1: " << g2oSww_new.translation() << endl;
                        //cout << "Rw2w1: " << g2oSww_new.rotation().toRotationMatrix() << endl;
                        //cout << "Angle Rw2w1: " << 180*phi/3.14 << endl;
                        //cout << "scale w2w1: " << g2oSww_new.scale() << endl;

                        if (fabs(phi(0))<0.008f && fabs(phi(1))<0.008f && fabs(phi(2))<0.349f)
                        {
                            if(mpCurrentKF->GetMap()->IsInertial())
                            {
                                // If inertial, force only yaw
                                if ((mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO) &&
                                        mpCurrentKF->GetMap()->GetIniertialBA2()) // TODO, maybe with GetIniertialBA1
                                {
                                    phi(0)=0;
                                    phi(1)=0;
                                    g2oSww_new = g2o::Sim3(ExpSO3(phi),g2oSww_new.translation(),1.0);
                                    mg2oLoopScw = g2oTwc.inverse()*g2oSww_new;
                                }
                            }

                            mvpLoopMapPoints = mvpLoopMPs;//*mvvpLoopMapPoints[nCurrentIndex];

#ifdef REGISTER_TIMES
                            std::chrono::steady_clock::time_point time_StartLoop = std::chrono::steady_clock::now();
#endif
                            CorrectLoop();
#ifdef REGISTER_TIMES
                            std::chrono::steady_clock::time_point time_EndLoop = std::chrono::steady_clock::now();
                            double timeLoop = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLoop - time_StartLoop).count();
                            vTimeLoopTotal_ms.push_back(timeLoop);
#endif
                        }
                        else
                        {
                            cout << "BAD LOOP!!!" << endl;
                        }
                    }
                    else
                    {
                        mvpLoopMapPoints = mvpLoopMPs;
#ifdef REGISTER_TIMES
                        std::chrono::steady_clock::time_point time_StartLoop = std::chrono::steady_clock::now();
#endif
                        CorrectLoop();

#ifdef REGISTER_TIMES
                        std::chrono::steady_clock::time_point time_EndLoop = std::chrono::steady_clock::now();
                        double timeLoop = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLoop - time_StartLoop).count();
                        vTimeLoopTotal_ms.push_back(timeLoop);
#endif
                    }

                    // Reset all variables
                    mpLoopLastCurrentKF->SetErase();
                    mpLoopMatchedKF->SetErase();
                    mnLoopNumCoincidences = 0;
                    mvpLoopMatchedMPs.clear();
                    mvpLoopMPs.clear();
                    mnLoopNumNotFound = 0;
                    mbLoopDetected = false;
                }

            }
            mpLastCurrentKF = mpCurrentKF;
        }
		// 查看是否有外部线程请求复位当前线程
        ResetIfRequested();

        // 查看外部线程是否有终止当前线程的请求,如果有的话就跳出这个线程的主函数的主循环
        if(CheckFinish()){
            // cout << "LC: Finish requested" << endl;
            break;
        }

        usleep(5000);
    }

    //ofstream f_stats;
    //f_stats.open("PlaceRecognition_stats" + mpLocalMapper->strSequence + ".txt");
    //f_stats << "# current_timestamp, matched_timestamp, [0:Loop, 1:Merge]" << endl;
    //f_stats << fixed;
    //for(int i=0; i< vdPR_CurrentTime.size(); ++i)
    //{
    //    f_stats  << 1e9*vdPR_CurrentTime[i] << "," << 1e9*vdPR_MatchedTime[i] << "," << vnPR_TypeRecogn[i] << endl;
    //}

    //f_stats.close();

    SetFinish();
}

// 将某个关键帧加入到回环检测的过程中,由局部建图线程调用
void LoopClosing::InsertKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexLoopQueue);
    // 注意：这里第0个关键帧不能够参与到回环检测的过程中,因为第0关键帧定义了整个地图的世界坐标系
    if(pKF->mnId!=0)
        mlpLoopKeyFrameQueue.push_back(pKF);
}

/*
 * 查看列表中是否有等待被插入的关键帧
 * @return 如果存在，返回true
 */
bool LoopClosing::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexLoopQueue);
    return(!mlpLoopKeyFrameQueue.empty());
}

bool LoopClosing::NewDetectCommonRegions()
{
    {
        // Step 1 从队列中取出一个关键帧,作为当前检测闭环关键帧
        unique_lock<mutex> lock(mMutexLoopQueue);
        // 从队列头开始取，也就是先取早进来的关键帧
        mpCurrentKF = mlpLoopKeyFrameQueue.front();
        // 取出关键帧后从队列里弹出该关键帧
        mlpLoopKeyFrameQueue.pop_front();
        // Avoid that a keyframe can be erased while it is being process by this thread
        // 设置当前关键帧不要在优化的过程中被删除
        mpCurrentKF->SetNotErase();
        mpCurrentKF->mbCurrentPlaceRecognition = true;

        mpLastMap = mpCurrentKF->GetMap();
    }

    if(mpLastMap->IsInertial() && !mpLastMap->GetIniertialBA1())
    {
        mpKeyFrameDB->add(mpCurrentKF);
        mpCurrentKF->SetErase();
        return false;
    }

    if(mpTracker->mSensor == System::STEREO && mpLastMap->GetAllKeyFrames().size() < 5) //12
    {
        mpKeyFrameDB->add(mpCurrentKF);
        mpCurrentKF->SetErase();
        return false;
    }

    if(mpLastMap->GetAllKeyFrames().size() < 12)
    {
        mpKeyFrameDB->add(mpCurrentKF);
        mpCurrentKF->SetErase();
        return false;
    }

    //Check the last candidates with geometric validation
    // Loop candidates
    bool bLoopDetectedInKF = false;
    bool bCheckSpatial = false;

    if(mnLoopNumCoincidences > 0)
    {
        bCheckSpatial = true;
        // Find from the last KF candidates
        cv::Mat mTcl = mpCurrentKF->GetPose() * mpLoopLastCurrentKF->GetPoseInverse();
        g2o::Sim3 gScl(Converter::toMatrix3d(mTcl.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTcl.rowRange(0, 3).col(3)),1.0);
        g2o::Sim3 gScw = gScl * mg2oLoopSlw;
        int numProjMatches = 0;
        vector<MapPoint*> vpMatchedMPs;
        bool bCommonRegion = DetectAndReffineSim3FromLastKF(mpCurrentKF, mpLoopMatchedKF, gScw, numProjMatches, mvpLoopMPs, vpMatchedMPs);
        if(bCommonRegion)
        {

            bLoopDetectedInKF = true;

            mnLoopNumCoincidences++;
            mpLoopLastCurrentKF->SetErase();
            mpLoopLastCurrentKF = mpCurrentKF;
            mg2oLoopSlw = gScw;
            mvpLoopMatchedMPs = vpMatchedMPs;


            mbLoopDetected = mnLoopNumCoincidences >= 3;
            mnLoopNumNotFound = 0;

            if(!mbLoopDetected)
            {
                //f_succes_pr << mpCurrentKF->mNameFile << " " << "8"<< endl;
                //f_succes_pr << "% Number of spatial consensous: " << std::to_string(mnLoopNumCoincidences) << endl;
                cout << "PR: Loop detected with Reffine Sim3" << endl;
            }
        }
        else
        {
            bLoopDetectedInKF = false;
            /*f_succes_pr << mpCurrentKF->mNameFile << " " << "8"<< endl;
            f_succes_pr << "% Number of spatial consensous: " << std::to_string(mnLoopNumCoincidences) << endl;*/

            mnLoopNumNotFound++;
            if(mnLoopNumNotFound >= 2)
            {
                /*for(int i=0; i<mvpLoopLastKF.size(); ++i)
                {
                    mvpLoopLastKF[i]->SetErase();
                    mvpLoopCandidateKF[i]->SetErase();
                    mvpLoopLastKF[i]->mbCurrentPlaceRecognition = true;
                }

                mvpLoopCandidateKF.clear();
                mvpLoopLastKF.clear();
                mvg2oSim3LoopTcw.clear();
                mvnLoopNumMatches.clear();
                mvvpLoopMapPoints.clear();
                mvvpLoopMatchedMapPoints.clear();*/

                mpLoopLastCurrentKF->SetErase();
                mpLoopMatchedKF->SetErase();
                //mg2oLoopScw;
                mnLoopNumCoincidences = 0;
                mvpLoopMatchedMPs.clear();
                mvpLoopMPs.clear();
                mnLoopNumNotFound = 0;
            }

        }
    }

    //Merge candidates
    bool bMergeDetectedInKF = false;
    if(mnMergeNumCoincidences > 0)
    {
        // Find from the last KF candidates

        cv::Mat mTcl = mpCurrentKF->GetPose() * mpMergeLastCurrentKF->GetPoseInverse();
        g2o::Sim3 gScl(Converter::toMatrix3d(mTcl.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTcl.rowRange(0, 3).col(3)),1.0);
        g2o::Sim3 gScw = gScl * mg2oMergeSlw;
        int numProjMatches = 0;
        vector<MapPoint*> vpMatchedMPs;
        bool bCommonRegion = DetectAndReffineSim3FromLastKF(mpCurrentKF, mpMergeMatchedKF, gScw, numProjMatches, mvpMergeMPs, vpMatchedMPs);
        if(bCommonRegion)
        {
            //cout << "BoW: Merge reffined Sim3 transformation suscelful" << endl;
            bMergeDetectedInKF = true;

            mnMergeNumCoincidences++;
            mpMergeLastCurrentKF->SetErase();
            mpMergeLastCurrentKF = mpCurrentKF;
            mg2oMergeSlw = gScw;
            mvpMergeMatchedMPs = vpMatchedMPs;

            mbMergeDetected = mnMergeNumCoincidences >= 3;
        }
        else
        {
            //cout << "BoW: Merge reffined Sim3 transformation failed" << endl;
            mbMergeDetected = false;
            bMergeDetectedInKF = false;

            mnMergeNumNotFound++;
            if(mnMergeNumNotFound >= 2)
            {
                /*cout << "+++++++Merge detected failed in KF" << endl;

                for(int i=0; i<mvpMergeLastKF.size(); ++i)
                {
                    mvpMergeLastKF[i]->SetErase();
                    mvpMergeCandidateKF[i]->SetErase();
                    mvpMergeLastKF[i]->mbCurrentPlaceRecognition = true;
                }

                mvpMergeCandidateKF.clear();
                mvpMergeLastKF.clear();
                mvg2oSim3MergeTcw.clear();
                mvnMergeNumMatches.clear();
                mvvpMergeMapPoints.clear();
                mvvpMergeMatchedMapPoints.clear();*/

                mpMergeLastCurrentKF->SetErase();
                mpMergeMatchedKF->SetErase();
                mnMergeNumCoincidences = 0;
                mvpMergeMatchedMPs.clear();
                mvpMergeMPs.clear();
                mnMergeNumNotFound = 0;
            }


        }
    }

    if(mbMergeDetected || mbLoopDetected)
    {
        //f_time_pr << "Geo" << " " << timeGeoKF_ms.count() << endl;
        mpKeyFrameDB->add(mpCurrentKF);
        return true;
    }

    //-------------
    //TODO: This is only necessary if we use a minimun score for pick the best candidates
    const vector<KeyFrame*> vpConnectedKeyFrames = mpCurrentKF->GetVectorCovisibleKeyFrames();
    const DBoW2::BowVector &CurrentBowVec = mpCurrentKF->mBowVec;
    /*float minScore = 1;
    for(size_t i=0; i<vpConnectedKeyFrames.size(); i++)
    {
        KeyFrame* pKF = vpConnectedKeyFrames[i];
        if(pKF->isBad())
            continue;
        const DBoW2::BowVector &BowVec = pKF->mBowVec;

        float score = mpORBVocabulary->score(CurrentBowVec, BowVec);

        if(score<minScore)
            minScore = score;
    }*/
    //-------------

    // Extract candidates from the bag of words
    vector<KeyFrame*> vpMergeBowCand, vpLoopBowCand;
    //cout << "LC: bMergeDetectedInKF: " << bMergeDetectedInKF << "   bLoopDetectedInKF: " << bLoopDetectedInKF << endl;
    if(!bMergeDetectedInKF || !bLoopDetectedInKF)
    {
#ifdef REGISTER_TIMES
        std::chrono::steady_clock::time_point time_StartDetectBoW = std::chrono::steady_clock::now();
#endif
        // Search in BoW
        mpKeyFrameDB->DetectNBestCandidates(mpCurrentKF, vpLoopBowCand, vpMergeBowCand,3);
#ifdef REGISTER_TIMES
        std::chrono::steady_clock::time_point time_EndDetectBoW = std::chrono::steady_clock::now();
        timeDetectBoW = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndDetectBoW - time_StartDetectBoW).count();
#endif
    }

    // Check the BoW candidates if the geometric candidate list is empty
    //Loop candidates
/*#ifdef COMPILEDWITHC11
    std::chrono::steady_clock::time_point timeStartGeoBoW = std::chrono::steady_clock::now();
#else
    std::chrono::monotonic_clock::time_point timeStartGeoBoW = std::chrono::monotonic_clock::now();
#endif*/

    if(!bLoopDetectedInKF && !vpLoopBowCand.empty())
    {
        mbLoopDetected = DetectCommonRegionsFromBoW(vpLoopBowCand, mpLoopMatchedKF, mpLoopLastCurrentKF, mg2oLoopSlw, mnLoopNumCoincidences, mvpLoopMPs, mvpLoopMatchedMPs);
    }
    // Merge candidates

    //cout << "LC: Find BoW candidates" << endl;

    if(!bMergeDetectedInKF && !vpMergeBowCand.empty())
    {
        mbMergeDetected = DetectCommonRegionsFromBoW(vpMergeBowCand, mpMergeMatchedKF, mpMergeLastCurrentKF, mg2oMergeSlw, mnMergeNumCoincidences, mvpMergeMPs, mvpMergeMatchedMPs);
    }

    //cout << "LC: add to KFDB" << endl;
    mpKeyFrameDB->add(mpCurrentKF);

    if(mbMergeDetected || mbLoopDetected)
    {
        return true;
    }

    //cout << "LC: erase KF" << endl;
    mpCurrentKF->SetErase();
    mpCurrentKF->mbCurrentPlaceRecognition = false;

    return false;
}

bool LoopClosing::DetectAndReffineSim3FromLastKF(KeyFrame* pCurrentKF, KeyFrame* pMatchedKF, g2o::Sim3 &gScw, int &nNumProjMatches,
                                                 std::vector<MapPoint*> &vpMPs, std::vector<MapPoint*> &vpMatchedMPs)
{
    set<MapPoint*> spAlreadyMatchedMPs;
    nNumProjMatches = FindMatchesByProjection(pCurrentKF, pMatchedKF, gScw, spAlreadyMatchedMPs, vpMPs, vpMatchedMPs);


    int nProjMatches = 30;
    int nProjOptMatches = 50;
    int nProjMatchesRep = 100;

    /*if(mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO)
    {
        nProjMatches = 50;
        nProjOptMatches = 50;
        nProjMatchesRep = 80;
    }*/

    if(nNumProjMatches >= nProjMatches)
    {
        cv::Mat mScw = Converter::toCvMat(gScw);
        cv::Mat mTwm = pMatchedKF->GetPoseInverse();
        g2o::Sim3 gSwm(Converter::toMatrix3d(mTwm.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTwm.rowRange(0, 3).col(3)),1.0);
        g2o::Sim3 gScm = gScw * gSwm;
        Eigen::Matrix<double, 7, 7> mHessian7x7;

        bool bFixedScale = mbFixScale;       // TODO CHECK; Solo para el monocular inertial
        if(mpTracker->mSensor==System::IMU_MONOCULAR && !pCurrentKF->GetMap()->GetIniertialBA2())
            bFixedScale=false;
        int numOptMatches = Optimizer::OptimizeSim3(mpCurrentKF, pMatchedKF, vpMatchedMPs, gScm, 10, bFixedScale, mHessian7x7, true);



        if(numOptMatches > nProjOptMatches)
        {
            g2o::Sim3 gScw_estimation(Converter::toMatrix3d(mScw.rowRange(0, 3).colRange(0, 3)),
                           Converter::toVector3d(mScw.rowRange(0, 3).col(3)),1.0);

            vector<MapPoint*> vpMatchedMP;
            vpMatchedMP.resize(mpCurrentKF->GetMapPointMatches().size(), static_cast<MapPoint*>(NULL));

            nNumProjMatches = FindMatchesByProjection(pCurrentKF, pMatchedKF, gScw_estimation, spAlreadyMatchedMPs, vpMPs, vpMatchedMPs);
            //cout << "REFFINE-SIM3: Projection with optimize Sim3 from last KF with " << nNumProjMatches << " matches" << endl;
            if(nNumProjMatches >= nProjMatchesRep)
            {
                gScw = gScw_estimation;
                return true;
            }
        }
    }
    return false;
}

bool LoopClosing::DetectCommonRegionsFromBoW(std::vector<KeyFrame*> &vpBowCand, KeyFrame* &pMatchedKF2, KeyFrame* &pLastCurrentKF, g2o::Sim3 &g2oScw,
                                             int &nNumCoincidences, std::vector<MapPoint*> &vpMPs, std::vector<MapPoint*> &vpMatchedMPs)
{
    int nBoWMatches = 20;
    int nBoWInliers = 15;
    int nSim3Inliers = 20;
    int nProjMatches = 50;
    int nProjOptMatches = 80;
    /*if(mpTracker->mSensor==System::IMU_MONOCULAR ||mpTracker->mSensor==System::IMU_STEREO)
    {
        nBoWMatches = 20;
        nBoWInliers = 15;
        nSim3Inliers = 20;
        nProjMatches = 35;
        nProjOptMatches = 50;
    }*/

    set<KeyFrame*> spConnectedKeyFrames = mpCurrentKF->GetConnectedKeyFrames();

    int nNumCovisibles = 5;

    ORBmatcher matcherBoW(0.9, true);
    ORBmatcher matcher(0.75, true);
    int nNumGuidedMatching = 0;

    // Varibles to select the best numbe
    KeyFrame* pBestMatchedKF;
    int nBestMatchesReproj = 0;
    int nBestNumCoindicendes = 0;
    g2o::Sim3 g2oBestScw;
    std::vector<MapPoint*> vpBestMapPoints;
    std::vector<MapPoint*> vpBestMatchedMapPoints;

    int numCandidates = vpBowCand.size();
    vector<int> vnStage(numCandidates, 0);
    vector<int> vnMatchesStage(numCandidates, 0);

    int index = 0;
    for(KeyFrame* pKFi : vpBowCand)
    {
        //cout << endl << "-------------------------------" << endl;
        if(!pKFi || pKFi->isBad())
            continue;


        // Current KF against KF with covisibles version
        std::vector<KeyFrame*> vpCovKFi = pKFi->GetBestCovisibilityKeyFrames(nNumCovisibles);
        vpCovKFi.push_back(vpCovKFi[0]);
        vpCovKFi[0] = pKFi;

        std::vector<std::vector<MapPoint*> > vvpMatchedMPs;
        vvpMatchedMPs.resize(vpCovKFi.size());
        std::set<MapPoint*> spMatchedMPi;
        int numBoWMatches = 0;

        KeyFrame* pMostBoWMatchesKF = pKFi;
        int nMostBoWNumMatches = 0;

        std::vector<MapPoint*> vpMatchedPoints = std::vector<MapPoint*>(mpCurrentKF->GetMapPointMatches().size(), static_cast<MapPoint*>(NULL));
        std::vector<KeyFrame*> vpKeyFrameMatchedMP = std::vector<KeyFrame*>(mpCurrentKF->GetMapPointMatches().size(), static_cast<KeyFrame*>(NULL));

        int nIndexMostBoWMatchesKF=0;
        for(int j=0; j<vpCovKFi.size(); ++j)
        {
            if(!vpCovKFi[j] || vpCovKFi[j]->isBad())
                continue;

            int num = matcherBoW.SearchByBoW(mpCurrentKF, vpCovKFi[j], vvpMatchedMPs[j]);
            //cout << "BoW: " << num << " putative matches with KF " << vpCovKFi[j]->mnId << endl;
            if (num > nMostBoWNumMatches)
            {
                nMostBoWNumMatches = num;
                nIndexMostBoWMatchesKF = j;
            }
        }

        bool bAbortByNearKF = false;
        for(int j=0; j<vpCovKFi.size(); ++j)
        {
            if(spConnectedKeyFrames.find(vpCovKFi[j]) != spConnectedKeyFrames.end())
            {
                bAbortByNearKF = true;
                //cout << "BoW: Candidate KF aborted by proximity" << endl;
                break;
            }

            //cout << "Matches: " << num << endl;
            for(int k=0; k < vvpMatchedMPs[j].size(); ++k)
            {
                MapPoint* pMPi_j = vvpMatchedMPs[j][k];
                if(!pMPi_j || pMPi_j->isBad())
                    continue;

                if(spMatchedMPi.find(pMPi_j) == spMatchedMPi.end())
                {
                    spMatchedMPi.insert(pMPi_j);
                    numBoWMatches++;

                    vpMatchedPoints[k]= pMPi_j;
                    vpKeyFrameMatchedMP[k] = vpCovKFi[j];
                }
            }
        }

        //cout <<"BoW: " << numBoWMatches << " independent putative matches" << endl;
        if(!bAbortByNearKF && numBoWMatches >= nBoWMatches) // TODO pick a good threshold
        {
            /*cout << "-------------------------------" << endl;
            cout << "Geometric validation with " << numBoWMatches << endl;
            cout << "KFc: " << mpCurrentKF->mnId << "; KFm: " << pMostBoWMatchesKF->mnId << endl;*/
            // Geometric validation

            bool bFixedScale = mbFixScale;
            if(mpTracker->mSensor==System::IMU_MONOCULAR && !mpCurrentKF->GetMap()->GetIniertialBA2())
                bFixedScale=false;

            Sim3Solver solver = Sim3Solver(mpCurrentKF, pMostBoWMatchesKF, vpMatchedPoints, bFixedScale, vpKeyFrameMatchedMP);
            solver.SetRansacParameters(0.99, nBoWInliers, 300); // at least 15 inliers

            bool bNoMore = false;
            vector<bool> vbInliers;
            int nInliers;
            bool bConverge = false;
            cv::Mat mTcm;
            while(!bConverge && !bNoMore)
            {
                mTcm = solver.iterate(20,bNoMore, vbInliers, nInliers, bConverge);
            }

            //cout << "Num inliers: " << nInliers << endl;
            if(bConverge)
            {
                //cout <<"BoW: " << nInliers << " inliers in Sim3Solver" << endl;

                // Match by reprojection
                //int nNumCovisibles = 5;
                vpCovKFi.clear();
                vpCovKFi = pMostBoWMatchesKF->GetBestCovisibilityKeyFrames(nNumCovisibles);
                int nInitialCov = vpCovKFi.size();
                vpCovKFi.push_back(pMostBoWMatchesKF);
                set<KeyFrame*> spCheckKFs(vpCovKFi.begin(), vpCovKFi.end());

                set<MapPoint*> spMapPoints;
                vector<MapPoint*> vpMapPoints;
                vector<KeyFrame*> vpKeyFrames;
                for(KeyFrame* pCovKFi : vpCovKFi)
                {
                    for(MapPoint* pCovMPij : pCovKFi->GetMapPointMatches())
                    {
                        if(!pCovMPij || pCovMPij->isBad())
                            continue;

                        if(spMapPoints.find(pCovMPij) == spMapPoints.end())
                        {
                            spMapPoints.insert(pCovMPij);
                            vpMapPoints.push_back(pCovMPij);
                            vpKeyFrames.push_back(pCovKFi);
                        }
                    }
                }
                //cout << "Point cloud: " << vpMapPoints.size() << endl;


                g2o::Sim3 gScm(Converter::toMatrix3d(solver.GetEstimatedRotation()),Converter::toVector3d(solver.GetEstimatedTranslation()),solver.GetEstimatedScale());
                g2o::Sim3 gSmw(Converter::toMatrix3d(pMostBoWMatchesKF->GetRotation()),Converter::toVector3d(pMostBoWMatchesKF->GetTranslation()),1.0);
                g2o::Sim3 gScw = gScm*gSmw; // Similarity matrix of current from the world position
                cv::Mat mScw = Converter::toCvMat(gScw);


                vector<MapPoint*> vpMatchedMP;
                vpMatchedMP.resize(mpCurrentKF->GetMapPointMatches().size(), static_cast<MapPoint*>(NULL));
                vector<KeyFrame*> vpMatchedKF;
                vpMatchedKF.resize(mpCurrentKF->GetMapPointMatches().size(), static_cast<KeyFrame*>(NULL));
                int numProjMatches = matcher.SearchByProjection(mpCurrentKF, mScw, vpMapPoints, vpKeyFrames, vpMatchedMP, vpMatchedKF, 8, 1.5);

                if(numProjMatches >= nProjMatches)
                {
                    // Optimize Sim3 transformation with every matches
                    Eigen::Matrix<double, 7, 7> mHessian7x7;

                    bool bFixedScale = mbFixScale;
                    if(mpTracker->mSensor==System::IMU_MONOCULAR && !mpCurrentKF->GetMap()->GetIniertialBA2())
                        bFixedScale=false;

                    int numOptMatches = Optimizer::OptimizeSim3(mpCurrentKF, pKFi, vpMatchedMP, gScm, 10, mbFixScale, mHessian7x7, true);
                    //cout <<"BoW: " << numOptMatches << " inliers in the Sim3 optimization" << endl;
                    //cout << "Inliers in Sim3 optimization: " << numOptMatches << endl;

                    if(numOptMatches >= nSim3Inliers)
                    {
                        //cout <<"BoW: " << numOptMatches << " inliers in Sim3 optimization" << endl;
                        g2o::Sim3 gSmw(Converter::toMatrix3d(pMostBoWMatchesKF->GetRotation()),Converter::toVector3d(pMostBoWMatchesKF->GetTranslation()),1.0);
                        g2o::Sim3 gScw = gScm*gSmw; // Similarity matrix of current from the world position
                        cv::Mat mScw = Converter::toCvMat(gScw);

                        vector<MapPoint*> vpMatchedMP;
                        vpMatchedMP.resize(mpCurrentKF->GetMapPointMatches().size(), static_cast<MapPoint*>(NULL));
                        int numProjOptMatches = matcher.SearchByProjection(mpCurrentKF, mScw, vpMapPoints, vpMatchedMP, 5, 1.0);
                        //cout <<"BoW: " << numProjOptMatches << " matches after of the Sim3 optimization" << endl;

                        if(numProjOptMatches >= nProjOptMatches)
                        {
                            int nNumKFs = 0;
                            //vpMatchedMPs = vpMatchedMP;
                            //vpMPs = vpMapPoints;
                            // Check the Sim3 transformation with the current KeyFrame covisibles
                            vector<KeyFrame*> vpCurrentCovKFs = mpCurrentKF->GetBestCovisibilityKeyFrames(nNumCovisibles);
                            //cout << "---" << endl;
                            //cout << "BoW: Geometrical validation" << endl;
                            int j = 0;
                            while(nNumKFs < 3 && j<vpCurrentCovKFs.size())
                            //for(int j=0; j<vpCurrentCovKFs.size(); ++j)
                            {
                                KeyFrame* pKFj = vpCurrentCovKFs[j];
                                cv::Mat mTjc = pKFj->GetPose() * mpCurrentKF->GetPoseInverse();
                                g2o::Sim3 gSjc(Converter::toMatrix3d(mTjc.rowRange(0, 3).colRange(0, 3)),Converter::toVector3d(mTjc.rowRange(0, 3).col(3)),1.0);
                                g2o::Sim3 gSjw = gSjc * gScw;
                                int numProjMatches_j = 0;
                                vector<MapPoint*> vpMatchedMPs_j;
                                bool bValid = DetectCommonRegionsFromLastKF(pKFj,pMostBoWMatchesKF, gSjw,numProjMatches_j, vpMapPoints, vpMatchedMPs_j);

                                if(bValid)
                                {
                                    nNumKFs++;
                                }

                                j++;
                            }

                            if(nNumKFs < 3)
                            {
                                vnStage[index] = 8;
                                vnMatchesStage[index] = nNumKFs;
                            }

                            if(nBestMatchesReproj < numProjOptMatches)
                            {
                                nBestMatchesReproj = numProjOptMatches;
                                nBestNumCoindicendes = nNumKFs;
                                pBestMatchedKF = pMostBoWMatchesKF;
                                g2oBestScw = gScw;
                                vpBestMapPoints = vpMapPoints;
                                vpBestMatchedMapPoints = vpMatchedMP;
                            }


                        }

                    }

                }
            }
        }
        index++;
    }

    if(nBestMatchesReproj > 0)
    {
        pLastCurrentKF = mpCurrentKF;
        nNumCoincidences = nBestNumCoindicendes;
        pMatchedKF2 = pBestMatchedKF;
        pMatchedKF2->SetNotErase();
        g2oScw = g2oBestScw;
        vpMPs = vpBestMapPoints;
        vpMatchedMPs = vpBestMatchedMapPoints;

        return nNumCoincidences >= 3;
    }
    else
    {
        int maxStage = -1;
        int maxMatched;
        for(int i=0; i<vnStage.size(); ++i)
        {
            if(vnStage[i] > maxStage)
            {
                maxStage = vnStage[i];
                maxMatched = vnMatchesStage[i];
            }
        }

//        f_succes_pr << mpCurrentKF->mNameFile << " " << std::to_string(maxStage) << endl;
//        f_succes_pr << "% NumCand: " << std::to_string(numCandidates) << "; matches: " << std::to_string(maxMatched) << endl;
    }
    return false;
}

bool LoopClosing::DetectCommonRegionsFromLastKF(KeyFrame* pCurrentKF, KeyFrame* pMatchedKF, g2o::Sim3 &gScw, int &nNumProjMatches,
                                                std::vector<MapPoint*> &vpMPs, std::vector<MapPoint*> &vpMatchedMPs)
{
    set<MapPoint*> spAlreadyMatchedMPs(vpMatchedMPs.begin(), vpMatchedMPs.end());
    nNumProjMatches = FindMatchesByProjection(pCurrentKF, pMatchedKF, gScw, spAlreadyMatchedMPs, vpMPs, vpMatchedMPs);
    //cout << "Projection from last KF with " << nNumProjMatches << " matches" << endl;

    int nProjMatches = 30;
    if(nNumProjMatches >= nProjMatches)
    {
        /*cout << "PR_From_LastKF: KF " << pCurrentKF->mnId << " has " << nNumProjMatches << " with KF " << pMatchedKF->mnId << endl;

        int max_x = -1, min_x = 1000000;
        int max_y = -1, min_y = 1000000;
        for(MapPoint* pMPi : vpMatchedMPs)
        {
            if(!pMPi || pMPi->isBad())
            {
                continue;
            }

            tuple<size_t,size_t> indexes = pMPi->GetIndexInKeyFrame(pMatchedKF);
            int index = get<0>(indexes);
            if(index >= 0)
            {
                int coord_x = pCurrentKF->mvKeysUn[index].pt.x;
                if(coord_x < min_x)
                {
                    min_x = coord_x;
                }
                if(coord_x > max_x)
                {
                    max_x = coord_x;
                }
                int coord_y = pCurrentKF->mvKeysUn[index].pt.y;
                if(coord_y < min_y)
                {
                    min_y = coord_y;
                }
                if(coord_y > max_y)
                {
                    max_y = coord_y;
                }
            }
        }*/
        //cout << "PR_From_LastKF: Coord in X -> " << min_x << ", " << max_x << "; and Y -> " << min_y << ", " << max_y << endl;
        //cout << "PR_From_LastKF: features area in X -> " << (max_x - min_x) << " and Y -> " << (max_y - min_y) << endl;


        return true;
    }

    return false;
}

int LoopClosing::FindMatchesByProjection(KeyFrame* pCurrentKF, KeyFrame* pMatchedKFw, g2o::Sim3 &g2oScw,
                                         set<MapPoint*> &spMatchedMPinOrigin, vector<MapPoint*> &vpMapPoints,
                                         vector<MapPoint*> &vpMatchedMapPoints)
{
    int nNumCovisibles = 5;
    vector<KeyFrame*> vpCovKFm = pMatchedKFw->GetBestCovisibilityKeyFrames(nNumCovisibles);
    int nInitialCov = vpCovKFm.size();
    vpCovKFm.push_back(pMatchedKFw);
    set<KeyFrame*> spCheckKFs(vpCovKFm.begin(), vpCovKFm.end());
    set<KeyFrame*> spCurrentCovisbles = pCurrentKF->GetConnectedKeyFrames();
    for(int i=0; i<nInitialCov; ++i)
    {
        vector<KeyFrame*> vpKFs = vpCovKFm[i]->GetBestCovisibilityKeyFrames(nNumCovisibles);
        int nInserted = 0;
        int j = 0;
        while(j < vpKFs.size() && nInserted < nNumCovisibles)
        {
            if(spCheckKFs.find(vpKFs[j]) == spCheckKFs.end() && spCurrentCovisbles.find(vpKFs[j]) == spCurrentCovisbles.end())
            {
                spCheckKFs.insert(vpKFs[j]);
                ++nInserted;
            }
            ++j;
        }
        vpCovKFm.insert(vpCovKFm.end(), vpKFs.begin(), vpKFs.end());
    }
    set<MapPoint*> spMapPoints;
    vpMapPoints.clear();
    vpMatchedMapPoints.clear();
    for(KeyFrame* pKFi : vpCovKFm)
    {
        for(MapPoint* pMPij : pKFi->GetMapPointMatches())
        {
            if(!pMPij || pMPij->isBad())
                continue;

            if(spMapPoints.find(pMPij) == spMapPoints.end())
            {
                spMapPoints.insert(pMPij);
                vpMapPoints.push_back(pMPij);
            }
        }
    }
    //cout << "Point cloud: " << vpMapPoints.size() << endl;

    cv::Mat mScw = Converter::toCvMat(g2oScw);

    ORBmatcher matcher(0.9, true);

    vpMatchedMapPoints.resize(pCurrentKF->GetMapPointMatches().size(), static_cast<MapPoint*>(NULL));
    int num_matches = matcher.SearchByProjection(pCurrentKF, mScw, vpMapPoints, vpMatchedMapPoints, 3, 1.5);

    return num_matches;
}

void LoopClosing::CorrectLoop()
{
    cout << "Loop detected!" << endl;

    // Send a stop signal to Local Mapping
    // Avoid new keyframes are inserted while correcting the loop
    // Step 0：结束局部地图线程、全局BA，为闭环矫正做准备
    // 请求局部地图停止，防止在回环矫正时局部地图线程中InsertKeyFrame函数插入新的关键帧
    mpLocalMapper->RequestStop();
    mpLocalMapper->EmptyQueue(); // Proccess keyframes in the queue


    // If a Global Bundle Adjustment is running, abort it
    cout << "Request GBA abort" << endl;
    if(isRunningGBA())
    {
        // 如果有全局BA在运行，终止掉，迎接新的全局BA
        unique_lock<mutex> lock(mMutexGBA);
        mbStopGBA = true;
        // 记录全局BA次数
        mnFullBAIdx++;

        if(mpThreadGBA)
        {
            cout << "GBA running... Abort!" << endl;
            mpThreadGBA->detach();
            delete mpThreadGBA;
        }
    }

    // Wait until Local Mapping has effectively stopped
    // 一直等到局部地图线程结束再继续
    while(!mpLocalMapper->isStopped())
    {
        usleep(1000);
    }

    // Ensure current keyframe is updated
	// Step 1：根据共视关系更新当前关键帧与其它关键帧之间的连接关系
    // 因为之前闭环检测、计算Sim3中改变了该关键帧的地图点，所以需要更新
    cout << "start updating connections" << endl;
    mpCurrentKF->UpdateConnections();

    // Retrive keyframes connected to the current keyframe and compute corrected Sim3 pose by propagation
    // Step 2：通过位姿传播，得到Sim3优化后，与当前帧相连的关键帧的位姿，以及它们的地图点
    // 当前帧与世界坐标系之间的Sim变换在ComputeSim3函数中已经确定并优化，
    // 通过相对位姿关系，可以确定这些相连的关键帧与世界坐标系之间的Sim3变换

    // 取出当前关键帧及其共视关键帧，称为“当前关键帧组”
    mvpCurrentConnectedKFs = mpCurrentKF->GetVectorCovisibleKeyFrames();
    mvpCurrentConnectedKFs.push_back(mpCurrentKF);

    // CorrectedSim3：存放闭环g2o优化后当前关键帧的共视关键帧的世界坐标系下Sim3 变换
    // NonCorrectedSim3：存放没有矫正的当前关键帧的共视关键帧的世界坐标系下Sim3 变换
    KeyFrameAndPose CorrectedSim3, NonCorrectedSim3;
    // 先将mpCurrentKF的Sim3变换存入，认为是准的，所以固定不动
    CorrectedSim3[mpCurrentKF]=mg2oLoopScw;
    // 当前关键帧到世界坐标系下的变换矩阵
    cv::Mat Twc = mpCurrentKF->GetPoseInverse();

    Map* pLoopMap = mpCurrentKF->GetMap();
    // 对地图点操作
    {
        // Get Map Mutex
        // 锁定地图点
        unique_lock<mutex> lock(pLoopMap->mMutexMapUpdate);

        const bool bImuInit = pLoopMap->isImuInitialized();
        // Step 2.1：通过mg2oLoopScw（认为是准的）来进行位姿传播，得到当前关键帧的共视关键帧的世界坐标系下Sim3 位姿（还没有修正）
        // 遍历"当前关键帧组""
        for(vector<KeyFrame*>::iterator vit=mvpCurrentConnectedKFs.begin(), vend=mvpCurrentConnectedKFs.end(); vit!=vend; vit++)
        {
            KeyFrame* pKFi = *vit;

            cv::Mat Tiw = pKFi->GetPose();

            if(pKFi!=mpCurrentKF)      //跳过当前关键帧，因为当前关键帧的位姿已经在前面优化过了，在这里是参考基准
            {
                // 得到当前关键帧 mpCurrentKF 到其共视关键帧 pKFi 的相对变换
                cv::Mat Tic = Tiw*Twc;
                cv::Mat Ric = Tic.rowRange(0,3).colRange(0,3);
                cv::Mat tic = Tic.rowRange(0,3).col(3);

                // g2oSic：当前关键帧 mpCurrentKF 到其共视关键帧 pKFi 的Sim3 相对变换
                // 这里是non-correct, 所以scale=1.0
                g2o::Sim3 g2oSic(Converter::toMatrix3d(Ric),Converter::toVector3d(tic),1.0);
                // 当前帧的位姿固定不动，其它的关键帧根据相对关系得到Sim3调整的位姿
                g2o::Sim3 g2oCorrectedSiw = g2oSic*mg2oLoopScw;
                //Pose corrected with the Sim3 of the loop closure
                // 存放闭环g2o优化后当前关键帧的共视关键帧的Sim3 位姿
                CorrectedSim3[pKFi]=g2oCorrectedSiw;
            }

            cv::Mat Riw = Tiw.rowRange(0,3).colRange(0,3);
            cv::Mat tiw = Tiw.rowRange(0,3).col(3);
            g2o::Sim3 g2oSiw(Converter::toMatrix3d(Riw),Converter::toVector3d(tiw),1.0);
            //Pose without correction
            NonCorrectedSim3[pKFi]=g2oSiw;
        }

        // Correct all MapPoints obsrved by current keyframe and neighbors, so that they align with the other side of the loop
		// Step 2.2：得到矫正的当前关键帧的共视关键帧位姿后，修正这些关键帧的地图点
        // 遍历待矫正的共视关键帧（不包括当前关键帧）
        for(KeyFrameAndPose::iterator mit=CorrectedSim3.begin(), mend=CorrectedSim3.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;
            g2o::Sim3 g2oCorrectedSiw = mit->second;
            g2o::Sim3 g2oCorrectedSwi = g2oCorrectedSiw.inverse();

            g2o::Sim3 g2oSiw =NonCorrectedSim3[pKFi];

            vector<MapPoint*> vpMPsi = pKFi->GetMapPointMatches();
            // 遍历待矫正共视关键帧中的每一个地图点
            for(size_t iMP=0, endMPi = vpMPsi.size(); iMP<endMPi; iMP++)
            {
                MapPoint* pMPi = vpMPsi[iMP];
                if(!pMPi)
                    continue;
                if(pMPi->isBad())
                    continue;
                if(pMPi->mnCorrectedByKF==mpCurrentKF->mnId) // 标记，防止重复矫正
                    continue;

                // 矫正过程本质上也是基于当前关键帧的优化后的位姿展开的
                // Project with non-corrected pose and project back with corrected pose
                // 将该未校正的eigP3Dw先从世界坐标系映射到未校正的pKFi相机坐标系，然后再反映射到校正后的世界坐标系下
                cv::Mat P3Dw = pMPi->GetWorldPos();
                // 地图点世界坐标系下坐标
                Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
                // map(P) 内部做了变换 R*P +t  
                // 下面变换是：eigP3Dw： world →g2oSiw→ i →g2oCorrectedSwi→ world
                Eigen::Matrix<double,3,1> eigCorrectedP3Dw = g2oCorrectedSwi.map(g2oSiw.map(eigP3Dw));

                cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);
                pMPi->SetWorldPos(cvCorrectedP3Dw);
                // 记录矫正该地图点的关键帧id，防止重复
                pMPi->mnCorrectedByKF = mpCurrentKF->mnId;
                // 记录该地图点所在的关键帧id
                pMPi->mnCorrectedReference = pKFi->mnId;
                // 因为地图点更新了，需要更新其平均观测方向以及观测距离范围
                pMPi->UpdateNormalAndDepth();
            }

            // Update keyframe pose with corrected Sim3. First transform Sim3 to SE3 (scale translation)
            // Step 2.3：将共视关键帧的Sim3转换为SE3，根据更新的Sim3，更新关键帧的位姿
            // 其实是现在已经有了更新后的关键帧组中关键帧的位姿,但是在上面的操作时只是暂时存储到了 KeyFrameAndPose 类型的变量中,还没有写回到关键帧对象中
            // 调用toRotationMatrix 可以自动归一化旋转矩阵
            Eigen::Matrix3d eigR = g2oCorrectedSiw.rotation().toRotationMatrix();
            Eigen::Vector3d eigt = g2oCorrectedSiw.translation();
            double s = g2oCorrectedSiw.scale();
            // cout << "scale for loop-closing: " << s << endl;
            // 平移向量中包含有尺度信息，还需要用尺度归一化
            eigt *=(1./s); //[R t/s;0 1]

            cv::Mat correctedTiw = Converter::toCvSE3(eigR,eigt);
            // 设置矫正后的新的pose
            pKFi->SetPose(correctedTiw);

            // Correct velocity according to orientation correction
            if(bImuInit)
            {
                Eigen::Matrix3d Rcor = eigR.transpose()*g2oSiw.rotation().toRotationMatrix();
                pKFi->SetVelocity(Converter::toCvMat(Rcor)*pKFi->GetVelocity());
            }

            // Make sure connections are updated
            // Step 2.4：根据共视关系更新当前帧与其它关键帧之间的连接
            // 地图点的位置改变了,可能会引起共视关系\权值的改变 
            pKFi->UpdateConnections();
        }
        // TODO Check this index increasement
        pLoopMap->IncreaseChangeIndex();


        // Start Loop Fusion
        // Update matched map points and replace if duplicated
        // Step 3：检查当前帧的地图点与经过闭环匹配后该帧的地图点是否存在冲突，对冲突的进行替换或填补
        // mvpCurrentMatchedPoints 是当前关键帧和闭环关键帧组的所有地图点进行投影得到的匹配点
        for(size_t i=0; i<mvpLoopMatchedMPs.size(); i++)
        {
            if(mvpLoopMatchedMPs[i])
            {
                //取出同一个索引对应的两种地图点，决定是否要替换
                // 匹配投影得到的地图点
                MapPoint* pLoopMP = mvpLoopMatchedMPs[i];
                // 原来的地图点
                MapPoint* pCurMP = mpCurrentKF->GetMapPoint(i);
                if(pCurMP)
                    // 如果有重复的MapPoint，则用匹配的地图点代替现有的
                    // 因为匹配的地图点是经过一系列操作后比较精确的，现有的地图点很可能有累计误差
                    pCurMP->Replace(pLoopMP);
                else
                {
                    // 如果当前帧没有该MapPoint，则直接添加
                    mpCurrentKF->AddMapPoint(pLoopMP,i);
                    pLoopMP->AddObservation(mpCurrentKF,i);
                    pLoopMP->ComputeDistinctiveDescriptors();
                }
            }
        }
    }

    // Project MapPoints observed in the neighborhood of the loop keyframe
    // into the current keyframe and neighbors using corrected poses.
    // Fuse duplications.
    // Step 4：将闭环相连关键帧组mvpLoopMapPoints 投影到当前关键帧组中，进行匹配，融合，新增或替换当前关键帧组中KF的地图点
    // 因为 闭环相连关键帧组mvpLoopMapPoints 在地图中时间比较久经历了多次优化，认为是准确的
    // 而当前关键帧组中的关键帧的地图点是最近新计算的，可能有累积误差
    // CorrectedSim3：存放矫正后当前关键帧的共视关键帧，及其世界坐标系下Sim3 变换
    SearchAndFuse(CorrectedSim3, mvpLoopMapPoints);
    //cout << "LC: end SearchAndFuse" << endl;


    // After the MapPoint fusion, new links in the covisibility graph will appear attaching both sides of the loop
    // Step 5：更新当前关键帧之间的共视相连关系，得到因闭环时MapPoints融合而新得到的连接关系
    // LoopConnections：存储因为闭环地图点调整而新生成的连接关系
    map<KeyFrame*, set<KeyFrame*> > LoopConnections;

    // Step 5.1：遍历当前帧相连关键帧组（一级相连）
    for(vector<KeyFrame*>::iterator vit=mvpCurrentConnectedKFs.begin(), vend=mvpCurrentConnectedKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        // Step 5.2：得到与当前帧相连关键帧的相连关键帧（二级相连）
        vector<KeyFrame*> vpPreviousNeighbors = pKFi->GetVectorCovisibleKeyFrames();

        // Update connections. Detect new links.
        // Step 5.3：更新一级相连关键帧的连接关系(会把当前关键帧添加进去,因为地图点已经更新和替换了)
        pKFi->UpdateConnections();
        // Step 5.4：取出该帧更新后的连接关系
        LoopConnections[pKFi]=pKFi->GetConnectedKeyFrames();
        // Step 5.5：从连接关系中去除闭环之前的二级连接关系，剩下的连接就是由闭环得到的连接关系
        for(vector<KeyFrame*>::iterator vit_prev=vpPreviousNeighbors.begin(), vend_prev=vpPreviousNeighbors.end(); vit_prev!=vend_prev; vit_prev++)
        {
            LoopConnections[pKFi].erase(*vit_prev);
        }
        // Step 5.6：从连接关系中去除闭环之前的一级连接关系，剩下的连接就是由闭环得到的连接关系
        for(vector<KeyFrame*>::iterator vit2=mvpCurrentConnectedKFs.begin(), vend2=mvpCurrentConnectedKFs.end(); vit2!=vend2; vit2++)
        {
            LoopConnections[pKFi].erase(*vit2);
        }
    }
    //cout << "LC: end updating covisibility graph" << endl;

    // Optimize graph
    //cout << "start opt essentialgraph" << endl;
    bool bFixedScale = mbFixScale;
    // TODO CHECK; Solo para el monocular inertial
    if(mpTracker->mSensor==System::IMU_MONOCULAR && !mpCurrentKF->GetMap()->GetIniertialBA2())
        bFixedScale=false;


    //cout << "Optimize essential graph" << endl;
    if(pLoopMap->IsInertial() && pLoopMap->isImuInitialized())
    {
        //cout << "With 4DoF" << endl;
        Optimizer::OptimizeEssentialGraph4DoF(pLoopMap, mpLoopMatchedKF, mpCurrentKF, NonCorrectedSim3, CorrectedSim3, LoopConnections);
    }
    else
    {
        //cout << "With 7DoF" << endl;
		// Step 6：进行EssentialGraph优化，LoopConnections是形成闭环后新生成的连接关系，不包括步骤7中当前帧与闭环匹配帧之间的连接关系
        Optimizer::OptimizeEssentialGraph(pLoopMap, mpLoopMatchedKF, mpCurrentKF, NonCorrectedSim3, CorrectedSim3, LoopConnections, bFixedScale);
    }


    //cout << "Optimize essential graph finished" << endl;
    //usleep(5*1000*1000);

    mpAtlas->InformNewBigChange();

    // Add loop edge
    // Step 7：添加当前帧与闭环匹配帧之间的边（这个连接关系不优化）
    // !这两句话应该放在OptimizeEssentialGraph之前，因为OptimizeEssentialGraph的步骤4.2中有优化
    mpLoopMatchedKF->AddLoopEdge(mpCurrentKF);
    mpCurrentKF->AddLoopEdge(mpLoopMatchedKF);

    // Launch a new thread to perform Global Bundle Adjustment (Only if few keyframes, if not it would take too much time)
    if(!pLoopMap->isImuInitialized() || (pLoopMap->KeyFramesInMap()<200 && mpAtlas->CountMaps()==1))
    {
	    // Step 8：新建一个线程用于全局BA优化
    	// OptimizeEssentialGraph只是优化了一些主要关键帧的位姿，这里进行全局BA可以全局优化所有位姿和MapPoints

        mbRunningGBA = true;
        mbFinishedGBA = false;
        mbStopGBA = false;

        mpThreadGBA = new thread(&LoopClosing::RunGlobalBundleAdjustment, this, pLoopMap, mpCurrentKF->mnId);
    }

    // Loop closed. Release Local Mapping.
    mpLocalMapper->Release();    

    mLastLoopKFid = mpCurrentKF->mnId; //TODO old varible, it is not use in the new algorithm
}

void LoopClosing::MergeLocal()
{
    Verbose::PrintMess("MERGE: Merge Visual detected!!!!", Verbose::VERBOSITY_NORMAL);

    int numTemporalKFs = 15; //TODO (set by parameter): Temporal KFs in the local window if the map is inertial.

    //Relationship to rebuild the essential graph, it is used two times, first in the local window and later in the rest of the map
    KeyFrame* pNewChild;
    KeyFrame* pNewParent;

    vector<KeyFrame*> vpLocalCurrentWindowKFs;
    vector<KeyFrame*> vpMergeConnectedKFs;

    // Flag that is true only when we stopped a running BA, in this case we need relaunch at the end of the merge
    bool bRelaunchBA = false;

    Verbose::PrintMess("MERGE: Check Full Bundle Adjustment", Verbose::VERBOSITY_DEBUG);
    // If a Global Bundle Adjustment is running, abort it
    if(isRunningGBA())
    {
        unique_lock<mutex> lock(mMutexGBA);
        mbStopGBA = true;

        mnFullBAIdx++;

        if(mpThreadGBA)
        {
            mpThreadGBA->detach();
            delete mpThreadGBA;
        }
        bRelaunchBA = true;
    }

    Verbose::PrintMess("MERGE: Request Stop Local Mapping", Verbose::VERBOSITY_DEBUG);
    mpLocalMapper->RequestStop();
    // Wait until Local Mapping has effectively stopped
    while(!mpLocalMapper->isStopped())
    {
        usleep(1000);
    }
    Verbose::PrintMess("MERGE: Local Map stopped", Verbose::VERBOSITY_DEBUG);

    mpLocalMapper->EmptyQueue();

    // Merge map will become in the new active map with the local window of KFs and MPs from the current map.
    // Later, the elements of the current map will be transform to the new active map reference, in order to keep real time tracking
    Map* pCurrentMap = mpCurrentKF->GetMap();
    Map* pMergeMap = mpMergeMatchedKF->GetMap();

    // Ensure current keyframe is updated
    mpCurrentKF->UpdateConnections();

    //Get the current KF and its neighbors(visual->covisibles; inertial->temporal+covisibles)
    set<KeyFrame*> spLocalWindowKFs;
    //Get MPs in the welding area from the current map
    set<MapPoint*> spLocalWindowMPs;
    if(pCurrentMap->IsInertial() && pMergeMap->IsInertial()) //TODO Check the correct initialization
    {
        KeyFrame* pKFi = mpCurrentKF;
        int nInserted = 0;
        while(pKFi && nInserted < numTemporalKFs)
        {
            spLocalWindowKFs.insert(pKFi);
            pKFi = mpCurrentKF->mPrevKF;
            nInserted++;

            set<MapPoint*> spMPi = pKFi->GetMapPoints();
            spLocalWindowMPs.insert(spMPi.begin(), spMPi.end());
        }

        pKFi = mpCurrentKF->mNextKF;
        while(pKFi)
        {
            spLocalWindowKFs.insert(pKFi);

            set<MapPoint*> spMPi = pKFi->GetMapPoints();
            spLocalWindowMPs.insert(spMPi.begin(), spMPi.end());
        }
    }
    else
    {
        spLocalWindowKFs.insert(mpCurrentKF);
    }

    vector<KeyFrame*> vpCovisibleKFs = mpCurrentKF->GetBestCovisibilityKeyFrames(numTemporalKFs);
    spLocalWindowKFs.insert(vpCovisibleKFs.begin(), vpCovisibleKFs.end());
    const int nMaxTries = 3;
    int nNumTries = 0;
    while(spLocalWindowKFs.size() < numTemporalKFs && nNumTries < nMaxTries)
    {
        vector<KeyFrame*> vpNewCovKFs;
        vpNewCovKFs.empty();
        for(KeyFrame* pKFi : spLocalWindowKFs)
        {
            vector<KeyFrame*> vpKFiCov = pKFi->GetBestCovisibilityKeyFrames(numTemporalKFs/2);
            for(KeyFrame* pKFcov : vpKFiCov)
            {
                if(pKFcov && !pKFcov->isBad() && spLocalWindowKFs.find(pKFcov) == spLocalWindowKFs.end())
                {
                    vpNewCovKFs.push_back(pKFcov);
                }

            }
        }

        spLocalWindowKFs.insert(vpNewCovKFs.begin(), vpNewCovKFs.end());
        nNumTries++;
    }

    for(KeyFrame* pKFi : spLocalWindowKFs)
    {
        if(!pKFi || pKFi->isBad())
            continue;

        set<MapPoint*> spMPs = pKFi->GetMapPoints();
        spLocalWindowMPs.insert(spMPs.begin(), spMPs.end());
    }

    set<KeyFrame*> spMergeConnectedKFs;
    if(pCurrentMap->IsInertial() && pMergeMap->IsInertial()) //TODO Check the correct initialization
    {
        KeyFrame* pKFi = mpMergeMatchedKF;
        int nInserted = 0;
        while(pKFi && nInserted < numTemporalKFs)
        {
            spMergeConnectedKFs.insert(pKFi);
            pKFi = mpCurrentKF->mPrevKF;
            nInserted++;
        }

        pKFi = mpMergeMatchedKF->mNextKF;
        while(pKFi)
        {
            spMergeConnectedKFs.insert(pKFi);
        }
    }
    else
    {
        spMergeConnectedKFs.insert(mpMergeMatchedKF);
    }
    vpCovisibleKFs = mpMergeMatchedKF->GetBestCovisibilityKeyFrames(numTemporalKFs);
    spMergeConnectedKFs.insert(vpCovisibleKFs.begin(), vpCovisibleKFs.end());
    nNumTries = 0;
    while(spMergeConnectedKFs.size() < numTemporalKFs && nNumTries < nMaxTries)
    {
        vector<KeyFrame*> vpNewCovKFs;
        for(KeyFrame* pKFi : spMergeConnectedKFs)
        {
            vector<KeyFrame*> vpKFiCov = pKFi->GetBestCovisibilityKeyFrames(numTemporalKFs/2);
            for(KeyFrame* pKFcov : vpKFiCov)
            {
                if(pKFcov && !pKFcov->isBad() && spMergeConnectedKFs.find(pKFcov) == spMergeConnectedKFs.end())
                {
                    vpNewCovKFs.push_back(pKFcov);
                }

            }
        }

        spMergeConnectedKFs.insert(vpNewCovKFs.begin(), vpNewCovKFs.end());
        nNumTries++;
    }

    set<MapPoint*> spMapPointMerge;
    for(KeyFrame* pKFi : spMergeConnectedKFs)
    {
        set<MapPoint*> vpMPs = pKFi->GetMapPoints();
        spMapPointMerge.insert(vpMPs.begin(),vpMPs.end());
    }

    vector<MapPoint*> vpCheckFuseMapPoint;
    vpCheckFuseMapPoint.reserve(spMapPointMerge.size());
    std::copy(spMapPointMerge.begin(), spMapPointMerge.end(), std::back_inserter(vpCheckFuseMapPoint));

    //
    cv::Mat Twc = mpCurrentKF->GetPoseInverse();

    cv::Mat Rwc = Twc.rowRange(0,3).colRange(0,3);
    cv::Mat twc = Twc.rowRange(0,3).col(3);
    g2o::Sim3 g2oNonCorrectedSwc(Converter::toMatrix3d(Rwc),Converter::toVector3d(twc),1.0);
    g2o::Sim3 g2oNonCorrectedScw = g2oNonCorrectedSwc.inverse();
    g2o::Sim3 g2oCorrectedScw = mg2oMergeScw; //TODO Check the transformation

    KeyFrameAndPose vCorrectedSim3, vNonCorrectedSim3;
    vCorrectedSim3[mpCurrentKF]=g2oCorrectedScw;
    vNonCorrectedSim3[mpCurrentKF]=g2oNonCorrectedScw;

    for(KeyFrame* pKFi : spLocalWindowKFs)
    {
        if(!pKFi || pKFi->isBad())
        {
            continue;
        }

        g2o::Sim3 g2oCorrectedSiw;

        if(pKFi!=mpCurrentKF)
        {
            cv::Mat Tiw = pKFi->GetPose();
            cv::Mat Riw = Tiw.rowRange(0,3).colRange(0,3);
            cv::Mat tiw = Tiw.rowRange(0,3).col(3);
            g2o::Sim3 g2oSiw(Converter::toMatrix3d(Riw),Converter::toVector3d(tiw),1.0);
            //Pose without correction
            vNonCorrectedSim3[pKFi]=g2oSiw;

            cv::Mat Tic = Tiw*Twc;
            cv::Mat Ric = Tic.rowRange(0,3).colRange(0,3);
            cv::Mat tic = Tic.rowRange(0,3).col(3);
            g2o::Sim3 g2oSic(Converter::toMatrix3d(Ric),Converter::toVector3d(tic),1.0);
            g2oCorrectedSiw = g2oSic*mg2oMergeScw;
            vCorrectedSim3[pKFi]=g2oCorrectedSiw;
        }
        else
        {
            g2oCorrectedSiw = g2oCorrectedScw;
        }
        pKFi->mTcwMerge  = pKFi->GetPose();

        // Update keyframe pose with corrected Sim3. First transform Sim3 to SE3 (scale translation)
        Eigen::Matrix3d eigR = g2oCorrectedSiw.rotation().toRotationMatrix();
        Eigen::Vector3d eigt = g2oCorrectedSiw.translation();
        double s = g2oCorrectedSiw.scale();

        pKFi->mfScale = s;
        eigt *=(1./s); //[R t/s;0 1]

        cv::Mat correctedTiw = Converter::toCvSE3(eigR,eigt);

        pKFi->mTcwMerge = correctedTiw;

        if(pCurrentMap->isImuInitialized())
        {
            Eigen::Matrix3d Rcor = eigR.transpose()*vNonCorrectedSim3[pKFi].rotation().toRotationMatrix();
            pKFi->mVwbMerge = Converter::toCvMat(Rcor)*pKFi->GetVelocity();
        }

    }

    for(MapPoint* pMPi : spLocalWindowMPs)
    {
        if(!pMPi || pMPi->isBad())
            continue;

        KeyFrame* pKFref = pMPi->GetReferenceKeyFrame();
        g2o::Sim3 g2oCorrectedSwi = vCorrectedSim3[pKFref].inverse();
        g2o::Sim3 g2oNonCorrectedSiw = vNonCorrectedSim3[pKFref];

        // Project with non-corrected pose and project back with corrected pose
        cv::Mat P3Dw = pMPi->GetWorldPos();
        Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
        Eigen::Matrix<double,3,1> eigCorrectedP3Dw = g2oCorrectedSwi.map(g2oNonCorrectedSiw.map(eigP3Dw));
        Eigen::Matrix3d eigR = g2oCorrectedSwi.rotation().toRotationMatrix();
        Eigen::Matrix3d Rcor = eigR * g2oNonCorrectedSiw.rotation().toRotationMatrix();

        cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);

        pMPi->mPosMerge = cvCorrectedP3Dw;
        pMPi->mNormalVectorMerge = Converter::toCvMat(Rcor) * pMPi->GetNormal();
    }

    {
        unique_lock<mutex> currentLock(pCurrentMap->mMutexMapUpdate); // We update the current map with the Merge information
        unique_lock<mutex> mergeLock(pMergeMap->mMutexMapUpdate); // We remove the Kfs and MPs in the merged area from the old map

        for(KeyFrame* pKFi : spLocalWindowKFs)
        {
            if(!pKFi || pKFi->isBad())
            {
                //cout << "Bad KF in correction" << endl;
                continue;
            }

            pKFi->mTcwBefMerge = pKFi->GetPose();
            pKFi->mTwcBefMerge = pKFi->GetPoseInverse();
            pKFi->SetPose(pKFi->mTcwMerge);

            // Make sure connections are updated
            pKFi->UpdateMap(pMergeMap);
            pKFi->mnMergeCorrectedForKF = mpCurrentKF->mnId;
            pMergeMap->AddKeyFrame(pKFi);
            pCurrentMap->EraseKeyFrame(pKFi);

            if(pCurrentMap->isImuInitialized())
            {
                pKFi->SetVelocity(pKFi->mVwbMerge);
            }
        }

        for(MapPoint* pMPi : spLocalWindowMPs)
        {
            if(!pMPi || pMPi->isBad())
                continue;

            pMPi->SetWorldPos(pMPi->mPosMerge);
            pMPi->SetNormalVector(pMPi->mNormalVectorMerge);
            pMPi->UpdateMap(pMergeMap);
            pMergeMap->AddMapPoint(pMPi);
            pCurrentMap->EraseMapPoint(pMPi);
        }

        mpAtlas->ChangeMap(pMergeMap);
        mpAtlas->SetMapBad(pCurrentMap);
        pMergeMap->IncreaseChangeIndex();
    }


    //Rebuild the essential graph in the local window
    pCurrentMap->GetOriginKF()->SetFirstConnection(false);
    pNewChild = mpCurrentKF->GetParent(); // Old parent, it will be the new child of this KF
    pNewParent = mpCurrentKF; // Old child, now it will be the parent of its own parent(we need eliminate this KF from children list in its old parent)
    mpCurrentKF->ChangeParent(mpMergeMatchedKF);
    while(pNewChild /*&& spLocalWindowKFs.find(pNewChild) != spLocalWindowKFs.end()*/)
    {
        pNewChild->EraseChild(pNewParent); // We remove the relation between the old parent and the new for avoid loop
        KeyFrame * pOldParent = pNewChild->GetParent();

        pNewChild->ChangeParent(pNewParent);
        //cout << "The new parent of KF " << pNewChild->mnId << " was " << pNewChild->GetParent()->mnId << endl;

        pNewParent = pNewChild;
        pNewChild = pOldParent;

    }

    //Update the connections between the local window
    mpMergeMatchedKF->UpdateConnections();
    //cout << "MERGE-VISUAL: Essential graph rebuilded" << endl;


    //std::copy(spMapPointCurrent.begin(), spMapPointCurrent.end(), std::back_inserter(vpCheckFuseMapPoint));
    vpMergeConnectedKFs = mpMergeMatchedKF->GetVectorCovisibleKeyFrames();
    vpMergeConnectedKFs.push_back(mpMergeMatchedKF);
    vpCheckFuseMapPoint.reserve(spMapPointMerge.size());
    std::copy(spMapPointMerge.begin(), spMapPointMerge.end(), std::back_inserter(vpCheckFuseMapPoint));

    // Project MapPoints observed in the neighborhood of the merge keyframe
    // into the current keyframe and neighbors using corrected poses.
    // Fuse duplications.
    SearchAndFuse(vCorrectedSim3, vpCheckFuseMapPoint);

    // Update connectivity
    for(KeyFrame* pKFi : spLocalWindowKFs)
    {
        if(!pKFi || pKFi->isBad())
            continue;

        pKFi->UpdateConnections();
    }
    for(KeyFrame* pKFi : spMergeConnectedKFs)
    {
        if(!pKFi || pKFi->isBad())
            continue;

        pKFi->UpdateConnections();
    }

    bool bStop = false;
    Verbose::PrintMess("MERGE: Start local BA ", Verbose::VERBOSITY_DEBUG);
    vpLocalCurrentWindowKFs.clear();
    vpMergeConnectedKFs.clear();
    std::copy(spLocalWindowKFs.begin(), spLocalWindowKFs.end(), std::back_inserter(vpLocalCurrentWindowKFs));
    std::copy(spMergeConnectedKFs.begin(), spMergeConnectedKFs.end(), std::back_inserter(vpMergeConnectedKFs));
    if (mpTracker->mSensor==System::IMU_MONOCULAR || mpTracker->mSensor==System::IMU_STEREO)
    {
        Optimizer::MergeInertialBA(mpLocalMapper->GetCurrKF(),mpMergeMatchedKF,&bStop, mpCurrentKF->GetMap(),vCorrectedSim3);
    }
    else
    {
        Optimizer::LocalBundleAdjustment(mpCurrentKF, vpLocalCurrentWindowKFs, vpMergeConnectedKFs,&bStop);
    }

    // Loop closed. Release Local Mapping.
    mpLocalMapper->Release();

    Verbose::PrintMess("MERGE: Finish the LBA", Verbose::VERBOSITY_DEBUG);


    ////
    //Update the non critical area from the current map to the merged map
    vector<KeyFrame*> vpCurrentMapKFs = pCurrentMap->GetAllKeyFrames();
    vector<MapPoint*> vpCurrentMapMPs = pCurrentMap->GetAllMapPoints();

    if(vpCurrentMapKFs.size() == 0)
    {
        Verbose::PrintMess("MERGE: There are not KFs outside of the welding area", Verbose::VERBOSITY_DEBUG);
    }
    else
    {
        Verbose::PrintMess("MERGE: Calculate the new position of the elements outside of the window", Verbose::VERBOSITY_DEBUG);
        //Apply the transformation
        {
            if(mpTracker->mSensor == System::MONOCULAR)
            {
                unique_lock<mutex> currentLock(pCurrentMap->mMutexMapUpdate); // We update the current map with the Merge information

                for(KeyFrame* pKFi : vpCurrentMapKFs)
                {
                    if(!pKFi || pKFi->isBad() || pKFi->GetMap() != pCurrentMap)
                    {
                        continue;
                    }

                    g2o::Sim3 g2oCorrectedSiw;

                    cv::Mat Tiw = pKFi->GetPose();
                    cv::Mat Riw = Tiw.rowRange(0,3).colRange(0,3);
                    cv::Mat tiw = Tiw.rowRange(0,3).col(3);
                    g2o::Sim3 g2oSiw(Converter::toMatrix3d(Riw),Converter::toVector3d(tiw),1.0);
                    //Pose without correction
                    vNonCorrectedSim3[pKFi]=g2oSiw;

                    cv::Mat Tic = Tiw*Twc;
                    cv::Mat Ric = Tic.rowRange(0,3).colRange(0,3);
                    cv::Mat tic = Tic.rowRange(0,3).col(3);
                    g2o::Sim3 g2oSim(Converter::toMatrix3d(Ric),Converter::toVector3d(tic),1.0);
                    g2oCorrectedSiw = g2oSim*mg2oMergeScw;
                    vCorrectedSim3[pKFi]=g2oCorrectedSiw;

                    // Update keyframe pose with corrected Sim3. First transform Sim3 to SE3 (scale translation)
                    Eigen::Matrix3d eigR = g2oCorrectedSiw.rotation().toRotationMatrix();
                    Eigen::Vector3d eigt = g2oCorrectedSiw.translation();
                    double s = g2oCorrectedSiw.scale();

                    pKFi->mfScale = s;
                    eigt *=(1./s); //[R t/s;0 1]

                    cv::Mat correctedTiw = Converter::toCvSE3(eigR,eigt);

                    pKFi->mTcwBefMerge = pKFi->GetPose();
                    pKFi->mTwcBefMerge = pKFi->GetPoseInverse();

                    pKFi->SetPose(correctedTiw);

                    if(pCurrentMap->isImuInitialized())
                    {
                        Eigen::Matrix3d Rcor = eigR.transpose()*vNonCorrectedSim3[pKFi].rotation().toRotationMatrix();
                        pKFi->SetVelocity(Converter::toCvMat(Rcor)*pKFi->GetVelocity()); // TODO: should add here scale s
                    }

                }
                for(MapPoint* pMPi : vpCurrentMapMPs)
                {
                    if(!pMPi || pMPi->isBad()|| pMPi->GetMap() != pCurrentMap)
                        continue;

                    KeyFrame* pKFref = pMPi->GetReferenceKeyFrame();
                    g2o::Sim3 g2oCorrectedSwi = vCorrectedSim3[pKFref].inverse();
                    g2o::Sim3 g2oNonCorrectedSiw = vNonCorrectedSim3[pKFref];

                    // Project with non-corrected pose and project back with corrected pose
                    cv::Mat P3Dw = pMPi->GetWorldPos();
                    Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
                    Eigen::Matrix<double,3,1> eigCorrectedP3Dw = g2oCorrectedSwi.map(g2oNonCorrectedSiw.map(eigP3Dw));

                    cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);
                    pMPi->SetWorldPos(cvCorrectedP3Dw);

                    pMPi->UpdateNormalAndDepth();
                }
            }
        }

        mpLocalMapper->RequestStop();
        // Wait until Local Mapping has effectively stopped
        while(!mpLocalMapper->isStopped())
        {
            usleep(1000);
        }

        // Optimize graph (and update the loop position for each element form the begining to the end)
        if(mpTracker->mSensor != System::MONOCULAR)
        {
            Optimizer::OptimizeEssentialGraph(mpCurrentKF, vpMergeConnectedKFs, vpLocalCurrentWindowKFs, vpCurrentMapKFs, vpCurrentMapMPs);
        }


        {
            // Get Merge Map Mutex
            unique_lock<mutex> currentLock(pCurrentMap->mMutexMapUpdate); // We update the current map with the Merge information
            unique_lock<mutex> mergeLock(pMergeMap->mMutexMapUpdate); // We remove the Kfs and MPs in the merged area from the old map

            for(KeyFrame* pKFi : vpCurrentMapKFs)
            {
                if(!pKFi || pKFi->isBad() || pKFi->GetMap() != pCurrentMap)
                {
                    continue;
                }

                // Make sure connections are updated
                pKFi->UpdateMap(pMergeMap);
                pMergeMap->AddKeyFrame(pKFi);
                pCurrentMap->EraseKeyFrame(pKFi);
            }

            for(MapPoint* pMPi : vpCurrentMapMPs)
            {
                if(!pMPi || pMPi->isBad())
                    continue;

                pMPi->UpdateMap(pMergeMap);
                pMergeMap->AddMapPoint(pMPi);
                pCurrentMap->EraseMapPoint(pMPi);
            }
        }
    }

    mpLocalMapper->Release();

    Verbose::PrintMess("MERGE:Completed!!!!!", Verbose::VERBOSITY_DEBUG);

    if(bRelaunchBA && (!pCurrentMap->isImuInitialized() || (pCurrentMap->KeyFramesInMap()<200 && mpAtlas->CountMaps()==1)))
    {
        // Launch a new thread to perform Global Bundle Adjustment
        Verbose::PrintMess("Relaunch Global BA", Verbose::VERBOSITY_DEBUG);
        mbRunningGBA = true;
        mbFinishedGBA = false;
        mbStopGBA = false;
        mpThreadGBA = new thread(&LoopClosing::RunGlobalBundleAdjustment,this, pMergeMap, mpCurrentKF->mnId);
    }

    mpMergeMatchedKF->AddMergeEdge(mpCurrentKF);
    mpCurrentKF->AddMergeEdge(mpMergeMatchedKF);

    pCurrentMap->IncreaseChangeIndex();
    pMergeMap->IncreaseChangeIndex();

    mpAtlas->RemoveBadMaps();

}

void LoopClosing::MergeLocal2()
{
    cout << "Merge detected!!!!" << endl;

    int numTemporalKFs = 11; //TODO (set by parameter): Temporal KFs in the local window if the map is inertial.

    //Relationship to rebuild the essential graph, it is used two times, first in the local window and later in the rest of the map
    KeyFrame* pNewChild;
    KeyFrame* pNewParent;

    vector<KeyFrame*> vpLocalCurrentWindowKFs;
    vector<KeyFrame*> vpMergeConnectedKFs;

    KeyFrameAndPose CorrectedSim3, NonCorrectedSim3;
    // NonCorrectedSim3[mpCurrentKF]=mg2oLoopScw;

    // Flag that is true only when we stopped a running BA, in this case we need relaunch at the end of the merge
    bool bRelaunchBA = false;

    cout << "Check Full Bundle Adjustment" << endl;
    // If a Global Bundle Adjustment is running, abort it
    if(isRunningGBA())
    {
        unique_lock<mutex> lock(mMutexGBA);
        mbStopGBA = true;

        mnFullBAIdx++;

        if(mpThreadGBA)
        {
            mpThreadGBA->detach();
            delete mpThreadGBA;
        }
        bRelaunchBA = true;
    }


    cout << "Request Stop Local Mapping" << endl;
    mpLocalMapper->RequestStop();
    // Wait until Local Mapping has effectively stopped
    while(!mpLocalMapper->isStopped())
    {
        usleep(1000);
    }
    cout << "Local Map stopped" << endl;

    Map* pCurrentMap = mpCurrentKF->GetMap();
    Map* pMergeMap = mpMergeMatchedKF->GetMap();

    {
        float s_on = mSold_new.scale();
        cv::Mat R_on = Converter::toCvMat(mSold_new.rotation().toRotationMatrix());
        cv::Mat t_on = Converter::toCvMat(mSold_new.translation());

        unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);

        mpLocalMapper->EmptyQueue();

        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        bool bScaleVel=false;
        if(s_on!=1)
            bScaleVel=true;
        mpAtlas->GetCurrentMap()->ApplyScaledRotation(R_on,s_on,bScaleVel,t_on);
        mpTracker->UpdateFrameIMU(s_on,mpCurrentKF->GetImuBias(),mpTracker->GetLastKeyFrame());

        std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    }

    const int numKFnew=pCurrentMap->KeyFramesInMap();

    if((mpTracker->mSensor==System::IMU_MONOCULAR || mpTracker->mSensor==System::IMU_STEREO)&& !pCurrentMap->GetIniertialBA2()){
        // Map is not completly initialized
        Eigen::Vector3d bg, ba;
        bg << 0., 0., 0.;
        ba << 0., 0., 0.;
        Optimizer::InertialOptimization(pCurrentMap,bg,ba);
        IMU::Bias b (ba[0],ba[1],ba[2],bg[0],bg[1],bg[2]);
        unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
        mpTracker->UpdateFrameIMU(1.0f,b,mpTracker->GetLastKeyFrame());

        // Set map initialized
        pCurrentMap->SetIniertialBA2();
        pCurrentMap->SetIniertialBA1();
        pCurrentMap->SetImuInitialized();

    }

    // Load KFs and MPs from merge map
    {
        // Get Merge Map Mutex (This section stops tracking!!)
        unique_lock<mutex> currentLock(pCurrentMap->mMutexMapUpdate); // We update the current map with the Merge information
        unique_lock<mutex> mergeLock(pMergeMap->mMutexMapUpdate); // We remove the Kfs and MPs in the merged area from the old map


        vector<KeyFrame*> vpMergeMapKFs = pMergeMap->GetAllKeyFrames();
        vector<MapPoint*> vpMergeMapMPs = pMergeMap->GetAllMapPoints();


        for(KeyFrame* pKFi : vpMergeMapKFs)
        {
            if(!pKFi || pKFi->isBad() || pKFi->GetMap() != pMergeMap)
            {
                continue;
            }

            // Make sure connections are updated
            pKFi->UpdateMap(pCurrentMap);
            pCurrentMap->AddKeyFrame(pKFi);
            pMergeMap->EraseKeyFrame(pKFi);
        }

        for(MapPoint* pMPi : vpMergeMapMPs)
        {
            if(!pMPi || pMPi->isBad() || pMPi->GetMap() != pMergeMap)
                continue;

            pMPi->UpdateMap(pCurrentMap);
            pCurrentMap->AddMapPoint(pMPi);
            pMergeMap->EraseMapPoint(pMPi);
        }

        // Save non corrected poses (already merged maps)
        vector<KeyFrame*> vpKFs = pCurrentMap->GetAllKeyFrames();
        for(KeyFrame* pKFi : vpKFs)
        {
            cv::Mat Tiw=pKFi->GetPose();
            cv::Mat Riw = Tiw.rowRange(0,3).colRange(0,3);
            cv::Mat tiw = Tiw.rowRange(0,3).col(3);
            g2o::Sim3 g2oSiw(Converter::toMatrix3d(Riw),Converter::toVector3d(tiw),1.0);
            NonCorrectedSim3[pKFi]=g2oSiw;
        }
    }

    pMergeMap->GetOriginKF()->SetFirstConnection(false);
    pNewChild = mpMergeMatchedKF->GetParent(); // Old parent, it will be the new child of this KF
    pNewParent = mpMergeMatchedKF; // Old child, now it will be the parent of its own parent(we need eliminate this KF from children list in its old parent)
    mpMergeMatchedKF->ChangeParent(mpCurrentKF);
    while(pNewChild)
    {
        pNewChild->EraseChild(pNewParent); // We remove the relation between the old parent and the new for avoid loop
        KeyFrame * pOldParent = pNewChild->GetParent();
        pNewChild->ChangeParent(pNewParent);
        pNewParent = pNewChild;
        pNewChild = pOldParent;

    }

    vector<MapPoint*> vpCheckFuseMapPoint; // MapPoint vector from current map to allow to fuse duplicated points with the old map (merge)
    vector<KeyFrame*> vpCurrentConnectedKFs;

    mvpMergeConnectedKFs.push_back(mpMergeMatchedKF);
    vector<KeyFrame*> aux = mpMergeMatchedKF->GetVectorCovisibleKeyFrames();
    mvpMergeConnectedKFs.insert(mvpMergeConnectedKFs.end(), aux.begin(), aux.end());
    if (mvpMergeConnectedKFs.size()>6)
        mvpMergeConnectedKFs.erase(mvpMergeConnectedKFs.begin()+6,mvpMergeConnectedKFs.end());

    mpCurrentKF->UpdateConnections();
    vpCurrentConnectedKFs.push_back(mpCurrentKF);
    aux = mpCurrentKF->GetVectorCovisibleKeyFrames();
    vpCurrentConnectedKFs.insert(vpCurrentConnectedKFs.end(), aux.begin(), aux.end());
    if (vpCurrentConnectedKFs.size()>6)
        vpCurrentConnectedKFs.erase(vpCurrentConnectedKFs.begin()+6,vpCurrentConnectedKFs.end());

    set<MapPoint*> spMapPointMerge;
    for(KeyFrame* pKFi : mvpMergeConnectedKFs)
    {
        set<MapPoint*> vpMPs = pKFi->GetMapPoints();
        spMapPointMerge.insert(vpMPs.begin(),vpMPs.end());
        if(spMapPointMerge.size()>1000)
            break;
    }

    vpCheckFuseMapPoint.reserve(spMapPointMerge.size());
    std::copy(spMapPointMerge.begin(), spMapPointMerge.end(), std::back_inserter(vpCheckFuseMapPoint));

    SearchAndFuse(vpCurrentConnectedKFs, vpCheckFuseMapPoint);

    for(KeyFrame* pKFi : vpCurrentConnectedKFs)
    {
        if(!pKFi || pKFi->isBad())
            continue;

        pKFi->UpdateConnections();
    }
    for(KeyFrame* pKFi : mvpMergeConnectedKFs)
    {
        if(!pKFi || pKFi->isBad())
            continue;

        pKFi->UpdateConnections();
    }

    if (numKFnew<10){
        mpLocalMapper->Release();
        return;
    }

    // Perform BA
    bool bStopFlag=false;
    KeyFrame* pCurrKF = mpTracker->GetLastKeyFrame();
    Optimizer::MergeInertialBA(pCurrKF, mpMergeMatchedKF, &bStopFlag, pCurrentMap,CorrectedSim3);

    // Release Local Mapping.
    mpLocalMapper->Release();


    return;
}

void LoopClosing::SearchAndFuse(const KeyFrameAndPose &CorrectedPosesMap, vector<MapPoint*> &vpMapPoints)
{
    ORBmatcher matcher(0.8);

    int total_replaces = 0;

    for(KeyFrameAndPose::const_iterator mit=CorrectedPosesMap.begin(), mend=CorrectedPosesMap.end(); mit!=mend;mit++)
    {
        int num_replaces = 0;
        KeyFrame* pKFi = mit->first;
        Map* pMap = pKFi->GetMap();

        g2o::Sim3 g2oScw = mit->second;
        cv::Mat cvScw = Converter::toCvMat(g2oScw);

        vector<MapPoint*> vpReplacePoints(vpMapPoints.size(),static_cast<MapPoint*>(NULL));
        int numFused = matcher.Fuse(pKFi,cvScw,vpMapPoints,4,vpReplacePoints);

        // Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);
        const int nLP = vpMapPoints.size();
        for(int i=0; i<nLP;i++)
        {
            MapPoint* pRep = vpReplacePoints[i];
            if(pRep)
            {


                num_replaces += 1;
                pRep->Replace(vpMapPoints[i]);

            }
        }

        total_replaces += num_replaces;
    }
}


void LoopClosing::SearchAndFuse(const vector<KeyFrame*> &vConectedKFs, vector<MapPoint*> &vpMapPoints)
{
    ORBmatcher matcher(0.8);

    int total_replaces = 0;

    for(auto mit=vConectedKFs.begin(), mend=vConectedKFs.end(); mit!=mend;mit++)
    {
        int num_replaces = 0;
        KeyFrame* pKF = (*mit);
        Map* pMap = pKF->GetMap();
        cv::Mat cvScw = pKF->GetPose();

        vector<MapPoint*> vpReplacePoints(vpMapPoints.size(),static_cast<MapPoint*>(NULL));
        matcher.Fuse(pKF,cvScw,vpMapPoints,4,vpReplacePoints);

        // Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);
        const int nLP = vpMapPoints.size();
        for(int i=0; i<nLP;i++)
        {
            MapPoint* pRep = vpReplacePoints[i];
            if(pRep)
            {
                num_replaces += 1;
                pRep->Replace(vpMapPoints[i]);
            }
        }
        total_replaces += num_replaces;
    }
}

// 由外部线程调用,请求复位当前线程

void LoopClosing::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        mbResetRequested = true;
    }

    while(1)
    {
        {
        unique_lock<mutex> lock2(mMutexReset);
        if(!mbResetRequested)
            break;
        }
        usleep(5000);
    }
}

void LoopClosing::RequestResetActiveMap(Map *pMap)
{
    {
        unique_lock<mutex> lock(mMutexReset);
        mbResetActiveMapRequested = true;
        mpMapToReset = pMap;
    }

    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetActiveMapRequested)
                break;
        }
        usleep(3000);
    }
}
// 当前线程调用,检查是否有外部线程请求复位当前线程,如果有的话就复位回环检测线程
void LoopClosing::ResetIfRequested()
{
    unique_lock<mutex> lock(mMutexReset);
    // 如果有来自于外部的线程的复位请求,那么就复位当前线程
    if(mbResetRequested)
    {
        cout << "Loop closer reset requested..." << endl;
        mlpLoopKeyFrameQueue.clear();   // 清空参与和进行回环检测的关键帧队列
        mLastLoopKFid=0;                // 上一次没有和任何关键帧形成闭环关系
        mbResetRequested=false;         // 复位请求标志复位
        mbResetActiveMapRequested = false;
    }
    else if(mbResetActiveMapRequested)
    {

        for (list<KeyFrame*>::const_iterator it=mlpLoopKeyFrameQueue.begin(); it != mlpLoopKeyFrameQueue.end();)
        {
            KeyFrame* pKFi = *it;
            if(pKFi->GetMap() == mpMapToReset)
            {
                it = mlpLoopKeyFrameQueue.erase(it);
            }
            else
                ++it;
        }

        mLastLoopKFid=mpAtlas->GetLastInitKFid(); //TODO old variable, it is not use in the new algorithm
        mbResetActiveMapRequested=false;

    }
}

void LoopClosing::RunGlobalBundleAdjustment(Map* pActiveMap, unsigned long nLoopKF)
{
    Verbose::PrintMess("Starting Global Bundle Adjustment", Verbose::VERBOSITY_NORMAL);

    const bool bImuInit = pActiveMap->isImuInitialized();

#ifdef REGISTER_TIMES
    std::chrono::steady_clock::time_point time_StartFGBA = std::chrono::steady_clock::now();
#endif

    if(!bImuInit)
        Optimizer::GlobalBundleAdjustemnt(pActiveMap,10,&mbStopGBA,nLoopKF,false);
    else
        Optimizer::FullInertialBA(pActiveMap,7,false,nLoopKF,&mbStopGBA);

#ifdef REGISTER_TIMES
    std::chrono::steady_clock::time_point time_StartMapUpdate = std::chrono::steady_clock::now();

    double timeFullGBA = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_StartMapUpdate - time_StartFGBA).count();
    vTimeFullGBA_ms.push_back(timeFullGBA);
#endif


    int idx =  mnFullBAIdx;
    // Optimizer::GlobalBundleAdjustemnt(mpMap,10,&mbStopGBA,nLoopKF,false);

    // Update all MapPoints and KeyFrames
    // Local Mapping was active during BA, that means that there might be new keyframes
    // not included in the Global BA and they are not consistent with the updated map.
    // We need to propagate the correction through the spanning tree
    {
        unique_lock<mutex> lock(mMutexGBA);
        if(idx!=mnFullBAIdx)
            return;

        if(!bImuInit && pActiveMap->isImuInitialized())
            return;

        if(!mbStopGBA)
        {
            Verbose::PrintMess("Global Bundle Adjustment finished", Verbose::VERBOSITY_NORMAL);
            Verbose::PrintMess("Updating map ...", Verbose::VERBOSITY_NORMAL);

            mpLocalMapper->RequestStop();
            // Wait until Local Mapping has effectively stopped

            while(!mpLocalMapper->isStopped() && !mpLocalMapper->isFinished())
            {
                usleep(1000);
            }

            // Get Map Mutex
            unique_lock<mutex> lock(pActiveMap->mMutexMapUpdate);
            // cout << "LC: Update Map Mutex adquired" << endl;

            //pActiveMap->PrintEssentialGraph();
            // Correct keyframes starting at map first keyframe
            list<KeyFrame*> lpKFtoCheck(pActiveMap->mvpKeyFrameOrigins.begin(),pActiveMap->mvpKeyFrameOrigins.end());

            while(!lpKFtoCheck.empty())
            {
                KeyFrame* pKF = lpKFtoCheck.front();
                const set<KeyFrame*> sChilds = pKF->GetChilds();
                //cout << "---Updating KF " << pKF->mnId << " with " << sChilds.size() << " childs" << endl;
                //cout << " KF mnBAGlobalForKF: " << pKF->mnBAGlobalForKF << endl;
                cv::Mat Twc = pKF->GetPoseInverse();
                //cout << "Twc: " << Twc << endl;
                //cout << "GBA: Correct KeyFrames" << endl;
                for(set<KeyFrame*>::const_iterator sit=sChilds.begin();sit!=sChilds.end();sit++)
                {
                    KeyFrame* pChild = *sit;
                    if(!pChild || pChild->isBad())
                        continue;

                    if(pChild->mnBAGlobalForKF!=nLoopKF)
                    {
                        //cout << "++++New child with flag " << pChild->mnBAGlobalForKF << "; LoopKF: " << nLoopKF << endl;
                        //cout << " child id: " << pChild->mnId << endl;
                        cv::Mat Tchildc = pChild->GetPose()*Twc;
                        //cout << "Child pose: " << Tchildc << endl;
                        //cout << "pKF->mTcwGBA: " << pKF->mTcwGBA << endl;
                        pChild->mTcwGBA = Tchildc*pKF->mTcwGBA;//*Tcorc*pKF->mTcwGBA;

                        cv::Mat Rcor = pChild->mTcwGBA.rowRange(0,3).colRange(0,3).t()*pChild->GetRotation();
                        if(!pChild->GetVelocity().empty()){
                            //cout << "Child velocity: " << pChild->GetVelocity() << endl;
                            pChild->mVwbGBA = Rcor*pChild->GetVelocity();
                        }
                        else
                            Verbose::PrintMess("Child velocity empty!! ", Verbose::VERBOSITY_NORMAL);


                        //cout << "Child bias: " << pChild->GetImuBias() << endl;
                        pChild->mBiasGBA = pChild->GetImuBias();


                        pChild->mnBAGlobalForKF=nLoopKF;

                    }
                    lpKFtoCheck.push_back(pChild);
                }

                //cout << "-------Update pose" << endl;
                pKF->mTcwBefGBA = pKF->GetPose();
                //cout << "pKF->mTcwBefGBA: " << pKF->mTcwBefGBA << endl;
                pKF->SetPose(pKF->mTcwGBA);
                /*cv::Mat Tco_cn = pKF->mTcwBefGBA * pKF->mTcwGBA.inv();
                cv::Vec3d trasl = Tco_cn.rowRange(0,3).col(3);
                double dist = cv::norm(trasl);
                cout << "GBA: KF " << pKF->mnId << " had been moved " << dist << " meters" << endl;
                double desvX = 0;
                double desvY = 0;
                double desvZ = 0;
                if(pKF->mbHasHessian)
                {
                    cv::Mat hessianInv = pKF->mHessianPose.inv();

                    double covX = hessianInv.at<double>(3,3);
                    desvX = std::sqrt(covX);
                    double covY = hessianInv.at<double>(4,4);
                    desvY = std::sqrt(covY);
                    double covZ = hessianInv.at<double>(5,5);
                    desvZ = std::sqrt(covZ);
                    pKF->mbHasHessian = false;
                }
                if(dist > 1)
                {
                    cout << "--To much distance correction: It has " << pKF->GetConnectedKeyFrames().size() << " connected KFs" << endl;
                    cout << "--It has " << pKF->GetCovisiblesByWeight(80).size() << " connected KF with 80 common matches or more" << endl;
                    cout << "--It has " << pKF->GetCovisiblesByWeight(50).size() << " connected KF with 50 common matches or more" << endl;
                    cout << "--It has " << pKF->GetCovisiblesByWeight(20).size() << " connected KF with 20 common matches or more" << endl;

                    cout << "--STD in meters(x, y, z): " << desvX << ", " << desvY << ", " << desvZ << endl;


                    string strNameFile = pKF->mNameFile;
                    cv::Mat imLeft = cv::imread(strNameFile, CV_LOAD_IMAGE_UNCHANGED);

                    cv::cvtColor(imLeft, imLeft, CV_GRAY2BGR);

                    vector<MapPoint*> vpMapPointsKF = pKF->GetMapPointMatches();
                    int num_MPs = 0;
                    for(int i=0; i<vpMapPointsKF.size(); ++i)
                    {
                        if(!vpMapPointsKF[i] || vpMapPointsKF[i]->isBad())
                        {
                            continue;
                        }
                        num_MPs += 1;
                        string strNumOBs = to_string(vpMapPointsKF[i]->Observations());
                        cv::circle(imLeft, pKF->mvKeys[i].pt, 2, cv::Scalar(0, 255, 0));
                        cv::putText(imLeft, strNumOBs, pKF->mvKeys[i].pt, CV_FONT_HERSHEY_DUPLEX, 1, cv::Scalar(255, 0, 0));
                    }
                    cout << "--It has " << num_MPs << " MPs matched in the map" << endl;

                    string namefile = "./test_GBA/GBA_" + to_string(nLoopKF) + "_KF" + to_string(pKF->mnId) +"_D" + to_string(dist) +".png";
                    cv::imwrite(namefile, imLeft);
                }*/


                if(pKF->bImu)
                {
                    //cout << "-------Update inertial values" << endl;
                    pKF->mVwbBefGBA = pKF->GetVelocity();
                    if (pKF->mVwbGBA.empty())
                        Verbose::PrintMess("pKF->mVwbGBA is empty", Verbose::VERBOSITY_NORMAL);

                    assert(!pKF->mVwbGBA.empty());
                    pKF->SetVelocity(pKF->mVwbGBA);
                    pKF->SetNewBias(pKF->mBiasGBA);                    
                }

                lpKFtoCheck.pop_front();
            }

            //cout << "GBA: Correct MapPoints" << endl;
            // Correct MapPoints
            const vector<MapPoint*> vpMPs = pActiveMap->GetAllMapPoints();

            // 遍历每一个地图点
            for(size_t i=0; i<vpMPs.size(); i++)
            {
                MapPoint* pMP = vpMPs[i];

                if(pMP->isBad())
                    continue;

                // NOTICE 并不是所有的地图点都会直接参与到全局BA优化中,但是大部分的地图点需要根据全局BA优化后的结果来重新纠正自己的位姿
                // 如果这个地图点直接参与到了全局BA优化的过程,那么就直接重新设置器位姿即可
                if(pMP->mnBAGlobalForKF==nLoopKF)
                {
                    // If optimized by Global BA, just update
                    pMP->SetWorldPos(pMP->mPosGBA);
                }
                else // 如故这个地图点并没有直接参与到全局BA优化的过程中,那么就使用器参考关键帧的新位姿来优化自己的位姿
                {
                    // Update according to the correction of its reference keyframe
                    KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();

                    // 说明这个关键帧，在前面的过程中也没有因为“当前关键帧”得到全局BA优化 
                    //? 可是,为什么会出现这种情况呢? 难道是因为这个地图点的参考关键帧设置成为了bad?
                    if(pRefKF->mnBAGlobalForKF!=nLoopKF)
                        continue;

                    if(pRefKF->mTcwBefGBA.empty())
                        continue;

                    // Map to non-corrected camera
                    cv::Mat Rcw = pRefKF->mTcwBefGBA.rowRange(0,3).colRange(0,3);
                    cv::Mat tcw = pRefKF->mTcwBefGBA.rowRange(0,3).col(3);
                    // 转换到其参考关键帧相机坐标系下的坐标
                    cv::Mat Xc = Rcw*pMP->GetWorldPos()+tcw;

                    // Backproject using corrected camera
                    // 然后使用已经纠正过的参考关键帧的位姿,再将该地图点变换到世界坐标系下
                    cv::Mat Twc = pRefKF->GetPoseInverse();
                    cv::Mat Rwc = Twc.rowRange(0,3).colRange(0,3);
                    cv::Mat twc = Twc.rowRange(0,3).col(3);

                    pMP->SetWorldPos(Rwc*Xc+twc);
                }
            }

            pActiveMap->InformNewBigChange();
            pActiveMap->IncreaseChangeIndex();

            // TODO Check this update
            // mpTracker->UpdateFrameIMU(1.0f, mpTracker->GetLastKeyFrame()->GetImuBias(), mpTracker->GetLastKeyFrame());

            mpLocalMapper->Release();

            Verbose::PrintMess("Map updated!", Verbose::VERBOSITY_NORMAL);
        }

        mbFinishedGBA = true;
        mbRunningGBA = false;
    }

#ifdef REGISTER_TIMES
    std::chrono::steady_clock::time_point time_EndMapUpdate = std::chrono::steady_clock::now();

    double timeMapUpdate = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMapUpdate - time_StartMapUpdate).count();
    vTimeMapUpdate_ms.push_back(timeMapUpdate);

    double timeGBA = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMapUpdate - time_StartFGBA).count();
    vTimeGBATotal_ms.push_back(timeGBA);
#endif
}

// 由外部线程调用,请求终止当前线程
void LoopClosing::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    // cout << "LC: Finish requested" << endl;
    mbFinishRequested = true;
}

// 当前线程调用,查看是否有外部线程请求当前线程
bool LoopClosing::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

// 有当前线程调用,执行完成该函数之后线程主函数退出,线程销毁
void LoopClosing::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;
}

// 由外部线程调用,判断当前回环检测线程是否已经正确终止了
bool LoopClosing::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}


} //namespace ORB_SLAM
