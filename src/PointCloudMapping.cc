/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2016  <copyright holder> <email>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "PointCloudMapping.h"
#include <KeyFrame.h>
#include <opencv2/highgui/highgui.hpp>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include "Converter.h"
#include "System.h"

#include<sys/time.h>
namespace ORB_SLAM3
{
// int currentloopcount = 0;
PointCloudMapping::PointCloudMapping(double resolution_, double meank_, double thresh_)
    : mabIsUpdating(false)
{
    this->resolution = resolution_;
    this->meank = meank_;
    this->thresh = thresh_;
    std::cout<<resolution<<" "<<meank<<" "<<thresh<<std::endl;
    statistical_filter = new pcl::StatisticalOutlierRemoval<pcl::PointXYZRGBA>(true);
    voxel = new pcl::VoxelGrid<pcl::PointXYZRGBA>();
    statistical_filter->setMeanK(meank);
    statistical_filter->setStddevMulThresh(thresh);
    voxel->setLeafSize(resolution, resolution, resolution);
    globalMap = pcl::PointCloud<pcl::PointXYZRGBA>::Ptr(new pcl::PointCloud<pcl::PointXYZRGBA>);

    viewerThread = make_shared<thread>(bind(&PointCloudMapping::viewer, this));
}

void PointCloudMapping::shutdown()
{
    {
        unique_lock<mutex> lck(shutDownMutex);
        shutDownFlag = true;
        keyFrameUpdated.notify_one();
    }
    viewerThread->join();
}

void PointCloudMapping::Clear()
{
    std::cout << "清除稠密地图" << std::endl;
    std::unique_lock<std::mutex> lck(mMutexGlobalMap);
    globalMap.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
}
void PointCloudMapping::insertKeyFrame(KeyFrame *kf, cv::Mat &color, cv::Mat &depth, int idk, vector<KeyFrame *> vpKFs)
{
    // cout << "receive a keyframe, 第" << kf->mnId << "个" << endl;
    // if (color.empty())
    //     return;
    // unique_lock<mutex> lck(keyframeMutex);
    // keyframes.push_back(kf);
    // currentvpKFs = vpKFs;
    // PointCloude pointcloude;
    // pointcloude.pcID = idk;
    // pointcloude.T = ORB_SLAM3::Converter::toSE3Quat(kf->GetPose());
    // pointcloude.pcE = generatePointCloud(kf, color, depth);
    // kf->mptrPointCloud = pointcloude.pcE;
    // pointcloud.push_back(pointcloude);
    // keyFrameUpdated.notify_one();
}

void PointCloudMapping::insertKeyFrame(KeyFrame *kf)
{
    // cout << "receive a keyframe, 第" << kf->mnId << "个" << endl;
    if (kf->imLeftRgb.empty())
        return;
    unique_lock<mutex> lck(keyframeMutex);
    mlNewKeyFrames.emplace_back(kf);
    if(mlNewKeyFrames.size() > 35)
        mlNewKeyFrames.pop_front();

}

void PointCloudMapping::generatePointCloud(KeyFrame *kf) //,Eigen::Isometry3d T
{
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr pPointCloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
    // point cloud is null ptr
    for (int m = 0; m < kf->imDepth.rows; m += 4)
    {
        for (int n = 0; n < kf->imDepth.cols; n += 4)
        {
            float d = kf->imDepth.ptr<float>(m)[n];
            if (d < 0.05 || d > 6)
                continue;
            pcl::PointXYZRGBA p;
            p.z = d;
            p.x = (n - kf->cx) * p.z / kf->fx;
            p.y = (m - kf->cy) * p.z / kf->fy;

            p.b = kf->imLeftRgb.ptr<uchar>(m)[n * 3];
            p.g = kf->imLeftRgb.ptr<uchar>(m)[n * 3 + 1];
            p.r = kf->imLeftRgb.ptr<uchar>(m)[n * 3 + 2];

            pPointCloud->points.push_back(p);
        }
    }
    pPointCloud->height = 1;
    pPointCloud->width = pPointCloud->points.size();
    pPointCloud->is_dense = true;
    kf->mptrPointCloud = pPointCloud;
}

void PointCloudMapping::viewer()
{
    pcl::visualization::CloudViewer viewer("viewer");
    // KeyFrame * pCurKF;
    while (1)
    {
        {
            unique_lock<mutex> lck_shutdown(shutDownMutex);
            if (shutDownFlag)
            {
                break;
            }
        }
        // std::cout<<"sssss"<<std::endl;
        if (bStop || mabIsUpdating)
        {
            //cout<<"loopbusy || bStop"<<endl;
            continue;
        }
        // std::cout<<"1111111"<<std::endl;
        int N;
        std::list<KeyFrame *> lNewKeyFrames;
        {
            unique_lock<mutex> lck(keyframeMutex);
            N = mlNewKeyFrames.size();
            lNewKeyFrames = mlNewKeyFrames;
            if(N == 0)
                continue;
            else
            {
                mlNewKeyFrames.clear();
            }

        }
        // timeval start, finish; //定义开始，结束变量
        //初始化
        // cout<<"待处理点云个数 = "<<N<<endl;
        double generatePointCloudTime = 0, transformPointCloudTime = 0; 
        for (auto pKF : lNewKeyFrames)
        {
            if (pKF->isBad())
                continue;
            // gettimeofday(&start,NULL);
            generatePointCloud(pKF);
            // gettimeofday(&finish,NULL);//初始化结束时间
            // generatePointCloudTime += finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec)/1000000.0;

            // gettimeofday(&start,NULL);
            pcl::PointCloud<pcl::PointXYZRGBA>::Ptr p(new pcl::PointCloud<pcl::PointXYZRGBA>);
            pcl::transformPointCloud(*(pKF->mptrPointCloud), *(p), Converter::toMatrix4d(pKF->GetPoseInverse()));

            {
                std::unique_lock<std::mutex> lck(mMutexGlobalMap);
                *globalMap += *p;
            }
            // gettimeofday(&finish,NULL);//初始化结束时间
            // transformPointCloudTime += finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec)/1000000.0;
        }
        // gettimeofday(&start,NULL);
        pcl::PointCloud<pcl::PointXYZRGBA>::Ptr tmp(new pcl::PointCloud<pcl::PointXYZRGBA>);

        // 去除孤立点这个比较耗时，用处也不是很大，可以去掉
        // statistical_filter->setInputCloud(globalMap);  
        // statistical_filter->filter(*tmp);

        voxel->setInputCloud(globalMap);
        voxel->filter(*globalMap);
        // gettimeofday(&finish,NULL);//初始化结束时间
        // double filter = finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec)/1000000.0;//转换浮点型
        // std::cout<<"filter: "<<filter<<std::endl;

        // std::cout<<"generatePointCloudTime: "<<generatePointCloudTime<<std::endl;
        // std::cout<<"transformPointCloudTime: "<<transformPointCloudTime<<std::endl;
        // gettimeofday(&start,NULL);
        viewer.showCloud(globalMap);  // 这个比较费时，建议不需要实时显示的可以屏蔽或改成几次显示一次
        // gettimeofday(&finish,NULL);//初始化结束时间
        // double duration = finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec)/1000000.0;//转换浮点型
        // std::cout<<"showCloud: "<<duration<<std::endl;
    }
}

// 保存地图的函数，需要的自行调用~
void PointCloudMapping::save()
{
    std::unique_lock<std::mutex> lck(mMutexGlobalMap);
    pcl::io::savePCDFile("result.pcd", *globalMap);
    cout << "globalMap save finished" << endl;
}
void PointCloudMapping::updatecloud(Map &curMap)
{
    std::unique_lock<std::mutex> lck(updateMutex);
    
    mabIsUpdating = true;
    currentvpKFs = curMap.GetAllKeyFrames();
    // loopbusy = true;
    cout << "开始点云更新" << endl;
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr tmpGlobalMap(new pcl::PointCloud<pcl::PointXYZRGBA>);
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr curPointCloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr tmpGlobalMapFilter(new pcl::PointCloud<pcl::PointXYZRGBA>());
    for (int i = 0; i < currentvpKFs.size(); i++)
    {
        if (!mabIsUpdating)
        {
            std::cout << "中断点云更新" <<std::endl;
            return;
        }
        if (!currentvpKFs[i]->isBad() && currentvpKFs[i]->mptrPointCloud)
        {
            
            pcl::transformPointCloud(
                *(currentvpKFs[i]->mptrPointCloud), *(curPointCloud),
                Converter::toMatrix4d(currentvpKFs[i]->GetPoseInverse()));
            *tmpGlobalMap += *curPointCloud;

            voxel->setInputCloud(tmpGlobalMap);
            voxel->filter(*tmpGlobalMapFilter);
            tmpGlobalMap->swap(*tmpGlobalMapFilter);
        }
    }
    cout << "点云更新完成" << endl;
    {
        std::unique_lock<std::mutex> lck(mMutexGlobalMap);
        globalMap = tmpGlobalMap;
    }
    mabIsUpdating = false;
}
}