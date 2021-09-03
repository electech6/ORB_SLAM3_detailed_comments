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

#include<iostream>
#include<algorithm>
#include<fstream>
#include<chrono>
#include <ctime>
#include <sstream>

#include<opencv2/core/core.hpp>

#include<System.h>
#include "ImuTypes.h"

using namespace std;

void LoadImages(const string &strImagePath, const string &strPathTimes,
                vector<string> &vstrImages, vector<double> &vTimeStamps);

void LoadIMU(const string &strImuPath, vector<double> &vTimeStamps, vector<cv::Point3f> &vAcc, vector<cv::Point3f> &vGyro);

double ttrack_tot = 0;
int main(int argc, char **argv)
{
    // 输出运行的序列数目
    const int num_seq = (argc-3)/3;
    cout << "num_seq = " << num_seq << endl;
    bool bFileName= ((argc % 3) == 1);

    string file_name;
    if (bFileName)
        file_name = string(argv[argc-1]);

    cout << "file name: " << file_name << endl;

    // 按照下面提示至少输入6个参数
    if(argc < 6)
    {
        cerr << endl << "Usage: ./mono_inertial_tum_vi path_to_vocabulary path_to_settings path_to_image_folder_1 path_to_times_file_1 path_to_imu_data_1 (path_to_image_folder_2 path_to_times_file_2 path_to_imu_data_2 ... path_to_image_folder_N path_to_times_file_N path_to_imu_data_N) (trajectory_file_name)" << endl;
        return 1;
    }


    // Load all sequences:
    // 准备加载所有序列的数据
    int seq;
    vector< vector<string> > vstrImageFilenames;    //图像文件名
    vector< vector<double> > vTimestampsCam;        //图像时间戳
    vector< vector<cv::Point3f> > vAcc, vGyro;      //加速度计，陀螺仪
    vector< vector<double> > vTimestampsImu;        //IMU时间戳
    vector<int> nImages;                            //图像序列    
    vector<int> nImu;
    vector<int> first_imu(num_seq,0);               //记录和第一帧图像时间戳最接近的imu时间戳索引

    vstrImageFilenames.resize(num_seq);
    vTimestampsCam.resize(num_seq);
    vAcc.resize(num_seq);
    vGyro.resize(num_seq);
    vTimestampsImu.resize(num_seq);
    nImages.resize(num_seq);
    nImu.resize(num_seq);

    int tot_images = 0;
    // 遍历每个序列
    for (seq = 0; seq<num_seq; seq++)
    {
        // Step 1 加载图像名和对应的图像时间戳
        cout << "Loading images for sequence " << seq << "...";
        LoadImages(string(argv[3*(seq+1)]), string(argv[3*(seq+1)+1]), vstrImageFilenames[seq], vTimestampsCam[seq]);
        cout << "LOADED!" << endl;

        // Step 2 加载IMU数据
        cout << "Loading IMU for sequence " << seq << "...";
        LoadIMU(string(argv[3*(seq+1)+2]), vTimestampsImu[seq], vAcc[seq], vGyro[seq]);
        cout << "LOADED!" << endl;

        nImages[seq] = vstrImageFilenames[seq].size();
        tot_images += nImages[seq];
        nImu[seq] = vTimestampsImu[seq].size();

        //检查是否存在有效数目的图像和imu数据
        if((nImages[seq]<=0)||(nImu[seq]<=0))
        {
            cerr << "ERROR: Failed to load images or IMU for sequence" << seq << endl;
            return 1;
        }

        // Find first imu to be considered, supposing imu measurements start first
        // Step 3 默认IMU数据早于图像数据记录，找到和第一帧图像时间戳最接近的imu时间戳索引，记录在first_imu[seq]中
        while(vTimestampsImu[seq][first_imu[seq]]<=vTimestampsCam[seq][0]){
            first_imu[seq]++;
        }
        // 因为上面退出while循环时IMU时间戳刚刚超过图像时间戳，所以这里需要再减一个索引    
        first_imu[seq]--; // first imu measurement to be considered

    }

    // Vector for tracking time statistics
    vector<float> vTimesTrack;
    vTimesTrack.resize(tot_images);

    cout << endl << "-------" << endl;
    cout.precision(17);

    /*cout << "Start processing sequence ..." << endl;
    cout << "Images in the sequence: " << nImages << endl;
    cout << "IMU data in the sequence: " << nImu << endl << endl;*/

    // Create SLAM system. It initializes all system threads and gets ready to process frames.
    // Step 4 SLAM系统的初始化，包括读取配置文件、字典，创建跟踪、局部建图、闭环、显示线程
    ORB_SLAM3::System SLAM(argv[1],argv[2],ORB_SLAM3::System::IMU_MONOCULAR, true, 0, file_name);

    //遍历所有数据
    int proccIm = 0;
    for (seq = 0; seq<num_seq; seq++)
    {
        // Main loop
        cv::Mat im;
        //存放imu数据容器,包含该加速度,角速度,时间戳
        vector<ORB_SLAM3::IMU::Point> vImuMeas;
        proccIm = 0;
        //直方图均衡化,直方图均衡化的思想就是这样的:
        //假设我有灰度级255的图像，但是都是属于［100，110］的灰度，图像对比度就很低，我应该尽可能拉到整个［0，255］
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8, 8));
        for(int ni=0; ni<nImages[seq]; ni++, proccIm++)
        {
            // Read image from file
            // Step 5 读取每一帧图像并转换为灰度图存储在im,seq表示第几个数据集,ni表示这个数据集的第几个数据
            im = cv::imread(vstrImageFilenames[seq][ni],cv::IMREAD_GRAYSCALE);

            // clahe
            //直方图均衡化
            clahe->apply(im,im);


            // 取出对应的图像时间戳
            double tframe = vTimestampsCam[seq][ni];

            if(im.empty())
            {
                cerr << endl << "Failed to load image at: "
                     <<  vstrImageFilenames[seq][ni] << endl;
                return 1;
            }


            // Load imu measurements from previous frame
            //清空imu测量
            vImuMeas.clear();

            if(ni>0)
            {
                // cout << "t_cam " << tframe << endl;
                // Step 6 把上一图像帧和当前图像帧之间的imu信息存储在vImuMeas里
                // 注意第一个图像帧没有对应的imu数据 
                // seq: 数据集序列索引；first_imu[seq]：当前数据集中和当前帧图像最接近的imu时间戳索引；ni：图像的索引
                while(vTimestampsImu[seq][first_imu[seq]]<=vTimestampsCam[seq][ni])
                {
                    // 存储IMU的加速度计信息vAcc、陀螺仪信息vGyro、时间戳信息vTimestampsImu
                    vImuMeas.push_back(ORB_SLAM3::IMU::Point(vAcc[seq][first_imu[seq]].x,vAcc[seq][first_imu[seq]].y,vAcc[seq][first_imu[seq]].z,
                                                             vGyro[seq][first_imu[seq]].x,vGyro[seq][first_imu[seq]].y,vGyro[seq][first_imu[seq]].z,
                                                             vTimestampsImu[seq][first_imu[seq]]));
                    // cout << "t_imu = " << fixed << vImuMeas.back().t << endl;
                    // 更新
                    first_imu[seq]++;
                }
            }

            // cout << "first imu: " << first_imu[seq] << endl;
            /*cout << "first imu time: " << fixed << vTimestampsImu[first_imu] << endl;
            cout << "size vImu: " << vImuMeas.size() << endl;*/
    #ifdef COMPILEDWITHC11
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    #else
            std::chrono::monotonic_clock::time_point t1 = std::chrono::monotonic_clock::now();
    #endif

            // Pass the image to the SLAM system
            // cout << "tframe = " << tframe << endl;
            // Step 7 跟踪线程作为主线程运行
            SLAM.TrackMonocular(im,tframe,vImuMeas); // TODO change to monocular_inertial

    #ifdef COMPILEDWITHC11
            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    #else
            std::chrono::monotonic_clock::time_point t2 = std::chrono::monotonic_clock::now();
    #endif

            double ttrack= std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
            ttrack_tot += ttrack;
            // std::cout << "ttrack: " << ttrack << std::endl;

            vTimesTrack[ni]=ttrack;

            // Wait to load the next frame
            // 等待读取下一帧
            double T=0;
            if(ni<nImages[seq]-1)
                T = vTimestampsCam[seq][ni+1]-tframe;
            else if(ni>0)
                T = tframe-vTimestampsCam[seq][ni-1];

            if(ttrack<T)
                usleep((T-ttrack)*1e6); // 1e6

        }
        if(seq < num_seq - 1)
        {
            cout << "Changing the dataset" << endl;
            // Step 8 更换数据集 
            SLAM.ChangeDataset();
        }

    }

    // cout << "ttrack_tot = " << ttrack_tot << std::endl;
    // Stop all threads
    // Step 9 关闭SLAM中所有线程
    SLAM.Shutdown();


    // Tracking time statistics

    // Save camera trajectory
    // Step 10 保存相机位姿（轨迹）
    if (bFileName)
    {
        const string kf_file =  "kf_" + string(argv[argc-1]) + ".txt";
        const string f_file =  "f_" + string(argv[argc-1]) + ".txt";
        SLAM.SaveTrajectoryEuRoC(f_file);
        SLAM.SaveKeyFrameTrajectoryEuRoC(kf_file);
    }
    else
    {
        SLAM.SaveTrajectoryEuRoC("CameraTrajectory.txt");
        SLAM.SaveKeyFrameTrajectoryEuRoC("KeyFrameTrajectory.txt");
    }

    sort(vTimesTrack.begin(),vTimesTrack.end());
    float totaltime = 0;
    for(int ni=0; ni<nImages[0]; ni++)
    {
        totaltime+=vTimesTrack[ni];
    }
    cout << "-------" << endl << endl;
    cout << "median tracking time: " << vTimesTrack[nImages[0]/2] << endl;
    cout << "mean tracking time: " << totaltime/proccIm << endl;

    /*const string kf_file =  "kf_" + ss.str() + ".txt";
    const string f_file =  "f_" + ss.str() + ".txt";

    SLAM.SaveTrajectoryEuRoC(f_file);
    SLAM.SaveKeyFrameTrajectoryEuRoC(kf_file);*/

    return 0;
}

/**
 * @brief 加载图像及其时间戳
 * 
 * @param[in] strImagePath          图像路径名称
 * @param[in] strPathTimes          图像时间戳文件名
 * @param[in] vstrImages            图像路径及名称
 * @param[in] vTimeStamps           图像时间戳和vstrImages一一对应
 */
void LoadImages(const string &strImagePath, const string &strPathTimes,
                vector<string> &vstrImages, vector<double> &vTimeStamps)
{
    ifstream fTimes;
    cout << strImagePath << endl;
    cout << strPathTimes << endl;
    fTimes.open(strPathTimes.c_str());
    vTimeStamps.reserve(5000);
    vstrImages.reserve(5000);
    while(!fTimes.eof())
    {
        string s;
        getline(fTimes,s);
        if(!s.empty())
        {
            stringstream ss;
            ss << s;
            vstrImages.push_back(strImagePath + "/" + ss.str() + ".png");
            double t;
            ss >> t;
            vTimeStamps.push_back(t/1e9);

        }
    }
}

/**
 * @brief 加载IMU数据
 * 
 * @param[in] strImuPath        IMU文件路径
 * @param[in] vTimeStamps       IMU时间戳
 * @param[in] vAcc              加速度计数据
 * @param[in] vGyro             陀螺仪数据
 */
void LoadIMU(const string &strImuPath, vector<double> &vTimeStamps, vector<cv::Point3f> &vAcc, vector<cv::Point3f> &vGyro)
{
    ifstream fImu;
    fImu.open(strImuPath.c_str());
    // 预申请大小为5000的vector，不够可以再扩展
    vTimeStamps.reserve(5000);
    vAcc.reserve(5000);
    vGyro.reserve(5000);

    while(!fImu.eof())
    {
        string s;
        getline(fImu,s);
        // 跳过注释或无效数据
        if (s[0] == '#')
            continue;

        if(!s.empty())
        {
            string item;
            size_t pos = 0;
            double data[7];
            int count = 0;
            while ((pos = s.find(',')) != string::npos) {
                item = s.substr(0, pos);
                data[count++] = stod(item);
                s.erase(0, pos + 1);
            }
            item = s.substr(0, pos);
            data[6] = stod(item);

            // 注意这里的时间戳除以10的9次方转化为秒单位
            vTimeStamps.push_back(data[0]/1e9);
            vAcc.push_back(cv::Point3f(data[4],data[5],data[6]));
            vGyro.push_back(cv::Point3f(data[1],data[2],data[3]));
        }
    }
}
